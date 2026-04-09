"""
spectral_unmixing.py|
--------------------
Unmixing espectral de mu_a recuperado por IAD.

Modelo: 5 parámetros libres
    theta = (HbT [µM], SO2 [0-1], r_vess [µm], f_v_mel [0-1], k_mel [nm⁻¹])

Cromóforos:
    - HbO2 / Hb   → Prahl 1998, OMLC (cm⁻¹/M, conv. log10)
    - Melanina (melanosoma):
      µa_melanosoma(λ) = 519 mm⁻¹ · exp(−k_mel · (λ − 500 nm))
      Modelo exponencial anclado en 500 nm (Zonios et al. 2008, JBO 13(1):014017).
      Reemplaza la ley de potencia de Jacques 2013, cuyo rango de validez
      (650–800 nm) es inadecuado para el visible corto (423–680 nm).

Empaquetamiento vascular: van Veen et al. 2002, Opt. Lett. 27(4):246-248.

Unidades internas: mm⁻¹ (igual que el output de IAD).
"""

import numpy as np
import os
import re
import argparse
import pandas as pd
from pathlib import Path
import unicodedata
from scipy.interpolate import interp1d
from scipy.optimize import least_squares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ─────────────────────────────────────────────
# 0.  CONFIGURACIÓN GLOBAL DEL AJUSTE FISIOLÓGICO
# ─────────────────────────────────────────────
# Bounds fisiológicos del ajuste.
# Modifica aquí los límites de cada variable antes de ejecutar el script.
# Orden del solver:
# [HbT (µM), SO2, r_vess (µm), f_v_melanosome, k_mel (nm⁻¹)]
HBT_BOUNDS_UM = (10.0, 100.0) #Cambridge Textbook + Rajaram et al.
SO2_BOUNDS = (0.5, 0.99) # Tseng et al. 2009 
R_VESS_BOUNDS_UM = (3.0, 15.0)  # Rajaram et al.
F_V_MELANOSOME_BOUNDS = (0.003, 0.08) #Tseng et al. 2009 citando Jacques:
# Constante de decaimiento exponencial de la melanina [nm⁻¹].
# Modelo: µa_mel(λ) = 51.9 mm⁻¹ · exp(−k_mel · (λ − 500))
# Rango fisiológico in vivo ~0.008–0.015 nm⁻¹ (Zonios et al. 2008, JBO).
K_MEL_BOUNDS_NM1 = (0.004, 0.022) #Zonios 2008

# Punto inicial del ajuste en el mismo orden.
X0 = np.array([40.0, 0.75, 10.0, 0.030, 0.011], dtype=float)

PARAMETER_NAMES = (
    "HbT (µM)",
    "SO2",
    "r_vess (µm)",
    "f_v_melanosome",
    "k_mel (nm⁻¹)",
)


def get_solver_bounds() -> tuple[np.ndarray, np.ndarray]:
    """
    Construye los bounds del solver a partir de la configuración global.
    """
    lower = np.array([
        HBT_BOUNDS_UM[0],
        SO2_BOUNDS[0],
        R_VESS_BOUNDS_UM[0],
        F_V_MELANOSOME_BOUNDS[0],
        K_MEL_BOUNDS_NM1[0],
    ], dtype=float)
    upper = np.array([
        HBT_BOUNDS_UM[1],
        SO2_BOUNDS[1],
        R_VESS_BOUNDS_UM[1],
        F_V_MELANOSOME_BOUNDS[1],
        K_MEL_BOUNDS_NM1[1],
    ], dtype=float)

    if np.any(lower >= upper):
        raise ValueError(
            "Bounds fisiológicos inválidos: cada límite inferior debe ser menor "
            "que su límite superior."
        )

    return lower, upper


# ─────────────────────────────────────────────
# 1.  ESPECTROS BASE DE CROMÓFOROS
# ─────────────────────────────────────────────

# --- 1a. Hemoglobina  (Prahl / OMLC) ---
# Columnas: lambda_nm, epsilon_HbO2, epsilon_Hb   [cm⁻¹/M, base log10]
# Seleccionamos sólo el rango visible relevante (400–700 nm)
_HB_RAW = np.array([
    [400,266232,223296],[402,284224,236188],[404,308716,253368],
    [406,354208,270548],[408,422320,287356],[410,466840,303956],
    [412,500200,321344],[414,524280,342596],[416,521880,363848],
    [418,515520,385680],[420,480360,407560],[422,431880,429880],
    [424,376236,461200],[426,326032,481840],[428,283112,500840],
    [430,246072,528600],[432,214120,552160],[434,165332,552160],
    [436,132820,547040],[438,119140,501560],[440,102580,413280],
    [442,92780,363240],[444,81444,282724],[446,76324,237224],
    [448,67044,173320],[450,62816,103292],[452,58864,62640],
    [454,53552,36170],[456,49496,30698.8],[458,47496,25886.4],
    [460,44480,23388.8],[462,41320,20891.2],[464,39807.2,19260.8],
    [466,37073.2,18142.4],[468,34870.8,17025.6],[470,33209.2,16156.4],
    [472,31620,15310],[474,30113.6,15048.4],[476,28850.8,14792.8],
    [478,27718,14657.2],[480,26629.2,14550],[482,25701.6,14881.2],
    [484,25180.4,15212.4],[486,24669.6,15543.6],[488,24174.8,15898],
    [490,23684.4,16684],[492,23086.8,17469.6],[494,22457.6,18255.6],
    [496,21850.4,19041.2],[498,21260,19891.2],[500,20932.8,20862],
    [502,20596.4,21832.8],[504,20418,22803.6],[506,19946,23774.4],
    [508,19996,24745.2],[510,20035.2,25773.6],[512,20150.4,26936.8],
    [514,20429.2,28100],[516,21001.6,29263.2],[518,22509.6,30426.4],
    [520,24202.4,31589.6],[522,26450.4,32851.2],[524,29269.2,34397.6],
    [526,32496.4,35944],[528,35990,37490],[530,39956.8,39036.4],
    [532,43876,40584],[534,46924,42088],[536,49752,43592],
    [538,51712,45092],[540,53236,46592],[542,53292,48148],
    [544,52096,49708],[546,49868,51268],[548,46660,52496],
    [550,43016,53412],[552,39675.2,54080],[554,36815.2,54520],
    [556,34476.8,54540],[558,33456,54164],[560,32613.2,53788],
    [562,32620,52276],[564,33915.6,50572],[566,36495.2,48828],
    [568,40172,46948],[570,44496,45072],[572,49172,43340],
    [574,53308,41716],[576,55540,40092],[578,54728,38467.6],
    [580,50104,37020],[582,43304,35676.4],[584,34639.6,34332.8],
    [586,26600.4,32851.6],[588,19763.2,31075.2],[590,14400.8,28324.4],
    [592,10468.4,25470],[594,7678.8,22574.8],[596,5683.6,19800],
    [598,4504.4,17058.4],[600,3200,14677.2],[602,2664,13622.4],
    [604,2128,12567.6],[606,1789.2,11513.2],[608,1647.6,10477.6],
    [610,1506,9443.6],[612,1364.4,8591.2],[614,1222.8,7762],
    [616,1110,7344.8],[618,1026,6927.2],[620,942,6509.6],
    [622,858,6193.2],[624,774,5906.8],[626,707.6,5620],
    [628,658.8,5366.8],[630,610,5148.8],[632,561.2,4930.8],
    [634,512.4,4730.8],[636,478.8,4602.4],[638,460.4,4473.6],
    [640,442,4345.2],[642,423.6,4216.8],[644,405.2,4088.4],
    [646,390.4,3965.08],[648,379.2,3857.6],[650,368,3750.12],
    [652,356.8,3642.64],[654,345.6,3535.16],[656,335.2,3427.68],
    [658,325.6,3320.2],[660,319.6,3226.56],[662,314,3140.28],
    [664,308.4,3053.96],[666,302.8,2967.68],[668,298,2881.4],
    [670,294,2795.12],[672,290,2708.84],[674,285.6,2627.64],
    [676,282,2554.4],[678,279.2,2481.16],[680,277.6,2407.92],
])
# Conversión cm⁻¹/M → mm⁻¹/M  (÷10) y base log10 → natural (×2.303)
# Resultado: eps [mm⁻¹/µM] = raw_eps [cm⁻¹/M] × 2.303 × 0.1 × 1e-6
_CONV_HB = 2.303 * 0.1 * 1e-6   # mm⁻¹ por (µM)

HB_LAM    = _HB_RAW[:, 0]
EPS_HbO2  = _HB_RAW[:, 1] * _CONV_HB   # mm⁻¹/µM
EPS_Hb    = _HB_RAW[:, 2] * _CONV_HB   # mm⁻¹/µM

# --- 1b. Melanina (modelo exponencial, Zonios et al. 2008 JBO) ---
# µa_melanosoma(λ) = 519 cm⁻¹ · exp(−k_mel · (λ − 500 nm))
# Anclada al valor de referencia de Jacques en 500 nm para mantener
# compatibilidad física. El parámetro libre k_mel reemplaza el exponente m.
# Referencia: Zonios et al. (2008) J. Biomed. Opt. 13(1):014017.
MELANIN_REF_LAMBDA_NM = 500.0
MELANOSOME_MUA_REF_CM1 = 519.0
MELANOSOME_MUA_REF_MM1 = MELANOSOME_MUA_REF_CM1 * 0.1
K_MEL_DEFAULT_NM1 = 0.011   # valor inicial típico in vivo


def melanosome_mu_a_exponential(lam_nm: np.ndarray,
                                k_mel: float = K_MEL_DEFAULT_NM1) -> np.ndarray:
    """
    µa del interior del melanosoma [mm⁻¹] usando decaimiento exponencial en λ.

    µa_melanosoma(λ) = 51.9 mm⁻¹ · exp(−k_mel · (λ − 500))

    Anclada en λ=500 nm al valor de referencia de Jacques (519 cm⁻¹).
    k_mel [nm⁻¹]: constante de decaimiento. Rango in vivo: ~0.008–0.015 nm⁻¹
    (Zonios et al. 2008, JBO 13(1):014017).
    """
    lam_nm = np.asarray(lam_nm, dtype=float)
    return MELANOSOME_MUA_REF_MM1 * np.exp(-k_mel * (lam_nm - MELANIN_REF_LAMBDA_NM))


# ─────────────────────────────────────────────
# 2.  INTERPOLACIÓN A GRILLA DE DATOS
# ─────────────────────────────────────────────

def build_chromophore_matrix(lam_data: np.ndarray) -> dict:
    """
    Interpola todos los espectros base a la grilla de longitudes de onda
    del espectro medido.

    Retorna dict con arrays de la misma longitud que lam_data,
    listos para usar en el modelo.
    """
    def interp(lam_src, eps_src, lam_dst, fill_value="extrapolate"):
        f = interp1d(lam_src, eps_src, kind='linear',
                     bounds_error=False, fill_value=fill_value)
        vals = f(lam_dst)
        vals = np.clip(vals, 0, None)   # espectros físicos no pueden ser negativos
        return vals

    return {
        "HbO2": interp(HB_LAM, EPS_HbO2, lam_data),
        "Hb": interp(HB_LAM, EPS_Hb, lam_data),
        "lam_grid": np.asarray(lam_data, dtype=float),
    }


# ─────────────────────────────────────────────
# 3.  MODELO FÍSICO
# ─────────────────────────────────────────────

def cpack(mu_blood: np.ndarray, r_um: float) -> np.ndarray:
    """
    Factor de empaquetamiento vascular (van Veen et al. 2002).

    mu_blood : µa intravascular de sangre completa [mm⁻¹]
    r_um     : radio efectivo del vaso [µm]

    Retorna Cpack ∈ (0, 1].
    """
    r_mm = r_um * 1e-3          # µm → mm
    x = 2.0 * mu_blood * r_mm  # argumento adimensional
    # Límite x→0: Cpack→1; límite x→∞: Cpack→1/x
    # Usamos la forma analítica con protección numérica en x≈0
    safe_x = np.where(x < 1e-8, 1e-8, x)
    return (1.0 - np.exp(-safe_x)) / safe_x


def mu_a_model(theta: np.ndarray, eps: dict) -> np.ndarray:
    """
    µa predicho por el modelo en toda la grilla espectral.

    theta = [HbT (µM), SO2 (frac), r_vess (µm), f_v_melanosome (fracción), k_mel (nm⁻¹)]
    eps   = dict de espectros base interpolados (output de build_chromophore_matrix)

    Melanina: modelo exponencial (Zonios et al. 2008, JBO 13(1):014017).
    """
    HbT, SO2, r_vess, f_v_mel, k_mel = theta

    # Absorción intravascular (sangre completa, 2326 µM de Hb)
    mu_blood = (SO2 * eps["HbO2"] + (1.0 - SO2) * eps["Hb"]) * 2326.0   # mm⁻¹

    # Factor de empaquetamiento
    Cp = cpack(mu_blood, r_vess)

    # Fracción volumétrica de sangre en tejido
    nu_blood = HbT / 2326.0

    # Contribución vascular al tejido
    mu_vascular = Cp * mu_blood * nu_blood

    # Contribución de melanina (modelo exponencial)
    mu_mel_base = melanosome_mu_a_exponential(eps["lam_grid"], k_mel)
    mu_mel = f_v_mel * mu_mel_base

    return mu_vascular + mu_mel


def residuals(theta: np.ndarray, lam: np.ndarray, mu_a_meas: np.ndarray,
              eps: dict) -> np.ndarray:
    return mu_a_meas - mu_a_model(theta, eps)


# ─────────────────────────────────────────────
# 4.  SOLVER
# ─────────────────────────────────────────────

def fit_spectrum(lam: np.ndarray, mu_a_meas: np.ndarray,
                 eps: dict, x0: np.ndarray = X0,
                 verbose: bool = False):
    """
    Ajusta el modelo a un espectro µa medido.

    Retorna:
        result   : objeto OptimizeResult de scipy
        theta_opt: array de parámetros óptimos
        mu_a_fit : espectro reconstruido
        rms      : RMS del residual [mm⁻¹]
    """
    lower, upper = get_solver_bounds()
    x0 = np.asarray(x0, dtype=float)
    if x0.shape != lower.shape:
        raise ValueError(
            f"x0 tiene forma {x0.shape}, pero se esperaban {len(lower)} parámetros."
        )
    x0_fit = np.clip(x0, lower, upper)

    result = least_squares(
        residuals,
        x0_fit,
        bounds=(lower, upper),
        method='trf',
        args=(lam, mu_a_meas, eps),
        ftol=1e-10, xtol=1e-10, gtol=1e-10,
        max_nfev=5000,
    )
    theta_opt = result.x
    mu_a_fit  = mu_a_model(theta_opt, eps)
    rms       = np.sqrt(np.mean(result.fun**2))

    if verbose:
        print("\n── Resultado del ajuste ──")
        for n, v in zip(PARAMETER_NAMES, theta_opt):
            print(f"  {n:22s}: {v:.4f}")
        print(f"  {'RMS residual':22s}: {rms:.5f} mm⁻¹")
        print(f"  {'Convergencia':22s}: {result.message}")

    return result, theta_opt, mu_a_fit, rms


# ─────────────────────────────────────────────
# 5.  DIAGNÓSTICO: PERFIL DE COSTO vs r_vess
# ─────────────────────────────────────────────

def degeneracy_profile(lam: np.ndarray, mu_a_meas: np.ndarray,
                       eps: dict, r_grid: np.ndarray = None) -> tuple:
    """
    Fija r_vess en una grilla y minimiza los demás parámetros.
    Devuelve (r_grid, rms_profile, HbT_profile).

    Permite visualizar la degeneración HbT–r_vess descrita en el doc.
    """
    if r_grid is None:
        r_grid = np.linspace(5, 25, 20)

    rms_profile = np.zeros_like(r_grid)
    HbT_profile = np.zeros_like(r_grid)

    lower_fixed, upper_fixed = get_solver_bounds()

    for i, r_fixed in enumerate(r_grid):
        lower_fixed[2] = r_fixed - 1e-6
        upper_fixed[2] = r_fixed + 1e-6
        x0_r = np.clip(X0.copy(), lower_fixed, upper_fixed)
        x0_r[2] = r_fixed
        res = least_squares(
            residuals, x0_r,
            bounds=(lower_fixed, upper_fixed),
            method='trf', args=(lam, mu_a_meas, eps),
            ftol=1e-10, xtol=1e-10,
        )
        rms_profile[i] = np.sqrt(np.mean(res.fun**2))
        HbT_profile[i] = res.x[0]

    return r_grid, rms_profile, HbT_profile


# ─────────────────────────────────────────────
# 6.  EJECUCIÓN SOBRE ESTRUCTURA DE SUJETO
# ─────────────────────────────────────────────

INPUT_CSV_NAME = "resumen_resultados_phan_sierra.csv"
OUTPUT_PREFIX = "chromophore_unmixing_5param"
OUTPUT_SUBDIR_NAME = "chromophore_unmixing_results"
MEASUREMENTS_DIR_CANDIDATES = ("mediciones", "Mediciones")


def _path_from_user_text(raw_path: str) -> Path:
    """
    Convierte rutas introducidas manualmente. Si se ejecuta desde WSL y el
    usuario pega una ruta Windows tipo C:\\..., se traduce a /mnt/c/...
    """
    cleaned = raw_path.strip().strip('"').strip("'")
    if os.name != "nt" and len(cleaned) >= 3 and cleaned[1] == ":" and cleaned[2] in ("\\", "/"):
        drive = cleaned[0].lower()
        rest = cleaned[3:].replace("\\", "/")
        return Path(f"/mnt/{drive}/{rest}")
    return Path(cleaned).expanduser()


def get_measurements_root(project_root: Path | None = None) -> Path:
    """
    Devuelve la carpeta que contiene los sujetos medidos dentro del proyecto.
    Se acepta tanto `mediciones/` como `Mediciones/` para no romper árboles
    existentes.
    """
    base_dir = (project_root or Path(__file__).resolve().parent).resolve()
    for folder_name in MEASUREMENTS_DIR_CANDIDATES:
        candidate = base_dir / folder_name
        if candidate.exists() and candidate.is_dir():
            return candidate
    return base_dir / MEASUREMENTS_DIR_CANDIDATES[0]


def slugify_label(text: str) -> str:
    """
    Normaliza etiquetas para usarlas en nombres de archivo sin perder la
    referencia anatómica.
    """
    normalized = unicodedata.normalize("NFKD", str(text))
    ascii_text = normalized.encode("ascii", "ignore").decode("ascii").lower()
    slug = re.sub(r"[^a-z0-9]+", "_", ascii_text).strip("_")
    return slug


def list_root_subfolders(root_dir: Path) -> list:
    """
    Lista subcarpetas directas de la raíz del proyecto candidatas a selección
    manual. Se omiten carpetas ocultas y artefactos obvios de Python/IDE.
    """
    folders = []
    for child in sorted(root_dir.iterdir()):
        if not child.is_dir():
            continue
        if child.name.startswith(".") or child.name.startswith("__"):
            continue
        folders.append(child)
    return folders


def resolve_subject_root(selected_dir: Path, site_dir: Path) -> Path:
    """
    Devuelve la carpeta raíz del sujeto. Si el usuario seleccionó directamente
    un sitio, se usa la carpeta padre como sujeto.
    """
    return site_dir.parent if selected_dir.resolve() == site_dir.resolve() else selected_dir


def select_subject_folder(subject_dir_arg: str = "",
                          root_dir: Path | None = None) -> Path:
    """
    Selecciona la carpeta de sujeto. Si no se pasa `--subject-dir`, muestra en
    terminal las subcarpetas directas de `mediciones/` y pide el índice.
    """
    project_root = (root_dir or Path(__file__).resolve().parent).resolve()
    measurements_root = get_measurements_root(project_root)

    if subject_dir_arg:
        subject_path = _path_from_user_text(subject_dir_arg)
        if subject_path.is_absolute() or subject_path.exists():
            return subject_path
        candidate = measurements_root / subject_path
        if candidate.exists():
            return candidate
        return subject_path

    base_dir = measurements_root
    candidates = list_root_subfolders(base_dir)
    if not candidates:
        raise SystemExit(
            f"No encontré subcarpetas seleccionables en la carpeta de mediciones: {base_dir}"
        )

    print(f"\nCarpeta de mediciones: {base_dir}")
    print("Carpetas disponibles:")
    for idx, folder in enumerate(candidates):
        print(f"  [{idx}] {folder.name}")

    try:
        while True:
            selected_text = input("Índice de la carpeta de sujeto: ").strip()
            if not selected_text:
                print("Introduce un índice válido.")
                continue
            if selected_text.lower() in {"q", "quit", "exit"}:
                raise SystemExit("Selección cancelada por el usuario.")

            try:
                selected_idx = int(selected_text)
            except ValueError:
                print("Entrada inválida. Debes escribir un número entero.")
                continue

            if 0 <= selected_idx < len(candidates):
                return candidates[selected_idx]

            print(f"Índice fuera de rango. Usa un valor entre 0 y {len(candidates) - 1}.")
    except Exception as exc:
        if isinstance(exc, EOFError):
            raise SystemExit(
                "No hay entrada interactiva disponible. Ejecuta con "
                "--subject-dir '<ruta_sujeto>'."
            ) from exc
        raise


def find_site_jobs(subject_dir: Path) -> list:
    """
    Busca subcarpetas directas de un sujeto que tengan salida IAD nominal:
        <sujeto>/<sitio>/IAD_results/resumen_resultados_phan_sierra.csv

    También acepta que el usuario seleccione directamente una carpeta de sitio.
    """
    direct_csv = subject_dir / "IAD_results" / INPUT_CSV_NAME
    if direct_csv.exists():
        return [(subject_dir, subject_dir / "IAD_results", direct_csv)]

    jobs = []
    skipped = []
    for site_dir in sorted(p for p in subject_dir.iterdir() if p.is_dir()):
        iad_dir = site_dir / "IAD_results"
        csv_path = iad_dir / INPUT_CSV_NAME
        if csv_path.exists():
            jobs.append((site_dir, iad_dir, csv_path))
        else:
            skipped.append(site_dir.name)

    if skipped:
        print("Subcarpetas sin salida IAD esperada; se omiten:")
        for name in skipped:
            print(f"  - {name}")

    return jobs


def select_nominal_iad_branch(df: pd.DataFrame, source_csv: Path) -> tuple:
    """
    Evita mezclar ramas de incertidumbre/calibración IAD. Si hay varias ramas,
    usa explícitamente la rama nominal cuando existe.
    """
    branch_cols = [c for c in ("escenario", "escenario_mr", "factor_musp") if c in df.columns]
    if not branch_cols:
        return df, "sin columnas de escenario; se usa todo el CSV"

    branches = df[branch_cols].drop_duplicates()
    if len(branches) <= 1:
        branch_text = ", ".join(f"{c}={branches.iloc[0][c]}" for c in branch_cols)
        return df, branch_text

    mask = pd.Series(True, index=df.index)
    if "escenario" in df.columns:
        mask &= df["escenario"].astype(str).str.lower().eq("nominal")
    if "escenario_mr" in df.columns:
        mask &= df["escenario_mr"].astype(str).str.lower().eq("nominal")
    if "factor_musp" in df.columns:
        factor = pd.to_numeric(df["factor_musp"], errors="coerce")
        mask &= np.isclose(factor, 1.0, rtol=0.0, atol=1e-12)

    if not mask.any():
        raise ValueError(
            f"{source_csv} contiene varias ramas IAD, pero no hay rama nominal "
            "(escenario=nominal, escenario_mr=nominal, factor_musp=1.0)."
        )

    return df.loc[mask].copy(), "rama nominal IAD"


def prepare_iad_spectrum(df: pd.DataFrame, source_csv: Path) -> tuple:
    """
    Extrae el espectro µa(lambda) recuperado por IAD sin mezclar regiones ni
    escenarios. Retorna df_loc, region_id, phan_ubicacion y nota de filtrado.
    """
    required = {"lambda_nm", "mu_a_mm-1"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{source_csv} no tiene columnas requeridas: {sorted(missing)}")

    df_loc = df.copy()

    if "returncode" in df_loc.columns:
        ok_return = pd.to_numeric(df_loc["returncode"], errors="coerce").eq(0)
        if ok_return.any():
            df_loc = df_loc.loc[ok_return].copy()

    if "txt_found" in df_loc.columns:
        ok_txt = df_loc["txt_found"].astype(str).str.lower().isin(("true", "1", "yes", "si", "sí"))
        if ok_txt.any():
            df_loc = df_loc.loc[ok_txt].copy()

    if df_loc.empty:
        raise ValueError(
            f"{source_csv} no tiene puntos IAD válidos después de filtrar returncode/txt_found."
        )

    df_loc, branch_note = select_nominal_iad_branch(df_loc, source_csv)

    if "region_id" in df_loc.columns:
        regions = df_loc["region_id"].dropna().astype(str).unique()
        if len(regions) != 1:
            raise ValueError(
                f"{source_csv} contiene {len(regions)} region_id distintos después "
                "del filtrado. No se mezclan regiones anatómicas en un solo ajuste."
            )
        region_id = regions[0]
    else:
        region_id = source_csv.parents[1].name

    if "phan_ubicacion" in df_loc.columns:
        phan_vals = df_loc["phan_ubicacion"].dropna().astype(str).unique()
        phan_ubicacion = phan_vals[0] if len(phan_vals) == 1 else ""
    else:
        phan_ubicacion = ""

    df_loc["lambda_nm"] = pd.to_numeric(df_loc["lambda_nm"], errors="coerce")
    df_loc["mu_a_mm-1"] = pd.to_numeric(df_loc["mu_a_mm-1"], errors="coerce")
    df_loc = df_loc.dropna(subset=["lambda_nm", "mu_a_mm-1"])
    df_loc = df_loc.sort_values("lambda_nm").reset_index(drop=True)

    if df_loc.empty:
        raise ValueError(f"{source_csv} no tiene puntos válidos de lambda_nm y mu_a_mm-1.")
    if df_loc["lambda_nm"].duplicated().any():
        raise ValueError(
            f"{source_csv} tiene longitudes de onda duplicadas después del filtrado; "
            "revisa que no se estén mezclando escenarios."
        )
    if len(df_loc) < len(X0) + 1:
        raise ValueError(
            f"{source_csv} tiene sólo {len(df_loc)} puntos válidos; insuficiente para "
            f"ajustar {len(X0)} parámetros libres."
        )

    return df_loc, region_id, phan_ubicacion, branch_note


def save_unmixing_plot(output_png: Path, site_label: str, lam: np.ndarray,
                       mu_a_obs: np.ndarray, mu_a_fit: np.ndarray,
                       contrib_vasc: np.ndarray, contrib_mel: np.ndarray,
                       r_grid: np.ndarray,
                       rms_prof: np.ndarray, HbT_prof: np.ndarray,
                       theta_opt: np.ndarray, rms: float) -> None:
    HbT, SO2, r_vess, f_v_mel, k_mel = theta_opt

    fig = plt.figure(figsize=(14, 10))
    gs = gridspec.GridSpec(2, 2, hspace=0.40, wspace=0.35)

    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(lam, mu_a_obs, 'k-', lw=1.8, label='µa medido (IAD)', alpha=0.85)
    ax1.plot(lam, mu_a_fit, 'r--', lw=2.0, label=f'Modelo ajustado (RMS = {rms:.4f} mm⁻¹)')
    ax1.stackplot(
        lam,
        contrib_vasc, contrib_mel,
        labels=['Hemoglobina (con Cpack)', 'Melanina simplificada'],
        colors=['#D04030', '#8B4513'], alpha=0.25,
    )
    ax1.set_xlabel('λ (nm)', fontsize=12)
    ax1.set_ylabel('µa (mm⁻¹)', fontsize=12)
    ax1.set_title(f'Unmixing espectral - {site_label}', fontsize=13)
    ax1.legend(fontsize=9, loc='upper right')
    ax1.grid(alpha=0.3)

    ax2 = fig.add_subplot(gs[1, 0])
    residual = mu_a_obs - mu_a_fit
    ax2.plot(lam, residual, 'b-', lw=1.2)
    ax2.axhline(0, color='k', lw=0.8, ls='--')
    ax2.fill_between(lam, residual, 0, alpha=0.15, color='blue')
    ax2.set_xlabel('λ (nm)', fontsize=11)
    ax2.set_ylabel('Residual (mm⁻¹)', fontsize=11)
    ax2.set_title('Residual del ajuste', fontsize=12)
    ax2.grid(alpha=0.3)

    ax3 = fig.add_subplot(gs[1, 1])
    ax3_twin = ax3.twinx()
    l1, = ax3.plot(r_grid, rms_prof * 1000, 'g-o', ms=4, lw=1.5, label='RMS (×10³)')
    l2, = ax3_twin.plot(r_grid, HbT_prof, 'r-s', ms=4, lw=1.5, label='HbT (µM)')
    ax3.axvline(theta_opt[2], color='k', ls='--', lw=1.0, alpha=0.6,
                label=f'r* = {theta_opt[2]:.1f} µm')
    ax3.set_xlabel('r_vess (µm)', fontsize=11)
    ax3.set_ylabel('RMS ×10³ (mm⁻¹)', fontsize=11, color='g')
    ax3_twin.set_ylabel('HbT (µM)', fontsize=11, color='r')
    ax3.set_title('Perfil degeneración HbT-r_vess', fontsize=12)
    ax3.legend(handles=[l1, l2], loc='upper left', fontsize=8)
    ax3.grid(alpha=0.3)

    param_text = (
        f"HbT = {HbT:.1f} µM\n"
        f"SO₂ = {SO2*100:.1f} %\n"
        f"r_vess = {r_vess:.1f} µm\n"
        f"f_v_mel = {f_v_mel:.4f}\n"
        f"k_mel = {k_mel:.4f} nm⁻¹"
    )
    ax1.text(
        0.01, 0.05, param_text,
        transform=ax1.transAxes, fontsize=9,
        verticalalignment='bottom',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.8),
    )

    fig.savefig(output_png, dpi=150, bbox_inches='tight')
    plt.close(fig)


def build_site_output_paths(output_dir: Path, site_dir: Path,
                            region_id: str, phan_ubicacion: str,
                            subject_name: str = "") -> dict:
    """
    Construye rutas de salida por sitio dentro de una carpeta común del sujeto.
    Se añade una referencia anatómica al nombre para evitar sobrescrituras.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    anatomical_label = phan_ubicacion or region_id or site_dir.name
    anatomical_slug = slugify_label(anatomical_label)
    site_slug = slugify_label(site_dir.name)

    stem_parts = [OUTPUT_PREFIX]
    if anatomical_slug:
        stem_parts.append(anatomical_slug)
    if site_slug and site_slug != anatomical_slug:
        stem_parts.append(site_slug)
    output_stem = "_".join(stem_parts)
    subject_slug = slugify_label(subject_name)
    png_stem = output_stem
    if subject_slug:
        png_stem = f"{output_stem}_{subject_slug}"

    return {
        "results_csv": output_dir / f"{output_stem}_results.csv",
        "summary_csv": output_dir / f"{output_stem}_fit_summary.csv",
        "degeneracy_csv": output_dir / f"{output_stem}_degeneracy_profile.csv",
        "output_png": output_dir / f"{png_stem}.png",
    }


def resolve_subject_output_dir(selected_dir: Path, site_dir: Path) -> Path:
    """
    Devuelve la carpeta de salida del sujeto. Si el usuario seleccionó
    directamente un sitio, se usa la carpeta padre como sujeto.
    """
    subject_root = resolve_subject_root(selected_dir, site_dir)
    return subject_root / OUTPUT_SUBDIR_NAME


def load_subject_name(subject_root: Path) -> str:
    """
    Lee el nombre del sujeto desde datos_sujeto.txt.
    Si no existe o no hay un nombre útil, retorna cadena vacía.
    """
    subject_file = subject_root / "datos_sujeto.txt"
    if not subject_file.exists():
        return ""

    try:
        for raw_line in subject_file.read_text(encoding="utf-8").splitlines():
            line = raw_line.strip()
            if not line or ":" not in line:
                continue
            key, value = line.split(":", 1)
            if key.strip().lower() != "nombre":
                continue

            subject_name = value.strip()
            if subject_name.lower() in {"", "no especificado", "n/a", "na", "none"}:
                return ""
            return subject_name
    except Exception:
        return ""

    return ""


def run_unmixing_for_site(site_dir: Path, csv_path: Path, output_dir: Path,
                          subject_name: str = "") -> dict:
    df = pd.read_csv(csv_path)
    df_loc, region_id, phan_ubicacion, branch_note = prepare_iad_spectrum(df, csv_path)

    lam = df_loc["lambda_nm"].to_numpy(dtype=float)
    mu_a_obs = df_loc["mu_a_mm-1"].to_numpy(dtype=float)

    site_label = f"{site_dir.name} [{region_id}]"
    if phan_ubicacion:
        site_label += f" / {phan_ubicacion}"

    print(f"\nSitio: {site_label}")
    print(f"  CSV: {csv_path}")
    print(f"  Rama usada: {branch_note}")
    print(f"  λ: {lam[0]:.1f} - {lam[-1]:.1f} nm ({len(lam)} puntos)")
    print(f"  µa: {mu_a_obs.min():.4f} - {mu_a_obs.max():.4f} mm⁻¹")

    eps = build_chromophore_matrix(lam)
    result, theta_opt, mu_a_fit, rms = fit_spectrum(lam, mu_a_obs, eps, verbose=True)
    r_grid, rms_prof, HbT_prof = degeneracy_profile(lam, mu_a_obs, eps)

    HbT, SO2, r_vess, f_v_mel, k_mel = theta_opt
    mu_blood_opt = (SO2 * eps["HbO2"] + (1.0 - SO2) * eps["Hb"]) * 2326.0
    Cp_opt = cpack(mu_blood_opt, r_vess)
    contrib_vasc = Cp_opt * mu_blood_opt * (HbT / 2326.0)
    mu_mel_base = melanosome_mu_a_exponential(eps["lam_grid"], k_mel)
    contrib_mel = f_v_mel * mu_mel_base
    residual = mu_a_obs - mu_a_fit

    output_paths = build_site_output_paths(
        output_dir, site_dir, region_id, phan_ubicacion, subject_name=subject_name
    )
    results_csv = output_paths["results_csv"]
    summary_csv = output_paths["summary_csv"]
    degeneracy_csv = output_paths["degeneracy_csv"]
    output_png = output_paths["output_png"]

    pd.DataFrame({
        "wavelength_nm": lam,
        "mu_a_measured_mm-1": mu_a_obs,
        "mu_a_model_mm-1": mu_a_fit,
        "mu_a_vascular_mm-1": contrib_vasc,
        "mu_a_melanin_simplified_mm-1": contrib_mel,
        "mu_a_eumelanin_mm-1": np.nan,
        "mu_a_pheomelanin_mm-1": np.nan,
        "residual_mm-1": residual,
        "region_id": region_id,
        "phan_ubicacion": phan_ubicacion,
        "source_csv": str(csv_path),
    }).to_csv(results_csv, index=False)

    pd.DataFrame([{
        "site_folder": site_dir.name,
        "region_id": region_id,
        "phan_ubicacion": phan_ubicacion,
        "source_csv": str(csv_path),
        "branch_used": branch_note,
        "HbT_uM": HbT,
        "SO2_fraction": SO2,
        "r_vess_um": r_vess,
        "f_v_melanosome": f_v_mel,
        "k_mel_nm-1": k_mel,
        "melanosome_mua_ref_cm-1_at_500nm": MELANOSOME_MUA_REF_CM1,
        "c_eumel_mg_per_mL": np.nan,
        "c_feomel_mg_per_mL": np.nan,
        "model_version": "melanina_exponencial_Zonios2008",
        "rms_residual_mm-1": rms,
        "success": bool(result.success),
        "status": result.status,
        "message": result.message,
    }]).to_csv(summary_csv, index=False)

    pd.DataFrame({
        "r_vess_um": r_grid,
        "rms_residual_mm-1": rms_prof,
        "HbT_uM": HbT_prof,
    }).to_csv(degeneracy_csv, index=False)

    save_unmixing_plot(
        output_png, site_label, lam, mu_a_obs, mu_a_fit,
        contrib_vasc, contrib_mel,
        r_grid, rms_prof, HbT_prof, theta_opt, rms,
    )

    print("  Salidas:")
    print(f"    {results_csv}")
    print(f"    {summary_csv}")
    print(f"    {degeneracy_csv}")
    print(f"    {output_png}")

    return {
        "site": site_dir.name,
        "region_id": region_id,
        "rms": rms,
        "success": bool(result.success),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Unmixing cromofórico por sitio anatómico usando la estructura "
            "mediciones/<sujeto>/<sitio>/IAD_results/resumen_resultados_phan_sierra.csv"
        )
    )
    parser.add_argument(
        "--subject-dir",
        type=str,
        default="",
        help=(
            "Ruta de carpeta de sujeto. Si no se indica, se listan las "
            "subcarpetas directas de mediciones/ para elegir por índice. Si la "
            "ruta es relativa, primero se intenta resolver dentro de mediciones/."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    subject_dir = select_subject_folder(args.subject_dir).resolve()
    if not subject_dir.exists() or not subject_dir.is_dir():
        raise SystemExit(f"La carpeta seleccionada no existe o no es directorio: {subject_dir}")

    jobs = find_site_jobs(subject_dir)
    if not jobs:
        raise SystemExit(
            f"No encontré subcarpetas con {INPUT_CSV_NAME} dentro de IAD_results en: {subject_dir}"
        )

    print(f"\nCarpeta de sujeto: {subject_dir}")
    print(f"Sitios a procesar: {len(jobs)}")

    completed = []
    failures = []
    for site_dir, iad_dir, csv_path in jobs:
        try:
            subject_root = resolve_subject_root(subject_dir, site_dir)
            subject_name = load_subject_name(subject_root)
            output_dir = resolve_subject_output_dir(subject_dir, site_dir)
            completed.append(
                run_unmixing_for_site(
                    site_dir, csv_path, output_dir, subject_name=subject_name
                )
            )
        except Exception as exc:
            failures.append((site_dir, exc))
            print(f"\nERROR en {site_dir.name}: {exc}")

    print("\nResumen de unmixing:")
    for item in completed:
        print(
            f"  OK {item['site']} [{item['region_id']}]: "
            f"RMS={item['rms']:.5f} mm⁻¹, success={item['success']}"
        )
    for site_dir, exc in failures:
        print(f"  FALLÓ {site_dir.name}: {exc}")

    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()