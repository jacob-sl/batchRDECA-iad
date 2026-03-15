from pathlib import Path
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# CONFIGURACIÓN
# ============================================================

# Ruta al CSV con tus datos M_R
MR_CSV_PATH = Path(__file__).parent / "M_R_data.csv"
# Ruta al archivo plantilla .rxt
RXT_TEMPLATE_PATH = Path(__file__).parent / "sample-F.rxt"
# Ruta al ejecutable de IAD
IAD_EXE_PATH = Path(__file__).parent / "IADSCOTT" / "iad.exe"
# Carpeta principal de trabajo
IAD_FOLDER_NAME = "IAD_run"
# Subcarpeta donde se guardará una corrida por longitud de onda
IAD_PER_LAMBDA_FOLDER_NAME = "por_lambda"

# ============================================================
# MODELO ÓPTICO FIJO
# ============================================================

# g fijo y n fijo, como acordamos.
USAR_G_FIJO = True
G_FIJO = 0.8
N_TEJIDO = 1.41  # actualmente solo documental; IAD lo toma del .rxt/header

# Modo rápido omite el "sanity check" que IAD hace con una simulación Montecarlo.
USAR_MODO_RAPIDO = True
# Mostrar o no el comando completo por cada lambda
MOSTRAR_COMANDOS = False

# ── TOGGLE PRINCIPAL, UN SOLO ESPECTRO O UNA TANDA CON MARCAS TEMPORALES ──
MODO_TEMPORAL = False   # False → M_R_data.csv (un espectro)
                        # True  → M_R_tiempo_data.csv (serie temporal)

# Solo se usa cuando MODO_TEMPORAL = True
CSV_TEMPORAL_NAME = "M_R_tiempo_data.csv"
CSV_SALIDA_TEMPORAL = "IAD_run/resumen_iad_temporal_phan_sierra.csv"
INICIO_MEDICION = 22
MAX_MEDICIONES = None
VENTANA_PROMEDIO = 5
WORKERS = os.cpu_count() - 1 or 4


# ============================================================
# FUNCIÓN PHAN-SIERRA para μs'(λ)
# ============================================================
# Este script descarta la versión heurística previa. La versión válida
# usa solo datos Palm reportados en el suplemento de Phan et al.
#
# Datos medidos en literatura:
#   λ (nm):        471, 526, 591, 621, 659, 691, 731, 851
#   μs' (mm^-1):   2.57, 2.18, 1.78, 1.72, 1.60, 1.54, 1.51, 1.45
#   sd(μs')       0.421,0.278,0.147,0.126,0.108,0.0984,0.0883,0.0813
#   μa (mm^-1):   0.317,0.215,0.092,0.0288,0.0204,0.0147,0.0111,0.00978
#
# Referencia anatómica:
# Palm = glabrous skin, aproximación razonable para tejido palmar/digital.

PHAN_LAMBDA_NM = np.array([471.0, 526.0, 591.0, 621.0, 659.0, 691.0, 731.0, 851.0], dtype=float)
PHAN_MUSP_MM = np.array([2.57, 2.18, 1.78, 1.72, 1.60, 1.54, 1.51, 1.45], dtype=float)
PHAN_MUSP_SD_MM = np.array([0.421, 0.278, 0.147, 0.126, 0.108, 0.0984, 0.0883, 0.0813], dtype=float)
PHAN_MUA_MM = np.array([0.317, 0.215, 0.092, 0.0288, 0.0204, 0.0147, 0.0111, 0.00978], dtype=float)

# Selector de modelo para μs'(λ)
PHAN_SIERRA_MODE = "pchip"   # opciones: "powerlaw" o "pchip"
PHAN_SIERRA_ALLOW_EXTRAPOLATION = False

# Referencia para el ajuste compacto por ley de potencia
PHAN_SIERRA_LAMBDA_REF_NM = 526.0

_x = np.log(PHAN_LAMBDA_NM / PHAN_SIERRA_LAMBDA_REF_NM)
_y = np.log(PHAN_MUSP_MM)
_slope, _intercept = np.polyfit(_x, _y, 1)

PHAN_SIERRA_A_MM = float(np.exp(_intercept))
PHAN_SIERRA_B = float(-_slope)

PHAN_MUA_INTERP_MODE = "pchip"
PHAN_MUA_ALLOW_EXTRAPOLATION = False
PHAN_MUA_NORM_LAMBDA_REF_NM = 526.0

# Sensibilidad basada en literatura: nominal y banda ±1σ de Phan.
ESCENARIOS_MUSP = ("menos_1sd", "nominal", "mas_1sd")

try:
    from scipy.interpolate import PchipInterpolator
    _phan_sierra_pchip = PchipInterpolator(PHAN_LAMBDA_NM, PHAN_MUSP_MM, extrapolate=True)
    _phan_musp_sd_pchip = PchipInterpolator(PHAN_LAMBDA_NM, PHAN_MUSP_SD_MM, extrapolate=True)
    _phan_mua_pchip = PchipInterpolator(PHAN_LAMBDA_NM, PHAN_MUA_MM, extrapolate=True)
except Exception:
    _phan_sierra_pchip = None
    _phan_musp_sd_pchip = None
    _phan_mua_pchip = None

# ============================================================
# FUNCIONES
# ============================================================

def leer_csv_mr(csv_path: Path) -> pd.DataFrame:
    """Lee el archivo M_R_data.csv y valida que tenga las columnas esperadas."""
    df = pd.read_csv(csv_path)

    columnas_esperadas = {"wavelength_nm", "reflectance"}
    if not columnas_esperadas.issubset(df.columns):
        raise ValueError(
            f"El CSV debe contener las columnas {columnas_esperadas}. "
            f"Columnas encontradas: {list(df.columns)}"
        )

    df = df[["wavelength_nm", "reflectance"]].copy()
    df = df.dropna()
    return df



def extraer_encabezado_rxt(template_path: Path) -> list[str]:
    """Lee el archivo plantilla .rxt y devuelve solo el encabezado."""
    with template_path.open("r", encoding="utf-8") as f:
        lines = f.readlines()

    header_lines = []
    encontro_tabla = False

    for line in lines:
        stripped = line.strip()
        if stripped.startswith("#lambda"):
            encontro_tabla = True
            break
        header_lines.append(line.rstrip("\n"))

    if not encontro_tabla:
        raise ValueError(
            "No se encontró la línea '#lambda' en la plantilla .rxt. "
            "Revisa que el archivo tenga el formato esperado."
        )

    return header_lines



def _scalar_si_corresponde(valor):
    """Devuelve float si el resultado es escalar; si no, devuelve ndarray."""
    if np.ndim(valor) == 0:
        return float(valor)
    return np.asarray(valor, dtype=float)



def _validar_rango_phan(lam: np.ndarray, allow_extrapolation: bool, etiqueta: str) -> None:
    """Bloquea extrapolación fuera del soporte de Phan si no está permitida."""
    if allow_extrapolation:
        return

    lam_min = float(PHAN_LAMBDA_NM.min())
    lam_max = float(PHAN_LAMBDA_NM.max())
    fuera = (lam < lam_min) | (lam > lam_max)
    if np.any(fuera):
        raise ValueError(
            f"{etiqueta} se solicitó fuera del rango soportado por Phan "
            f"[{lam_min:.1f}, {lam_max:.1f}] nm. "
            "Se bloquea la extrapolación por consistencia con la literatura."
        )



def _interp_1d_phan(
    lambda_nm: float | np.ndarray,
    valores: np.ndarray,
    pchip_obj,
    *,
    mode: str = "pchip",
    allow_extrapolation: bool = False,
    etiqueta: str,
):
    """Interpola datos de Phan con PCHIP o lineal, con extrapolación controlada."""
    lam = np.asarray(lambda_nm, dtype=float)
    _validar_rango_phan(lam, allow_extrapolation=allow_extrapolation, etiqueta=etiqueta)

    if mode == "pchip":
        if pchip_obj is not None:
            salida = pchip_obj(lam)
        else:
            salida = np.interp(lam, PHAN_LAMBDA_NM, valores)
    elif mode == "linear":
        salida = np.interp(lam, PHAN_LAMBDA_NM, valores)
    else:
        raise ValueError(f"Modo desconocido para {etiqueta}: {mode}")

    return _scalar_si_corresponde(salida)



def phan_sierra_mu_sp_mm(
    lambda_nm: float | np.ndarray,
    mode: str = PHAN_SIERRA_MODE,
    allow_extrapolation: bool = PHAN_SIERRA_ALLOW_EXTRAPOLATION,
):
    """
    Devuelve μs'(λ) en mm^-1 usando la función Phan-Sierra.

    Parámetros
    ----------
    lambda_nm : float o array-like
        Longitud(es) de onda en nm.
    mode : str
        'powerlaw' -> ajuste ley de potencia
        'pchip'    -> spline monótona sobre puntos reales de Palm
    """
    lam = np.asarray(lambda_nm, dtype=float)
    _validar_rango_phan(lam, allow_extrapolation=allow_extrapolation, etiqueta="μs'(λ) de Phan-Sierra")

    if mode == "powerlaw":
        mu_sp = PHAN_SIERRA_A_MM * (lam / PHAN_SIERRA_LAMBDA_REF_NM) ** (-PHAN_SIERRA_B)
    elif mode == "pchip":
        mu_sp = _interp_1d_phan(
            lam,
            PHAN_MUSP_MM,
            _phan_sierra_pchip,
            mode="pchip",
            allow_extrapolation=allow_extrapolation,
            etiqueta="μs'(λ) de Phan-Sierra",
        )
    else:
        raise ValueError(f"Modo desconocido para Phan-Sierra: {mode}")

    return _scalar_si_corresponde(mu_sp)



def phan_musp_sd_mm(
    lambda_nm: float | np.ndarray,
    mode: str = "pchip",
    allow_extrapolation: bool = False,
):
    """Interpola la desviación estándar reportada por Phan para μs'."""
    return _interp_1d_phan(
        lambda_nm,
        PHAN_MUSP_SD_MM,
        _phan_musp_sd_pchip,
        mode=mode,
        allow_extrapolation=allow_extrapolation,
        etiqueta="sd de μs' de Phan",
    )



def phan_mu_a_mm(
    lambda_nm: float | np.ndarray,
    mode: str = PHAN_MUA_INTERP_MODE,
    allow_extrapolation: bool = PHAN_MUA_ALLOW_EXTRAPOLATION,
):
    """Interpola μa Palm reportado por Phan a las lambdas deseadas."""
    return _interp_1d_phan(
        lambda_nm,
        PHAN_MUA_MM,
        _phan_mua_pchip,
        mode=mode,
        allow_extrapolation=allow_extrapolation,
        etiqueta="μa Palm de Phan",
    )



def mu_sp_escenario_phan_mm(
    lambda_nm: float,
    escenario: str,
    mode: str = PHAN_SIERRA_MODE,
) -> tuple[float, float, float, float]:
    """
    Devuelve μs' de entrada para cada escenario basado en literatura:
    nominal y banda ±1σ de Phan.
    """
    mu_nominal = float(phan_sierra_mu_sp_mm(lambda_nm, mode=mode))
    mu_sd = float(phan_musp_sd_mm(lambda_nm, mode="pchip", allow_extrapolation=False))

    if escenario == "nominal":
        mu_input = mu_nominal
    elif escenario == "menos_1sd":
        mu_input = max(1e-9, mu_nominal - mu_sd)
    elif escenario == "mas_1sd":
        mu_input = mu_nominal + mu_sd
    else:
        raise ValueError(f"Escenario desconocido para μs': {escenario}")

    factor = mu_input / mu_nominal if mu_nominal > 0 else float("nan")
    return mu_input, mu_nominal, mu_sd, factor



def tabla_phan_sierra() -> pd.DataFrame:
    """Separa datos medidos de Phan y aproximaciones usadas en el modelo."""
    ajuste_powerlaw = np.asarray(phan_sierra_mu_sp_mm(PHAN_LAMBDA_NM, mode="powerlaw"), dtype=float)
    ajuste_pchip = np.asarray(phan_sierra_mu_sp_mm(PHAN_LAMBDA_NM, mode="pchip"), dtype=float)
    return pd.DataFrame(
        {
            "lambda_nm": PHAN_LAMBDA_NM,
            "mu_s_prime_phan_mm-1": PHAN_MUSP_MM,
            "mu_s_prime_sd_phan_mm-1": PHAN_MUSP_SD_MM,
            "mu_a_phan_mm-1": PHAN_MUA_MM,
            "mu_s_prime_phan_sierra_powerlaw_mm-1": ajuste_powerlaw,
            "mu_s_prime_phan_sierra_pchip_mm-1": ajuste_pchip,
        }
    )



def leer_resumen_iad_desde_csv(csv_path: Path) -> pd.DataFrame:
    """
    Lee el CSV de salida de IAD.

    Prioriza columnas con nombre; si no existen, usa:
    - lambda en la tercera columna
    - μa en la antepenúltima columna
    """
    df = pd.read_csv(csv_path)
    columnas = list(df.columns)

    if "lambda_nm" not in df.columns:
        if len(columnas) < 3:
            raise ValueError("El CSV no tiene una tercera columna para lambda.")
        df["lambda_nm"] = pd.to_numeric(df.iloc[:, 2], errors="coerce")

    if "mu_a_mm-1" not in df.columns:
        if len(columnas) < 3:
            raise ValueError("El CSV no tiene suficientes columnas para extraer μa antepenúltimo.")
        df["mu_a_mm-1"] = pd.to_numeric(df.iloc[:, -3], errors="coerce")

    if "escenario" not in df.columns:
        df["escenario"] = "nominal"

    if "returncode" in df.columns:
        df = df[pd.to_numeric(df["returncode"], errors="coerce").fillna(-1) == 0]

    if "txt_found" in df.columns:
        df = df[df["txt_found"].astype(str).str.lower().isin(["true", "1"])]

    df["lambda_nm"] = pd.to_numeric(df["lambda_nm"], errors="coerce")
    df["mu_a_mm-1"] = pd.to_numeric(df["mu_a_mm-1"], errors="coerce")
    return df.dropna(subset=["lambda_nm", "mu_a_mm-1"]).copy()



def _interpolar_en_x(x: np.ndarray, y: np.ndarray, x_ref: float, etiqueta: str) -> float:
    """Interpola linealmente un valor de referencia dentro del soporte de x."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size == 0:
        raise ValueError(f"No hay datos para interpolar {etiqueta}.")
    if x_ref < float(np.min(x)) or x_ref > float(np.max(x)):
        raise ValueError(
            f"{etiqueta}={x_ref:.1f} nm cae fuera del rango disponible "
            f"[{float(np.min(x)):.1f}, {float(np.max(x)):.1f}] nm."
        )
    orden = np.argsort(x)
    return float(np.interp(x_ref, x[orden], y[orden]))



def construir_comparacion_mu_a(sub_df: pd.DataFrame) -> pd.DataFrame:
    """Construye tabla comparativa μa_IAD vs μa_Phan sobre las mismas lambdas."""
    df = sub_df.copy().sort_values("lambda_nm")
    lambdas = df["lambda_nm"].to_numpy(dtype=float)
    mu_a_iad = df["mu_a_mm-1"].to_numpy(dtype=float)
    mu_a_phan = np.asarray(phan_mu_a_mm(lambdas), dtype=float)

    ref_iad = _interpolar_en_x(lambdas, mu_a_iad, PHAN_MUA_NORM_LAMBDA_REF_NM, "λ_ref de normalización para μa_IAD")
    ref_phan = float(phan_mu_a_mm(PHAN_MUA_NORM_LAMBDA_REF_NM))
    residuo = mu_a_iad - mu_a_phan

    df["mu_a_iad_mm-1"] = mu_a_iad
    df["mu_a_phan_mm-1"] = mu_a_phan
    df["error_abs_mm-1"] = np.abs(residuo)
    df["error_rel"] = np.abs(residuo) / np.maximum(np.abs(mu_a_phan), 1e-12)
    df["error_rel_pct"] = 100.0 * df["error_rel"]
    df["mu_a_iad_norm"] = mu_a_iad / max(abs(ref_iad), 1e-12)
    df["mu_a_phan_norm"] = mu_a_phan / max(abs(ref_phan), 1e-12)
    df["error_abs_norm"] = np.abs(df["mu_a_iad_norm"] - df["mu_a_phan_norm"])
    return df



def calcular_metricas_mu_a(df_comp: pd.DataFrame) -> dict:
    """Calcula métricas de error en escala absoluta y normalizada."""
    y_true = df_comp["mu_a_phan_mm-1"].to_numpy(dtype=float)
    y_pred = df_comp["mu_a_iad_mm-1"].to_numpy(dtype=float)
    y_true_norm = df_comp["mu_a_phan_norm"].to_numpy(dtype=float)
    y_pred_norm = df_comp["mu_a_iad_norm"].to_numpy(dtype=float)

    def _r2(y_ref: np.ndarray, y_hat: np.ndarray) -> float:
        ss_res = float(np.sum((y_ref - y_hat) ** 2))
        ss_tot = float(np.sum((y_ref - np.mean(y_ref)) ** 2))
        if ss_tot <= 0:
            return float("nan")
        return 1.0 - ss_res / ss_tot

    return {
        "n_puntos": int(len(df_comp)),
        "lambda_min_nm": float(df_comp["lambda_nm"].min()),
        "lambda_max_nm": float(df_comp["lambda_nm"].max()),
        "musp_mode": PHAN_SIERRA_MODE,
        "mua_interp_mode": PHAN_MUA_INTERP_MODE,
        "lambda_ref_norm_nm": PHAN_MUA_NORM_LAMBDA_REF_NM,
        "rmse_mm-1": float(np.sqrt(np.mean((y_pred - y_true) ** 2))),
        "mae_mm-1": float(np.mean(np.abs(y_pred - y_true))),
        "error_rel_medio_pct": float(100.0 * np.mean(np.abs(y_pred - y_true) / np.maximum(np.abs(y_true), 1e-12))),
        "r2": _r2(y_true, y_pred),
        "rmse_norm": float(np.sqrt(np.mean((y_pred_norm - y_true_norm) ** 2))),
        "mae_norm": float(np.mean(np.abs(y_pred_norm - y_true_norm))),
        "r2_norm": _r2(y_true_norm, y_pred_norm),
    }



def g_ma_et_al(lambda_nm: float) -> float:
    """Polinomio de Ma et al. (solo por si se quiere desactivar g fijo)."""
    L = lambda_nm
    g = (
        -5.603
        + 3.61e-2 * L
        - 8.17e-5 * L ** 2
        + 9.51e-8 * L ** 3
        - 5.92e-11 * L ** 4
        + 1.83e-14 * L ** 5
        - 2.11e-18 * L ** 6
    )
    return float(max(0.5, min(0.99, g)))



def construir_rxt_una_lambda(
    wavelength: float,
    reflectance: float,
    header_lines: list[str],
    output_rxt_path: Path,
) -> None:
    """Construye un archivo .rxt nuevo para una sola longitud de onda."""
    with output_rxt_path.open("w", encoding="utf-8", newline="\n") as f:
        for line in header_lines:
            f.write(line + "\n")

        f.write("\n")
        f.write("#lambda\tM_R\n")
        f.write(f"{wavelength:.10f}\t{reflectance:.10f}\n")



def ejecutar_iad(
    iad_exe_path: Path,
    rxt_input_path: Path,
    iad_args: list[str],
    working_dir: Path,
):
    """Ejecuta iad.exe desde Python."""
    cmd = [str(iad_exe_path), *iad_args, rxt_input_path.name]

    if MOSTRAR_COMANDOS:
        print("\nComando que se ejecutará:")
        print(" ".join(cmd))
        print()

    result = subprocess.run(
        cmd,
        cwd=working_dir,
        text=True,
        shell=False,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    return result



def buscar_txt_salida(iad_dir: Path, rxt_output_path: Path) -> Path | None:
    """Intenta encontrar el archivo .txt generado por IAD."""
    txt_esperado = iad_dir / f"{rxt_output_path.stem}.txt"
    if txt_esperado.exists():
        return txt_esperado

    candidatos = sorted(iad_dir.glob(f"{rxt_output_path.stem}*.txt"))
    if candidatos:
        return candidatos[0]

    return None



def extraer_resultado_iad(txt_path: Path):
    """Extrae la fila numérica del archivo .txt de salida de IAD."""
    with txt_path.open("r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue

        partes = stripped.split()

        try:
            float(partes[0])
        except Exception:
            continue

        if len(partes) >= 8:
            return {
                "lambda_nm": float(partes[0]),
                "M_R_measured": float(partes[1]),
                "M_R_fit": float(partes[2]),
                "M_T_measured": float(partes[3]),
                "M_T_fit": float(partes[4]),
                "mu_a_mm-1": float(partes[5]),
                "mu_s_prime_mm-1": float(partes[6]),
                "g": float(partes[7]),
            }

    return None



def leer_csv_mr_temporal(csv_path: Path) -> dict:
    """
    Lee M_R_tiempo_data.csv y agrupa por medicion.
    Columnas esperadas: medicion, tiempo, lambda, reflectancia
    """
    df = pd.read_csv(csv_path)

    columnas_esperadas = {"medicion", "tiempo", "lambda", "reflectancia"}
    if not columnas_esperadas.issubset(df.columns):
        raise ValueError(
            f"El CSV temporal debe contener {columnas_esperadas}. "
            f"Columnas encontradas: {list(df.columns)}"
        )

    mediciones = {}
    for medicion_id, grupo in df.groupby("medicion", sort=True):
        tiempo = float(grupo["tiempo"].iloc[0])
        espectro = sorted(zip(grupo["lambda"], grupo["reflectancia"]), key=lambda x: x[0])
        mediciones[int(medicion_id)] = {"tiempo": tiempo, "espectro": list(espectro)}

    return mediciones



def _progreso(completadas: int, total: int, t_inicio: float, ancho_barra: int = 30):
    """Imprime barra de progreso con % y ETA en una sola línea."""
    frac = completadas / total
    llenas = int(ancho_barra * frac)
    barra = "█" * llenas + "░" * (ancho_barra - llenas)

    elapsed = time.time() - t_inicio
    if completadas > 0:
        eta_s = elapsed / completadas * (total - completadas)
        if eta_s >= 3600:
            eta_str = f"{eta_s / 3600:.1f}h"
        elif eta_s >= 60:
            eta_str = f"{eta_s / 60:.0f}m {eta_s % 60:.0f}s"
        else:
            eta_str = f"{eta_s:.0f}s"
    else:
        eta_str = "---"

    linea = f"\r  [{barra}] {frac:6.1%}  ({completadas}/{total})  ETA: {eta_str}   "
    sys.stdout.write(linea)
    sys.stdout.flush()
    if completadas == total:
        sys.stdout.write("\n")



def limpiar_carpeta_por_lambda(carpeta: Path):
    """Limpia archivos viejos dentro de la carpeta por_lambda."""
    if not carpeta.exists():
        return

    for archivo in carpeta.iterdir():
        if archivo.is_file():
            archivo.unlink()



def graficar_phan_sierra(ruta_png: Path, mostrar: bool = False):
    """Grafica datos Palm de Phan y las dos versiones de Phan-Sierra."""
    lambda_plot = np.linspace(float(PHAN_LAMBDA_NM.min()), float(PHAN_LAMBDA_NM.max()), 400)
    mu_powerlaw = np.asarray(phan_sierra_mu_sp_mm(lambda_plot, mode="powerlaw"), dtype=float)
    mu_pchip = np.asarray(phan_sierra_mu_sp_mm(lambda_plot, mode="pchip"), dtype=float)

    plt.figure(figsize=(8, 5))
    plt.plot(lambda_plot, mu_powerlaw, label="Phan-Sierra (powerlaw)", lw=2.0, alpha=0.9)
    plt.plot(lambda_plot, mu_pchip, label="Phan-Sierra (pchip)", lw=2.2, alpha=0.9)
    plt.errorbar(
        PHAN_LAMBDA_NM,
        PHAN_MUSP_MM,
        yerr=PHAN_MUSP_SD_MM,
        fmt="o",
        ms=6,
        capsize=4,
        label="Palm de Phan et al. (media ± sd)",
    )
    plt.xlabel("Longitud de onda (nm)")
    plt.ylabel("μs' (1/mm)")
    plt.title("Función Phan-Sierra para μs'(λ)")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    if mostrar:
        plt.show()
    plt.savefig(ruta_png, dpi=300, bbox_inches="tight")
    plt.close()



def graficar_mu_sp_escenarios(df_mr: pd.DataFrame, ruta_png: Path):
    """Grafica μs' nominal y banda ±1σ de Phan sobre las lambdas medidas."""
    lambdas = np.sort(df_mr["wavelength_nm"].astype(float).unique())
    mu_nom = np.asarray(phan_sierra_mu_sp_mm(lambdas), dtype=float)
    mu_sd = np.asarray(phan_musp_sd_mm(lambdas), dtype=float)

    plt.figure(figsize=(10, 6))
    plt.plot(lambdas, mu_nom, linewidth=2.0, label=f"Nominal ({PHAN_SIERRA_MODE})")
    plt.plot(lambdas, np.maximum(mu_nom - mu_sd, 1e-9), linewidth=1.6, linestyle="--", label="Nominal - 1sd Phan")
    plt.plot(lambdas, mu_nom + mu_sd, linewidth=1.6, linestyle="--", label="Nominal + 1sd Phan")
    plt.fill_between(
        lambdas,
        np.maximum(mu_nom - mu_sd, 1e-9),
        mu_nom + mu_sd,
        alpha=0.15,
        label="Banda ±1sd de Phan",
    )

    plt.xlabel("Longitud de onda (nm)")
    plt.ylabel("μs' de entrada (1/mm)")
    plt.title("Sensibilidad de μs'(λ) basada en Phan")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(ruta_png, dpi=300, bbox_inches="tight")
    plt.close()



def graficar_mu_a_sensibilidad(df_resultados: pd.DataFrame, ruta_png: Path):
    """Grafica μa recuperado para los escenarios nominal y ±1σ."""
    df_plot = df_resultados.dropna(subset=["lambda_nm", "mu_a_mm-1", "escenario"]).copy()
    if df_plot.empty:
        return

    plt.figure(figsize=(10, 6))
    for escenario in ["menos_1sd", "nominal", "mas_1sd"]:
        sub = df_plot[df_plot["escenario"] == escenario].sort_values("lambda_nm")
        if not sub.empty:
            plt.plot(sub["lambda_nm"], sub["mu_a_mm-1"], marker="o", markersize=3.5, linewidth=1.5, label=escenario)

    plt.xlabel("Longitud de onda (nm)")
    plt.ylabel("μa recuperado (1/mm)")
    plt.title("Sensibilidad de μa a cambios en μs'(λ)")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(ruta_png, dpi=300, bbox_inches="tight")
    plt.close()



def graficar_comparacion_mu_a(df_comp: pd.DataFrame, ruta_png: Path, escenario: str):
    """Grafica μa_IAD contra μa_Phan interpolado a las mismas lambdas."""
    df_plot = df_comp.sort_values("lambda_nm")

    plt.figure(figsize=(10, 6))
    plt.plot(df_plot["lambda_nm"], df_plot["mu_a_iad_mm-1"], marker="o", markersize=3.5, linewidth=1.7, label="μa IAD")
    plt.plot(df_plot["lambda_nm"], df_plot["mu_a_phan_mm-1"], linewidth=2.1, label="μa Palm de Phan interpolado")
    plt.xlabel("Longitud de onda (nm)")
    plt.ylabel("μa (1/mm)")
    plt.title(f"μa_IAD vs μa_Phan ({escenario})")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(ruta_png, dpi=300, bbox_inches="tight")
    plt.close()



def graficar_comparacion_mu_a_normalizada(df_comp: pd.DataFrame, ruta_png: Path, escenario: str):
    """Grafica comparación normalizada de μa para separar forma y escala."""
    df_plot = df_comp.sort_values("lambda_nm")

    plt.figure(figsize=(10, 6))
    plt.plot(df_plot["lambda_nm"], df_plot["mu_a_iad_norm"], marker="o", markersize=3.5, linewidth=1.7, label="μa IAD normalizado")
    plt.plot(df_plot["lambda_nm"], df_plot["mu_a_phan_norm"], linewidth=2.1, label="μa Phan normalizado")
    plt.xlabel("Longitud de onda (nm)")
    plt.ylabel(f"μa / μa({PHAN_MUA_NORM_LAMBDA_REF_NM:.0f} nm)")
    plt.title(f"Comparación normalizada de μa ({escenario})")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(ruta_png, dpi=300, bbox_inches="tight")
    plt.close()



def generar_comparaciones_mu_a_desde_csv(csv_path: Path, iad_dir: Path) -> Path | None:
    """Lee el resumen de salida y genera comparación μa_IAD vs μa_Phan."""
    df_resumen = leer_resumen_iad_desde_csv(csv_path)
    if df_resumen.empty:
        print(f"[Comparacion μa] Sin filas válidas en {csv_path}")
        return None

    if "medicion" in df_resumen.columns:
        print("[Comparacion μa] CSV temporal detectado: la comparación automática se omite en este script.")
        return None

    metricas = []
    for escenario, sub in df_resumen.groupby("escenario", sort=True):
        df_comp = construir_comparacion_mu_a(sub)
        nombre_seguro = str(escenario).replace(" ", "_")
        ruta_comp = iad_dir / f"comparacion_mu_a_phan_{nombre_seguro}.csv"
        df_comp.to_csv(ruta_comp, index=False)
        graficar_comparacion_mu_a(df_comp, iad_dir / f"grafica_mu_a_iad_vs_phan_{nombre_seguro}.png", str(escenario))
        graficar_comparacion_mu_a_normalizada(
            df_comp,
            iad_dir / f"grafica_mu_a_normalizada_{nombre_seguro}.png",
            str(escenario),
        )

        fila_metricas = calcular_metricas_mu_a(df_comp)
        fila_metricas["escenario"] = escenario
        metricas.append(fila_metricas)

    if not metricas:
        return None

    ruta_metricas = iad_dir / "metricas_mu_a_phan_por_escenario.csv"
    pd.DataFrame(metricas).sort_values("escenario").to_csv(ruta_metricas, index=False)
    return ruta_metricas



def _correr_iad_una_lambda(
    wavelength: float,
    reflectance: float,
    header_lines: list,
    per_lambda_dir: Path,
    iad_exe_path: Path,
    escenario: str,
    factor_musp: float,
    med_id: int = 0,
) -> dict:
    """Ejecuta IAD para una sola (λ, reflectancia) y devuelve el dict de resultado."""
    reflectance = max(reflectance, 1e-4)
    mu_sp_mm, mu_sp_nominal, mu_sp_sd, factor_aplicado = mu_sp_escenario_phan_mm(
        wavelength,
        escenario=escenario,
        mode=PHAN_SIERRA_MODE,
    )
    g_valor = G_FIJO if USAR_G_FIJO else g_ma_et_al(wavelength)

    if med_id != 0:
        base_name = f"{escenario}_med{med_id:06d}_lambda_{wavelength:.2f}".replace(".", "p")
    else:
        base_name = f"{escenario}_lambda_{wavelength:.2f}".replace(".", "p")
    rxt_output_path = per_lambda_dir / f"{base_name}.rxt"

    construir_rxt_una_lambda(
        wavelength=wavelength,
        reflectance=reflectance,
        header_lines=header_lines,
        output_rxt_path=rxt_output_path,
    )

    iad_args = ["-g", f"{g_valor:.6f}", "-j", f"{mu_sp_mm:.6f}"]
    if USAR_MODO_RAPIDO:
        iad_args.extend(["-M", "0"])

    result = ejecutar_iad(
        iad_exe_path=iad_exe_path,
        rxt_input_path=rxt_output_path,
        iad_args=iad_args,
        working_dir=per_lambda_dir,
    )

    txt_output_path = buscar_txt_salida(per_lambda_dir, rxt_output_path)

    fila = {
        "escenario": escenario,
        "factor_musp": factor_aplicado,
        "lambda_nm": wavelength,
        "reflectance_input": reflectance,
        "mu_s_prime_nominal_mm-1": mu_sp_nominal,
        "mu_s_prime_sd_phan_mm-1": mu_sp_sd,
        "mu_s_prime_input_mm-1": mu_sp_mm,
        "g_input": g_valor,
        "n_fijo": N_TEJIDO,
        "returncode": result.returncode,
        "txt_found": txt_output_path is not None,
    }

    if txt_output_path is not None:
        extraido = extraer_resultado_iad(txt_output_path)
        if extraido is not None:
            fila.update(extraido)
        txt_output_path.unlink(missing_ok=True)

    rxt_output_path.unlink(missing_ok=True)
    return fila



def main():
    # ------------------------------------------------------------
    # 1. Validar rutas comunes
    # ------------------------------------------------------------
    if not RXT_TEMPLATE_PATH.exists():
        raise FileNotFoundError(f"No existe la plantilla RXT: {RXT_TEMPLATE_PATH}")
    if not IAD_EXE_PATH.exists():
        raise FileNotFoundError(f"No existe iad.exe: {IAD_EXE_PATH}")

    # ------------------------------------------------------------
    # 2. Crear carpetas de trabajo
    # ------------------------------------------------------------
    base_dir = Path(__file__).parent
    iad_dir = base_dir / IAD_FOLDER_NAME
    iad_dir.mkdir(exist_ok=True)

    per_lambda_dir = iad_dir / IAD_PER_LAMBDA_FOLDER_NAME
    per_lambda_dir.mkdir(exist_ok=True)

    limpiar_carpeta_por_lambda(per_lambda_dir)

    # ------------------------------------------------------------
    # 3. Leer encabezado de plantilla .rxt
    # ------------------------------------------------------------
    header_lines = extraer_encabezado_rxt(RXT_TEMPLATE_PATH)

    # ------------------------------------------------------------
    # 4. Guardar documentación del ajuste Phan-Sierra
    # ------------------------------------------------------------
    print(f"[Phan-Sierra] modo = {PHAN_SIERRA_MODE}")
    print(f"[Phan-Sierra] A = {PHAN_SIERRA_A_MM:.8f} mm^-1")
    print(f"[Phan-Sierra] B = {PHAN_SIERRA_B:.8f}")
    print(f"[Phan-Sierra] lambda_ref = {PHAN_SIERRA_LAMBDA_REF_NM:.1f} nm")
    print(f"[Phan-Sierra] extrapolacion permitida = {PHAN_SIERRA_ALLOW_EXTRAPOLATION}")
    print(f"[Phan-Sierra] escenarios = {', '.join(ESCENARIOS_MUSP)}")

    ruta_tabla_phan_sierra = iad_dir / "tabla_phan_sierra.csv"
    tabla_phan_sierra().to_csv(ruta_tabla_phan_sierra, index=False)
    graficar_phan_sierra(iad_dir / "grafica_phan_sierra.png")

    # ============================================================
    # MODO NORMAL — un solo espectro (M_R_data.csv)
    # ============================================================
    if not MODO_TEMPORAL:
        if not MR_CSV_PATH.exists():
            raise FileNotFoundError(f"No existe el archivo CSV: {MR_CSV_PATH}")

        df_mr = leer_csv_mr(MR_CSV_PATH)
        graficar_mu_sp_escenarios(df_mr, iad_dir / "grafica_mu_sp_escenarios.png")

        resultados = []
        total = len(df_mr) * len(ESCENARIOS_MUSP)
        t_inicio = time.time()
        contador = 0

        for escenario in ESCENARIOS_MUSP:
            for _, row in df_mr.iterrows():
                wavelength = float(row["wavelength_nm"])
                reflectance = float(row["reflectance"])

                fila = _correr_iad_una_lambda(
                    wavelength,
                    reflectance,
                    header_lines,
                    per_lambda_dir,
                    IAD_EXE_PATH,
                    escenario=escenario,
                    factor_musp=1.0,
                )
                resultados.append(fila)
                contador += 1
                _progreso(contador, total, t_inicio)

        df_resultados = pd.DataFrame(resultados).sort_values(["escenario", "lambda_nm"])
        ruta_resumen = iad_dir / "resumen_resultados_phan_sierra.csv"
        df_resultados.to_csv(ruta_resumen, index=False)

        graficar_mu_a_sensibilidad(df_resultados, iad_dir / "grafica_mu_a_sensibilidad.png")
        ruta_metricas = generar_comparaciones_mu_a_desde_csv(ruta_resumen, iad_dir)

        print("\nEjecución terminada.")
        print(f"Tabla Phan-Sierra guardada en: {ruta_tabla_phan_sierra}")
        print(f"Resumen guardado en: {ruta_resumen}")
        if ruta_metricas is not None:
            print(f"Métricas μa vs Phan guardadas en: {ruta_metricas}")
        print(f"Gráficas guardadas en: {iad_dir}")

    # ============================================================
    # MODO TEMPORAL — serie temporal (M_R_tiempo_data.csv)
    # ============================================================
    else:
        csv_temporal_path = base_dir / CSV_TEMPORAL_NAME
        if not csv_temporal_path.exists():
            raise FileNotFoundError(f"No existe el CSV temporal: {csv_temporal_path}")

        mediciones = leer_csv_mr_temporal(csv_temporal_path)

        ids_medicion = sorted(mediciones.keys())
        ids_medicion = [m for m in ids_medicion if m >= INICIO_MEDICION]
        if MAX_MEDICIONES is not None:
            ids_medicion = ids_medicion[:MAX_MEDICIONES]

        if VENTANA_PROMEDIO > 1:
            grupos = [ids_medicion[i:i + VENTANA_PROMEDIO] for i in range(0, len(ids_medicion), VENTANA_PROMEDIO)]
            mediciones_prom = {}
            for grupo in grupos:
                rep_id = grupo[0]
                tiempo_prom = sum(mediciones[m]["tiempo"] for m in grupo) / len(grupo)
                wls = [wl for wl, _ in mediciones[grupo[0]]["espectro"]]
                n = len(grupo)
                refs_prom = [sum(mediciones[m]["espectro"][j][1] for m in grupo) / n for j in range(len(wls))]
                mediciones_prom[rep_id] = {
                    "tiempo": tiempo_prom,
                    "espectro": list(zip(wls, refs_prom)),
                }
            ids_medicion = [g[0] for g in grupos]
            mediciones = mediciones_prom

        total_med = len(ids_medicion)

        tareas = [
            (escenario, med_id, mediciones[med_id]["tiempo"], float(wl), float(ref))
            for escenario in ESCENARIOS_MUSP
            for med_id in ids_medicion
            for wl, ref in mediciones[med_id]["espectro"]
        ]
        total_tareas = len(tareas)

        print(
            f"Mediciones a procesar: {total_med}  |  "
            f"Tareas IAD totales: {total_tareas}  |  Workers: {WORKERS}"
        )

        def _tarea(escenario, med_id, tiempo, wl, ref):
            fila = _correr_iad_una_lambda(
                wl,
                ref,
                header_lines,
                per_lambda_dir,
                IAD_EXE_PATH,
                escenario=escenario,
                factor_musp=1.0,
                med_id=med_id,
            )
            fila["medicion"] = med_id
            fila["tiempo"] = tiempo
            return fila

        resultados = []
        t_inicio = time.time()

        if WORKERS <= 1:
            for i, tarea in enumerate(tareas, 1):
                resultados.append(_tarea(*tarea))
                _progreso(i, total_tareas, t_inicio)
        else:
            with ThreadPoolExecutor(max_workers=WORKERS) as executor:
                futuros = {executor.submit(_tarea, *t): t for t in tareas}
                completadas = 0
                for fut in as_completed(futuros):
                    resultados.append(fut.result())
                    completadas += 1
                    _progreso(completadas, total_tareas, t_inicio)

        df_resultados = pd.DataFrame(resultados)
        cols_primeras = [
            "escenario",
            "factor_musp",
            "medicion",
            "tiempo",
            "lambda_nm",
            "reflectance_input",
            "mu_a_mm-1",
            "mu_s_prime_mm-1",
            "g",
        ]
        cols_resto = [c for c in df_resultados.columns if c not in cols_primeras]
        df_resultados = df_resultados[cols_primeras + cols_resto]
        df_resultados = df_resultados.sort_values(["escenario", "medicion", "lambda_nm"])

        ruta_resumen = base_dir / CSV_SALIDA_TEMPORAL
        ruta_resumen.parent.mkdir(exist_ok=True)
        df_resultados.to_csv(ruta_resumen, index=False)

        print("\nEjecución temporal terminada.")
        print(f"Mediciones procesadas: {total_med}")
        print(f"Resumen temporal guardado en: {ruta_resumen}")

    print(f"Carpeta principal IAD: {iad_dir}")
    print(f"Archivos por lambda: {per_lambda_dir}")


if __name__ == "__main__":
    main()
