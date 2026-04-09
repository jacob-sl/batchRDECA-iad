from pathlib import Path
import json
import os
import subprocess
import sys
import time
import unicodedata
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")   # backend no interactivo: evita errores de Tk en hilos worker
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# ============================================================
# CONFIGURACIÓN DEL SISTEMA
# ============================================================
#None = todos los sujetos o mediciones; int = número máximo a procesar.
INICIO_MEDICION = None #Espectro a partir del cual se comenzará el procesado.
MAX_MEDICIONES = None #En caso de que no se quiera procesar todas las mediciones, se puede limitar el número máximo. 
VENTANA_PROMEDIO = 50 # Se promedian estos espectros, reduce resolución temporal pero mejora Signal to Noise Ratio. Solo para modo temporal.

WORKERS = os.cpu_count() - 1 or 4 #Hilos a usar para paralelización del IAD

# ============================================================
# CONFIGURACIÓN DELMODELO ÓPTICO
# ============================================================

# ── TOGGLE PRINCIPAL ──
MODO = "single"   # "single"    → un espectro (M_R_data.csv)
                     # "temporal"  → serie temporal (M_R_tiempo_data.csv)
                     # "batch_all" → procesa TODOS los sujetos en Mediciones/

# Truncado opcional del M_R_data.csv single antes de correr IAD.
# Rango inclusivo en nm: LAMBDA_MIN <= wavelength_nm <= LAMBDA_MAX.
TRUNCAR_MR_SINGLE_POR_LAMBDA = True
MR_SINGLE_LAMBDA_MIN_NM = 430.0
MR_SINGLE_LAMBDA_MAX_NM = 680.0

# g fijo y N fijo.
USAR_G_FIJO = True
G_FIJO = 0.8
N_TEJIDO = 1.41  # Solo genera el formato; IAD toma la refracción del .rxt/header
PHAN_UBICACION = "Palm"
PHAN_UBICACION_DEFAULT = PHAN_UBICACION
# Opciones: "Forehead", "Cheek", "Ventral Forearm", "Palm", "Back",
#           "Upper Arm", "Dorsal Forearm", "Neck", "Shin", "Chest"
# Modo rápido omite el "sanity check" que IAD hace con una simulación Montecarlo.
USAR_MODO_RAPIDO = True
# Mostrar o no el comando completo por cada lambda
MOSTRAR_COMANDOS = False

# Modelo de μs'(λ): Phan-Sierra
PHAN_SIERRA_MODE = "pchip"   # opciones: "powerlaw" o "pchip"
PHAN_SIERRA_ALLOW_EXTRAPOLATION = True
PHAN_MUSP_SD_ALLOW_EXTRAPOLATION = True
PHAN_SIERRA_LAMBDA_REF_NM = 526.0  # λ_ref del ajuste por ley de potencia

# Referencia de μa(λ) de Phan usada solo en postanálisis
# (métricas y overlay visual). No entra como input a IAD.
PHAN_MUA_INTERP_MODE = "pchip"
PHAN_MUA_ALLOW_EXTRAPOLATION = True
PHAN_MUA_NORM_LAMBDA_REF_NM = 526.0  # λ_ref para normalizar μa_IAD y μa_Phan
# Sensibilidad basada en literatura: escenarios de μs' a computar.
# Ejemplos:
# ESCENARIOS_MUSP = ("nominal",)
# ESCENARIOS_MUSP = ("menos_1sd", "nominal", "mas_1sd")
ESCENARIOS_MUSP_VALIDOS = ("menos_1sd", "nominal", "mas_1sd")
ESCENARIOS_MUSP = ("nominal",)
# Escenarios de incertidumbre de medición (σ_M_R por scan).
# Requiere que M_R_data.csv tenga columna reflectance_std (generada por single_adq.py).
# Si no existe esa columna, solo "nominal" es válido.
ESCENARIOS_MR_STD_VALIDOS = ("menos_1sd", "nominal", "mas_1sd")
ESCENARIOS_MR_STD = ("nominal",)
GENERAR_GRAFICA_PHAN_SIERRA_PNG = True
GENERAR_GRAFICA_PHAN_SIERRA_PDF = True
GENERAR_DASHBOARD_MUA_PNG = True
GENERAR_DASHBOARD_MUA_PDF = True

# ============================================================
# GESTIÓN DE CARPETAS
# ============================================================
# Nombres de archivos CSV dentro de series/
CSV_SINGLE_NAME = "M_R_data.csv"
CSV_TEMPORAL_NAME = "M_R_tiempo_data.csv"
# Carpeta de mediciones
MEDICIONES_DIR = Path(__file__).parent / "Mediciones"
# Selección automática del último sujeto (True) o selector interactivo (False)
AUTO_ULTIMO_SUJETO = False
# Ruta al archivo plantilla .rxt
RXT_TEMPLATE_PATH = Path(__file__).parent / "sample-F.rxt"
# Ruta al ejecutable de IAD
IAD_EXE_PATH = Path(__file__).parent / "IADSCOTT" / "iad.exe"
# Carpeta principal de trabajo
IAD_FOLDER_NAME = "IAD_run"
# Subcarpeta donde se guardará una corrida por longitud de onda
IAD_PER_LAMBDA_FOLDER_NAME = "por_lambda"

# Mapeo controlado entre carpetas de adquisicion single y datos Phan-Sierra.
# "pulgar_izquierdo" se mapea a Palm por decision experimental: piel glabra.
REGIONES_ANATOMICAS_PHAN = {
    "pulgar_izquierdo": "Palm",
    "palma": "Palm",
    "antebrazo_dorsal": "Dorsal Forearm",
    "antebrazo_ventral": "Ventral Forearm",
}

# ============================================================
# FUNCIÓN PHAN-SIERRA para μs'(λ)
# ============================================================
# Datos de μs' y μa medidos por SFDI en 10 ubicaciones anatómicas.
# Fuente: Phan et al., JBO 26(2) 026001, Supplementary Tables 2 & 3.
# λ (nm): 471, 526, 591, 621, 659, 691, 731, 851
# Nota: para el ajuste por ley de potencia se excluye el último punto (851 nm).
# Se conserva en las tablas y en la interpolación PCHIP como referencia literaria.

PHAN_DATOS = {
    "Forehead":        {"musp": [2.12, 1.75, 1.56, 1.53, 1.46, 1.46, 1.51, 1.65],
                        "musp_sd": [1.47, 0.812, 0.518, 0.42, 0.346, 0.29, 0.242, 0.174],
                        "mua": [0.848, 0.533, 0.238, 0.128, 0.097, 0.0707, 0.0513, 0.0282]},
    "Cheek":           {"musp": [1.82, 1.63, 1.46, 1.43, 1.36, 1.37, 1.42, 1.53],
                        "musp_sd": [0.987, 0.711, 0.491, 0.405, 0.34, 0.286, 0.241, 0.173],
                        "mua": [0.688, 0.47, 0.211, 0.116, 0.0871, 0.0637, 0.0459, 0.022]},
    "Ventral Forearm": {"musp": [2.25, 1.93, 1.63, 1.54, 1.45, 1.41, 1.42, 1.46],
                        "musp_sd": [1.0, 0.71, 0.415, 0.333, 0.269, 0.216, 0.169, 0.115],
                        "mua": [0.644, 0.424, 0.201, 0.111, 0.0841, 0.0619, 0.0451, 0.0255]},
    "Palm":            {"musp": [2.57, 2.18, 1.78, 1.72, 1.60, 1.54, 1.51, 1.45],
                        "musp_sd": [0.421, 0.278, 0.147, 0.126, 0.108, 0.0984, 0.0883, 0.0813],
                        "mua": [0.317, 0.215, 0.092, 0.0288, 0.0204, 0.0147, 0.0111, 0.00978]},
    "Back":            {"musp": [1.64, 1.47, 1.36, 1.35, 1.30, 1.30, 1.34, 1.42],
                        "musp_sd": [0.99, 0.736, 0.499, 0.425, 0.363, 0.306, 0.246, 0.154],
                        "mua": [0.674, 0.469, 0.242, 0.149, 0.113, 0.0828, 0.0587, 0.0236]},
    "Upper Arm":       {"musp": [2.23, 1.90, 1.58, 1.48, 1.39, 1.36, 1.37, 1.41],
                        "musp_sd": [1.05, 0.739, 0.418, 0.326, 0.263, 0.209, 0.167, 0.12],
                        "mua": [0.578, 0.378, 0.173, 0.0957, 0.0717, 0.0531, 0.0388, 0.0213]},
    "Dorsal Forearm":  {"musp": [1.70, 1.51, 1.37, 1.34, 1.29, 1.28, 1.31, 1.38],
                        "musp_sd": [0.989, 0.726, 0.456, 0.366, 0.298, 0.239, 0.188, 0.12],
                        "mua": [0.752, 0.504, 0.249, 0.15, 0.116, 0.0855, 0.0614, 0.0305]},
    "Neck":            {"musp": [1.74, 1.60, 1.42, 1.38, 1.30, 1.27, 1.29, 1.31],
                        "musp_sd": [0.876, 0.693, 0.447, 0.374, 0.314, 0.246, 0.196, 0.117],
                        "mua": [0.549, 0.363, 0.18, 0.111, 0.0807, 0.0581, 0.0401, 0.0162]},
    "Shin":            {"musp": [1.85, 1.59, 1.40, 1.34, 1.27, 1.24, 1.25, 1.30],
                        "musp_sd": [0.945, 0.683, 0.419, 0.33, 0.273, 0.222, 0.179, 0.141],
                        "mua": [0.591, 0.391, 0.204, 0.123, 0.0942, 0.0699, 0.0511, 0.0275]},
    "Chest":           {"musp": [1.96, 1.71, 1.45, 1.38, 1.29, 1.25, 1.25, 1.28],
                        "musp_sd": [1.06, 0.773, 0.44, 0.344, 0.28, 0.232, 0.188, 0.141],
                        "mua": [0.472, 0.322, 0.158, 0.0883, 0.0649, 0.047, 0.033, 0.0163]},
}

PHAN_LAMBDA_NM  = np.array([471.0, 526.0, 591.0, 621.0, 659.0, 691.0, 731.0, 851.0], dtype=float)
PHAN_MUSP_MM    = np.array([], dtype=float)
PHAN_MUSP_SD_MM = np.array([], dtype=float)
PHAN_MUA_MM     = np.array([], dtype=float)
PHAN_POWERLAW_LAMBDA_NM = np.array([], dtype=float)
PHAN_POWERLAW_MUSP_MM = np.array([], dtype=float)
PHAN_SIERRA_A_MM = float("nan")
PHAN_SIERRA_B = float("nan")
_phan_sierra_pchip = None
_phan_musp_sd_pchip = None
_phan_mua_pchip = None


def configurar_phan_ubicacion(ubicacion: str) -> None:
    """Actualiza el contexto Phan-Sierra global para una ubicacion anatomica."""
    global PHAN_UBICACION
    global PHAN_MUSP_MM, PHAN_MUSP_SD_MM, PHAN_MUA_MM
    global PHAN_POWERLAW_LAMBDA_NM, PHAN_POWERLAW_MUSP_MM
    global PHAN_SIERRA_A_MM, PHAN_SIERRA_B
    global _phan_sierra_pchip, _phan_musp_sd_pchip, _phan_mua_pchip

    if ubicacion not in PHAN_DATOS:
        raise ValueError(
            f"PHAN_UBICACION='{ubicacion}' no valida. "
            f"Opciones: {list(PHAN_DATOS)}"
        )

    PHAN_UBICACION = ubicacion
    _ub = PHAN_DATOS[ubicacion]
    PHAN_MUSP_MM = np.array(_ub["musp"], dtype=float)
    PHAN_MUSP_SD_MM = np.array(_ub["musp_sd"], dtype=float)
    PHAN_MUA_MM = np.array(_ub["mua"], dtype=float)
    PHAN_POWERLAW_LAMBDA_NM = PHAN_LAMBDA_NM[:-1]
    PHAN_POWERLAW_MUSP_MM = PHAN_MUSP_MM[:-1]

    _x = np.log(PHAN_POWERLAW_LAMBDA_NM / PHAN_SIERRA_LAMBDA_REF_NM)
    _y = np.log(PHAN_POWERLAW_MUSP_MM)
    _slope, _intercept = np.polyfit(_x, _y, 1)
    PHAN_SIERRA_A_MM = float(np.exp(_intercept))
    PHAN_SIERRA_B = float(-_slope)

    try:
        from scipy.interpolate import PchipInterpolator
        _phan_sierra_pchip = PchipInterpolator(PHAN_LAMBDA_NM, PHAN_MUSP_MM, extrapolate=True)
        _phan_musp_sd_pchip = PchipInterpolator(PHAN_LAMBDA_NM, PHAN_MUSP_SD_MM, extrapolate=True)
        _phan_mua_pchip = PchipInterpolator(PHAN_LAMBDA_NM, PHAN_MUA_MM, extrapolate=True)
    except Exception:
        _phan_sierra_pchip = None
        _phan_musp_sd_pchip = None
        _phan_mua_pchip = None

if isinstance(ESCENARIOS_MUSP, str):
    ESCENARIOS_MUSP = (ESCENARIOS_MUSP,)   # "nominal" → ("nominal",)

escenarios_invalidos = [e for e in ESCENARIOS_MUSP if e not in ESCENARIOS_MUSP_VALIDOS]
if escenarios_invalidos:
    raise ValueError(
        "ESCENARIOS_MUSP contiene valores no válidos: "
        f"{escenarios_invalidos}. "
        f"Opciones: {ESCENARIOS_MUSP_VALIDOS}"
    )
if not ESCENARIOS_MUSP:
    raise ValueError("ESCENARIOS_MUSP no puede estar vacío.")

if isinstance(ESCENARIOS_MR_STD, str):
    ESCENARIOS_MR_STD = (ESCENARIOS_MR_STD,)
_inv_mr = [e for e in ESCENARIOS_MR_STD if e not in ESCENARIOS_MR_STD_VALIDOS]
if _inv_mr:
    raise ValueError(
        f"ESCENARIOS_MR_STD inválidos: {_inv_mr}. "
        f"Opciones: {ESCENARIOS_MR_STD_VALIDOS}"
    )
if not ESCENARIOS_MR_STD:
    raise ValueError("ESCENARIOS_MR_STD no puede estar vacío.")

configurar_phan_ubicacion(PHAN_UBICACION_DEFAULT)

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

    cols = ["wavelength_nm", "reflectance"]
    if "reflectance_std" in df.columns:
        cols.append("reflectance_std")
    else:
        df["reflectance_std"] = 0.0
        cols.append("reflectance_std")

    df = df[cols].copy()
    df = df.dropna()
    df = truncar_df_mr_single_por_lambda(df, csv_path)
    return df


def truncar_df_mr_single_por_lambda(df: pd.DataFrame, csv_path: Path) -> pd.DataFrame:
    """Aplica truncado opcional por longitud de onda al M_R_data.csv single."""
    if not TRUNCAR_MR_SINGLE_POR_LAMBDA:
        return df

    lambda_min = float(MR_SINGLE_LAMBDA_MIN_NM)
    lambda_max = float(MR_SINGLE_LAMBDA_MAX_NM)
    if lambda_min >= lambda_max:
        raise ValueError(
            "Rango invalido para truncar M_R single: "
            f"MR_SINGLE_LAMBDA_MIN_NM={lambda_min} debe ser menor que "
            f"MR_SINGLE_LAMBDA_MAX_NM={lambda_max}."
        )

    n_original = len(df)
    df_truncado = df[
        (df["wavelength_nm"] >= lambda_min)
        & (df["wavelength_nm"] <= lambda_max)
    ].copy()

    if df_truncado.empty:
        lambda_csv_min = float(df["wavelength_nm"].min()) if not df.empty else float("nan")
        lambda_csv_max = float(df["wavelength_nm"].max()) if not df.empty else float("nan")
        raise ValueError(
            f"El truncado [{lambda_min:.2f}, {lambda_max:.2f}] nm dejó sin puntos a {csv_path}. "
            f"Rango disponible en CSV: [{lambda_csv_min:.2f}, {lambda_csv_max:.2f}] nm."
        )

    lambda_trunc_min = float(df_truncado["wavelength_nm"].min())
    lambda_trunc_max = float(df_truncado["wavelength_nm"].max())
    if not (lambda_trunc_min <= PHAN_MUA_NORM_LAMBDA_REF_NM <= lambda_trunc_max):
        raise ValueError(
            f"El truncado efectivo [{lambda_trunc_min:.2f}, {lambda_trunc_max:.2f}] nm "
            f"excluye PHAN_MUA_NORM_LAMBDA_REF_NM={PHAN_MUA_NORM_LAMBDA_REF_NM:.2f} nm, "
            "que se usa para normalizar las métricas de μa. "
            "Ajusta el rango de truncado o PHAN_MUA_NORM_LAMBDA_REF_NM."
        )

    print(
        f"  [M_R] Truncado por lambda: {lambda_min:.2f}-{lambda_max:.2f} nm "
        f"({len(df_truncado)}/{n_original} puntos) en {csv_path.name}"
    )
    return df_truncado


def normalizar_region_id(texto: str) -> str:
    """Normaliza nombres de carpeta a ids estables para mapear anatomia."""
    texto_ascii = unicodedata.normalize("NFKD", texto).encode("ascii", "ignore").decode("ascii")
    texto_limpio = "".join(ch.lower() if ch.isalnum() else "_" for ch in texto_ascii)
    return "_".join(parte for parte in texto_limpio.split("_") if parte)


def phan_ubicacion_para_carpeta_medicion(carpeta_medicion: Path) -> str:
    """
    Resuelve la ubicacion Phan-Sierra para una carpeta anatomica single.

    Primero usa metadata_medicion.json si existe; si no, usa el nombre de carpeta.
    """
    ruta_metadata = carpeta_medicion / "metadata_medicion.json"
    if ruta_metadata.exists():
        with ruta_metadata.open("r", encoding="utf-8") as f:
            metadata = json.load(f)
        phan_ubicacion = metadata.get("phan_ubicacion")
        if phan_ubicacion in PHAN_DATOS:
            return phan_ubicacion
        raise ValueError(
            f"{ruta_metadata} contiene phan_ubicacion='{phan_ubicacion}', "
            f"pero no existe en PHAN_DATOS. Opciones: {list(PHAN_DATOS)}"
        )

    region_id = normalizar_region_id(carpeta_medicion.name)
    if region_id in REGIONES_ANATOMICAS_PHAN:
        return REGIONES_ANATOMICAS_PHAN[region_id]

    raise ValueError(
        f"No se pudo mapear la carpeta anatomica '{carpeta_medicion.name}' a Phan-Sierra. "
        f"Usa una de estas carpetas: {list(REGIONES_ANATOMICAS_PHAN)}. "
        "Nota: 'antebrazo' solo es ambiguo; usa 'antebrazo_dorsal' o 'antebrazo_ventral'."
    )


def listar_mediciones_single_sujeto(carpeta_sujeto: Path, strict: bool = True) -> list[dict]:
    """
    Devuelve las mediciones single de un sujeto.

    Formato viejo: sujeto/series/M_R_data.csv.
    Formato nuevo: sujeto/<region_anatomica>/M_R_data.csv.
    """
    csv_legacy = carpeta_sujeto / "series" / CSV_SINGLE_NAME
    if csv_legacy.exists():
        return [
            {
                "csv_path": csv_legacy,
                "output_dir": csv_legacy.parent / "IAD_results",
                "region_id": "series",
                "region_label": "series",
                "phan_ubicacion": PHAN_UBICACION_DEFAULT,
                "formato": "legacy",
            }
        ]

    mediciones = []
    for subcarpeta in sorted(carpeta_sujeto.iterdir()):
        if not subcarpeta.is_dir():
            continue
        csv_path = subcarpeta / CSV_SINGLE_NAME
        if not csv_path.exists():
            continue

        try:
            phan_ubicacion = phan_ubicacion_para_carpeta_medicion(subcarpeta)
        except ValueError:
            if strict:
                raise
            continue
        mediciones.append(
            {
                "csv_path": csv_path,
                "output_dir": subcarpeta / "IAD_results",
                "region_id": normalizar_region_id(subcarpeta.name),
                "region_label": subcarpeta.name,
                "phan_ubicacion": phan_ubicacion,
                "formato": "anatomico",
            }
        )

    return mediciones



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
    allow_extrapolation: bool = PHAN_MUSP_SD_ALLOW_EXTRAPOLATION,
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
    """Interpola μa reportado por Phan para la ubicacion activa."""
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

    if escenario == "nominal":
        mu_input = mu_nominal
        mu_sd = float("nan")
    elif escenario == "menos_1sd":
        mu_sd = float(phan_musp_sd_mm(lambda_nm, mode="pchip"))
        mu_input = max(1e-9, mu_nominal - mu_sd)
    elif escenario == "mas_1sd":
        mu_sd = float(phan_musp_sd_mm(lambda_nm, mode="pchip"))
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
            # Se conserva como referencia literaria en la tabla exportada.
            # No se usa como input del modelo IAD.
            "mu_a_phan_mm-1": PHAN_MUA_MM,
            "mu_s_prime_phan_sierra_powerlaw_mm-1": ajuste_powerlaw,
            "mu_s_prime_phan_sierra_pchip_mm-1": ajuste_pchip,
        }
    )



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
        "phan_ubicacion": PHAN_UBICACION,
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


def _aplicar_estilo_publicacion(ax):
    """Estilo sobrio para figuras de óptica biomédica."""
    ax.set_facecolor("white")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.tick_params(axis="both", which="major", direction="out", length=4, width=0.8, color="#222222")
    ax.grid(True, color="#d7d7d7", linewidth=0.6, alpha=0.5)



def graficar_phan_sierra(ruta_png: Path, mostrar: bool = False):
    """Grafica datos Phan y las dos versiones de Phan-Sierra para la ubicacion activa."""
    lambda_plot = np.linspace(float(PHAN_LAMBDA_NM.min()), float(PHAN_LAMBDA_NM.max()), 400)
    mu_powerlaw = np.asarray(phan_sierra_mu_sp_mm(lambda_plot, mode="powerlaw"), dtype=float)
    mu_pchip = np.asarray(phan_sierra_mu_sp_mm(lambda_plot, mode="pchip"), dtype=float)

    fig, ax = plt.subplots(figsize=(8, 5))
    _aplicar_estilo_publicacion(ax)
    ax.plot(lambda_plot, mu_powerlaw, label="Phan-Sierra (powerlaw)", lw=2.0, color="#111111")
    ax.plot(lambda_plot, mu_pchip, label="Phan-Sierra (pchip)", lw=2.0, ls="--", color="#0072B2")
    ax.errorbar(
        PHAN_LAMBDA_NM,
        PHAN_MUSP_MM,
        yerr=PHAN_MUSP_SD_MM,
        fmt="o",
        ms=5,
        mfc="white",
        mec="#111111",
        mew=1.0,
        ecolor="#7f7f7f",
        elinewidth=0.9,
        capsize=2.5,
        capthick=0.9,
        label=f"{PHAN_UBICACION} de Phan et al. (media ± sd)",
    )
    ax.set_xlabel("Longitud de onda (nm)")
    ax.set_ylabel("μs' (1/mm)")
    ax.set_title("Función Phan-Sierra para μs'(λ)")
    ax.legend(frameon=False, loc="best", handlelength=2.6)
    fig.tight_layout()
    if mostrar:
        plt.show()
    if GENERAR_GRAFICA_PHAN_SIERRA_PNG:
        plt.savefig(ruta_png, dpi=300, bbox_inches="tight")
    if GENERAR_GRAFICA_PHAN_SIERRA_PDF:
        plt.savefig(ruta_png.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()


def exportar_referencia_phan_sierra(output_dir: Path) -> None:
    """Guarda la referencia Phan-Sierra activa junto a los resultados IAD."""
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"[Phan-Sierra] ubicacion = {PHAN_UBICACION}")
    print(f"[Phan-Sierra] modo = {PHAN_SIERRA_MODE}")
    print(f"[Phan-Sierra] A = {PHAN_SIERRA_A_MM:.8f} mm^-1")
    print(f"[Phan-Sierra] B = {PHAN_SIERRA_B:.8f}")
    print(f"[Phan-Sierra] lambda_ref = {PHAN_SIERRA_LAMBDA_REF_NM:.1f} nm")
    print(f"[Phan-Sierra] extrapolacion permitida = {PHAN_SIERRA_ALLOW_EXTRAPOLATION}")
    print(f"[Phan-Sierra] extrapolacion sd permitida = {PHAN_MUSP_SD_ALLOW_EXTRAPOLATION}")
    print(f"[Phan-Sierra] escenarios = {', '.join(ESCENARIOS_MUSP)}")

    tabla_phan_sierra().to_csv(output_dir / "tabla_phan_sierra.csv", index=False)
    graficar_phan_sierra(output_dir / "grafica_phan_sierra.png")



def graficar_dashboard_mu_a(df_resultados: pd.DataFrame, ruta_png: Path):
    """
    Gráfica de sensibilidad de μa recuperado.
    - Líneas separadas por ESCENARIOS_MUSP (variabilidad de μs' de Phan).
    - Banda fill_between por ESCENARIOS_MR_STD (variabilidad de medición σ_M_R).
    """
    df_plot = df_resultados.dropna(subset=["lambda_nm", "mu_a_mm-1", "escenario"]).copy()
    if "escenario_mr" not in df_plot.columns:
        df_plot["escenario_mr"] = "nominal"
    if df_plot.empty:
        return

    fig, ax = plt.subplots(figsize=(10, 6))
    _aplicar_estilo_publicacion(ax)

    estilos = {
        "nominal":   {"lw": 2.0, "ls": "-",  "marker": "o",  "ms": 3.2, "color": "#111111"},
        "menos_1sd": {"lw": 1.6, "ls": "--", "marker": None, "ms": 0,   "color": "#0072B2"},
        "mas_1sd":   {"lw": 1.6, "ls": "-.", "marker": None, "ms": 0,   "color": "#D55E00"},
    }
    for escenario in ESCENARIOS_MUSP:
        sub = df_plot[df_plot["escenario"] == escenario]
        if sub.empty:
            continue
        est = estilos.get(escenario, estilos["nominal"])

        # Curva principal: escenario_mr == "nominal"
        sub_nom = sub[sub["escenario_mr"] == "nominal"].sort_values("lambda_nm")
        if not sub_nom.empty:
            ax.plot(
                sub_nom["lambda_nm"], sub_nom["mu_a_mm-1"],
                marker=est["marker"], markersize=est["ms"],
                linewidth=est["lw"], linestyle=est["ls"],
                color=est["color"],
                label=escenario,
            )

        # Banda fill_between de variabilidad de medición (σ_M_R)
        if len(ESCENARIOS_MR_STD) > 1:
            sub_m = sub[sub["escenario_mr"] == "menos_1sd"].sort_values("lambda_nm")
            sub_p = sub[sub["escenario_mr"] == "mas_1sd"].sort_values("lambda_nm")
            if not sub_m.empty and not sub_p.empty:
                ax.fill_between(
                    sub_m["lambda_nm"].values,
                    sub_m["mu_a_mm-1"].values,
                    sub_p["mu_a_mm-1"].values,
                    alpha=0.18, color=est["color"],
                    label=f"{escenario} ±1σ_MR",
                )

    # Referencia μa de Phan: solo puntos medidos dentro del rango procesado.
    lambda_min = float(df_plot["lambda_nm"].min())
    lambda_max = float(df_plot["lambda_nm"].max())
    mascara_ref = (PHAN_LAMBDA_NM >= lambda_min) & (PHAN_LAMBDA_NM <= lambda_max)
    if np.any(mascara_ref):
        ax.plot(
            PHAN_LAMBDA_NM[mascara_ref],
            PHAN_MUA_MM[mascara_ref],
            linestyle="None",
            marker="o",
            markersize=5.5,
            markerfacecolor="white",
            markeredgewidth=1.0,
            markeredgecolor="#111111",
            label=f"μa Phan ({PHAN_UBICACION})",
        )

    ax.set_xlabel("Longitud de onda (nm)")
    ax.set_ylabel("μa recuperado (1/mm)")
    ax.set_title(f"Sensibilidad de μa — ref: {PHAN_UBICACION} (Phan et al.)")
    ax.legend(frameon=False, loc="best", handlelength=2.6)

    plt.tight_layout()
    if GENERAR_DASHBOARD_MUA_PNG:
        plt.savefig(ruta_png, dpi=200, bbox_inches="tight")
    if GENERAR_DASHBOARD_MUA_PDF:
        plt.savefig(ruta_png.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()



def generar_metricas_mu_a(df_resultados: pd.DataFrame, output_dir: Path) -> Path | None:
    """Computa métricas de error μa_IAD vs μa_Phan por escenario. Solo CSV, sin gráficas."""
    # Esta etapa se mantiene como benchmark/postanálisis y no altera la inversión IAD.
    df_resumen = df_resultados.dropna(subset=["lambda_nm", "mu_a_mm-1"]).copy()
    if df_resumen.empty:
        return None

    if "escenario" not in df_resumen.columns:
        df_resumen["escenario"] = "nominal"
    if "escenario_mr" not in df_resumen.columns:
        df_resumen["escenario_mr"] = "nominal"

    metricas = []
    for (escenario, escenario_mr), sub in df_resumen.groupby(
        ["escenario", "escenario_mr"], sort=True
    ):
        df_comp = construir_comparacion_mu_a(sub)
        fila_metricas = calcular_metricas_mu_a(df_comp)
        fila_metricas["escenario"] = escenario
        fila_metricas["escenario_mr"] = escenario_mr
        metricas.append(fila_metricas)

    if not metricas:
        return None

    ruta_metricas = output_dir / "metricas_mu_a_phan_por_escenario.csv"
    pd.DataFrame(metricas).sort_values(["escenario", "escenario_mr"]).to_csv(ruta_metricas, index=False)
    return ruta_metricas



def _correr_iad_una_lambda(
    wavelength: float,
    reflectance: float,
    header_lines: list,
    per_lambda_dir: Path,
    iad_exe_path: Path,
    escenario: str,
    med_id: int = 0,
    task_id: int = 0,
    escenario_mr: str = "nominal",
    reflectance_std: float = 0.0,
) -> dict:
    """Ejecuta IAD para una sola (λ, reflectancia) y devuelve el dict de resultado."""
    reflectance = max(reflectance, 1e-4)
    # Aplicar escenario de incertidumbre de medición
    if   escenario_mr == "mas_1sd":   reflectance = reflectance + reflectance_std
    elif escenario_mr == "menos_1sd": reflectance = max(1e-4, reflectance - reflectance_std)
    # "nominal": sin cambio

    mu_sp_mm, mu_sp_nominal, mu_sp_sd, factor_aplicado = mu_sp_escenario_phan_mm(
        wavelength,
        escenario=escenario,
        mode=PHAN_SIERRA_MODE,
    )
    g_valor = G_FIJO if USAR_G_FIJO else g_ma_et_al(wavelength)

    if task_id != 0:
        base_name = f"{escenario}_mr{escenario_mr}_task{task_id:06d}_lambda_{wavelength:.2f}".replace(".", "p")
    elif med_id != 0:
        base_name = f"{escenario}_mr{escenario_mr}_med{med_id:06d}_lambda_{wavelength:.2f}".replace(".", "p")
    else:
        base_name = f"{escenario}_mr{escenario_mr}_lambda_{wavelength:.2f}".replace(".", "p")
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
        "escenario_mr": escenario_mr,
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



def seleccionar_sujeto(temporal: bool) -> Path:
    """
    Busca carpetas de sujeto en MEDICIONES_DIR y devuelve la ruta al CSV
    correspondiente (single o temporal).

    Si AUTO_ULTIMO_SUJETO es True, selecciona el último automáticamente.
    Si no, muestra un menú interactivo.
    """
    if not MEDICIONES_DIR.exists():
        raise FileNotFoundError(f"No existe la carpeta de mediciones: {MEDICIONES_DIR}")

    carpetas = sorted([
        d for d in MEDICIONES_DIR.iterdir()
        if d.is_dir() and d.name.startswith("sujeto_")
    ])

    # Filtrar: para temporal, solo carpetas que terminen en _temporal
    # Para single, solo carpetas que NO terminen en _temporal
    if temporal:
        carpetas = [d for d in carpetas if d.name.endswith("_temporal")]
    else:
        carpetas = [d for d in carpetas if not d.name.endswith("_temporal")]

    if not carpetas:
        tipo = "temporales" if temporal else "single"
        raise FileNotFoundError(
            f"No se encontraron carpetas de sujeto ({tipo}) en {MEDICIONES_DIR}"
        )

    csv_name = CSV_TEMPORAL_NAME if temporal else CSV_SINGLE_NAME

    # Verificar cuáles tienen el CSV esperado
    carpetas_validas = []
    for c in carpetas:
        csv_path = c / "series" / csv_name
        if csv_path.exists():
            carpetas_validas.append(c)

    if not carpetas_validas:
        raise FileNotFoundError(
            f"Ninguna carpeta de sujeto contiene series/{csv_name}"
        )

    if AUTO_ULTIMO_SUJETO:
        elegida = carpetas_validas[-1]
        print(f"[Auto] Usando último sujeto: {elegida.name}")
    else:
        print(f"\n{'='*50}")
        print(f"  Sujetos disponibles ({'temporal' if temporal else 'single'})")
        print(f"{'='*50}")
        for i, c in enumerate(carpetas_validas, 1):
            print(f"  {i:3d}) {c.name}")
        print(f"{'='*50}")
        print(f"  Enter = último ({carpetas_validas[-1].name})")

        seleccion = input("  Teclea el número del elemento de la lista: ").strip()
        if seleccion == "":
            elegida = carpetas_validas[-1]
        else:
            try:
                idx = int(seleccion) - 1
                if 0 <= idx < len(carpetas_validas):
                    elegida = carpetas_validas[idx]
                else:
                    raise ValueError
            except ValueError:
                raise ValueError(
                    f"Selección inválida: '{seleccion}'. "
                    f"Usa un número entre 1 y {len(carpetas_validas)}."
                )

    csv_path = elegida / "series" / csv_name
    print(f"[Sujeto] {elegida.name} → {csv_path}")
    return csv_path


def seleccionar_sujeto_single() -> tuple[Path, list[dict]]:
    """
    Selecciona una carpeta de sujeto single y devuelve todas sus mediciones.

    Soporta el formato viejo sujeto/series/M_R_data.csv y el formato nuevo
    sujeto/<region_anatomica>/M_R_data.csv.
    """
    if not MEDICIONES_DIR.exists():
        raise FileNotFoundError(f"No existe la carpeta de mediciones: {MEDICIONES_DIR}")

    carpetas = sorted([
        d for d in MEDICIONES_DIR.iterdir()
        if d.is_dir() and d.name.startswith("sujeto_") and not d.name.endswith("_temporal")
    ])

    mediciones_por_sujeto = {}
    for carpeta in carpetas:
        mediciones = listar_mediciones_single_sujeto(carpeta, strict=False)
        if mediciones:
            mediciones_por_sujeto[carpeta] = mediciones

    carpetas_validas = list(mediciones_por_sujeto)
    if not carpetas_validas:
        raise FileNotFoundError(
            f"No se encontraron sujetos single con {CSV_SINGLE_NAME} en {MEDICIONES_DIR}"
        )

    if AUTO_ULTIMO_SUJETO:
        elegida = carpetas_validas[-1]
        print(f"[Auto] Usando ultimo sujeto single: {elegida.name}")
    else:
        print(f"\n{'='*50}")
        print("  Sujetos disponibles (single)")
        print(f"{'='*50}")
        for i, carpeta in enumerate(carpetas_validas, 1):
            n_mediciones = len(mediciones_por_sujeto[carpeta])
            print(f"  {i:3d}) {carpeta.name} ({n_mediciones} medicion(es))")
        print(f"{'='*50}")
        print(f"  Enter = ultimo ({carpetas_validas[-1].name})")

        seleccion = input("  Teclea el numero del elemento de la lista: ").strip()
        if seleccion == "":
            elegida = carpetas_validas[-1]
        else:
            try:
                idx = int(seleccion) - 1
                if 0 <= idx < len(carpetas_validas):
                    elegida = carpetas_validas[idx]
                else:
                    raise ValueError
            except ValueError:
                raise ValueError(
                    f"Seleccion invalida: '{seleccion}'. "
                    f"Usa un numero entre 1 y {len(carpetas_validas)}."
                )

    mediciones = listar_mediciones_single_sujeto(elegida, strict=True)
    print(f"[Sujeto single] {elegida.name}")
    for medicion in mediciones:
        print(
            f"  - {medicion['region_label']} -> {medicion['phan_ubicacion']} "
            f"({medicion['csv_path']})"
        )
    return elegida, mediciones



def listar_todos_sujetos() -> list[tuple[Path, bool]]:
    """
    Escanea MEDICIONES_DIR y retorna lista de (carpeta, es_temporal)
    para todas las carpetas de sujeto que tengan el CSV esperado.
    """
    if not MEDICIONES_DIR.exists():
        raise FileNotFoundError(f"No existe la carpeta de mediciones: {MEDICIONES_DIR}")

    resultados = []
    for d in sorted(MEDICIONES_DIR.iterdir()):
        if not d.is_dir() or not d.name.startswith("sujeto_"):
            continue

        es_temporal = d.name.endswith("_temporal")
        if es_temporal:
            csv_path = d / "series" / CSV_TEMPORAL_NAME
            if csv_path.exists():
                resultados.append((d, es_temporal))
        elif listar_mediciones_single_sujeto(d, strict=False):
            resultados.append((d, es_temporal))

    if not resultados:
        raise FileNotFoundError(
            f"No se encontraron carpetas de sujeto con CSVs válidos en {MEDICIONES_DIR}"
        )

    return resultados



def procesar_sujeto_single(
    csv_path: Path,
    output_dir: Path,
    header_lines: list[str],
    per_lambda_dir: Path,
    phan_ubicacion: str | None = None,
    region_id: str | None = None,
) -> None:
    """Procesa un sujeto single: IAD por escenarios → CSV + métricas + dashboard."""
    if phan_ubicacion is not None:
        configurar_phan_ubicacion(phan_ubicacion)

    output_dir.mkdir(parents=True, exist_ok=True)
    exportar_referencia_phan_sierra(output_dir)
    limpiar_carpeta_por_lambda(per_lambda_dir)

    df_mr = leer_csv_mr(csv_path)

    # Validar backward compatibility: si reflectance_std es todo ceros y se piden variaciones
    if (df_mr["reflectance_std"] == 0).all() and ESCENARIOS_MR_STD != ("nominal",):
        raise ValueError(
            "ESCENARIOS_MR_STD pide variaciones ±1σ pero reflectance_std es 0 en todo el CSV. "
            "Regenera M_R_data.csv con single_adq.py actualizado, "
            "o usa ESCENARIOS_MR_STD = ('nominal',)."
        )

    resultados = []
    tareas = [
        (idx, escenario, escenario_mr,
         float(row["wavelength_nm"]), float(row["reflectance"]), float(row["reflectance_std"]))
        for escenario    in ESCENARIOS_MUSP
        for escenario_mr in ESCENARIOS_MR_STD
        for idx, (_, row) in enumerate(df_mr.iterrows(), 1)
    ]
    total = len(tareas)
    t_inicio = time.time()

    print(
        f"  Mediciones: {len(df_mr)}  |  "
        f"Escenarios μs': {len(ESCENARIOS_MUSP)}  |  "
        f"Escenarios M_R: {len(ESCENARIOS_MR_STD)}  |  "
        f"Tareas IAD: {total}  |  Workers: {WORKERS}"
    )

    def _tarea(task_id, escenario, escenario_mr, wavelength, reflectance, reflectance_std):
        return _correr_iad_una_lambda(
            wavelength,
            reflectance,
            header_lines,
            per_lambda_dir,
            IAD_EXE_PATH,
            escenario=escenario,
            task_id=task_id,
            escenario_mr=escenario_mr,
            reflectance_std=reflectance_std,
        )

    if WORKERS <= 1:
        for i, tarea in enumerate(tareas, 1):
            resultados.append(_tarea(*tarea))
            _progreso(i, total, t_inicio)
    else:
        with ThreadPoolExecutor(max_workers=WORKERS) as executor:
            futuros = {executor.submit(_tarea, *t): t for t in tareas}
            completadas = 0
            for fut in as_completed(futuros):
                resultados.append(fut.result())
                completadas += 1
                _progreso(completadas, total, t_inicio)

    df_resultados = pd.DataFrame(resultados).sort_values(["escenario", "escenario_mr", "lambda_nm"])
    if region_id is not None:
        df_resultados.insert(0, "region_id", region_id)
        df_resultados.insert(1, "phan_ubicacion", PHAN_UBICACION)
    else:
        df_resultados.insert(0, "phan_ubicacion", PHAN_UBICACION)
    ruta_resumen = output_dir / "resumen_resultados_phan_sierra.csv"
    df_resultados.to_csv(ruta_resumen, index=False)

    ruta_metricas = generar_metricas_mu_a(df_resultados, output_dir)
    graficar_dashboard_mu_a(df_resultados, output_dir / "dashboard_mu_a.png")

    print(f"\n  Resumen: {ruta_resumen}")
    if ruta_metricas is not None:
        print(f"  Métricas: {ruta_metricas}")
    print(f"  Dashboard: {output_dir / 'dashboard_mu_a.png'}")



def procesar_sujeto_temporal(
    csv_path: Path,
    output_dir: Path,
    header_lines: list[str],
    per_lambda_dir: Path,
) -> None:
    """Procesa un sujeto temporal: IAD paralelo → CSV de resumen."""
    limpiar_carpeta_por_lambda(per_lambda_dir)

    mediciones = leer_csv_mr_temporal(csv_path)

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
        f"  Mediciones: {total_med}  |  "
        f"Tareas IAD: {total_tareas}  |  Workers: {WORKERS}"
    )

    def _tarea(escenario, med_id, tiempo, wl, ref):
        fila = _correr_iad_una_lambda(
            wl,
            ref,
            header_lines,
            per_lambda_dir,
            IAD_EXE_PATH,
            escenario=escenario,
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

    ruta_resumen = output_dir / "resumen_iad_temporal_phan_sierra.csv"
    df_resultados.to_csv(ruta_resumen, index=False)

    print(f"\n  Resumen temporal: {ruta_resumen}")



def main():
    # ------------------------------------------------------------
    # 1. Validar rutas comunes
    # ------------------------------------------------------------
    if not RXT_TEMPLATE_PATH.exists():
        raise FileNotFoundError(f"No existe la plantilla RXT: {RXT_TEMPLATE_PATH}")
    if not IAD_EXE_PATH.exists():
        raise FileNotFoundError(f"No existe iad.exe: {IAD_EXE_PATH}")

    # ------------------------------------------------------------
    # 2. Crear carpetas de trabajo (scratch space)
    # ------------------------------------------------------------
    base_dir = Path(__file__).parent
    iad_dir = base_dir / IAD_FOLDER_NAME
    iad_dir.mkdir(exist_ok=True)

    per_lambda_dir = iad_dir / IAD_PER_LAMBDA_FOLDER_NAME
    per_lambda_dir.mkdir(exist_ok=True)

    # ------------------------------------------------------------
    # 3. Leer encabezado de plantilla .rxt
    # ------------------------------------------------------------
    header_lines = extraer_encabezado_rxt(RXT_TEMPLATE_PATH)

    # ------------------------------------------------------------
    # 4. Configuracion Phan-Sierra
    # ------------------------------------------------------------
    print(f"[Phan-Sierra] ubicacion por defecto = {PHAN_UBICACION_DEFAULT}")
    print(f"[Phan-Sierra] modo = {PHAN_SIERRA_MODE}")
    print(f"[Phan-Sierra] escenarios = {', '.join(ESCENARIOS_MUSP)}")

    # ============================================================
    # MODO SINGLE — un solo sujeto, un espectro
    # ============================================================
    if MODO == "single":
        carpeta_sujeto, mediciones = seleccionar_sujeto_single()
        for i, medicion in enumerate(mediciones, 1):
            output_dir = medicion["output_dir"]
            print(f"\n[Single {i}/{len(mediciones)}] {carpeta_sujeto.name} / {medicion['region_label']}")
            print(f"  CSV: {medicion['csv_path']}")
            print(f"  Phan-Sierra: {medicion['phan_ubicacion']}")
            print(f"  Resultados: {output_dir}")
            procesar_sujeto_single(
                medicion["csv_path"],
                output_dir,
                header_lines,
                per_lambda_dir,
                phan_ubicacion=medicion["phan_ubicacion"],
                region_id=medicion["region_id"],
            )

    # ============================================================
    # MODO TEMPORAL — un solo sujeto, serie temporal
    # ============================================================
    elif MODO == "temporal":
        csv_path = seleccionar_sujeto(temporal=True)
        output_dir = csv_path.parent / "IAD_results"
        output_dir.mkdir(exist_ok=True)
        exportar_referencia_phan_sierra(iad_dir)
        print(f"\n[Temporal] Procesando → {output_dir}")
        procesar_sujeto_temporal(csv_path, output_dir, header_lines, per_lambda_dir)

    # ============================================================
    # MODO BATCH_ALL — todos los sujetos secuencialmente
    # ============================================================
    elif MODO == "batch_all":
        sujetos = listar_todos_sujetos()
        print(f"\n[Batch] {len(sujetos)} sujetos encontrados en {MEDICIONES_DIR}")

        for i, (carpeta, es_temporal) in enumerate(sujetos, 1):
            tipo = "temporal" if es_temporal else "single"
            print(f"\n{'='*60}")
            print(f"[Batch {i}/{len(sujetos)}] {carpeta.name} ({tipo})")
            print(f"{'='*60}")

            try:
                if es_temporal:
                    csv_path = carpeta / "series" / CSV_TEMPORAL_NAME
                    output_dir = carpeta / "series" / "IAD_results"
                    output_dir.mkdir(exist_ok=True)
                    exportar_referencia_phan_sierra(iad_dir)
                    procesar_sujeto_temporal(csv_path, output_dir, header_lines, per_lambda_dir)
                else:
                    mediciones = listar_mediciones_single_sujeto(carpeta)
                    for j, medicion in enumerate(mediciones, 1):
                        print(
                            f"\n  [Single {j}/{len(mediciones)}] "
                            f"{medicion['region_label']} -> {medicion['phan_ubicacion']}"
                        )
                        procesar_sujeto_single(
                            medicion["csv_path"],
                            medicion["output_dir"],
                            header_lines,
                            per_lambda_dir,
                            phan_ubicacion=medicion["phan_ubicacion"],
                            region_id=medicion["region_id"],
                        )
            except Exception as e:
                print(f"  [ERROR] {carpeta.name}: {e}")
                continue

        print(f"\n[Batch] Completado. Resultados en cada carpeta de sujeto.")

    else:
        raise ValueError(f"MODO desconocido: '{MODO}'. Usa 'single', 'temporal' o 'batch_all'.")

    print(f"\nTrabajo IAD temporal/scratch en: {iad_dir}")


if __name__ == "__main__":
    main()
