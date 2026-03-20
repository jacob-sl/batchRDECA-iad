from pathlib import Path
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
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

# g se calcula por longitud de onda usando el polinomio de Ma et al.
# ajustado por Laura Sánchez (ec. 64 de su tesis).
# Si quieres forzar un g fijo para comparación, cambia USAR_G_FIJO = True
# y ajusta el valor en G_FIJO.
USAR_G_FIJO = False
G_FIJO      = 0.8   # solo se usa si USAR_G_FIJO = True

# Modo rápido omite el "sanity check" que IAD hace con una simulación Montecarlo.
USAR_MODO_RAPIDO = True
# Mostrar o no el comando completo por cada lambda (Estético)
MOSTRAR_COMANDOS = False 
# ── TOGGLE PRINCIPAL, UN SOLO ESPECTRO O UNA TANDA CON MARCAS TEMPORALES ──────────────────────────────────────────────────────────
MODO_TEMPORAL = False   # False → M_R_data.csv (un espectro)
                        # True  → M_R_tiempo_data.csv (serie temporal)

# Solo se usa cuando MODO_TEMPORAL = True
CSV_TEMPORAL_NAME = "M_R_tiempo_data.csv"
CSV_SALIDA_TEMPORAL = "IAD_run/resumen_iad_temporal.csv"
INICIO_MEDICION = 22    # Descartar las primeras N mediciones (señal basura)
MAX_MEDICIONES  = 10    # None = todas; (CANTIDAD DE ESPECTROS A PROCESAR, después de INICIO)
VENTANA_PROMEDIO = 1    # 1 = sin promediar; N > 1 = promediar N mediciones consecutivas
WORKERS        = os.cpu_count()-1 or 4   # detecta núcleos automáticamente


# ============================================================
# FUNCIONES
# ============================================================

def leer_csv_mr(csv_path: Path) -> pd.DataFrame:
    """
    Lee el archivo M_R_data.csv y valida que tenga las columnas esperadas.
    """
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
    """
    Lee el archivo plantilla .rxt y devuelve solo el encabezado,
    es decir, todo hasta antes de la tabla espectral.
    """
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


def mu_sp_laura_jacques_mm(lambda_nm: float) -> float:
    """
    Modelo tipo Laura/Jacques para mu_s'(lambda).

    OJO:
    Aquí se usan los parámetros que hemos venido manejando:
        mu_s'(lambda) = a * (lambda / 500)^(-b)

    con:
        a = 46 cm^-1
        b = 1.421

    Como IAD trabaja en 1/mm, convertimos:
        46 cm^-1 = 4.6 mm^-1

    Nota del cambio de parámetros:
    En la tesis de Laura Sánchez, laura eligió un valor de a = 46 cm^-1, y b= 1.421, sin embargo, esos valores tomados de la revisión de jaques es para la piel, pero nosotros estamos usando un dedo, pulpa altamente vascularizada. Tenemos que buscar una nueva fuente.

    """
    A_CM = 46.0
    B_EXP = 1.421

    mu_sp_cm = A_CM * (lambda_nm / 500.0) ** (-B_EXP)
    mu_sp_mm = mu_sp_cm / 10.0
    return mu_sp_mm


def g_ma_et_al(lambda_nm: float) -> float:
    """
    Factor de anisotropía g(lambda) para piel humana.

    Polinomio de 6° orden ajustado a los datos de Ma et al. (2005),
    tal como aparece en la ec. 64 de la tesis de Laura Sánchez y en la
    revisión de Jacques (2013).

    Válido aproximadamente en 325–1557 nm.
    Produce valores en el rango 0.75–0.84 para 500–600 nm,
    coherentes con las mediciones directas en piel visible.

    Parámetros del polinomio (λ en nm):
        g = -5.603
            + 3.61e-2  · λ
            - 8.17e-5  · λ²
            + 9.51e-8  · λ³
            - 5.92e-11 · λ⁴
            + 1.83e-14 · λ⁵
            - 2.11e-18 · λ⁶
    """
    L = lambda_nm
    g = (
        -5.603
        + 3.61e-2  * L
        - 8.17e-5  * L**2
        + 9.51e-8  * L**3
        - 5.92e-11 * L**4
        + 1.83e-14 * L**5
        - 2.11e-18 * L**6
    )
    # Sanidad: g debe estar en (0, 1) para tejido biológico
    return float(max(0.5, min(0.99, g)))


def construir_rxt_una_lambda(
    wavelength: float,
    reflectance: float,
    header_lines: list[str],
    output_rxt_path: Path
) -> None:
    """
    Construye un archivo .rxt nuevo para una sola longitud de onda.
    """
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
    working_dir: Path
):
    """
    Ejecuta iad.exe desde Python.
    """
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
    """
    Intenta encontrar el archivo .txt generado por IAD.
    """
    txt_esperado = iad_dir / f"{rxt_output_path.stem}.txt"
    if txt_esperado.exists():
        return txt_esperado

    candidatos = sorted(iad_dir.glob(f"{rxt_output_path.stem}*.txt"))
    if candidatos:
        return candidatos[0]

    return None


def extraer_resultado_iad(txt_path: Path):
    """
    Extrae la fila numérica del archivo .txt de salida de IAD.
    Se asume que el archivo contiene una sola lambda.
    """
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

        # Esperamos:
        # lambda M_R_measured M_R_fit M_T_measured M_T_fit mu_a mu_s' g ...
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

    Retorna: dict ordenado { medicion_id → {'tiempo': float,
                                             'espectro': [(lambda, reflectancia), ...]} }
    El espectro de cada medicion viene ordenado por longitud de onda.
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
        espectro = sorted(zip(grupo["lambda"], grupo["reflectancia"]),
                          key=lambda x: x[0])
        mediciones[int(medicion_id)] = {"tiempo": tiempo, "espectro": list(espectro)}

    return mediciones


def _progreso(completadas: int, total: int, t_inicio: float, ancho_barra: int = 30):
    """Imprime barra de progreso con % y ETA en una sola línea (sobreescribe)."""
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
    """
    Limpia archivos viejos dentro de la carpeta por_lambda para no mezclar corridas.
    """
    if not carpeta.exists():
        return

    for archivo in carpeta.iterdir():
        if archivo.is_file():
            archivo.unlink()


def graficar_mu_a_final(df_resultados: pd.DataFrame, ruta_png: Path):
    """
    Reconstruye y guarda la gráfica final mu_a vs lambda.
    """
    df_plot = df_resultados.dropna(subset=["lambda_nm", "mu_a_mm-1"]).copy()
    df_plot = df_plot.sort_values("lambda_nm")

    plt.figure(figsize=(10, 6))
    plt.plot(df_plot["lambda_nm"], df_plot["mu_a_mm-1"], marker='o', markersize=4, linewidth=1.6)
    plt.xlabel("Longitud de onda (nm)")
    plt.ylabel("mu_a (1/mm)")
    plt.title("Coeficiente de absorción reconstruido con IAD")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(ruta_png, dpi=300, bbox_inches="tight")
    plt.close()


def _correr_iad_una_lambda(
    wavelength: float,
    reflectance: float,
    header_lines: list,
    per_lambda_dir: Path,
    iad_exe_path: Path,
    med_id: int = 0,
) -> dict:
    """Ejecuta IAD para una sola (lambda, reflectancia) y devuelve el dict de resultado."""
    reflectance = max(reflectance, 1e-4)
    mu_sp_mm = mu_sp_laura_jacques_mm(wavelength)
    g_valor = G_FIJO if USAR_G_FIJO else g_ma_et_al(wavelength)

    if med_id != 0:
        base_name = f"med{med_id:06d}_lambda_{wavelength:.2f}".replace(".", "p")
    else:
        base_name = f"lambda_{wavelength:.2f}".replace(".", "p")
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
        "lambda_nm": wavelength,
        "reflectance_input": reflectance,
        "mu_s_prime_input_mm-1": mu_sp_mm,
        "g_input": g_valor,
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
    # 3. Leer encabezado de plantilla .rxt (común a ambos modos)
    # ------------------------------------------------------------
    header_lines = extraer_encabezado_rxt(RXT_TEMPLATE_PATH)

    # ============================================================
    # MODO NORMAL — un solo espectro (M_R_data.csv)
    # ============================================================
    if not MODO_TEMPORAL:
        if not MR_CSV_PATH.exists():
            raise FileNotFoundError(f"No existe el archivo CSV: {MR_CSV_PATH}")

        df_mr = leer_csv_mr(MR_CSV_PATH)
        resultados = []
        total = len(df_mr)
        t_inicio = time.time()

        for i, (_, row) in enumerate(df_mr.iterrows(), 1):
            wavelength = float(row["wavelength_nm"])
            reflectance = float(row["reflectance"])

            fila = _correr_iad_una_lambda(
                wavelength, reflectance, header_lines, per_lambda_dir, IAD_EXE_PATH
            )
            resultados.append(fila)
            _progreso(i, total, t_inicio)

        df_resultados = pd.DataFrame(resultados).sort_values("lambda_nm")
        ruta_resumen = iad_dir / "resumen_resultados_laura_jacques.csv"
        df_resultados.to_csv(ruta_resumen, index=False)

        ruta_grafica = iad_dir / "grafica_mu_a_final.png"
        graficar_mu_a_final(df_resultados, ruta_grafica)

        print("\nEjecución terminada.")
        print(f"Resumen guardado en: {ruta_resumen}")
        print(f"Gráfica final guardada en: {ruta_grafica}")

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
            grupos = [ids_medicion[i:i+VENTANA_PROMEDIO]
                      for i in range(0, len(ids_medicion), VENTANA_PROMEDIO)]
            mediciones_prom = {}
            for grupo in grupos:
                rep_id = grupo[0]
                tiempo_prom = sum(mediciones[m]["tiempo"] for m in grupo) / len(grupo)
                wls = [wl for wl, _ in mediciones[grupo[0]]["espectro"]]
                n = len(grupo)
                refs_prom = [
                    sum(mediciones[m]["espectro"][j][1] for m in grupo) / n
                    for j in range(len(wls))
                ]
                mediciones_prom[rep_id] = {
                    "tiempo": tiempo_prom,
                    "espectro": list(zip(wls, refs_prom)),
                }
            ids_medicion = [g[0] for g in grupos]
            mediciones = mediciones_prom

        total_med = len(ids_medicion)

        # Construir lista plana de tareas
        tareas = [
            (med_id, mediciones[med_id]["tiempo"], float(wl), float(ref))
            for med_id in ids_medicion
            for wl, ref in mediciones[med_id]["espectro"]
        ]
        total_tareas = len(tareas)

        print(f"Mediciones a procesar: {total_med}  |  Tareas IAD totales: {total_tareas}  |  Workers: {WORKERS}")

        def _tarea(med_id, tiempo, wl, ref):
            fila = _correr_iad_una_lambda(
                wl, ref, header_lines, per_lambda_dir, IAD_EXE_PATH, med_id
            )
            fila["medicion"] = med_id
            fila["tiempo"] = tiempo
            return fila

        resultados = []
        t_inicio = time.time()

        if WORKERS <= 1:
            # ── Secuencial ──
            for i, (med_id, tiempo, wl, ref) in enumerate(tareas, 1):
                resultados.append(_tarea(med_id, tiempo, wl, ref))
                _progreso(i, total_tareas, t_inicio)
        else:
            # ── Paralelo ──
            with ThreadPoolExecutor(max_workers=WORKERS) as executor:
                futuros = {
                    executor.submit(_tarea, *t): t
                    for t in tareas
                }
                completadas = 0
                for fut in as_completed(futuros):
                    resultados.append(fut.result())
                    completadas += 1
                    _progreso(completadas, total_tareas, t_inicio)

        df_resultados = pd.DataFrame(resultados)
        # Reordenar columnas para que medicion y tiempo queden primero
        cols_primeras = ["medicion", "tiempo", "lambda_nm", "reflectance_input",
                         "mu_a_mm-1", "mu_s_prime_mm-1", "g"]
        cols_resto = [c for c in df_resultados.columns if c not in cols_primeras]
        df_resultados = df_resultados[cols_primeras + cols_resto]
        df_resultados = df_resultados.sort_values(["medicion", "lambda_nm"])

        ruta_resumen = base_dir / CSV_SALIDA_TEMPORAL
        ruta_resumen.parent.mkdir(exist_ok=True)
        df_resultados.to_csv(ruta_resumen, index=False)

        print(f"\nEjecución temporal terminada.")
        print(f"Mediciones procesadas: {total_med}")
        print(f"Resumen temporal guardado en: {ruta_resumen}")

    print(f"Carpeta principal IAD: {iad_dir}")
    print(f"Archivos por lambda: {per_lambda_dir}")


if __name__ == "__main__":
    main()