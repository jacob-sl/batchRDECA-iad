# Importación de librerías necesarias
import atexit
import json
import os
import re
import sys
import clr
from datetime import datetime
from System import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
from scipy.signal import butter, filtfilt

# MUST HAVE pythonNET installed!
sys.path.append(r"C:\Program Files (x86)\Microsoft.NET\Primary Interop Assemblies\\")
clr.AddReference("Thorlabs.ccs.interop64")

# Import the NET reference
import Thorlabs.ccs.interop64

# Conexión con el espectrofotómetro Thorlabs CCS200
spec = Thorlabs.ccs.interop64.TLCCS(
    "USB0::0x1313::0x8089::M00496376::RAW",
    Boolean(True),
    Boolean(True)
)
print("Espectrofotómetro conectado exitosamente")


def _cleanup_spec():
    try:
        spec.Dispose()
        print("Espectrofotómetro desconectado (cleanup automático)")
    except Exception:
        pass


atexit.register(_cleanup_spec)

# ============================================================
# ========== PARÁMETROS AJUSTABLES ============================
# ============================================================

# Tiempo de integración inicial en segundos
TIEMPO_INTEGRACION_INICIAL = 0.1
TIEMPO_INTEGRACION = TIEMPO_INTEGRACION_INICIAL

# Número de mediciones a promediar para R_0 y R_1
NUM_MEDICIONES_PROMEDIO = 10

# Tiempo de espera entre mediciones individuales para R_0 y R_1 (segundos)
TIEMPO_ESPERA = 0.01

# Número de espectros temporales a adquirir para R_M
NUM_ESPECTROS_TEMPORAL = 500

# Reflectancia del estándar
R_STD = 0.99

# ========== PARÁMETROS DE OPTIMIZACIÓN DEL TIEMPO DE INTEGRACIÓN ==========
# Para CCS200 en este flujo: el máximo es 1.0 (valores normalizados)
PORCENTAJE_SATURACION = 0.95
VALOR_MAXIMO_DETECTOR = 1.0
UMBRAL_SATURACION = VALOR_MAXIMO_DETECTOR * PORCENTAJE_SATURACION  # 0.95

# ========== CONFIGURACIÓN DEL FILTRO BUTTERWORTH ==========
FRECUENCIA_CORTE_BUTTER = 0.1  # Frecuencia de corte normalizada (0-1, Nyquist=1)
ORDEN_FILTRO_BUTTER = 6

# ========== RANGO DE LONGITUDES DE ONDA (TRUNCADO) ==========
LAMBDA_MIN = 500  # Longitud de onda mínima en nm
LAMBDA_MAX = 590  # Longitud de onda máxima en nm

# Número de muestras deseadas después del diezmado
MUESTRAS_OBJETIVO = 200 

# ============================================================

print("Configuración:")
print(f"  - Tiempo de integración inicial: {TIEMPO_INTEGRACION_INICIAL*1000:.2f} ms")
print(f"  - Mediciones por promedio (R_0 y R_1): {NUM_MEDICIONES_PROMEDIO}")
print(f"  - Tiempo entre mediciones (R_0 y R_1): {TIEMPO_ESPERA} s")
print(f"  - Espectros temporales a adquirir (R_M): {NUM_ESPECTROS_TEMPORAL}")
print(f"  - Butterworth: orden {ORDEN_FILTRO_BUTTER}, fc={FRECUENCIA_CORTE_BUTTER}")
print(f"  - Rango truncado: {LAMBDA_MIN} - {LAMBDA_MAX} nm")
print(f"  - Muestras objetivo tras diezmado: {MUESTRAS_OBJETIVO}")
print("Configuración de saturación:")
print(f"  - Valor máximo del detector: {VALOR_MAXIMO_DETECTOR}")
print(f"  - Porcentaje objetivo: {PORCENTAJE_SATURACION*100:.0f}%")
print(f"  - Umbral de saturación: {UMBRAL_SATURACION:.3f}")
print(
    f"  - Duración mínima teórica de la serie temporal: "
    f"{NUM_ESPECTROS_TEMPORAL * TIEMPO_INTEGRACION_INICIAL:.2f} s "
    f"(sin contar overhead del sistema)"
)


# ============================================================
# ========== FUNCIONES =======================================
# ============================================================

def tomar_medicion(spec, tiempo_integracion):
    """
    Toma una sola medición del espectrofotómetro con el tiempo de integración dado.

    Parámetros:
    -----------
    spec : objeto TLCCS
        Objeto del espectrofotómetro Thorlabs
    tiempo_integracion : float
        Tiempo de integración en segundos

    Retorna:
    --------
    intensities : numpy.ndarray
        Array con las intensidades medidas (3648 puntos)
    """
    spec.setIntegrationTime(Double(tiempo_integracion))
    spec.startScan()

    scan = Array.CreateInstance(Double, 3648)
    spec.getScanData(scan)

    intensities = np.fromiter(scan, dtype=np.float64, count=3648)
    return intensities


def tomar_serie_adquisiciones(spec, tiempo_integracion, num_mediciones, tiempo_espera=0.1):
    """
    Toma múltiples mediciones y las promedia.
    Se usa para R_0 y R_1.

    Parámetros:
    -----------
    spec : objeto TLCCS
        Objeto del espectrofotómetro
    tiempo_integracion : float
        Tiempo de integración en segundos
    num_mediciones : int
        Número de mediciones a promediar
    tiempo_espera : float
        Tiempo de espera entre mediciones (segundos)

    Retorna:
    --------
    intensidad_promedio : numpy.ndarray
        Array con las intensidades promediadas
    """
    print(f"Tomando {num_mediciones} mediciones para promediar...")

    mediciones = []

    for i in range(num_mediciones):
        intensity = tomar_medicion(spec, tiempo_integracion)
        mediciones.append(intensity)

        print(f"  Medición {i+1}/{num_mediciones} completada")

        if i < num_mediciones - 1:
            time.sleep(tiempo_espera)

    intensidad_promedio = np.mean(mediciones, axis=0)

    print("✓ Serie completada. Promedio calculado.\n")
    return intensidad_promedio


def optimizar_tiempo_integracion(
    spec,
    tiempo_inicial,
    umbral_saturacion_absoluto=UMBRAL_SATURACION,
    factor_incremento=1.15,
    max_iteraciones=20
):
    """
    Optimiza el tiempo de integración aumentándolo gradualmente
    hasta encontrar el máximo tiempo usable sin saturación.

    Parámetros:
    -----------
    spec : objeto TLCCS
        Objeto del espectrofotómetro
    tiempo_inicial : float
        Tiempo de integración inicial en segundos
    umbral_saturacion_absoluto : float
        Umbral absoluto de saturación
    factor_incremento : float
        Factor por el cual se incrementa el tiempo en cada iteración
    max_iteraciones : int
        Número máximo de intentos

    Retorna:
    --------
    tiempo_optimo : float
        Tiempo de integración óptimo encontrado
    intensidad_optima : numpy.ndarray
        Intensidad medida con el tiempo óptimo
    """
    print("=" * 60)
    print("OPTIMIZACIÓN DEL TIEMPO DE INTEGRACIÓN")
    print("=" * 60)
    print(f"Tiempo inicial: {tiempo_inicial*1000:.2f} ms")
    print(f"Umbral de saturación absoluto: {umbral_saturacion_absoluto:.3f}")
    print(f"Factor de incremento: {factor_incremento} (aumenta {(factor_incremento-1)*100:.0f}% por iteración)")
    print("Estrategia: aumentar el tiempo hasta hallar el máximo usable sin saturar")
    print()

    tiempo_actual = tiempo_inicial
    tiempo_optimo = None
    intensidad_optima = None
    encontro_valor_valido = False

    for iteracion in range(max_iteraciones):
        print(f"Iteración {iteracion + 1}/{max_iteraciones}:")
        print(f"  Probando con tiempo de integración: {tiempo_actual*1000:.2f} ms")

        intensidad = tomar_medicion(spec, tiempo_actual)
        max_intensidad = np.max(intensidad)

        print(f"  Intensidad máxima: {max_intensidad:.6f}")

        if max_intensidad < umbral_saturacion_absoluto:
            margen = umbral_saturacion_absoluto - max_intensidad
            porcentaje_usado = (max_intensidad / umbral_saturacion_absoluto) * 100

            print("  ✓ Tiempo válido. Guardando como candidato óptimo.")
            print(f"  ✓ Usando {porcentaje_usado:.1f}% del umbral (margen: {margen:.6f})")

            tiempo_optimo = tiempo_actual
            intensidad_optima = intensidad
            encontro_valor_valido = True

            tiempo_actual *= factor_incremento
            print(f"  → Intentando aumentar a {tiempo_actual*1000:.2f} ms\n")
        else:
            exceso = max_intensidad - umbral_saturacion_absoluto
            print(f"  ⚠ ¡SATURACIÓN DETECTADA! ({max_intensidad:.6f} > {umbral_saturacion_absoluto:.6f})")
            print(f"  ⚠ Exceso: {exceso:.6f}")

            if encontro_valor_valido:
                print("  → El tiempo anterior era el óptimo.\n")
            else:
                tiempo_actual /= factor_incremento
                if tiempo_actual < 0.001:
                    print("  ⚠ Tiempo mínimo alcanzado (1 ms).")
                    tiempo_optimo = 0.001
                    intensidad_optima = tomar_medicion(spec, 0.001)
                    break
                print(f"  ⚠ Saturación detectada. Reduciendo a {tiempo_actual*1000:.2f} ms\n")
                continue
            break
    else:
        print("⚠ Alcanzado el número máximo de iteraciones.")
        if encontro_valor_valido:
            print(f"  Usando último tiempo válido: {tiempo_optimo*1000:.2f} ms")
        else:
            print(f"  No se encontró valor óptimo. Usando tiempo inicial: {tiempo_inicial*1000:.2f} ms")
            tiempo_optimo = tiempo_inicial
            intensidad_optima = tomar_medicion(spec, tiempo_inicial)

    print("\n" + "=" * 60)
    if encontro_valor_valido:
        print(f"TIEMPO DE INTEGRACIÓN ÓPTIMO: {tiempo_optimo*1000:.2f} ms")
        print(f"Intensidad máxima alcanzada: {np.max(intensidad_optima):.6f}")
        print(f"Porcentaje del umbral: {(np.max(intensidad_optima)/umbral_saturacion_absoluto)*100:.1f}%")
        print("Este es el máximo tiempo sin saturar el detector")
    else:
        print(f"⚠ NO SE PUDO OPTIMIZAR - Usando tiempo inicial: {tiempo_optimo*1000:.2f} ms")
        print(f"Intensidad con tiempo inicial: {np.max(intensidad_optima):.6f}")
        print("⚠ ADVERTENCIA: El tiempo inicial ya causa saturación. Reduce TIEMPO_INTEGRACION_INICIAL")
    print("=" * 60)

    return tiempo_optimo, intensidad_optima


def truncar_a_rango(intensidad_data, wavelengths_arr=None):
    """
    Trunca datos de intensidad al rango LAMBDA_MIN-LAMBDA_MAX.

    Parámetros:
    -----------
    intensidad_data : numpy.ndarray
        Array con intensidades medidas (3648 puntos)
    wavelengths_arr : numpy.ndarray, opcional
        Array de longitudes de onda. Si no se proporciona, usa el global `wavelengths`.

    Retorna:
    --------
    wavelengths_trunc : numpy.ndarray
        Longitudes de onda truncadas
    intensidad_trunc : numpy.ndarray
        Intensidades truncadas
    """
    if wavelengths_arr is None:
        wavelengths_arr = wavelengths
    mask_rango = (wavelengths_arr >= LAMBDA_MIN) & (wavelengths_arr <= LAMBDA_MAX)
    return wavelengths_arr[mask_rango], intensidad_data[mask_rango]


# Pre-cálculo de coeficientes Butterworth (constantes para todo el script)
_BUTTER_B, _BUTTER_A = butter(ORDEN_FILTRO_BUTTER, FRECUENCIA_CORTE_BUTTER, btype='low')


def aplicar_butterworth(datos, fc=None, orden=None):
    """
    Aplica filtro Butterworth pasabajas.
    Si existen NaN, los interpola antes del filtrado y los restaura después.

    Parámetros:
    -----------
    datos : numpy.ndarray
        Datos a filtrar
    fc : float
        Frecuencia de corte normalizada
    orden : int
        Orden del filtro

    Retorna:
    --------
    datos_filtrados : numpy.ndarray
        Datos filtrados
    """
    if fc is None:
        fc = FRECUENCIA_CORTE_BUTTER
    if orden is None:
        orden = ORDEN_FILTRO_BUTTER

    mask_nan = np.isnan(datos)
    if mask_nan.all():
        return datos.copy()

    datos_interp = datos.copy()
    if mask_nan.any():
        indices = np.arange(len(datos))
        datos_interp[mask_nan] = np.interp(
            indices[mask_nan],
            indices[~mask_nan],
            datos[~mask_nan]
        )

    # Reutilizar coeficientes pre-calculados si son los parámetros por defecto
    if fc == FRECUENCIA_CORTE_BUTTER and orden == ORDEN_FILTRO_BUTTER:
        b, a = _BUTTER_B, _BUTTER_A
    else:
        b, a = butter(orden, fc, btype='low')

    datos_filtrados = filtfilt(b, a, datos_interp)
    datos_filtrados[mask_nan] = np.nan

    return datos_filtrados


def adquirir_serie_temporal_rm(spec, tiempo_integracion, num_espectros):
    """
    Adquiere una serie temporal de espectros R_M sin promedio.

    Cada espectro recibe:
    - número de medición
    - tiempo relativo en segundos desde el inicio de la serie temporal

    Parámetros:
    -----------
    spec : objeto TLCCS
        Objeto del espectrofotómetro
    tiempo_integracion : float
        Tiempo de integración en segundos
    num_espectros : int
        Número de espectros a adquirir

    Retorna:
    --------
    tiempos_relativos : numpy.ndarray
        Tiempos relativos de cada espectro
    espectros_raw : list[numpy.ndarray]
        Lista con todos los espectros crudos adquiridos
    tiempo_total : float
        Tiempo total transcurrido de la serie temporal
    """
    print("=" * 50)
    print("SERIE TEMPORAL R_M")
    print("=" * 50)
    print(f"Se adquirirán {num_espectros} espectros temporales.")
    print(f"Tiempo de integración usado: {tiempo_integracion*1000:.2f} ms")
    print()

    tiempos_relativos = []
    espectros_raw = []

    # Configurar el tiempo de integración una sola vez antes del loop para
    # evitar que el CCS200 realice un ciclo interno en cada iteración,
    # lo que duplicaba el tiempo aparente por medición.
    spec.setIntegrationTime(Double(tiempo_integracion))

    inicio_adquisicion = time.perf_counter()

    for medicion in range(1, num_espectros + 1):
        # Registrar el timestamp al INICIO del scan (antes de la exposición)
        tiempo_relativo = time.perf_counter() - inicio_adquisicion
        spec.startScan()
        scan = Array.CreateInstance(Double, 3648)
        spec.getScanData(scan)
        intensidad_raw = np.fromiter(scan, dtype=np.float64, count=3648)

        tiempos_relativos.append(tiempo_relativo)
        espectros_raw.append(intensidad_raw)

        print(f"  Medición {medicion}/{num_espectros} completada")

    tiempo_total = time.perf_counter() - inicio_adquisicion

    print()
    print("✓ Serie temporal completada.")
    print(f"✓ Espectros adquiridos: {len(espectros_raw)}")
    print(f"✓ Tiempo total transcurrido: {tiempo_total:.3f} s\n")

    return np.array(tiempos_relativos, dtype=np.float64), espectros_raw, tiempo_total


def crear_carpeta_sujeto(base_dir):
    """
    Crea una carpeta con formato sujeto_XXX_dd_mm_YYYY y subcarpeta /series.
    """
    os.makedirs(base_dir, exist_ok=True)

    patron = re.compile(r"^sujeto_(\d{3})_\d{2}_\d{2}_\d{4}_temporal$")
    ids_existentes = []

    for nombre in os.listdir(base_dir):
        ruta = os.path.join(base_dir, nombre)
        if not os.path.isdir(ruta):
            continue
        match = patron.match(nombre)
        if match:
            ids_existentes.append(int(match.group(1)))

    siguiente_id = max(ids_existentes) + 1 if ids_existentes else 1
    fecha = datetime.now().strftime("%d_%m_%Y")

    while True:
        nombre_sujeto = f"sujeto_{siguiente_id:03d}_{fecha}_temporal"
        ruta_sujeto = os.path.join(base_dir, nombre_sujeto)
        if not os.path.exists(ruta_sujeto):
            break
        siguiente_id += 1

    ruta_series = os.path.join(ruta_sujeto, "series")
    os.makedirs(ruta_series, exist_ok=False)

    return ruta_sujeto, ruta_series, siguiente_id


def solicitar_dato_sujeto(prompt, default="No especificado"):
    valor = input(prompt).strip()
    return valor if valor else default


def confirmar_por_tecla(mensaje, tecla="o"):
    """
    Confirmación por tecla para acciones opcionales.
    """
    entrada = input(f"{mensaje} (presiona '{tecla}' + Enter para continuar): ").strip().lower()
    return entrada == tecla.lower()


def confirmar_simple(mensaje, default=True):
    """
    Confirmación simple tipo sí/no.
    """
    if default:
        entrada = input(f"{mensaje} [S/n]: ").strip().lower()
        return entrada in ("", "s", "si", "sí", "y", "yes")

    entrada = input(f"{mensaje} [s/N]: ").strip().lower()
    return entrada in ("s", "si", "sí", "y", "yes")


def guardar_y_mostrar_figura(fig, ruta_series, nombre_archivo, dpi=300):
    ruta_png = os.path.join(ruta_series, f"{nombre_archivo}.png")
    fig.savefig(ruta_png, dpi=dpi, bbox_inches="tight")
    print(f"Gráfica guardada en: {ruta_png}")
    plt.show()
    plt.close(fig)
    return ruta_png


def guardar_datos_sujeto(ruta_sujeto, sujeto_id, nombre, edad, municipio_nacimiento, fitzpatrick, afecciones):
    ruta_txt = os.path.join(ruta_sujeto, "datos_sujeto.txt")
    contenido = [
        f"sujeto_id: {sujeto_id:03d}",
        f"fecha_registro: {datetime.now().strftime('%d/%m/%Y %H:%M:%S')}",
        f"nombre: {nombre}",
        f"edad: {edad}",
        f"municipio_nacimiento: {municipio_nacimiento}",
        f"fitzpatrick: {fitzpatrick}",
        f"afecciones_conocidas: {afecciones}"
    ]

    with open(ruta_txt, "w", encoding="utf-8") as f:
        f.write("\n".join(contenido))

    return ruta_txt


def rutas_calibracion(directorio_calibraciones):
    ruta_npz = os.path.join(directorio_calibraciones, "calibracion_actual.npz")
    ruta_meta = os.path.join(directorio_calibraciones, "calibracion_actual.json")
    return ruta_npz, ruta_meta


def calibracion_es_compatible(wavelengths_guardadas, wavelengths_actuales, r0, r1):
    if len(wavelengths_guardadas) != len(wavelengths_actuales):
        return False
    if len(r0) != len(wavelengths_actuales) or len(r1) != len(wavelengths_actuales):
        return False
    return np.allclose(wavelengths_guardadas, wavelengths_actuales, atol=1e-6)


def cargar_calibracion(ruta_npz, ruta_meta):
    if not (os.path.exists(ruta_npz) and os.path.exists(ruta_meta)):
        return None

    try:
        with open(ruta_meta, "r", encoding="utf-8") as f:
            metadata = json.load(f)

        with np.load(ruta_npz) as data:
            wavelengths_guardadas = data["wavelengths"]
            r0 = data["R_0"]
            r1 = data["R_1"]
    except Exception as e:
        print(f"⚠ No se pudo cargar la calibración guardada: {e}")
        return None

    return {
        "metadata": metadata,
        "wavelengths": wavelengths_guardadas,
        "R_0": r0,
        "R_1": r1
    }


def guardar_calibracion(ruta_npz, ruta_meta, wavelengths_data, r0, r1, tiempo_integracion, num_mediciones):
    os.makedirs(os.path.dirname(ruta_npz), exist_ok=True)

    np.savez_compressed(
        ruta_npz,
        wavelengths=wavelengths_data,
        R_0=r0,
        R_1=r1
    )

    metadata = {
        "fecha_calibracion": datetime.now().isoformat(timespec="seconds"),
        "tiempo_integracion_s": float(tiempo_integracion),
        "num_mediciones_promedio": int(num_mediciones),
        "lambda_min_nm": float(LAMBDA_MIN),
        "lambda_max_nm": float(LAMBDA_MAX),
        "nota_uso": "R_0 y R_1 se reutilizan entre sujetos. Recalibrar solo si cambia el sistema o el montaje."
    }

    with open(ruta_meta, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)

    return metadata


def preparar_referencias_procesadas(r0, r1):
    """
    Prepara R_0 y R_1 para:
    - exportación en CSV
    - cálculo posterior de M_R temporal

    Retorna:
    --------
    dict con:
        wavelengths_trunc
        wavelengths_diezmado
        r0_butter
        r1_butter
        r0_diezmado
        r1_diezmado
        factor_diezmado
        denominador
    """
    wavelengths_trunc, r0_trunc = truncar_a_rango(r0)
    _, r1_trunc = truncar_a_rango(r1)

    r0_butter = aplicar_butterworth(r0_trunc)
    r1_butter = aplicar_butterworth(r1_trunc)

    factor_diezmado = max(1, len(wavelengths_trunc) // MUESTRAS_OBJETIVO)

    wavelengths_diezmado = wavelengths_trunc[::factor_diezmado]
    r0_diezmado = r0_butter[::factor_diezmado]
    r1_diezmado = r1_butter[::factor_diezmado]

    denominador = r1_butter - r0_butter
    denominador = np.where(denominador == 0, np.nan, denominador)

    return {
        "wavelengths_trunc": wavelengths_trunc,
        "wavelengths_diezmado": wavelengths_diezmado,
        "r0_butter": r0_butter,
        "r1_butter": r1_butter,
        "r0_diezmado": r0_diezmado,
        "r1_diezmado": r1_diezmado,
        "factor_diezmado": factor_diezmado,
        "denominador": denominador
    }


def guardar_rm_trueraw_tiempo_csv(ruta_csv, tiempos_relativos, espectros_raw, wavelengths_completas):
    """
    Guarda la serie temporal cruda real sin truncado, sin filtro y sin diezmado.

    Formato:
    medicion,tiempo,lambda,intensidad
    """
    rows = []
    for medicion, (tiempo_relativo, espectro_raw) in enumerate(
        zip(tiempos_relativos, espectros_raw), start=1
    ):
        n = len(wavelengths_completas)
        med_col = np.full(n, medicion)
        t_col = np.full(n, tiempo_relativo)
        rows.append(np.column_stack([med_col, t_col, wavelengths_completas, espectro_raw]))

    data = np.vstack(rows)
    df = pd.DataFrame(data, columns=["medicion", "tiempo", "lambda", "intensidad"])
    df["medicion"] = df["medicion"].astype(int)
    df.to_csv(ruta_csv, index=False)

    return ruta_csv


def guardar_mr_tiempo_csv(ruta_csv, tiempos_relativos, espectros_raw, refs_proc):
    """
    Guarda la serie temporal de M_R ya:
    - truncada
    - filtrada
    - diezmada

    Formato:
    medicion,tiempo,lambda,reflectancia
    """
    wavelengths_diezmado = refs_proc["wavelengths_diezmado"]
    r0_butter = refs_proc["r0_butter"]
    denominador = refs_proc["denominador"]
    factor_diezmado = refs_proc["factor_diezmado"]

    rows = []
    for medicion, (tiempo_relativo, espectro_raw) in enumerate(
        zip(tiempos_relativos, espectros_raw), start=1
    ):
        _, rm_trunc = truncar_a_rango(espectro_raw)
        rm_butter = aplicar_butterworth(rm_trunc)

        numerador = rm_butter - r0_butter
        mr_values = R_STD * (numerador / denominador)
        mr_diezmado = mr_values[::factor_diezmado]

        n = len(wavelengths_diezmado)
        med_col = np.full(n, medicion)
        t_col = np.full(n, tiempo_relativo)
        rows.append(np.column_stack([med_col, t_col, wavelengths_diezmado, mr_diezmado]))

    data = np.vstack(rows)
    df = pd.DataFrame(data, columns=["medicion", "tiempo", "lambda", "reflectancia"])
    df["medicion"] = df["medicion"].astype(int)
    df.to_csv(ruta_csv, index=False)

    return ruta_csv


print("Funciones definidas correctamente")


# ============================================================
# ========== OBTENCIÓN DE LONGITUDES DE ONDA =================
# ============================================================

ref16 = Int16(0)
nullable_double = Nullable[Double](0)
wvdata = Array.CreateInstance(Double, 3648)
spec.getWavelengthData(ref16, wvdata, nullable_double, nullable_double)
wavelengths = np.asarray(list(wvdata))

print(f"Longitudes de onda obtenidas: {len(wavelengths)} puntos")
print(f"Rango: {wavelengths[0]:.2f} nm - {wavelengths[-1]:.2f} nm")


# ============================================================
# ========== REGISTRO DE SUJETO ==============================
# ============================================================

directorio_raiz = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Mediciones")
directorio_calibraciones = os.path.join(directorio_raiz, "calibraciones")
ruta_calibracion_npz, ruta_calibracion_meta = rutas_calibracion(directorio_calibraciones)
ruta_sujeto, ruta_series, sujeto_id = crear_carpeta_sujeto(directorio_raiz)

print("\nRegistro de nuevo sujeto")
nombre_sujeto = solicitar_dato_sujeto("Nombre del sujeto: ")
edad_sujeto = solicitar_dato_sujeto("Edad del sujeto (años): ")
municipio_nacimiento_sujeto = solicitar_dato_sujeto("Municipio de nacimiento: ")
fitzpatrick_sujeto = solicitar_dato_sujeto("Fototipo Fitzpatrick (I-VI): ")
afecciones_sujeto = solicitar_dato_sujeto("Afecciones conocidas: ")

ruta_datos_sujeto = guardar_datos_sujeto(
    ruta_sujeto=ruta_sujeto,
    sujeto_id=sujeto_id,
    nombre=nombre_sujeto,
    edad=edad_sujeto,
    municipio_nacimiento=municipio_nacimiento_sujeto,
    fitzpatrick=fitzpatrick_sujeto,
    afecciones=afecciones_sujeto
)

print(f"Sujeto registrado en: {ruta_sujeto}")
print(f"  - TXT de sujeto: {ruta_datos_sujeto}")
print(f"  - Series: {ruta_series}")

rutas_graficas = []


# ============================================================
# ========== GESTIÓN DE CALIBRACIÓN ==========================
# ============================================================

usar_calibracion_guardada = False
calibracion = cargar_calibracion(ruta_calibracion_npz, ruta_calibracion_meta)

if calibracion is not None:
    metadata_cal = calibracion["metadata"]
    compatible = calibracion_es_compatible(
        calibracion["wavelengths"],
        wavelengths,
        calibracion["R_0"],
        calibracion["R_1"]
    )

    if compatible:
        fecha_cal = metadata_cal.get("fecha_calibracion", "")
        tiempo_cal = float(metadata_cal.get("tiempo_integracion_s", TIEMPO_INTEGRACION_INICIAL))

        print("\nCalibración encontrada en disco:")
        print(f"  - Fecha: {fecha_cal if fecha_cal else 'desconocida'}")
        print(f"  - Tiempo integración guardado: {tiempo_cal*1000:.2f} ms")
        print("  - Uso recomendado: reutilizar R_0/R_1 entre sujetos y recalibrar solo ante cambios de sistema o montaje.")

        forzar_recalibracion = confirmar_simple(
            "¿Hubo cambios en sistema o montaje que requieran nueva calibración (R_0/R_1)?",
            default=False
        )
        usar_calibracion_guardada = not forzar_recalibracion

        if usar_calibracion_guardada:
            print("  - Se reutilizará la calibración guardada.")
        else:
            print("  - Se medirá una nueva calibración (R_0 y R_1).")
    else:
        print("\n⚠ La calibración guardada no es compatible con la configuración actual.")
        print("   Se medirá una calibración nueva (R_0 y R_1).")
else:
    print("\nNo se encontró calibración guardada. Se medirá una calibración nueva (R_0 y R_1).")

if usar_calibracion_guardada:
    metadata_cal = calibracion["metadata"]
    R_0 = calibracion["R_0"]
    R_1 = calibracion["R_1"]
    TIEMPO_INTEGRACION = float(metadata_cal.get("tiempo_integracion_s", TIEMPO_INTEGRACION_INICIAL))
    optimizacion_realizada = False

    print("\n✓ Calibración cargada desde disco.")
    print(f"  - Tiempo de integración usado: {TIEMPO_INTEGRACION*1000:.2f} ms")
    print(f"  - R_0 cargado. Shape: {R_0.shape}")
    print(f"  - R_1 cargado. Shape: {R_1.shape}")
else:
    if confirmar_por_tecla("¿Deseas ejecutar la optimización de tiempo de integración?", tecla="o"):
        TIEMPO_INTEGRACION, intensidad_opt = optimizar_tiempo_integracion(
            spec,
            TIEMPO_INTEGRACION_INICIAL,
            umbral_saturacion_absoluto=UMBRAL_SATURACION,
            factor_incremento=1.10,
            max_iteraciones=30
        )
        optimizacion_realizada = True
    else:
        TIEMPO_INTEGRACION = TIEMPO_INTEGRACION_INICIAL
        intensidad_opt = tomar_medicion(spec, TIEMPO_INTEGRACION)
        optimizacion_realizada = False
        print(f"\nOptimización omitida. Se usará el tiempo inicial: {TIEMPO_INTEGRACION*1000:.2f} ms")

    # TOMA DE R_0
    if not confirmar_simple("¿Iniciar serie R_0?"):
        raise SystemExit("Proceso cancelado por el usuario antes de R_0.")

    print("=" * 50)
    print("SERIE R_0")
    print("=" * 50)

    R_0 = tomar_serie_adquisiciones(
        spec,
        TIEMPO_INTEGRACION,
        NUM_MEDICIONES_PROMEDIO,
        TIEMPO_ESPERA
    )

    print(f"R_0 guardado. Shape: {R_0.shape}")
    print(f"Rango de intensidades: {R_0.min():.6f} - {R_0.max():.6f}")
    print(f"Tiempo de integración usado para R_0: {TIEMPO_INTEGRACION*1000:.2f} ms")

    # TOMA DE R_1
    if not confirmar_simple("¿Iniciar serie R_1?"):
        raise SystemExit("Proceso cancelado por el usuario antes de R_1.")

    print("=" * 50)
    print("SERIE R_1")
    print("=" * 50)

    R_1 = tomar_serie_adquisiciones(
        spec,
        TIEMPO_INTEGRACION,
        NUM_MEDICIONES_PROMEDIO,
        TIEMPO_ESPERA
    )

    print(f"R_1 guardado. Shape: {R_1.shape}")
    print(f"Rango de intensidades: {R_1.min():.6f} - {R_1.max():.6f}")

    meta_guardada = guardar_calibracion(
        ruta_calibracion_npz,
        ruta_calibracion_meta,
        wavelengths,
        R_0,
        R_1,
        TIEMPO_INTEGRACION,
        NUM_MEDICIONES_PROMEDIO
    )

    print("\n✓ Calibración guardada en disco.")
    print(f"  - Archivo datos: {ruta_calibracion_npz}")
    print(f"  - Archivo metadata: {ruta_calibracion_meta}")
    print(f"  - Fecha calibración: {meta_guardada['fecha_calibracion']}")


# ============================================================
# ========== GRÁFICA SOLO SI HUBO OPTIMIZACIÓN ===============
# ============================================================

if not usar_calibracion_guardada and optimizacion_realizada:
    mask_rango_opt = (wavelengths >= LAMBDA_MIN) & (wavelengths <= LAMBDA_MAX)
    wl_opt = wavelengths[mask_rango_opt]
    int_opt = intensidad_opt[mask_rango_opt]

    print(f"\nEstadísticas de la medición de referencia ({LAMBDA_MIN}-{LAMBDA_MAX} nm):")
    print(f"  - Intensidad máxima: {np.max(int_opt):.6f}")
    print(f"  - Intensidad promedio: {np.mean(int_opt):.6f}")
    print(f"  - Intensidad mínima: {np.min(int_opt):.6f}")
    print(f"  - Mejora respecto al inicial: {(TIEMPO_INTEGRACION/TIEMPO_INTEGRACION_INICIAL):.2f}x")

    fig_espectro = plt.figure(figsize=(10, 5))
    plt.plot(wl_opt, int_opt, linewidth=1.5, color='blue')
    plt.axhline(
        y=UMBRAL_SATURACION,
        color='red',
        linestyle='--',
        label=f'Umbral saturación: {UMBRAL_SATURACION:.3f} ({PORCENTAJE_SATURACION*100:.0f}%)'
    )
    plt.axhline(
        y=VALOR_MAXIMO_DETECTOR,
        color='orange',
        linestyle=':',
        alpha=0.5,
        label=f'Máximo detector: {VALOR_MAXIMO_DETECTOR}'
    )
    plt.xlabel('Longitud de onda (nm)')
    plt.ylabel('Intensidad')
    plt.title(f'Espectro con tiempo de integración optimizado ({TIEMPO_INTEGRACION*1000:.2f} ms)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    rutas_graficas.append(
        guardar_y_mostrar_figura(fig_espectro, ruta_series, "espectro_tiempo_integracion")
    )


# ============================================================
# ========== PREPARACIÓN DE REFERENCIAS PROCESADAS ===========
# ============================================================

refs_proc = preparar_referencias_procesadas(R_0, R_1)

print(f"\nButterworth aplicado a referencias (orden={ORDEN_FILTRO_BUTTER}, fc={FRECUENCIA_CORTE_BUTTER})")
print(f"Diezmado aplicado con factor = {refs_proc['factor_diezmado']}")
print(f"Puntos truncados originales: {len(refs_proc['wavelengths_trunc'])}")
print(f"Puntos finales por espectro para M_R: {len(refs_proc['wavelengths_diezmado'])}")

# Exportar R_0 y R_1 como hasta ahora (procesados, truncados y diezmados)
df_R0_export = pd.DataFrame({
    'wavelength_nm': refs_proc['wavelengths_diezmado'],
    'intensity': refs_proc['r0_diezmado']
})
df_R1_export = pd.DataFrame({
    'wavelength_nm': refs_proc['wavelengths_diezmado'],
    'intensity': refs_proc['r1_diezmado']
})

ruta_r0_csv = os.path.join(ruta_series, 'R_0_data.csv')
ruta_r1_csv = os.path.join(ruta_series, 'R_1_data.csv')

df_R0_export.to_csv(ruta_r0_csv, index=False)
df_R1_export.to_csv(ruta_r1_csv, index=False)

print("\nReferencias exportadas:")
print(f"  - R_0_data.csv: {ruta_r0_csv}")
print(f"  - R_1_data.csv: {ruta_r1_csv}")


# ============================================================
# ========== SERIE TEMPORAL DE R_M ===========================
# ============================================================

if not confirmar_simple("¿Iniciar serie temporal R_M?"):
    raise SystemExit("Proceso cancelado por el usuario antes de la serie temporal R_M.")

tiempos_relativos_rm, espectros_rm_raw, tiempo_total_temporal = adquirir_serie_temporal_rm(
    spec,
    TIEMPO_INTEGRACION,
    NUM_ESPECTROS_TEMPORAL
)

# Exportar R_M verdadero crudo
ruta_rm_trueraw_csv = os.path.join(ruta_series, "R_M_trueraw_tiempo_data.csv")
guardar_rm_trueraw_tiempo_csv(
    ruta_rm_trueraw_csv,
    tiempos_relativos_rm,
    espectros_rm_raw,
    wavelengths
)

print(f"Archivo crudo verdadero guardado en: {ruta_rm_trueraw_csv}")

# Exportar M_R temporal ya truncado, filtrado y diezmado
ruta_mr_tiempo_csv = os.path.join(ruta_series, "M_R_tiempo_data.csv")
guardar_mr_tiempo_csv(
    ruta_mr_tiempo_csv,
    tiempos_relativos_rm,
    espectros_rm_raw,
    refs_proc
)

print(f"Archivo M_R temporal guardado en: {ruta_mr_tiempo_csv}")


# ============================================================
# ========== RESUMEN FINAL ===================================
# ============================================================

print("\n" + "=" * 60)
print("RESUMEN FINAL DE LA ADQUISICIÓN")
print("=" * 60)
print(f"Sujeto registrado en: {ruta_sujeto}")
print(f"  - TXT de sujeto: {ruta_datos_sujeto}")
print(f"  - Carpeta series: {ruta_series}")
print()
print("Archivos exportados:")
print(f"  - R_0_data.csv: {ruta_r0_csv}")
print(f"  - R_1_data.csv: {ruta_r1_csv}")
print(f"  - R_M_trueraw_tiempo_data.csv: {ruta_rm_trueraw_csv}")
print(f"  - M_R_tiempo_data.csv: {ruta_mr_tiempo_csv}")
print()
print("Resumen de la serie temporal:")
print(f"  - Número total de espectros adquiridos: {len(espectros_rm_raw)}")
print(f"  - Tiempo total transcurrido: {tiempo_total_temporal:.3f} s")
print(f"  - Tiempo de integración usado: {TIEMPO_INTEGRACION*1000:.2f} ms")
print(f"  - Primera marca temporal: {tiempos_relativos_rm[0]:.6f} s")
print(f"  - Última marca temporal: {tiempos_relativos_rm[-1]:.6f} s")
print(f"  - Longitudes de onda crudas por espectro: {len(wavelengths)}")
print(f"  - Longitudes de onda en M_R por espectro (tras truncado + diezmado): {len(refs_proc['wavelengths_diezmado'])}")

if rutas_graficas:
    print("\nGráficas guardadas:")
    for ruta_grafica in rutas_graficas:
        print(f"  - {ruta_grafica}")

print("=" * 60)


# ============================================================
# ========== DESCONEXIÓN =====================================
# ============================================================

atexit.unregister(_cleanup_spec)
spec.Dispose()
print("Espectrofotómetro desconectado correctamente")