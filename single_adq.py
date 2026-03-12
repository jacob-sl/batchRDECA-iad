# Importación de librerías necesarias
import atexit, json, os, re, sys, clr
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
spec = Thorlabs.ccs.interop64.TLCCS("USB0::0x1313::0x8089::M00496376::RAW", Boolean(True), Boolean(True))
print("Espectrofotómetro conectado exitosamente")

def _cleanup_spec():
    try:
        spec.Dispose()
        print("Espectrofotómetro desconectado (cleanup automático)")
    except Exception:
        pass

atexit.register(_cleanup_spec)

# ========== PARÁMETROS AJUSTABLES ==========
# Tiempo de integración en segundos (ej: 0.02 = 20ms, 0.1 = 100ms)
TIEMPO_INTEGRACION_INICIAL = 0.15
TIEMPO_INTEGRACION = TIEMPO_INTEGRACION_INICIAL
# Número de mediciones a promediar por cada adquisición
NUM_MEDICIONES_PROMEDIO = 5
# Tiempo de espera entre mediciones individuales (segundos)
TIEMPO_ESPERA = 0.01 
# ========== PARÁMETROS DE OPTIMIZACIÓN DEL TIEMPO DE INTEGRACIÓN ==========
# Umbral absoluto de saturación
# Para CCS200: el máximo es 1.0 (valores normalizados)
PORCENTAJE_SATURACION = 0.98  # 98% del máximo
VALOR_MAXIMO_DETECTOR = 1.0  # Valor máximo del detector (normalizado)
UMBRAL_SATURACION = VALOR_MAXIMO_DETECTOR * PORCENTAJE_SATURACION  # = 0.98

# ========== CONFIGURACIÓN DEL FILTRO BUTTERWORTH ==========
FRECUENCIA_CORTE_BUTTER = 0.1  # Frecuencia de corte normalizada (0-1, Nyquist=1)
ORDEN_FILTRO_BUTTER = 6         # Orden del filtro Butterworth

# ========== RANGO DE LONGITUDES DE ONDA (TRUNCADO) ==========
LAMBDA_MIN = 500  # Longitud de onda mínima en nm
LAMBDA_MAX = 625  # Longitud de onda máxima en nm
muestras_objetivo = 60  # Número de muestras deseadas después del diezmado
# ============================================================
print(f"Configuración:")
print(f"  - Tiempo de integración: {TIEMPO_INTEGRACION_INICIAL*1000} ms")
print(f"  - Mediciones por promedio: {NUM_MEDICIONES_PROMEDIO}")
print(f"  - Tiempo entre mediciones: {TIEMPO_ESPERA} s")
print(f"  - Butterworth: orden {ORDEN_FILTRO_BUTTER}, fc={FRECUENCIA_CORTE_BUTTER}")
print(f"Configuración de saturación:")
print(f"  - Valor máximo del detector: {VALOR_MAXIMO_DETECTOR}")
print(f"  - Porcentaje objetivo: {PORCENTAJE_SATURACION*100:.0f}%")
print(f"  - Umbral de saturación: {UMBRAL_SATURACION:.3f}")

#Definicion de funciones.
def tomar_medicion(spec, tiempo_integracion):
    """
    Toma una medición del espectrofotómetro con el tiempo de integración especificado.
    
    Parámetros:
    -----------
    spec : objeto TLCCS
        Objeto del espectrofotómetro Thorlabs
    tiempo_integracion : float
        Tiempo de integración en segundos
    
    Retorna:
    --------
    intensities : numpy.ndarray
        Array con las intensidades medidas
    
    Nota:
    -----
    OPTIMIZADO: No obtiene wavelengths en cada medición (son fijos).
    Usa la variable global 'wavelengths' obtenida una sola vez.
    """
    # Configurar tiempo de integración
    spec.setIntegrationTime(Double(tiempo_integracion))
    
    # Iniciar escaneo
    spec.startScan()
    
    # Obtener datos del escaneo
    scan = Array.CreateInstance(Double, 3648)
    spec.getScanData(scan)
    intensities = np.fromiter(scan, dtype=np.float64, count=3648)
    
    return intensities

def tomar_serie_adquisiciones(spec, tiempo_integracion, num_mediciones, tiempo_espera=0.1):
    """
    Toma múltiples mediciones y las promedia.
    
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
    
    Nota:
    -----
    OPTIMIZADO: Usa wavelengths global, no lo devuelve.
    """
    print(f"Tomando {num_mediciones} mediciones para promediar...")

    # Configurar el tiempo de integración una sola vez antes del loop
    spec.setIntegrationTime(Double(tiempo_integracion))

    mediciones = []

    for i in range(num_mediciones):
        spec.startScan()
        scan = Array.CreateInstance(Double, 3648)
        spec.getScanData(scan)
        intensity = np.fromiter(scan, dtype=np.float64, count=3648)
        mediciones.append(intensity)

        print(f"  Medición {i+1}/{num_mediciones} completada")

        if i < num_mediciones - 1:  # No esperar después de la última medición
            time.sleep(tiempo_espera)

    # Promediar todas las mediciones
    intensidad_promedio = np.mean(mediciones, axis=0)
    
    print(f"✓ Serie completada. Promedio calculado.\n")
    
    return intensidad_promedio

def optimizar_tiempo_integracion(spec, tiempo_inicial, umbral_saturacion_absoluto=60000, factor_incremento=1.15, max_iteraciones=20):
    """
    Optimiza el tiempo de integración AUMENTÁNDOLO gradualmente hasta el máximo sin saturar.
    
    Parámetros:
    -----------
    spec : objeto TLCCS
        Objeto del espectrofotómetro
    tiempo_inicial : float
        Tiempo de integración inicial en segundos (punto de partida conservador)
    umbral_saturacion_absoluto : float
        Umbral ABSOLUTO de saturación. Si la intensidad máxima supera este valor, usa el tiempo anterior.
        Ejemplo: 60000 para CCS200 (el máximo típico es ~65535 para 16 bits)
    factor_incremento : float
        Factor por el cual se incrementa el tiempo en cada iteración (ej: 1.15 = aumenta 15%)
    max_iteraciones : int
        Número máximo de intentos
    
    Retorna:
    --------
    tiempo_optimo : float
        Tiempo de integración óptimo encontrado (máximo sin saturar)
    intensidad_optima : numpy.ndarray
        Array de intensidad con el tiempo óptimo
    
    Nota:
    -----
    OPTIMIZADO: Usa wavelengths global, no lo devuelve.
    """
    print("=" * 60)
    print("🔧 OPTIMIZACIÓN DEL TIEMPO DE INTEGRACIÓN")
    print("=" * 60)
    print(f"Tiempo inicial: {tiempo_inicial*1000:.2f} ms")
    print(f"Umbral de saturación ABSOLUTO: {umbral_saturacion_absoluto:.0f}")
    print(f"Factor de incremento: {factor_incremento} (aumenta {(factor_incremento-1)*100:.0f}% por iteración)")
    print(f"Estrategia: Aumentar tiempo hasta encontrar el MÁXIMO sin saturar")
    print()
    
    tiempo_actual = tiempo_inicial
    tiempo_optimo = None
    intensidad_optima = None
    encontro_valor_valido = False
    
    for iteracion in range(max_iteraciones):
        # Tomar medición con el tiempo actual (solo intensidad)
        print(f"Iteración {iteracion + 1}/{max_iteraciones}:")
        print(f"  Probando con tiempo de integración: {tiempo_actual*1000:.2f} ms")
        
        intensidad = tomar_medicion(spec, tiempo_actual)
        
        # Encontrar el valor máximo
        max_intensidad = np.max(intensidad)
        
        print(f"  Intensidad máxima: {max_intensidad:.6f}")
        
        # Verificar si está por debajo del umbral ABSOLUTO
        if max_intensidad < umbral_saturacion_absoluto:
            margen = umbral_saturacion_absoluto - max_intensidad
            porcentaje_usado = (max_intensidad / umbral_saturacion_absoluto) * 100
            
            print(f"  ✓ Tiempo válido. Guardando como candidato óptimo.")
            print(f"  ✓ Usando {porcentaje_usado:.1f}% del umbral (margen: {margen:.6f})")
            
            # Guardar este como el mejor hasta ahora
            tiempo_optimo = tiempo_actual
            intensidad_optima = intensidad
            encontro_valor_valido = True
            
            # Intentar aumentar más
            tiempo_actual *= factor_incremento
            print(f"  → Intentando aumentar a {tiempo_actual*1000:.2f} ms\n")
        else:
            exceso = max_intensidad - umbral_saturacion_absoluto
            print(f"  ⚠ ¡SATURACIÓN DETECTADA! ({max_intensidad:.6f} > {umbral_saturacion_absoluto:.6f})")
            print(f"  ⚠ Exceso: {exceso:.6f}")
            
            if encontro_valor_valido:
                print(f"  → El tiempo anterior era el óptimo.\n")
            else:
                # Reducir el tiempo en lugar de quedarse con datos saturados
                tiempo_actual /= factor_incremento
                if tiempo_actual < 0.001:  # Mínimo 1 ms
                    print(f"  ⚠ Tiempo mínimo alcanzado (1 ms).")
                    tiempo_optimo = 0.001
                    intensidad_optima = tomar_medicion(spec, 0.001)
                    break
                print(f"  ⚠ Saturación detectada. Reduciendo a {tiempo_actual*1000:.2f} ms\n")
                continue
            break
    else:
        print(f"⚠ Alcanzado el número máximo de iteraciones.")
        if encontro_valor_valido:
            print(f"  Usando último tiempo válido: {tiempo_optimo*1000:.2f} ms")
        else:
            print(f"  No se encontró valor óptimo. Usando tiempo inicial: {tiempo_inicial*1000:.2f} ms")
            tiempo_optimo = tiempo_inicial
            # Tomar una medición final con el tiempo inicial
            intensidad_optima = tomar_medicion(spec, tiempo_inicial)
    
    print("\n" + "=" * 60)
    if encontro_valor_valido:
        print(f"TIEMPO DE INTEGRACIÓN ÓPTIMO: {tiempo_optimo*1000:.2f} ms")
        print(f"Intensidad máxima alcanzada: {np.max(intensidad_optima):.6f}")
        print(f"Porcentaje del umbral: {(np.max(intensidad_optima)/umbral_saturacion_absoluto)*100:.1f}%")
        print(f"Este es el MÁXIMO tiempo sin saturar el detector")
    else:
        print(f"⚠ NO SE PUDO OPTIMIZAR - Usando tiempo inicial: {tiempo_optimo*1000:.2f} ms")
        print(f"Intensidad con tiempo inicial: {np.max(intensidad_optima):.6f}")
        print(f"⚠ ADVERTENCIA: El tiempo inicial ya causa saturación. Reduce TIEMPO_INTEGRACION_INICIAL")
    print("=" * 60)
    
    return tiempo_optimo, intensidad_optima

def truncar_a_rango(intensidad_data):
    """
    Trunca datos de intensidad al rango LAMBDA_MIN-LAMBDA_MAX nm.
    
    Parámetros:
    -----------
    intensidad_data : numpy.ndarray
        Array con las intensidades medidas (3648 puntos)
    
    Retorna:
    --------
    wavelengths_trunc : numpy.ndarray
        Longitudes de onda truncadas
    intensidad_trunc : numpy.ndarray
        Intensidades truncadas
    """
    mask_rango = (wavelengths >= LAMBDA_MIN) & (wavelengths <= LAMBDA_MAX)
    return wavelengths[mask_rango], intensidad_data[mask_rango]

def aplicar_butterworth(datos, fc=None, orden=None):
    """
    Aplica filtro Butterworth pasabajas.
    Maneja NaN interpolando antes del filtrado y restaurándolos después.
    
    Parámetros:
    -----------
    datos : numpy.ndarray
        Array con los datos a filtrar
    fc : float
        Frecuencia de corte normalizada (0-1, Nyquist=1)
    orden : int
        Orden del filtro
    
    Retorna:
    --------
    datos_filtrados : numpy.ndarray
        Datos filtrados (NaN restaurados en posiciones originales)
    """
    if fc is None:
        fc = FRECUENCIA_CORTE_BUTTER
    if orden is None:
        orden = ORDEN_FILTRO_BUTTER

    # Manejar NaN: interpolar, filtrar, restaurar
    mask_nan = np.isnan(datos)
    if mask_nan.all():
        return datos.copy()

    datos_interp = datos.copy()
    if mask_nan.any():
        indices = np.arange(len(datos))
        datos_interp[mask_nan] = np.interp(
            indices[mask_nan], indices[~mask_nan], datos[~mask_nan]
        )

    b, a = butter(orden, fc, btype='low')
    datos_filtrados = filtfilt(b, a, datos_interp)

    # Restaurar NaN en posiciones originales
    datos_filtrados[mask_nan] = np.nan
    return datos_filtrados

def crear_carpeta_sujeto(base_dir):
    """
    Crea una carpeta con formato sujeto_XXX_dd_mm_YYYY y subcarpeta /series.
    """
    os.makedirs(base_dir, exist_ok=True)

    patron = re.compile(r"^s[ui]jeto_(\d{3})_\d{2}_\d{2}_\d{4}$")
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
        nombre_sujeto = f"sujeto_{siguiente_id:03d}_{fecha}"
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
    if len(r0) != len(wavelengths_actuales) or len(r1) != len(wavelengths_actuales):
        return False
    if len(wavelengths_guardadas) != len(wavelengths_actuales):
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

print("Funciones definidas correctamente")

# Obtener longitudes de onda (solo necesitamos hacerlo una vez)
ref16 = Int16(0)
nullable_double = Nullable[Double](0)
wvdata = Array.CreateInstance(Double, 3648)
spec.getWavelengthData(ref16, wvdata, nullable_double, nullable_double)
wavelengths = np.asarray(list(wvdata))

print(f"Longitudes de onda obtenidas: {len(wavelengths)} puntos")
print(f"Rango: {wavelengths[0]:.2f} nm - {wavelengths[-1]:.2f} nm")

# Registro de sujeto antes de realizar cualquier medición
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

ruta_salida_iad = os.path.join(ruta_sujeto, "IAD_salida.txt")
with open(ruta_salida_iad, "w", encoding="utf-8") as f:
    f.write("Salida del IAD pendiente de generar.\n")

print(f"Sujeto registrado en: {ruta_sujeto}")
print(f"  - TXT de sujeto: {ruta_datos_sujeto}")
print(f"  - Salida IAD: {ruta_salida_iad}")
print(f"  - Series: {ruta_series}")
rutas_graficas = []

# Gestión de calibración persistente (R_0, R_1)
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
    intensidad_opt = tomar_medicion(spec, TIEMPO_INTEGRACION)
    optimizacion_realizada = False

    print("\n✓ Calibración cargada desde disco.")
    print(f"  - Tiempo de integración usado: {TIEMPO_INTEGRACION*1000:.2f} ms")
    print(f"  - R_0 cargado. Shape: {R_0.shape}")
    print(f"  - R_1 cargado. Shape: {R_1.shape}")
else:
    # Confirmación por tecla para la optimización de tiempo de integración
    if confirmar_por_tecla("¿Deseas ejecutar la optimización de tiempo de integración?", tecla="o"):
        TIEMPO_INTEGRACION, intensidad_opt = optimizar_tiempo_integracion(
            spec,
            TIEMPO_INTEGRACION_INICIAL,
            umbral_saturacion_absoluto=UMBRAL_SATURACION,  # Umbral fijo absoluto
            factor_incremento=1.10,                         # Aumenta 10% en cada iteración
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

# Mostrar resultado solo si se realizó optimización
if optimizacion_realizada:
    mask_rango_opt = (wavelengths >= LAMBDA_MIN) & (wavelengths <= LAMBDA_MAX)
    wl_opt = wavelengths[mask_rango_opt]
    int_opt = intensidad_opt[mask_rango_opt]

    print(f"\nEstadísticas de la medición de referencia ({LAMBDA_MIN}-{LAMBDA_MAX} nm):")
    print(f"  - Intensidad máxima: {np.max(int_opt):.2f}")
    print(f"  - Intensidad promedio: {np.mean(int_opt):.2f}")
    print(f"  - Intensidad mínima: {np.min(int_opt):.2f}")
    print(f"  - Mejora respecto al inicial: {(TIEMPO_INTEGRACION/TIEMPO_INTEGRACION_INICIAL):.2f}x")

    # Graficar el espectro optimizado (truncado al rango LAMBDA_MIN-LAMBDA_MAX)
    fig_espectro = plt.figure(figsize=(10, 5))
    plt.plot(wl_opt, int_opt, linewidth=1.5, color='blue')
    plt.axhline(y=UMBRAL_SATURACION, color='red', linestyle='--', label=f'Umbral saturación: {UMBRAL_SATURACION:.0f} ({PORCENTAJE_SATURACION*100:.0f}%)')
    plt.axhline(y=VALOR_MAXIMO_DETECTOR, color='orange', linestyle=':', alpha=0.5, label=f'Máximo detector: {VALOR_MAXIMO_DETECTOR}')
    plt.xlabel('Longitud de onda (nm)')
    plt.ylabel('Intensidad')
    plt.title(f'Espectro con tiempo de integración MÁXIMO optimizado ({TIEMPO_INTEGRACION*1000:.2f} ms)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    rutas_graficas.append(
        guardar_y_mostrar_figura(fig_espectro, ruta_series, "espectro_tiempo_integracion")
    )

# TOMA DE R_M
if not confirmar_simple("¿Iniciar serie R_M?"):
    raise SystemExit("Proceso cancelado por el usuario antes de R_M.")

print("=" * 50)
print("SERIE R_M")
print("=" * 50)

R_M = tomar_serie_adquisiciones(
    spec, 
    TIEMPO_INTEGRACION, 
    NUM_MEDICIONES_PROMEDIO, 
    TIEMPO_ESPERA
)

print(f"R_M guardado. Shape: {R_M.shape}")
print(f"Rango de intensidades: {R_M.min():.6f} - {R_M.max():.6f}")

# ========== PROCESAMIENTO DE SEÑAL ==========
# 1. Truncar al rango de interés
wl_trunc, R_0_trunc = truncar_a_rango(R_0)
_, R_1_trunc = truncar_a_rango(R_1)
_, R_M_trunc = truncar_a_rango(R_M)

# 2. Aplicar filtro Butterworth ANTES de calcular reflectancia
R_0_butter = aplicar_butterworth(R_0_trunc)
R_1_butter = aplicar_butterworth(R_1_trunc)
R_M_butter = aplicar_butterworth(R_M_trunc)

print(f"Butterworth aplicado (orden={ORDEN_FILTRO_BUTTER}, fc={FRECUENCIA_CORTE_BUTTER})")

# 3. Cálculo de Reflectancia Medida (M_R)
# Fórmula: M_R = R_std * (R_M - R_0) / (R_1 - R_0)
R_STD = 0.99  # Reflectancia del estándar

numerador = R_M_butter - R_0_butter
denominador = R_1_butter - R_0_butter

# Evitar división por cero
denominador = np.where(denominador == 0, np.nan, denominador)

M_R_values = R_STD * (numerador / denominador)

# Crear DataFrames
df_R0 = pd.DataFrame({'wavelength_nm': wl_trunc, 'intensity': R_0_butter})
df_R1 = pd.DataFrame({'wavelength_nm': wl_trunc, 'intensity': R_1_butter})
df_RM = pd.DataFrame({'wavelength_nm': wl_trunc, 'intensity': R_M_butter})
df_MR = pd.DataFrame({'wavelength_nm': wl_trunc, 'reflectance': M_R_values})

# Mostrar estadísticas
print("=" * 50)
print("REFLECTANCIA MEDIDA (M_R)")
print("=" * 50)
print(f"Reflectancia del estándar: {R_STD}")
print(f"Rango de longitudes de onda: {df_MR['wavelength_nm'].min():.2f} - {df_MR['wavelength_nm'].max():.2f} nm")
print(f"Reflectancia mínima: {np.nanmin(M_R_values):.4f}")
print(f"Reflectancia máxima: {np.nanmax(M_R_values):.4f}")
print(f"Reflectancia promedio: {np.nanmean(M_R_values):.4f}")
print(df_MR.head())

# Graficar la reflectancia medida (solo M_R)
fig_mr = plt.figure(figsize=(12, 6))
plt.plot(df_MR['wavelength_nm'], df_MR['reflectance'], linewidth=1.5, color='green', label='Reflectancia Medida (M_R)')
plt.xlabel('Longitud de onda (nm)')
plt.ylabel('Reflectancia')
plt.ylim(0.05, 0.3)
plt.grid(True, alpha=0.3)
plt.legend()
plt.title('Reflectancia Medida del Tejido (M_R)')
plt.tight_layout()
rutas_graficas.append(
    guardar_y_mostrar_figura(fig_mr, ruta_series, "reflectancia_medida")
)

# 4. Diezmado para reducir puntos (para IAD)
# El Butterworth ya actúa como filtro anti-aliasing para el factor de diezmado usado
factor_diezmado = max(1, len(wl_trunc) // muestras_objetivo)

df_R0_diezmado = df_R0.iloc[::factor_diezmado].reset_index(drop=True)
df_R1_diezmado = df_R1.iloc[::factor_diezmado].reset_index(drop=True)
df_RM_diezmado = df_RM.iloc[::factor_diezmado].reset_index(drop=True)
df_MR_diezmado = df_MR.iloc[::factor_diezmado].reset_index(drop=True)

print(f"\nDiezmado aplicado con factor = {factor_diezmado}")
print(f"Puntos originales: {len(wl_trunc)} -> Puntos finales: {len(df_R0_diezmado)}")

# 5. Exportación de datos
df_R0_diezmado.to_csv(os.path.join(ruta_series, 'R_0_data.csv'), index=False)
df_R1_diezmado.to_csv(os.path.join(ruta_series, 'R_1_data.csv'), index=False)
df_RM_diezmado.to_csv(os.path.join(ruta_series, 'R_M_data.csv'), index=False)
df_MR_diezmado.to_csv(os.path.join(ruta_series, 'M_R_data.csv'), index=False)

print(f"\nDatos guardados en: {ruta_sujeto}")
print(f"  - TXT de sujeto: {ruta_datos_sujeto}")
print(f"  - Salida IAD: {ruta_salida_iad}")
print(f"  - Series: {ruta_series}")

# Graficar M_R diezmado
fig_mr_diezmado = plt.figure(figsize=(12, 6))
plt.plot(df_MR_diezmado['wavelength_nm'], df_MR_diezmado['reflectance'], linewidth=1.5, color='green', label='M_R diezmado', marker='o', markersize=3)
plt.xlabel('Longitud de onda (nm)')
plt.ylabel('Reflectancia')
plt.ylim(0.05, 0.3)
plt.grid(True, alpha=0.3)
plt.legend()
plt.title(f'Reflectancia Medida diezmada ({len(df_MR_diezmado)} puntos)')
plt.tight_layout()
rutas_graficas.append(
    guardar_y_mostrar_figura(fig_mr_diezmado, ruta_series, "MR_diezmado")
)

print("  - Gráficas guardadas:")
for ruta_grafica in rutas_graficas:
    print(f"    * {ruta_grafica}")

# Desconexión (atexit maneja el caso de error/excepción)
atexit.unregister(_cleanup_spec)
spec.Dispose()
print("Espectrofotómetro desconectado correctamente")
