# demasiado-iad

Pipeline de adquisición espectral y reconstrucción de propiedades ópticas de tejido biológico mediante el método **Inverse Adding-Doubling (IAD)** de Scott Prahl.

El sistema conecta un espectrofotómetro **Thorlabs CCS200** con el ejecutable `iad.exe` para obtener los coeficientes de absorción (`mu_a`), scattering reducido (`mu_s'`) y factor de anisotropía (`g`) de muestras de piel humana in vivo.

---

## Estructura del repositorio

```
EmpaquetadoIADFull/
├── single_adq.py          # Adquisición de un solo espectro (R_0, R_1, R_M)
├── ad_temp.py             # Adquisición temporal (serie de espectros R_M con marcas de tiempo)
├── batch_IAD.py           # Ejecución batch del IAD sobre los datos adquiridos
├── viz_mr_temporal.py     # Visualización 3D y animación de reflectancia M_R temporal
├── viz_iad_temporal.py    # Visualización 3D y animación de parámetros IAD temporales
├── sample-F.rxt           # Plantilla .rxt con la geometría de la esfera integradora
├── IADSCOTT/              # Ejecutable iad.exe y libiad.dll de Scott Prahl
│   ├── iad.exe
│   ├── libiad.dll
│   └── ...
├── Mediciones/            # Datos de sujetos (generado automáticamente)
│   ├── calibraciones/     # Calibración persistente (R_0, R_1)
│   └── sujeto_XXX_.../    # Carpeta por sujeto con series y datos
├── IAD_run/               # Salida del batch IAD (generado automáticamente)
│   ├── por_lambda/        # Archivos .rxt temporales por longitud de onda
│   └── resumen_*.csv      # Resultados consolidados
├── M_R_data.csv           # Reflectancia de un solo espectro (entrada para batch IAD, modo normal)
├── M_R_tiempo_data.csv    # Reflectancia temporal (entrada para batch IAD, modo temporal)
├── TLCCS.c / TLCCS.fp     # Referencia del driver Thorlabs CCS (no se ejecutan directamente)
└── .gitignore
```

> Las carpetas `IADSCOTT/`, `IAD_run/` y `Mediciones/` están en `.gitignore` porque contienen binarios, datos generados o datos sensibles de sujetos.

---

## Requisitos

### Hardware
- Espectrofotómetro **Thorlabs CCS200** conectado por USB (solo para `single_adq.py` y `ad_temp.py`).

### Software
- **Windows** (los scripts de adquisición usan .NET interop para el driver Thorlabs).
- **Python 3.10+** con las siguientes dependencias:
  - `numpy`, `pandas`, `matplotlib`, `scipy`
  - `pythonnet` (`pip install pythonnet`) — necesario para `single_adq.py` y `ad_temp.py`.
- **Thorlabs CCS200 driver** instalado (deja su DLL interop en `C:\Program Files (x86)\Microsoft.NET\Primary Interop Assemblies\`).
- **iad.exe** de Scott Prahl dentro de `IADSCOTT/` (necesario para `batch_IAD.py`).
- Los scripts usan como plantilla el sample-F.rtx del directorio raíz para los parámetros necesarios, como el diámetro del haz, estándares de reflectancia, etc. M_R_data y M_R_tiempo_data.csv aún se necesitan colocar manualmente, extrayendolos de la carpeta de mediciones sujeto a sujeto.

> `batch_IAD.py`, `viz_mr_temporal.py` y `viz_iad_temporal.py` **no requieren** el espectrofotómetro ni pythonnet. Solo necesitan los CSV de datos.

---

## Flujo de trabajo general

El pipeline tiene **3 etapas**, cada una con su script. Las etapas 1A y 1B son alternativas entre sí (una u otra, según si quieres un solo espectro o una serie temporal):

```
┌────────────────────────────────────────────────────────────────────┐
│  ETAPA 1A: Adquisición simple          (single_adq.py)            │
│    Espectrofotómetro → R_0, R_1, R_M → M_R_data.csv              │
│                                                                    │
│  ETAPA 1B: Adquisición temporal        (ad_temp.py)               │
│    Espectrofotómetro → R_0, R_1, serie R_M → M_R_tiempo_data.csv │
└──────────────────────────────┬─────────────────────────────────────┘
                               │
                               ▼
┌────────────────────────────────────────────────────────────────────┐
│  ETAPA 2: Reconstrucción IAD           (batch_IAD.py)             │
│    M_R_data.csv ó M_R_tiempo_data.csv → iad.exe → mu_a, mu_s', g │
│    Salida: IAD_run/resumen_*.csv                                  │
└──────────────────────────────┬─────────────────────────────────────┘
                               │
                               ▼
┌────────────────────────────────────────────────────────────────────┐
│  ETAPA 3: Visualización                                           │
│    viz_mr_temporal.py   → Waterfall 3D + animación de M_R         │
│    viz_iad_temporal.py  → Waterfall 3D + animación de mu_a/mu_s'/g│
└────────────────────────────────────────────────────────────────────┘
```

---

## Scripts principales

### 1. `single_adq.py` — Adquisición de un solo espectro

**Qué hace:** Conecta con el espectrofotómetro Thorlabs CCS200, adquiere tres señales (R_0, R_1, R_M), calcula la reflectancia medida M_R y la exporta en CSV. Registra datos del sujeto y gestiona la calibración de forma persistente.

**Cuándo usarlo:** Cuando solo necesitas **un espectro estático** del tejido (sin evolución temporal).

**Flujo interactivo paso a paso:**

1. **Registro del sujeto** — El script pide nombre, edad, municipio de nacimiento, fototipo Fitzpatrick y afecciones conocidas. Crea automáticamente una carpeta `Mediciones/sujeto_XXX_dd_mm_YYYY/series/`.

2. **Calibración (R_0 y R_1)** — Si existe una calibración previa compatible en disco, pregunta si reutilizarla. Si no, pide realizar dos mediciones de referencia:
   - **R_0**: Señal de fondo (sin muestra, sin referencia). Se promedian `NUM_MEDICIONES_PROMEDIO` scans.
   - **R_1**: Señal del estándar de reflectancia (R_std = 0.99). Se promedian igual.

3. **Optimización del tiempo de integración** (opcional) — Aumenta gradualmente el tiempo de integración hasta justo antes de saturar el detector, maximizando la relación señal/ruido.

4. **Medición R_M** — Señal con la muestra de tejido sobre el puerto. Se promedian `NUM_MEDICIONES_PROMEDIO` scans.

5. **Procesamiento** — Trunca al rango de longitudes de onda configurado (`LAMBDA_MIN`–`LAMBDA_MAX`), aplica filtro Butterworth pasabajas, calcula M_R con la fórmula:
   ```
   M_R = R_std * (R_M - R_0) / (R_1 - R_0)
   ```
   y aplica diezmado para reducir el número de puntos espectrales.

6. **Exportación** — Guarda en la carpeta del sujeto:
   - `R_0_data.csv`, `R_1_data.csv`, `R_M_data.csv` — Señales procesadas.
   - `M_R_data.csv` — Reflectancia medida lista para IAD.
   - `datos_sujeto.txt` — Información del sujeto.
   - Gráficas PNG de los espectros.

**Parámetros ajustables (al inicio del archivo):**

| Parámetro | Descripción | Valor ejemplo |
|---|---|---|
| `TIEMPO_INTEGRACION_INICIAL` | Tiempo de integración de arranque (segundos) | `0.15` |
| `NUM_MEDICIONES_PROMEDIO` | Scans a promediar por cada señal | `5` |
| `TIEMPO_ESPERA` | Pausa entre scans consecutivos (segundos) | `0.01` |
| `PORCENTAJE_SATURACION` | Fracción del máximo del detector usada como umbral | `0.98` |
| `FRECUENCIA_CORTE_BUTTER` | Frecuencia de corte del filtro Butterworth (0–1, normalizada) | `0.1` |
| `ORDEN_FILTRO_BUTTER` | Orden del filtro Butterworth | `6` |
| `LAMBDA_MIN` / `LAMBDA_MAX` | Rango espectral de interés (nm) | `500` / `625` |
| `muestras_objetivo` | Puntos espectrales tras diezmado | `60` |

---

### 2. `ad_temp.py` — Adquisición temporal (serie de espectros)

**Qué hace:** Igual que `single_adq.py` para la calibración y registro de sujeto, pero en lugar de un solo R_M promediado, adquiere una **serie temporal** de N espectros consecutivos sin promediar, cada uno con su marca de tiempo. Exporta tanto los datos crudos como la reflectancia M_R procesada por espectro.

**Cuándo usarlo:** Cuando necesitas observar cómo cambian las propiedades ópticas del tejido **a lo largo del tiempo** (por ejemplo, durante la aplicación de un agente químico, calentamiento, etc.).

**Diferencias clave con `single_adq.py`:**
- R_M no se promedia: cada scan individual es un espectro de la serie temporal.
- Genera dos CSVs de salida adicionales:
  - `M_R_tiempo_data.csv` — Reflectancia procesada (truncada, filtrada, diezmada) con columnas `medicion, tiempo, lambda, reflectancia`.
  - `R_M_trueraw_tiempo_data.csv` — Intensidades crudas completas sin ningún procesamiento (3648 puntos por espectro, todas las longitudes de onda).
- Las marcas de tiempo se registran al inicio de cada scan con `time.perf_counter()`.
- La carpeta del sujeto se nombra con sufijo `_temporal`.

**Parámetros adicionales respecto a `single_adq.py`:**

| Parámetro | Descripción | Valor ejemplo |
|---|---|---|
| `NUM_ESPECTROS_TEMPORAL` | Número total de espectros en la serie temporal | `500` |
| `MUESTRAS_OBJETIVO` | Puntos espectrales tras diezmado (por espectro) | `200` |
| `R_STD` | Reflectancia del estándar de calibración | `0.99` |

**Formato de `M_R_tiempo_data.csv`:**
```csv
medicion,tiempo,lambda,reflectancia
1,0.000012,500.23,0.1842
1,0.000012,501.05,0.1851
...
2,0.1816,500.23,0.1839
...
```
Cada fila es un punto (lambda, reflectancia) de una medición dada. El campo `tiempo` es el tiempo relativo en segundos desde el inicio de la serie. El campo `medicion` es el índice del espectro (1-based).

---

### 3. `batch_IAD.py` — Ejecución batch del Inverse Adding-Doubling

**Qué hace:** Toma los datos de reflectancia (CSV) generados por los scripts de adquisición y ejecuta `iad.exe` de Scott Prahl **una vez por cada longitud de onda** para reconstruir las propiedades ópticas del tejido: `mu_a` (absorción), `mu_s'` (scattering reducido) y `g` (anisotropía).

**Cuándo usarlo:** Después de haber adquirido datos con `single_adq.py` o `ad_temp.py`. Este script **no requiere** el espectrofotómetro; solo necesita los CSVs y el ejecutable `iad.exe`.

**Tiene dos modos de operación, controlados por `MODO_TEMPORAL`:**

#### Modo normal (`MODO_TEMPORAL = False`)
- Lee `M_R_data.csv` (un solo espectro).
- Ejecuta IAD secuencialmente, una vez por lambda.
- Produce `IAD_run/resumen_resultados_laura_jacques.csv` y una gráfica de `mu_a` vs lambda.

#### Modo temporal (`MODO_TEMPORAL = True`)
- Lee `M_R_tiempo_data.csv` (serie temporal de espectros).
- Ejecuta IAD para cada combinación (medicion, lambda), opcionalmente en **paralelo** con múltiples workers.
- Permite submuestreo: procesar solo cada N espectros (`PASO_MEDICION`) o limitar el total (`MAX_MEDICIONES`).
- Produce `IAD_run/resumen_iad_temporal.csv` con columnas `medicion, tiempo, lambda_nm, mu_a_mm-1, mu_s_prime_mm-1, g, ...`.

**Cómo funciona internamente (por cada longitud de onda):**

1. Genera un archivo `.rxt` temporal con la geometría de la esfera integradora (tomada de `sample-F.rxt`) y el dato de reflectancia para esa lambda.
2. Calcula el scattering reducido inicial `mu_s'` usando el modelo Laura/Jacques: `mu_s'(λ) = 46 * (λ/500)^(-1.421)` cm⁻¹.
3. Calcula `g(λ)` con el polinomio de 6to orden de Ma et al. (o usa un `G_FIJO` configurable).
4. Invoca `iad.exe -g <g> -j <mu_s'> [-M 0] archivo.rxt`.
5. Parsea el archivo `.txt` de salida y extrae los resultados.
6. Limpia los archivos temporales.

**Parámetros ajustables:**

| Parámetro | Descripción | Valor ejemplo |
|---|---|---|
| `MODO_TEMPORAL` | `False` = un espectro, `True` = serie temporal | `True` |
| `USAR_G_FIJO` / `G_FIJO` | Si `True`, usa un g constante en lugar del polinomio | `True` / `0.8` |
| `USAR_MODO_RAPIDO` | Si `True`, pasa `-M 0` a IAD (omite verificación Monte Carlo) | `True` |
| `MAX_MEDICIONES` | Límite de espectros a procesar (`None` = todos) | `10` |
| `PASO_MEDICION` | Submuestreo: 1 = todos, 10 = 1 de cada 10 | `1` |
| `WORKERS` | Hilos para ejecución paralela (modo temporal) | `os.cpu_count()-2` |
| `MOSTRAR_COMANDOS` | Imprime el comando iad.exe de cada ejecución | `False` |

**Archivo `sample-F.rxt` (plantilla):**
Define la geometría física del sistema óptico: índice de refracción de la muestra y los slides, espesores, diámetro del haz, propiedades de la(s) esfera(s) integradora(s) y reflectividad del estándar de calibración. El script lo usa como encabezado y le inyecta los datos de reflectancia por lambda.

---

### 4. `viz_mr_temporal.py` — Visualización de reflectancia temporal

**Qué hace:** Lee `M_R_tiempo_data.csv` y genera dos visualizaciones interactivas:
1. **Waterfall 3D** — Gráfica 3D donde cada línea es un espectro de reflectancia a un instante de tiempo, con gradiente de color rojo de claro a oscuro a medida que avanza el tiempo.
2. **Animación con controles** — Gráfica 2D animada con slider para navegar entre espectros y botones Play/Pause.

**Cuándo usarlo:** Después de `ad_temp.py`, para inspeccionar visualmente cómo cambia el espectro de reflectancia en el tiempo **antes** de correr el IAD.

**Parámetros ajustables:**

| Parámetro | Descripción | Valor ejemplo |
|---|---|---|
| `LAMBDA_VIS_MIN` / `LAMBDA_VIS_MAX` | Rango de lambdas a visualizar (nm) | `530` / `580` |
| `ESPECTRO_INICIO` / `ESPECTRO_FIN` | Rango de espectros por índice (1-based) | `1` / `200` |
| `PROMEDIO_ESPECTROS` | Espectros consecutivos a promediar en el waterfall | `1` |

> Este script espera que `M_R_tiempo_data.csv` esté en la misma carpeta que el script. Cópialo desde la carpeta del sujeto o crea un enlace simbólico.

---

### 5. `viz_iad_temporal.py` — Visualización de parámetros IAD temporales

**Qué hace:** Lee `IAD_run/resumen_iad_temporal.csv` (la salida del batch IAD en modo temporal) y genera las mismas dos visualizaciones que `viz_mr_temporal.py`, pero para los **parámetros ópticos reconstruidos** en vez de la reflectancia.

**Cuándo usarlo:** Después de `batch_IAD.py` con `MODO_TEMPORAL = True`, para ver la evolución temporal de `mu_a`, `mu_s'` o `g`.

**Parámetros ajustables:**

| Parámetro | Descripción | Valor ejemplo |
|---|---|---|
| `PARAM` | Columna a graficar: `"mu_a_mm-1"`, `"mu_s_prime_mm-1"`, o `"g"` | `"mu_a_mm-1"` |
| `PARAM_LABEL` | Etiqueta para los ejes | `"mu_a (1/mm)"` |
| `LAMBDA_VIS_MIN` / `LAMBDA_VIS_MAX` | Rango de lambdas a visualizar | `530` / `580` |
| `ESPECTRO_INICIO` / `ESPECTRO_FIN` | Rango de espectros por índice | `1` / `None` |
| `PROMEDIO_ESPECTROS` | Espectros a promediar en el waterfall | `1` |

---

## Guía rápida: ejemplo de uso completo

### Caso A: Un solo espectro

```bash
# 1. Adquirir (requiere espectrofotómetro conectado)
python single_adq.py
# → Sigue las instrucciones interactivas (registro, calibración, R_0, R_1, R_M)
# → Genera Mediciones/sujeto_XXX_.../series/M_R_data.csv

# 2. Copiar el M_R_data.csv a la raíz del proyecto
copy Mediciones\sujeto_001_...\series\M_R_data.csv .

# 3. Configurar batch_IAD.py: MODO_TEMPORAL = False, luego ejecutar
python batch_IAD.py
# → Genera IAD_run/resumen_resultados_laura_jacques.csv
# → Genera IAD_run/grafica_mu_a_final.png
```

### Caso B: Serie temporal

```bash
# 1. Adquirir (requiere espectrofotómetro conectado)
python ad_temp.py
# → Sigue las instrucciones interactivas
# → Genera Mediciones/sujeto_XXX_..._temporal/series/M_R_tiempo_data.csv

# 2. Copiar el M_R_tiempo_data.csv a la raíz del proyecto
copy Mediciones\sujeto_001_..._temporal\series\M_R_tiempo_data.csv .

# 3. (Opcional) Visualizar reflectancia antes del IAD
python viz_mr_temporal.py

# 4. Configurar batch_IAD.py: MODO_TEMPORAL = True, ajustar MAX_MEDICIONES y WORKERS
python batch_IAD.py
# → Genera IAD_run/resumen_iad_temporal.csv

# 5. Visualizar resultados IAD
python viz_iad_temporal.py
```

---

## Calibración persistente

Ambos scripts de adquisición (`single_adq.py` y `ad_temp.py`) guardan la calibración (R_0, R_1) en `Mediciones/calibraciones/`:
- `calibracion_actual.npz` — Arrays de R_0 y R_1.
- `calibracion_actual.json` — Metadata (fecha, tiempo de integración, rango espectral).

En ejecuciones posteriores, si la calibración es compatible (mismas longitudes de onda), el script ofrece **reutilizarla** sin tener que medir R_0 y R_1 de nuevo. Solo es necesario recalibrar si cambia el montaje óptico o el sistema.

---

## Notas sobre el procesamiento de señal

- **Filtro Butterworth** — Se aplica un pasabajas de orden 6 a las señales crudas (R_0, R_1, R_M) antes de calcular la reflectancia, para suavizar ruido espectral de alta frecuencia.
- **Diezmado** — Tras filtrar, se reduce el número de puntos espectrales (el Butterworth actúa como filtro anti-aliasing). Esto es importante porque IAD se ejecuta una vez por lambda, y reducir puntos acorta dramáticamente el tiempo de cómputo.
- **Modelo de scattering** — `batch_IAD.py` estima el scattering reducido inicial `mu_s'(λ)` con un modelo potencial tipo Jacques/Laura: `mu_s' = 46 * (λ/500)^(-1.421)` cm⁻¹, basado en valores publicados para piel humana. IAD usa esto como punto de partida para la inversión.
- **Factor de anisotropía** — Puede calcularse por lambda con el polinomio de Ma et al. (2005) o fijarse a un valor constante (`G_FIJO`).
