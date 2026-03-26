# Recuperación de Espectros de coeficiente de absorcion en tanda

Pipeline de adquisicion espectral y reconstruccion de propiedades opticas de tejido biologico mediante el metodo **Inverse Adding-Doubling (IAD)** de Scott Prahl, y tomando como referencia el trabajo de B. Quistian y L. Sanchez, optimizado para las necesidades de L. Aguilar y su servidor.

El sistema conecta un espectrofotometro **Thorlabs CCS200** con el ejecutable `iad.exe` para obtener espectros de reflectancia. Posteriormente, se realiza un procesado utilizando el programa IAD de Scott Prahl, con el objetivo ultimo de recuperar los coeficientes de absorcion (`mu_a`) de muestras de piel humana in vivo, dados el coeficiente de scattering reducido (`mu_s'`) y factor de anisotropia (`g`). Estos scripts permiten hacer la recuperacion en el espectro completo para una sola medicion y para una serie de las mismas, asi como realizar una reconstruccion temporal para observar cambios en los `mu_a` debido a fenomenos dinamicos.

La recuperacion del coeficiente de absorcion del tejido no es un tema trivial. Al contar solo con datos experimentales de medicion de reflectancia, estamos forzados a elegir priors que reflejen las condiciones pertinentes para el problema inverso.

---

## Modelo de Scattering: Funcion Phan-Sierra

El repositorio incluye dos scripts de reconstruccion IAD:

- **`batch_IAD.py`** — Usa la funcion de Jacques/Laura para estimar `mu_s'(lambda)`: `mu_s' = 46 * (lambda/500)^(-1.421)` cm^-1. Esta aproximacion es generica y no siempre es valida para todas las regiones anatomicas.

- **`batch_IAD_funcion_sierra.py`** — Usa datos mas recientes de **Phan et al. (2021)** provenientes de *"Characterizing reduced scattering coefficient of normal human skin across different anatomic locations and Fitzpatrick skin types using spatial frequency domain imaging"* (JBO 26(2) 026001). A esta funcion le llamamos **Phan-Sierra**.

### Ubicaciones anatomicas disponibles

El script `batch_IAD_funcion_sierra.py` incluye datos de `mu_s'`, su desviacion estandar, y `mu_a` de referencia para **10 ubicaciones anatomicas** del paper de Phan et al., medidas por SFDI a 8 longitudes de onda (471, 526, 591, 621, 659, 691, 731, 851 nm):

| Ubicacion | Descripcion |
|---|---|
| `Forehead` | Frente |
| `Cheek` | Mejilla |
| `Ventral Forearm` | Antebrazo ventral |
| `Palm` | Palma de la mano |
| `Back` | Espalda |
| `Upper Arm` | Brazo superior |
| `Dorsal Forearm` | Antebrazo dorsal |
| `Neck` | Cuello |
| `Shin` | Espinilla |
| `Chest` | Pecho |

Para seleccionar la ubicacion, se modifica la linea de configuracion al inicio del script:

```python
PHAN_UBICACION = "Palm"  # Cambiar a cualquiera de las 10 opciones
```

### Variantes del modelo

La funcion Phan-Sierra esta implementada en dos variantes, controladas por `PHAN_SIERRA_MODE`:

- **`"powerlaw"`** — Ajuste de ley de potencia: `mu_s'(lambda) = A * (lambda/lambda_ref)^(-B)`, donde A y B se obtienen por regresion log-log sobre los 8 puntos de Phan.
- **`"pchip"`** — Interpolacion PCHIP (spline monotona) directa sobre los puntos de Phan. Mas fiel a los datos pero no extrapola bien fuera del rango 471–851 nm.

### Escenarios de sensibilidad (modo single)

En modo `single`, el script genera tres escenarios de `mu_s'` basados en la varianza reportada por Phan:

- **`nominal`** — Media de `mu_s'` de la ubicacion seleccionada.
- **`menos_1sd`** — Media - 1 desviacion estandar.
- **`mas_1sd`** — Media + 1 desviacion estandar.

Esto permite evaluar la sensibilidad de `mu_a` recuperado a la incertidumbre en `mu_s'`. En modo temporal, **solo se computa el escenario nominal** para evitar triplicar el volumen de datos y tiempo de computo.

### Referencia de mu_a en el dashboard

El dashboard de sensibilidad (`dashboard_mu_a.png`) incluye ademas una curva de referencia de `mu_a` de Phan et al. para la ubicacion seleccionada, graficada como linea negra punteada. Esto permite comparar el `mu_a` recuperado por IAD contra el valor de referencia publicado.

---

## Estructura del repositorio

```
EmpaquetadoIADFull/
├── single_adq.py                  # Adquisicion de un solo espectro (R_0, R_1, R_M)
├── ad_temp.py                     # Adquisicion temporal (serie de espectros R_M con marcas de tiempo)
├── batch_IAD.py                   # Batch IAD con modelo Jacques/Laura (legacy)
├── batch_IAD_funcion_sierra.py    # Batch IAD con modelo Phan-Sierra (principal)
├── viz_mr_temporal.py             # Visualizacion 3D y animacion de reflectancia M_R temporal
├── viz_iad_temporal.py            # Visualizacion 3D y animacion de parametros IAD temporales
├── iad_gui.py                     # Esbozo de GUI (no operativo)
├── sample-F.rxt                   # Plantilla .rxt con la geometria de la esfera integradora
├── IADSCOTT/                      # Ejecutable iad.exe y libiad.dll de Scott Prahl
│   ├── iad.exe
│   ├── libiad.dll
│   └── ...
├── Mediciones/                    # Datos de sujetos (generado automaticamente)
│   ├── calibraciones/             # Calibracion persistente (R_0, R_1)
│   ├── sujeto_XXX_fecha/         # Carpeta de sujeto (modo single)
│   │   └── series/
│   │       ├── M_R_data.csv
│   │       └── IAD_results/      # Resultados IAD por sujeto
│   └── sujeto_XXX_fecha_temporal/ # Carpeta de sujeto (modo temporal)
│       └── series/
│           ├── M_R_tiempo_data.csv
│           └── IAD_results/
├── IAD_run/                       # Espacio de trabajo temporal del IAD
│   ├── por_lambda/                # Archivos .rxt temporales por longitud de onda
│   ├── grafica_phan_sierra.png    # Visualizacion de la funcion Phan-Sierra
│   └── tabla_phan_sierra.csv      # Tabla numerica de la funcion Phan-Sierra
├── TLCCS.c / TLCCS.fp            # Referencia del driver Thorlabs CCS (no se ejecutan)
└── .gitignore
```

> Las carpetas `IADSCOTT/`, `IAD_run/` y `Mediciones/` estan en `.gitignore` porque contienen binarios, datos generados o datos sensibles de sujetos. Los contenidos del .zip descargable del repositorio de Scott Prahl deben colocarse dentro de `IADSCOTT/`. El .rxt que viene por defecto en ese zip se puede ignorar, ya que la plantilla utilizada para el procesado en tanda es `sample-F.rxt` en el directorio raiz del proyecto.

---

## Requisitos

### Hardware
- Espectrofotometro **Thorlabs CCS200** conectado por USB (solo para `single_adq.py` y `ad_temp.py`).

### Software
- **Windows** (los scripts de adquisicion usan .NET interop para el driver Thorlabs).
- **Python 3.10+** con las siguientes dependencias:
  - `numpy`, `pandas`, `matplotlib`, `scipy`
  - `pythonnet` (`pip install pythonnet`) — necesario para `single_adq.py` y `ad_temp.py`.
- **Thorlabs CCS200 driver** instalado (deja su DLL interop en `C:\Program Files (x86)\Microsoft.NET\Primary Interop Assemblies\`).
- **iad.exe** de Scott Prahl dentro de `IADSCOTT/`.

> `batch_IAD.py`, `batch_IAD_funcion_sierra.py`, `viz_mr_temporal.py` y `viz_iad_temporal.py` **no requieren** el espectrofotometro ni pythonnet. Solo necesitan los CSV de datos y `iad.exe`.

---

## Flujo de trabajo general

El pipeline tiene **3 etapas**. Las etapas 1A y 1B son alternativas entre si (una u otra, segun si necesitas un solo espectro o una serie temporal):

```
┌──────────────────────────────────────────────────────────────────────┐
│  ETAPA 1A: Adquisicion simple          (single_adq.py)              │
│    Espectrofotometro → R_0, R_1, R_M → M_R_data.csv                │
│                                                                      │
│  ETAPA 1B: Adquisicion temporal        (ad_temp.py)                 │
│    Espectrofotometro → R_0, R_1, serie R_M → M_R_tiempo_data.csv   │
└──────────────────────────────┬───────────────────────────────────────┘
                               │
                               ▼
┌──────────────────────────────────────────────────────────────────────┐
│  ETAPA 2: Reconstruccion IAD                                        │
│    batch_IAD_funcion_sierra.py (recomendado, Phan-Sierra)           │
│    batch_IAD.py (legacy, Jacques/Laura)                             │
│                                                                      │
│    Selector interactivo de sujeto → iad.exe → mu_a, mu_s', g       │
│    Salida: sujeto_XXX/series/IAD_results/                           │
└──────────────────────────────┬───────────────────────────────────────┘
                               │
                               ▼
┌──────────────────────────────────────────────────────────────────────┐
│  ETAPA 3: Visualizacion                                             │
│    viz_mr_temporal.py   → Waterfall 3D + animacion de M_R           │
│    viz_iad_temporal.py  → Waterfall 3D + animacion de mu_a/mu_s'/g  │
└──────────────────────────────────────────────────────────────────────┘
```

---

## Scripts principales

### 1. `single_adq.py` — Adquisicion de un solo espectro

**Que hace:** Conecta con el espectrofotometro Thorlabs CCS200, adquiere tres senales (R_0, R_1, R_M), calcula la reflectancia medida M_R y la exporta en CSV. Registra datos del sujeto y gestiona la calibracion de forma persistente.

**Cuando usarlo:** Cuando solo necesitas **un espectro estatico** del tejido (sin evolucion temporal).

**Flujo interactivo paso a paso:**

1. **Registro del sujeto** — El script pide nombre, edad, municipio de nacimiento, fototipo Fitzpatrick y afecciones conocidas. Crea automaticamente una carpeta `Mediciones/sujeto_XXX_dd_mm_YYYY/series/`.

2. **Calibracion (R_0 y R_1)** — Si existe una calibracion previa compatible en disco, pregunta si reutilizarla. Si no, pide realizar dos mediciones de referencia:
   - **R_0**: Senal de fondo (sin muestra, sin referencia). Se promedian `NUM_MEDICIONES_PROMEDIO` scans.
   - **R_1**: Senal del estandar de reflectancia (R_std = 0.99). Se promedian igual.

3. **Optimizacion del tiempo de integracion** (opcional) — Aumenta gradualmente el tiempo de integracion hasta justo antes de saturar el detector, maximizando la relacion senal/ruido.

4. **Medicion R_M** — Senal con la muestra de tejido sobre el puerto. Se promedian `NUM_MEDICIONES_PROMEDIO` scans.

5. **Procesamiento** — Trunca al rango de longitudes de onda configurado (`LAMBDA_MIN`–`LAMBDA_MAX`), aplica filtro Butterworth pasabajas, calcula M_R con la formula:
   ```
   M_R = R_std * (R_M - R_0) / (R_1 - R_0)
   ```
   y aplica diezmado para reducir el numero de puntos espectrales.

6. **Exportacion** — Guarda en la carpeta del sujeto:
   - `R_0_data.csv`, `R_1_data.csv`, `R_M_data.csv` — Senales procesadas.
   - `M_R_data.csv` — Reflectancia medida lista para IAD.
   - `datos_sujeto.txt` — Informacion del sujeto.
   - Graficas PNG de los espectros.

**Parametros ajustables (al inicio del archivo):**

| Parametro | Descripcion | Valor ejemplo |
|---|---|---|
| `TIEMPO_INTEGRACION_INICIAL` | Tiempo de integracion de arranque (segundos) | `0.15` |
| `NUM_MEDICIONES_PROMEDIO` | Scans a promediar por cada senal | `5` |
| `TIEMPO_ESPERA` | Pausa entre scans consecutivos (segundos) | `0.01` |
| `PORCENTAJE_SATURACION` | Fraccion del maximo del detector usada como umbral | `0.98` |
| `FRECUENCIA_CORTE_BUTTER` | Frecuencia de corte del filtro Butterworth (0–1, normalizada) | `0.1` |
| `ORDEN_FILTRO_BUTTER` | Orden del filtro Butterworth | `6` |
| `LAMBDA_MIN` / `LAMBDA_MAX` | Rango espectral de interes (nm) | `500` / `625` |
| `muestras_objetivo` | Puntos espectrales tras diezmado | `60` |

---

### 2. `ad_temp.py` — Adquisicion temporal (serie de espectros)

**Que hace:** Igual que `single_adq.py` para la calibracion y registro de sujeto, pero en lugar de un solo R_M promediado, adquiere una **serie temporal** de N espectros consecutivos sin promediar, cada uno con su marca de tiempo. Exporta tanto los datos crudos como la reflectancia M_R procesada por espectro.

**Cuando usarlo:** Cuando necesitas observar como cambian las propiedades opticas del tejido **a lo largo del tiempo** (por ejemplo, durante la aplicacion de un agente quimico, calentamiento, etc.).

**Diferencias clave con `single_adq.py`:**
- R_M no se promedia: cada scan individual es un espectro de la serie temporal.
- Genera dos CSVs de salida adicionales:
  - `M_R_tiempo_data.csv` — Reflectancia procesada (truncada, filtrada, diezmada) con columnas `medicion, tiempo, lambda, reflectancia`.
  - `R_M_trueraw_tiempo_data.csv` — Intensidades crudas completas sin ningun procesamiento (3648 puntos por espectro, todas las longitudes de onda).
- Las marcas de tiempo se registran al inicio de cada scan con `time.perf_counter()`.
- La carpeta del sujeto se nombra con sufijo `_temporal`.

**Parametros adicionales respecto a `single_adq.py`:**

| Parametro | Descripcion | Valor ejemplo |
|---|---|---|
| `NUM_ESPECTROS_TEMPORAL` | Numero total de espectros en la serie temporal | `500` |
| `MUESTRAS_OBJETIVO` | Puntos espectrales tras diezmado (por espectro) | `200` |
| `R_STD` | Reflectancia del estandar de calibracion | `0.99` |

**Formato de `M_R_tiempo_data.csv`:**
```csv
medicion,tiempo,lambda,reflectancia
1,0.000012,500.23,0.1842
1,0.000012,501.05,0.1851
...
2,0.1816,500.23,0.1839
...
```
Cada fila es un punto (lambda, reflectancia) de una medicion dada. El campo `tiempo` es el tiempo relativo en segundos desde el inicio de la serie. El campo `medicion` es el indice del espectro (1-based).

---

### 3. `batch_IAD_funcion_sierra.py` — Reconstruccion IAD con modelo Phan-Sierra (principal)

**Que hace:** Toma los datos de reflectancia generados por los scripts de adquisicion y ejecuta `iad.exe` de Scott Prahl para reconstruir las propiedades opticas del tejido: `mu_a` (absorcion), `mu_s'` (scattering reducido) y `g` (anisotropia). Utiliza datos de referencia de Phan et al. (2021) para 10 ubicaciones anatomicas como priors del `mu_s'`.

**Este es el script de reconstruccion recomendado.** El script legacy `batch_IAD.py` usa el modelo Jacques/Laura que no diferencia entre ubicaciones anatomicas.

**Tiene tres modos de operacion, controlados por `MODO`:**

#### Modo `"single"` — Un solo sujeto, un espectro
- Muestra un selector interactivo de sujetos disponibles en `Mediciones/`.
- Ejecuta IAD en 3 escenarios de `mu_s'`: nominal, -1sigma, +1sigma.
- Genera dentro de la carpeta del sujeto (`series/IAD_results/`):
  - `resumen_resultados_phan_sierra.csv` — Resultados completos.
  - `dashboard_mu_a.png` — Grafica de sensibilidad de `mu_a` con curva de referencia Phan.
  - Metricas de comparacion entre escenarios.

#### Modo `"temporal"` — Un solo sujeto, serie temporal
- Selector interactivo de sujetos temporales.
- Ejecuta IAD **solo para el escenario nominal** (sin +/-1sigma para evitar triplicar computo).
- Ejecucion **paralela** con `ThreadPoolExecutor` y multiples workers.
- Soporta promediado de espectros consecutivos (`VENTANA_PROMEDIO`) y filtrado por indice de medicion (`INICIO_MEDICION`, `MAX_MEDICIONES`).
- Genera `resumen_iad_temporal_phan_sierra.csv` con columnas `medicion, tiempo, lambda_nm, mu_a_mm-1, mu_s_prime_mm-1, g, ...`.

#### Modo `"batch_all"` — Todos los sujetos
- Recorre automaticamente todas las carpetas de sujeto en `Mediciones/`.
- Detecta si cada sujeto es single o temporal y procesa acorde.

**Como funciona internamente (por cada longitud de onda):**

1. Genera un archivo `.rxt` temporal con la geometria de la esfera integradora (tomada de `sample-F.rxt`) y el dato de reflectancia para esa lambda.
2. Calcula `mu_s'` usando la funcion Phan-Sierra (powerlaw o pchip) para la ubicacion anatomica seleccionada.
3. Usa `g` fijo (configurable) o el polinomio de Ma et al.
4. Invoca `iad.exe -g <g> -j <mu_s'> [-M 0] archivo.rxt`.
5. Parsea el archivo `.txt` de salida y extrae los resultados.
6. Limpia los archivos temporales.

**Parametros ajustables:**

| Parametro | Descripcion | Valor ejemplo |
|---|---|---|
| `MODO` | `"single"`, `"temporal"` o `"batch_all"` | `"temporal"` |
| `PHAN_UBICACION` | Ubicacion anatomica de referencia (ver tabla arriba) | `"Palm"` |
| `PHAN_SIERRA_MODE` | Modelo de interpolacion: `"powerlaw"` o `"pchip"` | `"pchip"` |
| `USAR_G_FIJO` / `G_FIJO` | Si `True`, usa un g constante en lugar del polinomio | `True` / `0.8` |
| `USAR_MODO_RAPIDO` | Si `True`, pasa `-M 0` a IAD (omite verificacion Monte Carlo) | `True` |
| `AUTO_ULTIMO_SUJETO` | Si `True`, selecciona el ultimo sujeto sin menu interactivo | `False` |
| `INICIO_MEDICION` | Indice de medicion desde el cual empezar (temporal) | `22` |
| `MAX_MEDICIONES` | Limite de espectros a procesar (`None` = todos) | `None` |
| `VENTANA_PROMEDIO` | Espectros consecutivos a promediar antes del IAD (temporal) | `50` |
| `WORKERS` | Hilos para ejecucion paralela (temporal) | `os.cpu_count()-1` |
| `MOSTRAR_COMANDOS` | Imprime el comando iad.exe de cada ejecucion | `False` |

**Archivo `sample-F.rxt` (plantilla):**
Define la geometria fisica del sistema optico: indice de refraccion de la muestra y los slides, espesores, diametro del haz, propiedades de la(s) esfera(s) integradora(s) y reflectividad del estandar de calibracion. El script lo usa como encabezado y le inyecta los datos de reflectancia por lambda.

---

### 4. `batch_IAD.py` — Reconstruccion IAD con modelo Jacques/Laura (legacy)

**Que hace:** Igual que `batch_IAD_funcion_sierra.py` pero estima `mu_s'(lambda)` con el modelo potencial de Jacques/Laura: `mu_s' = 46 * (lambda/500)^(-1.421)` cm^-1. No diferencia entre ubicaciones anatomicas.

**Cuando usarlo:** Solo si se necesita comparar contra el modelo clasico. Para nuevas mediciones, usar `batch_IAD_funcion_sierra.py`.

Tiene dos modos: normal (`MODO_TEMPORAL = False`) y temporal (`MODO_TEMPORAL = True`).

---

### 5. `viz_mr_temporal.py` — Visualizacion de reflectancia temporal

**Que hace:** Lee `M_R_tiempo_data.csv` de un sujeto temporal (con selector interactivo) y genera dos visualizaciones interactivas:
1. **Waterfall 3D** — Grafica 3D donde cada linea es un espectro de reflectancia a un instante de tiempo, con gradiente de color lila a rojo a medida que avanza el tiempo.
2. **Animacion con controles** — Grafica 2D animada con slider para navegar entre espectros y botones Play/Pause.

**Cuando usarlo:** Despues de `ad_temp.py`, para inspeccionar visualmente como cambia el espectro de reflectancia en el tiempo **antes** de correr el IAD.

**Parametros ajustables:**

| Parametro | Descripcion | Valor ejemplo |
|---|---|---|
| `LAMBDA_VIS_MIN` / `LAMBDA_VIS_MAX` | Rango de lambdas a visualizar (nm) | `520` / `580` |
| `ESPECTRO_INICIO` / `ESPECTRO_FIN` | Rango de espectros por indice (1-based) | `20` / `50` |
| `PROMEDIO_ESPECTROS` | Espectros consecutivos a promediar en el waterfall | `1` |

---

### 6. `viz_iad_temporal.py` — Visualizacion de parametros IAD temporales

**Que hace:** Lee `resumen_iad_temporal_phan_sierra.csv` de un sujeto temporal (con selector interactivo) y genera las mismas dos visualizaciones que `viz_mr_temporal.py`, pero para los **parametros opticos reconstruidos** en vez de la reflectancia.

**Cuando usarlo:** Despues de `batch_IAD_funcion_sierra.py` con `MODO = "temporal"`, para ver la evolucion temporal de `mu_a`, `mu_s'` o `g`.

**Parametros ajustables:**

| Parametro | Descripcion | Valor ejemplo |
|---|---|---|
| `PARAM` | Columna a graficar: `"mu_a_mm-1"`, `"mu_s_prime_mm-1"`, o `"g"` | `"mu_a_mm-1"` |
| `PARAM_LABEL` | Etiqueta para los ejes | `"mu_a (1/mm)"` |
| `LAMBDA_VIS_MIN` / `LAMBDA_VIS_MAX` | Rango de lambdas a visualizar | `520` / `580` |
| `ESPECTRO_INICIO` / `ESPECTRO_FIN` | Rango de espectros por indice | `None` / `None` |
| `PROMEDIO_ESPECTROS` | Espectros a promediar en el waterfall | `1` |

---

## Guia rapida: ejemplo de uso completo

### Caso A: Un solo espectro (single)

```bash
# 1. Adquirir (requiere espectrofotometro conectado)
python single_adq.py
# → Sigue las instrucciones interactivas (registro, calibracion, R_0, R_1, R_M)
# → Genera Mediciones/sujeto_XXX_fecha/series/M_R_data.csv

# 2. Configurar batch_IAD_funcion_sierra.py:
#    MODO = "single"
#    PHAN_UBICACION = "Palm"  (o la ubicacion anatomica correspondiente)

# 3. Ejecutar
python batch_IAD_funcion_sierra.py
# → Selector interactivo de sujeto
# → Genera en Mediciones/sujeto_XXX/series/IAD_results/:
#      resumen_resultados_phan_sierra.csv
#      dashboard_mu_a.png (con curva de referencia Phan)
```

### Caso B: Serie temporal

```bash
# 1. Adquirir (requiere espectrofotometro conectado)
python ad_temp.py
# → Sigue las instrucciones interactivas
# → Genera Mediciones/sujeto_XXX_fecha_temporal/series/M_R_tiempo_data.csv

# 2. (Opcional) Visualizar reflectancia antes del IAD
python viz_mr_temporal.py

# 3. Configurar batch_IAD_funcion_sierra.py:
#    MODO = "temporal"
#    PHAN_UBICACION = "Ventral Forearm"  (o la que corresponda)
#    Ajustar INICIO_MEDICION, MAX_MEDICIONES, VENTANA_PROMEDIO, WORKERS

# 4. Ejecutar
python batch_IAD_funcion_sierra.py
# → Selector interactivo de sujeto temporal
# → Genera resumen_iad_temporal_phan_sierra.csv

# 5. Visualizar resultados IAD
python viz_iad_temporal.py
```

### Caso C: Procesar todos los sujetos de golpe

```bash
# Configurar MODO = "batch_all" y ejecutar
python batch_IAD_funcion_sierra.py
# → Procesa automaticamente todos los sujetos en Mediciones/
```

---

## Calibracion persistente

Ambos scripts de adquisicion (`single_adq.py` y `ad_temp.py`) guardan la calibracion (R_0, R_1) en `Mediciones/calibraciones/`:
- `calibracion_actual.npz` — Arrays de R_0 y R_1.
- `calibracion_actual.json` — Metadata (fecha, tiempo de integracion, rango espectral).

En ejecuciones posteriores, si la calibracion es compatible (mismas longitudes de onda), el script ofrece **reutilizarla** sin tener que medir R_0 y R_1 de nuevo. Solo es necesario recalibrar si cambia el montaje optico o el sistema.

---

## Notas sobre el procesamiento de senal

- **Filtro Butterworth** — Se aplica un pasabajas de orden 6 a las senales crudas (R_0, R_1, R_M) antes de calcular la reflectancia, para suavizar ruido espectral de alta frecuencia.
- **Diezmado** — Tras filtrar, se reduce el numero de puntos espectrales (el Butterworth actua como filtro anti-aliasing). Esto es importante porque IAD se ejecuta una vez por lambda, y reducir puntos acorta dramaticamente el tiempo de computo.
- **Modelo de scattering (Phan-Sierra)** — `batch_IAD_funcion_sierra.py` estima `mu_s'(lambda)` usando datos especificos por ubicacion anatomica de Phan et al. (2021). Incluye tanto la interpolacion PCHIP como el ajuste powerlaw. El script tambien incorpora los datos de `mu_a` de referencia de Phan para comparacion en el dashboard.
- **Modelo de scattering (Jacques/Laura)** — `batch_IAD.py` usa el modelo generico: `mu_s' = 46 * (lambda/500)^(-1.421)` cm^-1.
- **Factor de anisotropia** — Puede calcularse por lambda con el polinomio de Ma et al. (2005) o fijarse a un valor constante (`G_FIJO = 0.8` por defecto).

---

## Referencias

- **Prahl, S. A.** — Inverse Adding-Doubling. Software y documentacion: [omlc.org/software/iad](https://omlc.org/software/iad/).
- **Phan, T. T. et al. (2021)** — *Characterizing reduced scattering coefficient of normal human skin across different anatomic locations and Fitzpatrick skin types using spatial frequency domain imaging.* Journal of Biomedical Optics, 26(2), 026001.
- **Jacques, S. L. (2013)** — *Optical properties of biological tissues: a review.* Physics in Medicine & Biology, 58(11), R37.
- **Ma, X. et al. (2005)** — Polinomio de 6to orden para el factor de anisotropia g(lambda) de piel humana.
