"""
Visualización 3D waterfall y animación temporal de parámetros
recuperados por IAD (mu_a, mu_s', g) a partir de la serie temporal.

Uso: generar primero resumen_iad_temporal.csv con scriptIADLauraClaude.py
     (MODO_TEMPORAL = True), luego ejecutar este script.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.widgets import Button, Slider
import matplotlib.animation as animation


# ============================================================
# CONFIGURACIÓN
# ============================================================

# Archivo CSV de salida del modo temporal del IAD
CSV_NAME = "IAD_run/resumen_iad_temporal.csv"

# Parámetro a graficar — columnas disponibles:
#   "mu_a_mm-1"        Coeficiente de absorción (mm⁻¹)
#   "mu_s_prime_mm-1"  Coeficiente de scattering reducido (mm⁻¹)
#   "g"                Factor de anisotropía (adimensional)
PARAM = "mu_a_mm-1"
PARAM_LABEL = "mu_a (1/mm)"

# Rango de longitudes de onda para visualización (None = sin truncar)
LAMBDA_VIS_MIN = 500   # ej. 500
LAMBDA_VIS_MAX = 600   # ej. 600

# Rango de espectros a visualizar por índice (1-based, None = sin límite)
ESPECTRO_INICIO = 1    # ej. 50  → empieza en el espectro 50
ESPECTRO_FIN    = None # ej. 200 → termina en el espectro 200

# Número de espectros consecutivos a promediar para el waterfall 3D
PROMEDIO_ESPECTROS = 1


# ============================================================
# CARGA DE DATOS
# ============================================================

csv_path = Path(__file__).parent / CSV_NAME

if not csv_path.exists():
    raise FileNotFoundError(
        f"No se encontró '{CSV_NAME}' en {csv_path.parent}.\n"
        "Ejecuta scriptIADLauraClaude.py con MODO_TEMPORAL = True primero."
    )

df = pd.read_csv(csv_path)

# Pivotar: filas = mediciones (tiempo), columnas = lambda
pivot = df.pivot_table(index=["medicion", "tiempo"], columns="lambda_nm", values=PARAM)
pivot = pivot.sort_index()

tiempos = np.array([t for _, t in pivot.index])
lambdas_full = pivot.columns.values.astype(float)
espectros_full = pivot.values  # shape: (n_mediciones, n_lambdas)

# Truncar rango de longitudes de onda si se configuró
lam_min = LAMBDA_VIS_MIN if LAMBDA_VIS_MIN is not None else lambdas_full[0]
lam_max = LAMBDA_VIS_MAX if LAMBDA_VIS_MAX is not None else lambdas_full[-1]
mask_lambda = (lambdas_full >= lam_min) & (lambdas_full <= lam_max)

lambdas = lambdas_full[mask_lambda]
espectros = espectros_full[:, mask_lambda]

# Seleccionar rango de espectros por índice
idx_ini = (ESPECTRO_INICIO - 1) if ESPECTRO_INICIO is not None else 0
idx_fin = ESPECTRO_FIN if ESPECTRO_FIN is not None else len(tiempos)

tiempos = tiempos[idx_ini:idx_fin]
espectros = espectros[idx_ini:idx_fin]

n_mediciones = espectros.shape[0]

if LAMBDA_VIS_MIN is not None or LAMBDA_VIS_MAX is not None:
    print(f"Visualización truncada a [{lam_min:.1f}, {lam_max:.1f}] nm  ({len(lambdas)} puntos)")
print(f"Espectros: {idx_ini + 1}–{idx_ini + n_mediciones} de {len(pivot)}  "
      f"(t = {tiempos[0]:.3f}–{tiempos[-1]:.3f} s)")


# ============================================================
# 1. GRÁFICA 3D WATERFALL (con promediado)
# ============================================================

n_grupos = n_mediciones // PROMEDIO_ESPECTROS
espectros_avg = np.array([
    espectros[i * PROMEDIO_ESPECTROS:(i + 1) * PROMEDIO_ESPECTROS].mean(axis=0)
    for i in range(n_grupos)
])
tiempos_avg = np.array([
    tiempos[i * PROMEDIO_ESPECTROS:(i + 1) * PROMEDIO_ESPECTROS].mean()
    for i in range(n_grupos)
])

fig_3d = plt.figure(figsize=(12, 8))
ax3d = fig_3d.add_subplot(111, projection="3d")

# Gradiente azul: claro → saturado
colores_azul_3d = [
    mcolors.to_rgba((0.15, 0.35 + 0.50 * t, 1.0))
    for t in np.linspace(0, 1, n_grupos)
]

for i in range(n_grupos):
    ax3d.plot(lambdas, np.full_like(lambdas, tiempos_avg[i]), espectros_avg[i],
              linewidth=0.8, color=colores_azul_3d[i])

ax3d.set_xlabel("Longitud de onda (nm)", labelpad=10)
ax3d.set_ylabel("Tiempo (s)", labelpad=10)
ax3d.set_zlabel(PARAM_LABEL, labelpad=10)
ax3d.set_title(f"{PARAM_LABEL} — Waterfall 3D (IAD temporal)", fontsize=14, pad=15)
ax3d.grid(True, alpha=0.3)
ax3d.view_init(elev=25, azim=-60)
plt.tight_layout()


# ============================================================
# 2. ANIMACIÓN CON CONTROLES PLAY/PAUSE Y SLIDER
# ============================================================

fig_anim, ax_anim = plt.subplots(figsize=(10, 6))
plt.subplots_adjust(bottom=0.25)

z_min = np.nanmin(espectros) * 0.95
z_max = np.nanmax(espectros) * 1.05

# Gradiente azul para la animación (espectros individuales, sin promediar)
colores_azul = [
    mcolors.to_rgba((0.15, 0.35 + 0.50 * t, 1.0))
    for t in np.linspace(0, 1, n_mediciones)
]

line_anim, = ax_anim.plot(lambdas, espectros[0], linewidth=1.5, color=colores_azul[0])
ax_anim.set_xlim(lambdas[0], lambdas[-1])
ax_anim.set_ylim(z_min, z_max)
ax_anim.set_xlabel("Longitud de onda (nm)")
ax_anim.set_ylabel(PARAM_LABEL)
title_anim = ax_anim.set_title(
    f"Medición 1/{n_mediciones}  —  t = {tiempos[0]:.3f} s"
)
ax_anim.grid(True, alpha=0.3)

# Slider de tiempo
ax_slider = plt.axes([0.15, 0.10, 0.55, 0.04])
slider = Slider(ax_slider, "Medición", 1, n_mediciones, valinit=1, valstep=1)

# Botones Play / Pause
ax_play = plt.axes([0.78, 0.10, 0.08, 0.04])
ax_pause = plt.axes([0.87, 0.10, 0.08, 0.04])
btn_play = Button(ax_play, "Play", color="lightblue", hovercolor="steelblue")
btn_pause = Button(ax_pause, "Pause", color="lightsalmon", hovercolor="red")


class AnimState:
    playing = False


anim_state = AnimState()


def update_frame(idx):
    """Actualiza la gráfica al frame idx (0-based)."""
    line_anim.set_ydata(espectros[idx])
    line_anim.set_color(colores_azul[idx])
    title_anim.set_text(
        f"Medición {idx + 1}/{n_mediciones}  —  t = {tiempos[idx]:.3f} s"
    )
    slider.set_val(idx + 1)
    fig_anim.canvas.draw_idle()


def on_slider(val):
    idx = int(val) - 1
    line_anim.set_ydata(espectros[idx])
    line_anim.set_color(colores_azul[idx])
    title_anim.set_text(
        f"Medición {idx + 1}/{n_mediciones}  —  t = {tiempos[idx]:.3f} s"
    )
    fig_anim.canvas.draw_idle()


slider.on_changed(on_slider)


def animate(frame):
    if not anim_state.playing:
        return (line_anim,)
    idx = frame % n_mediciones
    update_frame(idx)
    return (line_anim,)


anim = animation.FuncAnimation(
    fig_anim, animate, frames=n_mediciones, interval=200, blit=False, repeat=True
)


def on_play(event):
    anim_state.playing = True


def on_pause(event):
    anim_state.playing = False


btn_play.on_clicked(on_play)
btn_pause.on_clicked(on_pause)


# ============================================================
# MOSTRAR TODO
# ============================================================

plt.show()
