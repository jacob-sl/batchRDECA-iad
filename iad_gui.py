"""
Pipeline IAD Completo - Interfaz Grafica Unificada
===================================================
Integra los 5 scripts del pipeline de propiedades opticas de tejido:
  1. Adquisicion individual y temporal (ad_temp.py / single_adq.py)
  2. Procesamiento IAD batch (batch_IAD.py)
  3. Visualizacion M_R individual y temporal (viz_mr_temporal.py)
  4. Visualizacion IAD individual y temporal (viz_iad_temporal.py)

Framework: tkinter + ttkbootstrap
"""

import atexit
import json
import os
import re
import subprocess
import sys
import threading
import time
import queue
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from scipy.signal import butter, filtfilt

import tkinter as tk

try:
    import ttkbootstrap as ttk
    from ttkbootstrap.constants import *
    from ttkbootstrap.scrolled import ScrolledText
    HAS_BOOTSTRAP = True
except ImportError:
    import tkinter.ttk as ttk
    from tkinter.scrolledtext import ScrolledText
    HAS_BOOTSTRAP = False

import tkinter.filedialog as filedialog
import tkinter.messagebox as messagebox


# ============================================================
# SECCION 1: CONSTANTES Y CONFIGURACION POR DEFECTO
# ============================================================

DEFAULTS = {
    "single": {
        "tiempo_integracion_inicial": 0.15,
        "num_mediciones_promedio": 5,
        "tiempo_espera": 0.01,
        "porcentaje_saturacion": 0.98,
        "valor_maximo_detector": 1.0,
        "frecuencia_corte_butter": 0.1,
        "orden_filtro_butter": 6,
        "lambda_min": 500,
        "lambda_max": 625,
        "muestras_objetivo": 60,
        "r_std": 0.99,
    },
    "temporal": {
        "tiempo_integracion_inicial": 0.1,
        "num_mediciones_promedio": 10,
        "tiempo_espera": 0.01,
        "num_espectros_temporal": 500,
        "porcentaje_saturacion": 0.95,
        "valor_maximo_detector": 1.0,
        "frecuencia_corte_butter": 0.1,
        "orden_filtro_butter": 6,
        "lambda_min": 500,
        "lambda_max": 590,
        "muestras_objetivo": 200,
        "r_std": 0.99,
    },
    "iad": {
        "usar_g_fijo": True,
        "g_fijo": 0.8,
        "usar_modo_rapido": True,
        "workers": max(1, (os.cpu_count() or 4) - 2),
        "max_mediciones": 10,
        "paso_medicion": 1,
    },
    "viz": {
        "lambda_vis_min": 530,
        "lambda_vis_max": 580,
        "espectro_inicio": 1,
        "espectro_fin": 200,
        "promedio_espectros": 1,
        "param_iad": "mu_a_mm-1",
    },
}

BASE_DIR = Path(__file__).parent
IAD_EXE_DEFAULT = BASE_DIR / "IADSCOTT" / "iad.exe"
RXT_TEMPLATE_DEFAULT = BASE_DIR / "sample-F.rxt"
MEDICIONES_DIR = BASE_DIR / "Mediciones"
CALIBRACIONES_DIR = MEDICIONES_DIR / "calibraciones"

# Rutas por defecto de los scripts originales
MR_SINGLE_CSV_DEFAULT = BASE_DIR / "M_R_data.csv"
MR_TEMPORAL_CSV_DEFAULT = BASE_DIR / "M_R_tiempo_data.csv"
IAD_OUTPUT_DIR_DEFAULT = BASE_DIR / "IAD_run"
IAD_TEMPORAL_CSV_DEFAULT = IAD_OUTPUT_DIR_DEFAULT / "resumen_iad_temporal.csv"
IAD_SINGLE_CSV_DEFAULT = IAD_OUTPUT_DIR_DEFAULT / "resumen_resultados_laura_jacques.csv"


# ============================================================
# SECCION 2: SPECTROMETER DRIVER
# ============================================================

class SpectrometerDriver:
    """Encapsula toda la interaccion con el espectrofotometro Thorlabs CCS200."""

    def __init__(self):
        self.spec = None
        self.wavelengths = None
        self.connected = False
        self._Double = None
        self._Array = None
        self._Int16 = None
        self._Nullable = None
        self._Boolean = None

    def connect(self):
        try:
            import clr
            sys.path.append(r"C:\Program Files (x86)\Microsoft.NET\Primary Interop Assemblies\\")
            clr.AddReference("Thorlabs.ccs.interop64")
            import Thorlabs.ccs.interop64
            from System import Double, Boolean, Array, Int16, Nullable

            self._Double = Double
            self._Array = Array
            self._Int16 = Int16
            self._Nullable = Nullable
            self._Boolean = Boolean

            self.spec = Thorlabs.ccs.interop64.TLCCS(
                "USB0::0x1313::0x8089::M00496376::RAW",
                Boolean(True), Boolean(True)
            )
            self.connected = True

            ref16 = Int16(0)
            nullable_double = Nullable[Double](0)
            wvdata = Array.CreateInstance(Double, 3648)
            self.spec.getWavelengthData(ref16, wvdata, nullable_double, nullable_double)
            self.wavelengths = np.asarray(list(wvdata))

            return True, f"Conectado. {len(self.wavelengths)} puntos ({self.wavelengths[0]:.1f}-{self.wavelengths[-1]:.1f} nm)"
        except ImportError:
            self.connected = False
            return False, "pythonnet no instalado. Solo procesamiento y visualizacion disponibles."
        except Exception as e:
            self.connected = False
            return False, f"No se pudo conectar: {e}"

    def disconnect(self):
        if self.spec is not None:
            try:
                self.spec.Dispose()
            except Exception:
                pass
            self.spec = None
            self.connected = False

    def take_single_measurement(self, integration_time):
        self.spec.setIntegrationTime(self._Double(integration_time))
        self.spec.startScan()
        scan = self._Array.CreateInstance(self._Double, 3648)
        self.spec.getScanData(scan)
        return np.fromiter(scan, dtype=np.float64, count=3648)

    def take_averaged_series(self, integration_time, n_measurements, wait_time,
                             progress_cb=None, log_cb=None):
        self.spec.setIntegrationTime(self._Double(integration_time))
        mediciones = []
        for i in range(n_measurements):
            self.spec.startScan()
            scan = self._Array.CreateInstance(self._Double, 3648)
            self.spec.getScanData(scan)
            intensity = np.fromiter(scan, dtype=np.float64, count=3648)
            mediciones.append(intensity)
            if progress_cb:
                progress_cb(i + 1, n_measurements)
            if log_cb:
                log_cb(f"  Medicion {i+1}/{n_measurements} completada")
            if i < n_measurements - 1:
                time.sleep(wait_time)
        return np.mean(mediciones, axis=0)

    def optimize_integration_time(self, initial_time, saturation_threshold,
                                  increment_factor=1.10, max_iterations=30,
                                  progress_cb=None, log_cb=None):
        tiempo_actual = initial_time
        tiempo_optimo = None
        intensidad_optima = None
        encontro_valido = False

        for iteracion in range(max_iterations):
            if log_cb:
                log_cb(f"Iter {iteracion+1}/{max_iterations}: probando {tiempo_actual*1000:.2f} ms")
            intensidad = self.take_single_measurement(tiempo_actual)
            max_int = np.max(intensidad)
            if log_cb:
                log_cb(f"  Max intensidad: {max_int:.6f}")
            if progress_cb:
                progress_cb(iteracion + 1, max_iterations)

            if max_int < saturation_threshold:
                porcentaje = (max_int / saturation_threshold) * 100
                if log_cb:
                    log_cb(f"  Valido ({porcentaje:.1f}% del umbral)")
                tiempo_optimo = tiempo_actual
                intensidad_optima = intensidad
                encontro_valido = True
                tiempo_actual *= increment_factor
            else:
                if log_cb:
                    log_cb(f"  SATURACION ({max_int:.6f} > {saturation_threshold:.6f})")
                if encontro_valido:
                    break
                else:
                    tiempo_actual /= increment_factor
                    if tiempo_actual < 0.001:
                        tiempo_optimo = 0.001
                        intensidad_optima = self.take_single_measurement(0.001)
                        break
                    continue

        if not encontro_valido and tiempo_optimo is None:
            tiempo_optimo = initial_time
            intensidad_optima = self.take_single_measurement(initial_time)

        if log_cb:
            log_cb(f"Tiempo optimo: {tiempo_optimo*1000:.2f} ms")
        return tiempo_optimo, intensidad_optima

    def acquire_temporal_series(self, integration_time, n_spectra,
                                progress_cb=None, log_cb=None):
        tiempos_relativos = []
        espectros_raw = []
        self.spec.setIntegrationTime(self._Double(integration_time))
        inicio = time.perf_counter()

        for medicion in range(1, n_spectra + 1):
            t_rel = time.perf_counter() - inicio
            self.spec.startScan()
            scan = self._Array.CreateInstance(self._Double, 3648)
            self.spec.getScanData(scan)
            intensidad = np.fromiter(scan, dtype=np.float64, count=3648)
            tiempos_relativos.append(t_rel)
            espectros_raw.append(intensidad)
            if progress_cb:
                progress_cb(medicion, n_spectra)
            if log_cb and medicion % 50 == 0:
                log_cb(f"  Espectro {medicion}/{n_spectra}")

        tiempo_total = time.perf_counter() - inicio
        if log_cb:
            log_cb(f"Serie completada. {len(espectros_raw)} espectros en {tiempo_total:.3f} s")
        return np.array(tiempos_relativos, dtype=np.float64), espectros_raw, tiempo_total


# ============================================================
# SECCION 3: SIGNAL PROCESSING
# ============================================================

class SignalProcessing:

    @staticmethod
    def truncate_to_range(wavelengths, intensity, lambda_min, lambda_max):
        mask = (wavelengths >= lambda_min) & (wavelengths <= lambda_max)
        return wavelengths[mask], intensity[mask]

    @staticmethod
    def apply_butterworth(data, fc, order):
        mask_nan = np.isnan(data)
        if mask_nan.all():
            return data.copy()
        datos_interp = data.copy()
        if mask_nan.any():
            indices = np.arange(len(data))
            datos_interp[mask_nan] = np.interp(
                indices[mask_nan], indices[~mask_nan], data[~mask_nan]
            )
        b, a = butter(order, fc, btype='low')
        datos_filtrados = filtfilt(b, a, datos_interp)
        datos_filtrados[mask_nan] = np.nan
        return datos_filtrados

    @staticmethod
    def compute_reflectance(r_m, r_0, r_1, r_std=0.99):
        numerador = r_m - r_0
        denominador = r_1 - r_0
        denominador = np.where(denominador == 0, np.nan, denominador)
        return r_std * (numerador / denominador)

    @staticmethod
    def decimate(data, target_samples):
        factor = max(1, len(data) // target_samples)
        return data[::factor], factor

    @staticmethod
    def prepare_processed_references(r0, r1, wavelengths, lambda_min, lambda_max,
                                     fc, order, target_samples, r_std=0.99):
        sp = SignalProcessing
        wl_trunc, r0_trunc = sp.truncate_to_range(wavelengths, r0, lambda_min, lambda_max)
        _, r1_trunc = sp.truncate_to_range(wavelengths, r1, lambda_min, lambda_max)
        r0_butter = sp.apply_butterworth(r0_trunc, fc, order)
        r1_butter = sp.apply_butterworth(r1_trunc, fc, order)
        factor = max(1, len(wl_trunc) // target_samples)
        wl_dec = wl_trunc[::factor]
        r0_dec = r0_butter[::factor]
        r1_dec = r1_butter[::factor]
        denominador = r1_butter - r0_butter
        denominador = np.where(denominador == 0, np.nan, denominador)
        return {
            "wavelengths_trunc": wl_trunc,
            "wavelengths_diezmado": wl_dec,
            "r0_butter": r0_butter,
            "r1_butter": r1_butter,
            "r0_diezmado": r0_dec,
            "r1_diezmado": r1_dec,
            "factor_diezmado": factor,
            "denominador": denominador,
        }


# ============================================================
# SECCION 4: IAD PROCESSOR
# ============================================================

class IADProcessor:

    def __init__(self, iad_exe_path, rxt_template_path):
        self.iad_exe = Path(iad_exe_path)
        self.rxt_template = Path(rxt_template_path)
        self.header_lines = None

    def validate_paths(self):
        msgs = []
        if not self.iad_exe.exists():
            msgs.append(f"iad.exe no encontrado: {self.iad_exe}")
        if not self.rxt_template.exists():
            msgs.append(f"Plantilla RXT no encontrada: {self.rxt_template}")
        return len(msgs) == 0, "\n".join(msgs)

    def load_template_header(self):
        with self.rxt_template.open("r", encoding="utf-8") as f:
            lines = f.readlines()
        header = []
        for line in lines:
            if line.strip().startswith("#lambda"):
                break
            header.append(line.rstrip("\n"))
        self.header_lines = header

    @staticmethod
    def mu_sp_laura_jacques_mm(lambda_nm):
        A_CM = 46.0
        B_EXP = 1.421
        return (A_CM * (lambda_nm / 500.0) ** (-B_EXP)) / 10.0

    @staticmethod
    def g_ma_et_al(lambda_nm):
        L = lambda_nm
        g = (-5.603 + 3.61e-2 * L - 8.17e-5 * L**2 + 9.51e-8 * L**3
             - 5.92e-11 * L**4 + 1.83e-14 * L**5 - 2.11e-18 * L**6)
        return float(max(0.5, min(0.99, g)))

    def _build_rxt(self, wavelength, reflectance, output_path):
        with output_path.open("w", encoding="utf-8", newline="\n") as f:
            for line in self.header_lines:
                f.write(line + "\n")
            f.write("\n#lambda\tM_R\n")
            f.write(f"{wavelength:.10f}\t{reflectance:.10f}\n")

    def _run_iad_one_lambda(self, wavelength, reflectance, per_lambda_dir,
                            use_fixed_g=True, g_fixed=0.8, fast_mode=True, med_id=0):
        mu_sp_mm = self.mu_sp_laura_jacques_mm(wavelength)
        g_valor = g_fixed if use_fixed_g else self.g_ma_et_al(wavelength)

        if med_id != 0:
            base_name = f"med{med_id:06d}_lambda_{wavelength:.2f}".replace(".", "p")
        else:
            base_name = f"lambda_{wavelength:.2f}".replace(".", "p")
        rxt_path = per_lambda_dir / f"{base_name}.rxt"

        self._build_rxt(wavelength, reflectance, rxt_path)

        iad_args = ["-g", f"{g_valor:.6f}", "-j", f"{mu_sp_mm:.6f}"]
        if fast_mode:
            iad_args.extend(["-M", "0"])

        cmd = [str(self.iad_exe), *iad_args, rxt_path.name]
        subprocess.run(cmd, cwd=per_lambda_dir, text=True, shell=False,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        txt_path = per_lambda_dir / f"{rxt_path.stem}.txt"
        if not txt_path.exists():
            candidatos = sorted(per_lambda_dir.glob(f"{rxt_path.stem}*.txt"))
            txt_path = candidatos[0] if candidatos else None

        fila = {
            "lambda_nm": wavelength,
            "reflectance_input": reflectance,
            "mu_s_prime_input_mm-1": mu_sp_mm,
            "g_input": g_valor,
        }

        if txt_path and txt_path.exists():
            extraido = self._parse_iad_output(txt_path)
            if extraido:
                fila.update(extraido)
            txt_path.unlink(missing_ok=True)
        rxt_path.unlink(missing_ok=True)
        return fila

    @staticmethod
    def _parse_iad_output(txt_path):
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
                    "mu_a_mm-1": float(partes[5]),
                    "mu_s_prime_mm-1": float(partes[6]),
                    "g": float(partes[7]),
                }
        return None

    @staticmethod
    def read_single_csv(csv_path):
        df = pd.read_csv(csv_path)
        if not {"wavelength_nm", "reflectance"}.issubset(df.columns):
            raise ValueError("CSV debe tener columnas: wavelength_nm, reflectance")
        return df[["wavelength_nm", "reflectance"]].dropna()

    @staticmethod
    def read_temporal_csv(csv_path):
        df = pd.read_csv(csv_path)
        if not {"medicion", "tiempo", "lambda", "reflectancia"}.issubset(df.columns):
            raise ValueError("CSV temporal debe tener: medicion, tiempo, lambda, reflectancia")
        mediciones = {}
        for med_id, grupo in df.groupby("medicion", sort=True):
            tiempo = float(grupo["tiempo"].iloc[0])
            espectro = sorted(zip(grupo["lambda"], grupo["reflectancia"]), key=lambda x: x[0])
            mediciones[int(med_id)] = {"tiempo": tiempo, "espectro": list(espectro)}
        return mediciones

    def run_single_spectrum(self, csv_path, output_dir, use_fixed_g, g_fixed,
                            fast_mode, progress_cb=None):
        self.load_template_header()
        df = self.read_single_csv(csv_path)
        per_lambda_dir = Path(output_dir) / "por_lambda"
        per_lambda_dir.mkdir(parents=True, exist_ok=True)
        for f in per_lambda_dir.iterdir():
            if f.is_file():
                f.unlink()
        resultados = []
        total = len(df)
        for i, (_, row) in enumerate(df.iterrows(), 1):
            fila = self._run_iad_one_lambda(
                float(row["wavelength_nm"]), float(row["reflectance"]),
                per_lambda_dir, use_fixed_g, g_fixed, fast_mode
            )
            resultados.append(fila)
            if progress_cb:
                progress_cb(i, total)
        df_res = pd.DataFrame(resultados).sort_values("lambda_nm")
        ruta_csv = Path(output_dir) / "resumen_resultados_laura_jacques.csv"
        df_res.to_csv(ruta_csv, index=False)
        return df_res, ruta_csv

    def run_temporal(self, csv_path, output_dir, use_fixed_g, g_fixed,
                     fast_mode, workers, max_mediciones, paso,
                     progress_cb=None):
        self.load_template_header()
        mediciones = self.read_temporal_csv(csv_path)
        per_lambda_dir = Path(output_dir) / "por_lambda"
        per_lambda_dir.mkdir(parents=True, exist_ok=True)
        for f in per_lambda_dir.iterdir():
            if f.is_file():
                f.unlink()

        ids = sorted(mediciones.keys())
        if max_mediciones:
            ids = ids[:max_mediciones]
        ids = ids[::paso]

        tareas = [
            (med_id, mediciones[med_id]["tiempo"], float(wl), float(ref))
            for med_id in ids
            for wl, ref in mediciones[med_id]["espectro"]
        ]
        total = len(tareas)

        def _tarea(med_id, tiempo, wl, ref):
            fila = self._run_iad_one_lambda(
                wl, ref, per_lambda_dir, use_fixed_g, g_fixed, fast_mode, med_id
            )
            fila["medicion"] = med_id
            fila["tiempo"] = tiempo
            return fila

        resultados = []
        completadas = 0

        if workers <= 1:
            for i, (mid, t, wl, ref) in enumerate(tareas, 1):
                resultados.append(_tarea(mid, t, wl, ref))
                if progress_cb:
                    progress_cb(i, total)
        else:
            with ThreadPoolExecutor(max_workers=workers) as executor:
                futuros = {executor.submit(_tarea, *t): t for t in tareas}
                for fut in as_completed(futuros):
                    resultados.append(fut.result())
                    completadas += 1
                    if progress_cb:
                        progress_cb(completadas, total)

        df_res = pd.DataFrame(resultados)
        cols_first = ["medicion", "tiempo", "lambda_nm", "reflectance_input",
                      "mu_a_mm-1", "mu_s_prime_mm-1", "g"]
        cols_rest = [c for c in df_res.columns if c not in cols_first]
        df_res = df_res[[c for c in cols_first if c in df_res.columns] + cols_rest]
        df_res = df_res.sort_values(["medicion", "lambda_nm"])
        ruta_csv = Path(output_dir) / "resumen_iad_temporal.csv"
        df_res.to_csv(ruta_csv, index=False)
        return df_res, ruta_csv


# ============================================================
# SECCION 5: VISUALIZATION ENGINE
# ============================================================

class VisualizationEngine:

    @staticmethod
    def plot_single_spectrum(fig, ax, wavelengths, values, ylabel="Reflectancia",
                             title="Espectro", color="green"):
        ax.clear()
        ax.plot(wavelengths, values, linewidth=1.5, color=color)
        ax.set_xlabel("Longitud de onda (nm)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()

    @staticmethod
    def plot_iad_results(fig, ax, df_results, param="mu_a_mm-1", ylabel="mu_a (1/mm)"):
        ax.clear()
        col = param if param in df_results.columns else "mu_a_mm-1"
        df_plot = df_results.dropna(subset=["lambda_nm", col]).sort_values("lambda_nm")
        ax.plot(df_plot["lambda_nm"], df_plot[col], marker='o', markersize=4,
                linewidth=1.6, color="steelblue")
        ax.set_xlabel("Longitud de onda (nm)")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{ylabel} reconstruido (IAD)")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()

    @staticmethod
    def create_waterfall_3d(fig, ax, lambdas, tiempos, espectros, ylabel, title,
                            color_base="red", avg_count=1):
        ax.clear()
        n = len(tiempos)
        n_grupos = n // max(1, avg_count)
        if n_grupos == 0:
            return

        espectros_avg = np.array([
            espectros[i * avg_count:(i + 1) * avg_count].mean(axis=0)
            for i in range(n_grupos)
        ])
        tiempos_avg = np.array([
            tiempos[i * avg_count:(i + 1) * avg_count].mean()
            for i in range(n_grupos)
        ])

        if color_base == "red":
            colores = [mcolors.to_rgba((1.0, 0.45 - 0.35 * t, 0.45 - 0.45 * t))
                       for t in np.linspace(0, 1, n_grupos)]
        else:
            colores = [mcolors.to_rgba((0.15, 0.35 + 0.50 * t, 1.0))
                       for t in np.linspace(0, 1, n_grupos)]

        for i in range(n_grupos):
            ax.plot(lambdas, np.full_like(lambdas, tiempos_avg[i]), espectros_avg[i],
                    linewidth=0.8, color=colores[i])

        ax.set_xlabel("Longitud de onda (nm)", labelpad=10)
        ax.set_ylabel("Tiempo (s)", labelpad=10)
        ax.set_zlabel(ylabel, labelpad=10)
        ax.set_title(title, fontsize=12, pad=15)
        ax.grid(True, alpha=0.3)
        ax.view_init(elev=25, azim=-60)
        fig.tight_layout()

    @staticmethod
    def setup_animation_2d(fig, ax, lambdas, tiempos, espectros, ylabel,
                           color_base="red"):
        ax.clear()
        n = len(tiempos)
        z_min = np.nanmin(espectros) * 0.95
        z_max = np.nanmax(espectros) * 1.05

        if color_base == "red":
            colores = [mcolors.to_rgba((1.0, 0.85 - 0.75 * t, 0.85 - 0.85 * t))
                       for t in np.linspace(0, 1, n)]
        else:
            colores = [mcolors.to_rgba((0.15, 0.35 + 0.50 * t, 1.0))
                       for t in np.linspace(0, 1, n)]

        line, = ax.plot(lambdas, espectros[0], linewidth=1.5, color=colores[0])
        ax.set_xlim(lambdas[0], lambdas[-1])
        ax.set_ylim(z_min, z_max)
        ax.set_xlabel("Longitud de onda (nm)")
        ax.set_ylabel(ylabel)
        title_obj = ax.set_title(f"Medicion 1/{n}  -  t = {tiempos[0]:.3f} s")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()

        return {
            "line": line,
            "title": title_obj,
            "colores": colores,
            "n": n,
            "lambdas": lambdas,
            "tiempos": tiempos,
            "espectros": espectros,
        }


# ============================================================
# SECCION 6: APLICACION GUI PRINCIPAL
# ============================================================

class IADApp:
    """Aplicacion principal con 7 pestanas."""

    def __init__(self):
        if HAS_BOOTSTRAP:
            self.root = ttk.Window(themename="darkly")
        else:
            self.root = tk.Tk()

        self.root.title("Pipeline IAD - Propiedades Opticas de Tejido")
        self.root.geometry("1300x920")
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)

        self.spectrometer = SpectrometerDriver()
        self.sp = SignalProcessing()
        self.task_queue = queue.Queue()

        self.R_0 = None
        self.R_1 = None
        self.calibration_wavelengths = None
        self.integration_time = None
        self.calibration_loaded = False

        self.param_vars = {}
        self._logs = {}
        self._pbars = {}
        self._canvases = {}
        self._figs = {}
        self._axes = {}

        # Estado de animaciones (compartido vmrt / viadt)
        self._anim_states = {}

        self._build_status_bar()
        self._build_notebook()
        self._poll_queue()

    def run(self):
        self.root.mainloop()

    # ── Barra de estado superior ──

    def _build_status_bar(self):
        frame = ttk.Frame(self.root)
        frame.pack(side="top", fill="x", padx=8, pady=(6, 2))

        self.hw_indicator = ttk.Label(frame, text="\u25cf", foreground="red",
                                      font=("Segoe UI", 14))
        self.hw_indicator.pack(side="left")

        self.status_label = ttk.Label(frame, text=" Espectrofotometro desconectado",
                                      font=("Segoe UI", 10))
        self.status_label.pack(side="left", padx=(2, 10))

        ttk.Button(frame, text="Conectar", command=self._connect_hw).pack(side="left", padx=2)
        ttk.Button(frame, text="Desconectar", command=self._disconnect_hw).pack(side="left", padx=2)

        ttk.Separator(self.root, orient="horizontal").pack(fill="x", padx=5)

    def _connect_hw(self):
        success, msg = self.spectrometer.connect()
        if success:
            self.hw_indicator.config(foreground="green")
            self.status_label.config(text=f" {msg}")
        else:
            self.hw_indicator.config(foreground="red")
            self.status_label.config(text=f" {msg}")

    def _disconnect_hw(self):
        self.spectrometer.disconnect()
        self.hw_indicator.config(foreground="red")
        self.status_label.config(text=" Espectrofotometro desconectado")

    # ── Notebook ──

    def _build_notebook(self):
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill="both", expand=True, padx=5, pady=5)

        self._build_tab_config()
        self._build_tab_acquisition()
        self._build_tab_viz_mr_single()
        self._build_tab_viz_mr_temporal()
        self._build_tab_iad()
        self._build_tab_viz_iad_single()
        self._build_tab_viz_iad_temporal()

    # ── Cola de tareas (threading seguro) ──

    def _poll_queue(self):
        try:
            while True:
                msg = self.task_queue.get_nowait()
                msg_type = msg[0]
                if msg_type == "progress":
                    _, tab, current, total = msg
                    self._update_progress(tab, current, total)
                elif msg_type == "log":
                    _, tab, text = msg
                    self._append_log(tab, text)
                elif msg_type == "done":
                    _, tab, result = msg
                    self._on_task_done(tab, result)
                elif msg_type == "error":
                    _, tab, err = msg
                    self._on_task_error(tab, err)
                elif msg_type == "plot":
                    _, tab, plot_fn = msg
                    plot_fn()
        except queue.Empty:
            pass
        self.root.after(50, self._poll_queue)

    def _run_in_thread(self, target, *args):
        t = threading.Thread(target=target, args=args, daemon=True)
        t.start()

    # ── Helpers genericos ──

    def _make_param_entry(self, parent, label, key, default, row, col=0, width=12):
        lbl = ttk.Label(parent, text=label)
        lbl.grid(row=row, column=col * 2, sticky="e", padx=(6, 2), pady=2)
        var = tk.StringVar(value=str(default))
        entry = ttk.Entry(parent, textvariable=var, width=width)
        entry.grid(row=row, column=col * 2 + 1, sticky="w", padx=(2, 6), pady=2)
        self.param_vars[key] = var
        return var

    def _get_float(self, key):
        try:
            return float(self.param_vars[key].get())
        except (ValueError, KeyError):
            return 0.0

    def _get_int(self, key):
        try:
            return int(float(self.param_vars[key].get()))
        except (ValueError, KeyError):
            return 0

    def _make_log(self, parent):
        log = ScrolledText(parent, height=8)
        log.pack(fill="both", expand=True, padx=5, pady=(2, 5))
        return log

    def _make_progress(self, parent):
        pbar = ttk.Progressbar(parent, mode="determinate", length=400)
        pbar.pack(fill="x", padx=5, pady=2)
        return pbar

    def _make_canvas(self, parent, figsize=(9, 4.5)):
        fig, ax = plt.subplots(figsize=figsize)
        canvas = FigureCanvasTkAgg(fig, master=parent)
        canvas.get_tk_widget().pack(fill="both", expand=True, padx=2, pady=2)
        toolbar = NavigationToolbar2Tk(canvas, parent)
        toolbar.update()
        return fig, ax, canvas

    def _append_log(self, tab, text):
        log = self._logs.get(tab)
        if log:
            log.insert("end", text + "\n")
            log.see("end")

    def _update_progress(self, tab, current, total):
        pbar = self._pbars.get(tab)
        if pbar and total > 0:
            pbar["value"] = (current / total) * 100

    def _on_task_done(self, tab, result):
        self._update_progress(tab, 100, 100)
        if tab == "acq":
            self._on_acq_task_done(result)
        elif tab == "iad":
            self._on_iad_task_done(result)

    def _on_task_error(self, tab, err):
        self._append_log(tab, f"ERROR: {err}")
        messagebox.showerror("Error", str(err))

    def _browse_dir(self, key):
        d = filedialog.askdirectory()
        if d:
            self.param_vars[key].set(d)

    def _browse_file(self, key, filetypes=None):
        if filetypes is None:
            filetypes = [("Todos", "*.*")]
        f = filedialog.askopenfilename(filetypes=filetypes)
        if f:
            self.param_vars[key].set(f)

    def _make_scrollable(self, parent):
        """Crea un frame scrollable dentro de parent. Retorna el frame interior."""
        canvas = tk.Canvas(parent, highlightthickness=0, borderwidth=0)
        vscroll = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        inner = ttk.Frame(canvas)

        inner.bind("<Configure>",
                   lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=inner, anchor="nw")
        canvas.configure(yscrollcommand=vscroll.set)

        vscroll.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)

        def _on_enter(event):
            canvas.bind_all("<MouseWheel>",
                            lambda e: canvas.yview_scroll(int(-1 * (e.delta / 120)), "units"))
        def _on_leave(event):
            canvas.unbind_all("<MouseWheel>")

        canvas.bind("<Enter>", _on_enter)
        canvas.bind("<Leave>", _on_leave)
        inner.bind("<Enter>", _on_enter)
        inner.bind("<Leave>", _on_leave)
        return inner

    def _add_confirm_button(self, parent):
        frame = ttk.Frame(parent)
        existing_children = parent.winfo_children()
        uses_grid = any(child.winfo_manager() == "grid" for child in existing_children)
        if uses_grid:
            grid_rows = [
                int(child.grid_info().get("row", 0))
                for child in existing_children
                if child.winfo_manager() == "grid"
            ]
            next_row = (max(grid_rows) + 1) if grid_rows else 0
            frame.grid(row=next_row, column=0, columnspan=4, sticky="ew",
                       padx=6, pady=(2, 4))
        else:
            frame.pack(fill="x", padx=6, pady=(2, 4))
        status_lbl = ttk.Label(frame, text="", foreground="green",
                               font=("Segoe UI", 8, "italic"))
        status_lbl.pack(side="left", padx=4)

        def _confirm():
            status_lbl.config(text="Cambios confirmados")
            frame.after(3000, lambda: status_lbl.config(text=""))

        ttk.Button(frame, text="Confirmar cambios",
                   command=_confirm).pack(side="right", padx=4)

    def _resolve_active_prefix(self):
        if hasattr(self, "acq_mode_var") and self.acq_mode_var.get() == "temporal":
            return "t"
        return "s"

    def _preview_subject_paths(self, temporal=False):
        base_raw = self.param_vars.get(
            "output_dir", tk.StringVar(value=str(MEDICIONES_DIR))
        ).get().strip()
        base = Path(base_raw) if base_raw else MEDICIONES_DIR
        suffix = "_temporal" if temporal else ""
        patron = re.compile(rf"^sujeto_(\d{{3}})_\d{{2}}_\d{{2}}_\d{{4}}{suffix}$")
        ids_existentes = []
        if base.exists():
            for nombre in os.listdir(base):
                m = patron.match(nombre)
                if m and (base / nombre).is_dir():
                    ids_existentes.append(int(m.group(1)))
        sig_id = max(ids_existentes) + 1 if ids_existentes else 1
        fecha = datetime.now().strftime("%d_%m_%Y")
        nombre = f"sujeto_{sig_id:03d}_{fecha}{suffix}"
        ruta_sujeto = base / nombre
        return ruta_sujeto, ruta_sujeto / "series", sig_id

    def _update_save_preview(self):
        if not hasattr(self, "acq_savepath_label"):
            return
        temporal = hasattr(self, "acq_mode_var") and self.acq_mode_var.get() == "temporal"
        _, ruta_series, sig_id = self._preview_subject_paths(temporal=temporal)
        archivo = "M_R_tiempo_data.csv" if temporal else "M_R_data.csv"
        self.acq_savepath_label.config(
            text=f"Sujeto siguiente: {sig_id:03d}\n{ruta_series / archivo}"
        )

    def _set_calibration_status(self, message, foreground="green"):
        detalles = []
        if self.integration_time:
            detalles.append(f"TI: {self.integration_time * 1000:.1f} ms")
        if self.R_0 is not None:
            detalles.append(f"R_0: {self.R_0.min():.4f}-{self.R_0.max():.4f}")
        if self.R_1 is not None:
            detalles.append(f"R_1: {self.R_1.min():.4f}-{self.R_1.max():.4f}")

        texto = message
        if detalles:
            texto += "  |  " + "  |  ".join(detalles)

        self.cal_status_label.config(text=texto, foreground=foreground)
        if hasattr(self, "acq_calibration_label"):
            self.acq_calibration_label.config(text=texto, foreground=foreground)

    def _refresh_calibration_plot(self):
        if "cfg_cal" not in self._axes:
            return

        fig = self._figs["cfg_cal"]
        ax = self._axes["cfg_cal"]
        canvas = self._canvases["cfg_cal"]
        ax.clear()
        ax.set_xlabel("Longitud de onda (nm)", fontsize=7)
        ax.set_ylabel("Intensidad", fontsize=7)
        ax.set_title("R_0 y R_1", fontsize=8)
        ax.tick_params(labelsize=6)
        ax.grid(True, alpha=0.3)

        wavelengths = self.calibration_wavelengths
        if wavelengths is None and self.spectrometer.wavelengths is not None:
            wavelengths = self.spectrometer.wavelengths

        if wavelengths is None or (self.R_0 is None and self.R_1 is None):
            ax.text(0.5, 0.5, "Sin calibracion cargada",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=8, color="gray")
            fig.tight_layout()
            canvas.draw_idle()
            return

        prefix = self._resolve_active_prefix()
        lmin = self._get_float(f"{prefix}_lmin")
        lmax = self._get_float(f"{prefix}_lmax")
        if lmax <= lmin:
            lmin, lmax = float(wavelengths[0]), float(wavelengths[-1])

        def _plot_curve(data, color, label):
            if data is None:
                return False
            n = min(len(wavelengths), len(data))
            wl = np.asarray(wavelengths[:n])
            vals = np.asarray(data[:n])
            mask = (wl >= lmin) & (wl <= lmax)
            if not np.any(mask):
                mask = np.ones_like(wl, dtype=bool)
            ax.plot(wl[mask], vals[mask], linewidth=1.2, color=color, label=label)
            return True

        drew = _plot_curve(self.R_0, "gray", "R_0")
        drew = _plot_curve(self.R_1, "orange", "R_1") or drew
        if drew:
            ax.legend(fontsize=7, loc="best")
        else:
            ax.text(0.5, 0.5, "Sin datos para graficar",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=8, color="gray")

        fig.tight_layout()
        canvas.draw_idle()
        if hasattr(self, "_acq_state"):
            self._update_acq_status()

    def _check_hw(self, tab="acq"):
        if not self.spectrometer.connected:
            messagebox.showwarning("Hardware",
                                   "Espectrofotometro no conectado.\n"
                                   "Usa el boton 'Conectar' en la barra superior.")
            return False
        return True

    # ================================================================
    # TAB 1: CONFIGURACION
    # ================================================================

    def _build_tab_config(self):
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="1. Configuracion")

        body = ttk.Panedwindow(tab, orient="horizontal")
        body.pack(fill="both", expand=True, padx=6, pady=6)

        # ── PANEL IZQUIERDO scrollable ──
        left_outer = ttk.Frame(body, width=520)
        body.add(left_outer, weight=1)
        left = self._make_scrollable(left_outer)

        # Datos del sujeto
        lf_suj = ttk.LabelFrame(left, text="  Datos del Sujeto  ")
        lf_suj.pack(fill="x", pady=(0, 4), padx=4)
        campos = [
            ("Nombre:", "suj_nombre", 30),
            ("Edad:", "suj_edad", 10),
            ("Municipio:", "suj_municipio", 22),
            ("Fitzpatrick (I-VI):", "suj_fitz", 8),
        ]
        for i, (lbl, key, w) in enumerate(campos):
            ttk.Label(lf_suj, text=lbl, anchor="e").grid(
                row=i, column=0, sticky="e", padx=(8, 4), pady=3)
            var = tk.StringVar()
            ttk.Entry(lf_suj, textvariable=var, width=w).grid(
                row=i, column=1, sticky="ew", padx=(0, 8), pady=3)
            self.param_vars[key] = var
        lf_suj.columnconfigure(1, weight=1)
        ttk.Label(lf_suj, text="Afecciones:", anchor="e").grid(
            row=4, column=0, sticky="ne", padx=(8, 4), pady=3)
        var_af = tk.StringVar()
        ttk.Entry(lf_suj, textvariable=var_af, width=30).grid(
            row=4, column=1, sticky="ew", padx=(0, 8), pady=3)
        self.param_vars["suj_afecciones"] = var_af
        self._add_confirm_button(lf_suj)

        # Carpeta de salida
        lf_out = ttk.LabelFrame(left, text="  Carpeta de Salida  ")
        lf_out.pack(fill="x", pady=(0, 4), padx=4)
        self.param_vars["output_dir"] = tk.StringVar(value=str(MEDICIONES_DIR))
        ttk.Entry(lf_out, textvariable=self.param_vars["output_dir"],
                  width=30).pack(side="left", padx=(8, 4), pady=6)
        ttk.Button(lf_out, text="...",
                   command=lambda: self._browse_dir("output_dir")).pack(side="left", pady=6)
        self.param_vars["output_dir"].trace_add("write", lambda *_: self._update_save_preview())
        self._add_confirm_button(lf_out)

        # ── BOTON GRANDE ──
        ttk.Separator(left, orient="horizontal").pack(fill="x", padx=4, pady=4)
        self.cfg_folder_status = ttk.Label(
            left, text="", foreground="gray",
            font=("Segoe UI", 8, "italic"), wraplength=390, justify="left")
        self.cfg_folder_status.pack(anchor="w", padx=6, pady=(0, 2))
        ttk.Button(left, text="Guardar y Crear Carpeta del Sujeto",
                   command=self._config_create_subject_folder).pack(
            fill="x", ipady=10, padx=6, pady=(0, 6))

        # ── CALIBRACION ──
        ttk.Separator(left, orient="horizontal").pack(fill="x", padx=4, pady=2)
        lf_cal = ttk.LabelFrame(left, text="  Calibracion (R_0 y R_1)  ")
        lf_cal.pack(fill="x", pady=(2, 4), padx=4)

        self.cal_status_label = ttk.Label(
            lf_cal, text="Sin calibracion cargada",
            foreground="gray", font=("Segoe UI", 8, "italic"), wraplength=390)
        self.cal_status_label.pack(anchor="w", padx=8, pady=(6, 2))

        meas_frame = ttk.Frame(lf_cal)
        meas_frame.pack(fill="x", padx=6, pady=(2, 0))
        ttk.Button(meas_frame, text="Optimizar TI",
                   command=self._cfg_optimize_ti).pack(fill="x", pady=1)
        ttk.Button(meas_frame, text="Medir R_0 (oscuridad)",
                   command=lambda: self._cfg_measure_ref("r0")).pack(fill="x", pady=1)
        ttk.Button(meas_frame, text="Medir R_1 (referencia blanca)",
                   command=lambda: self._cfg_measure_ref("r1")).pack(fill="x", pady=1)

        ttk.Separator(lf_cal, orient="horizontal").pack(fill="x", padx=6, pady=4)
        load_frame = ttk.Frame(lf_cal)
        load_frame.pack(fill="x", padx=6, pady=(0, 2))
        ttk.Button(load_frame, text="Cargar ultima guardada",
                   command=self._load_calibration_from_disk).pack(fill="x", pady=1)
        ttk.Button(load_frame, text="Cargar desde archivo NPZ...",
                   command=self._load_calibration_from_file).pack(fill="x", pady=1)
        self._add_confirm_button(lf_cal)

        # Canvas R_0 / R_1
        fig_cal, ax_cal = plt.subplots(figsize=(4.2, 2.0))
        canvas_cal = FigureCanvasTkAgg(fig_cal, master=lf_cal)
        canvas_cal.get_tk_widget().pack(fill="both", expand=True, padx=4, pady=(4, 2))
        ax_cal.set_xlabel("Longitud de onda (nm)", fontsize=7)
        ax_cal.set_ylabel("Intensidad", fontsize=7)
        ax_cal.set_title("R_0 y R_1", fontsize=8)
        ax_cal.tick_params(labelsize=6)
        ax_cal.grid(True, alpha=0.3)
        fig_cal.tight_layout()
        canvas_cal.draw()
        self._figs["cfg_cal"] = fig_cal
        self._axes["cfg_cal"] = ax_cal
        self._canvases["cfg_cal"] = canvas_cal

        self._pbars["cfg"] = ttk.Progressbar(lf_cal, mode="determinate")
        self._pbars["cfg"].pack(fill="x", padx=6, pady=2)
        self._logs["cfg"] = ScrolledText(lf_cal, height=4)
        self._logs["cfg"].pack(fill="x", padx=6, pady=(0, 6))

        # ── PANEL DERECHO scrollable: parametros tecnicos ──
        right_outer = ttk.Frame(body)
        body.add(right_outer, weight=1)
        inner = self._make_scrollable(right_outer)

        # Parametros adquisicion individual
        lf = ttk.LabelFrame(inner, text="  Adquisicion Individual  ")
        lf.pack(fill="x", padx=6, pady=(4, 4))
        ds = DEFAULTS["single"]
        self._make_param_entry(lf, "T. integracion inicial (s):", "s_ti", ds["tiempo_integracion_inicial"], 0)
        self._make_param_entry(lf, "Mediciones promedio:", "s_nmed", ds["num_mediciones_promedio"], 0, col=1)
        self._make_param_entry(lf, "T. espera entre med. (s):", "s_tespera", ds["tiempo_espera"], 1)
        self._make_param_entry(lf, "% Saturacion:", "s_sat", ds["porcentaje_saturacion"], 1, col=1)
        self._make_param_entry(lf, "Lambda min (nm):", "s_lmin", ds["lambda_min"], 2)
        self._make_param_entry(lf, "Lambda max (nm):", "s_lmax", ds["lambda_max"], 2, col=1)
        self._make_param_entry(lf, "Orden Butterworth:", "s_border", ds["orden_filtro_butter"], 3)
        self._make_param_entry(lf, "Frecuencia corte:", "s_bfc", ds["frecuencia_corte_butter"], 3, col=1)
        self._make_param_entry(lf, "Muestras objetivo:", "s_mobj", ds["muestras_objetivo"], 4)
        self._make_param_entry(lf, "R_std:", "s_rstd", ds["r_std"], 4, col=1)
        self._add_confirm_button(lf)

        # Parametros adquisicion temporal
        lf = ttk.LabelFrame(inner, text="  Adquisicion Temporal  ")
        lf.pack(fill="x", padx=6, pady=4)
        dt = DEFAULTS["temporal"]
        self._make_param_entry(lf, "T. integracion inicial (s):", "t_ti", dt["tiempo_integracion_inicial"], 0)
        self._make_param_entry(lf, "Mediciones promedio:", "t_nmed", dt["num_mediciones_promedio"], 0, col=1)
        self._make_param_entry(lf, "T. espera entre med. (s):", "t_tespera", dt["tiempo_espera"], 1)
        self._make_param_entry(lf, "% Saturacion:", "t_sat", dt["porcentaje_saturacion"], 1, col=1)
        self._make_param_entry(lf, "Lambda min (nm):", "t_lmin", dt["lambda_min"], 2)
        self._make_param_entry(lf, "Lambda max (nm):", "t_lmax", dt["lambda_max"], 2, col=1)
        self._make_param_entry(lf, "Orden Butterworth:", "t_border", dt["orden_filtro_butter"], 3)
        self._make_param_entry(lf, "Frecuencia corte:", "t_bfc", dt["frecuencia_corte_butter"], 3, col=1)
        self._make_param_entry(lf, "Muestras objetivo:", "t_mobj", dt["muestras_objetivo"], 4)
        self._make_param_entry(lf, "R_std:", "t_rstd", dt["r_std"], 4, col=1)
        self._make_param_entry(lf, "Num. espectros temporales:", "t_nesp", dt["num_espectros_temporal"], 5)
        self._add_confirm_button(lf)

        # Parametros procesamiento IAD
        lf = ttk.LabelFrame(inner, text="  Procesamiento IAD  ")
        lf.pack(fill="x", padx=6, pady=(4, 8))

        def _file_row(parent, label, key, default, filetypes=None, is_dir=False):
            r = ttk.Frame(parent)
            r.pack(fill="x", padx=6, pady=2)
            ttk.Label(r, text=label, width=16, anchor="e").pack(side="left")
            var = tk.StringVar(value=str(default))
            ttk.Entry(r, textvariable=var, width=42).pack(side="left", padx=4)
            if is_dir:
                ttk.Button(r, text="...", width=3,
                           command=lambda k=key: self._browse_dir(k)).pack(side="left")
            else:
                ttk.Button(r, text="...", width=3,
                           command=lambda k=key, ft=filetypes: self._browse_file(k, ft)).pack(side="left")
            self.param_vars[key] = var

        _file_row(lf, "iad.exe:", "iad_exe", IAD_EXE_DEFAULT, [("EXE", "*.exe")])
        _file_row(lf, "Plantilla RXT:", "iad_rxt", RXT_TEMPLATE_DEFAULT, [("RXT", "*.rxt")])
        _file_row(lf, "Carpeta salida:", "iad_outdir", IAD_OUTPUT_DIR_DEFAULT, is_dir=True)

        di = DEFAULTS["iad"]
        pg = ttk.Frame(lf)
        pg.pack(fill="x", padx=6, pady=4)
        self.iad_gfijo_var = tk.BooleanVar(value=di["usar_g_fijo"])
        ttk.Checkbutton(pg, text="Usar g fijo",
                        variable=self.iad_gfijo_var).grid(row=0, column=0, sticky="w", padx=4, pady=2)
        self._make_param_entry(pg, "g fijo:", "iad_g", di["g_fijo"], 0, col=1)
        self.iad_fast_var = tk.BooleanVar(value=di["usar_modo_rapido"])
        ttk.Checkbutton(pg, text="Modo rapido (-M 0)",
                        variable=self.iad_fast_var).grid(row=1, column=0, sticky="w", padx=4, pady=2)
        self._make_param_entry(pg, "Workers:", "iad_workers", di["workers"], 1, col=1)
        self._make_param_entry(pg, "Max mediciones:", "iad_maxmed", di["max_mediciones"], 2)
        self._make_param_entry(pg, "Paso medicion:", "iad_paso", di["paso_medicion"], 2, col=1)
        self._add_confirm_button(lf)

    def _config_create_subject_folder(self):
        """Crea la carpeta del sujeto con sus datos y muestra confirmacion."""
        nombre = self.param_vars.get("suj_nombre", tk.StringVar()).get().strip()
        base = Path(self.param_vars["output_dir"].get())

        # Determinar siguiente ID
        patron = re.compile(r"^sujeto_(\d{3})_\d{2}_\d{2}_\d{4}$")
        ids_existentes = []
        if base.exists():
            for n in os.listdir(base):
                m = patron.match(n)
                if m and (base / n).is_dir():
                    ids_existentes.append(int(m.group(1)))
        sig_id = max(ids_existentes) + 1 if ids_existentes else 1
        fecha = datetime.now().strftime("%d_%m_%Y")
        nombre_carpeta = f"sujeto_{sig_id:03d}_{fecha}"
        ruta_sujeto = base / nombre_carpeta

        try:
            ruta_series = ruta_sujeto / "series"
            ruta_series.mkdir(parents=True, exist_ok=True)

            datos = [
                f"sujeto_id: {sig_id:03d}",
                f"fecha_registro: {datetime.now().strftime('%d/%m/%Y %H:%M:%S')}",
                f"nombre: {nombre or 'No especificado'}",
                f"edad: {self.param_vars.get('suj_edad', tk.StringVar()).get() or 'No especificado'}",
                f"municipio: {self.param_vars.get('suj_municipio', tk.StringVar()).get() or 'No especificado'}",
                f"fitzpatrick: {self.param_vars.get('suj_fitz', tk.StringVar()).get() or 'No especificado'}",
                f"afecciones: {self.param_vars.get('suj_afecciones', tk.StringVar()).get() or 'No especificado'}",
            ]
            (ruta_sujeto / "datos_sujeto.txt").write_text("\n".join(datos), encoding="utf-8")

            self.cfg_folder_status.config(
                text=f"Carpeta creada: {nombre_carpeta}", foreground="green")
            self._update_save_preview()

            messagebox.showinfo(
                "Carpeta del Sujeto Creada",
                f"ID del sujeto: {sig_id:03d}\n"
                f"Carpeta creada en:\n{ruta_sujeto}\n\n"
                f"Puedes continuar con la adquisicion."
            )
        except Exception as e:
            self.cfg_folder_status.config(
                text=f"Error al crear carpeta: {e}", foreground="red")
            messagebox.showerror("Error", str(e))

    # ── Calibracion ──

    def _load_calibration_from_disk(self):
        npz_path = CALIBRACIONES_DIR / "calibracion_actual.npz"
        meta_path = CALIBRACIONES_DIR / "calibracion_actual.json"
        if not npz_path.exists() or not meta_path.exists():
            messagebox.showinfo("Info", "No se encontro calibracion guardada en disco.\n"
                                f"Buscando en: {CALIBRACIONES_DIR}")
            return
        try:
            with open(meta_path, "r", encoding="utf-8") as f:
                meta = json.load(f)
            with np.load(npz_path) as data:
                self.calibration_wavelengths = data["wavelengths"] if "wavelengths" in data else None
                self.R_0 = data["R_0"]
                self.R_1 = data["R_1"]
            fecha = meta.get("fecha_calibracion", "desconocida")
            ti = float(meta.get("tiempo_integracion_s", 0.1))
            self.integration_time = ti
            self.calibration_loaded = True
            self._set_calibration_status(f"Calibracion cargada desde disco ({fecha})")
            self._refresh_calibration_plot()
        except Exception as e:
            messagebox.showerror("Error", f"Error cargando calibracion: {e}")

    def _load_calibration_from_file(self):
        f = filedialog.askopenfilename(
            filetypes=[("NPZ", "*.npz")],
            initialdir=str(CALIBRACIONES_DIR) if CALIBRACIONES_DIR.exists() else str(BASE_DIR))
        if not f:
            return
        try:
            with np.load(f) as data:
                self.calibration_wavelengths = data["wavelengths"] if "wavelengths" in data else None
                self.R_0 = data["R_0"]
                self.R_1 = data["R_1"]
            self.calibration_loaded = True
            self._set_calibration_status("Calibracion cargada desde archivo")
            self._refresh_calibration_plot()
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def _cfg_optimize_ti(self):
        if not self._check_hw("cfg"):
            return

        prefix = self._resolve_active_prefix()
        self._append_log("cfg", "Optimizando tiempo de integracion para calibracion...")
        self._update_progress("cfg", 0, 100)

        def worker():
            try:
                ti_init = self._get_float(f"{prefix}_ti")
                umbral = self._get_float(f"{prefix}_sat")
                ti_opt, _ = self.spectrometer.optimize_integration_time(
                    ti_init, umbral,
                    progress_cb=lambda c, t: self.task_queue.put(("progress", "cfg", c, t)),
                    log_cb=lambda msg: self.task_queue.put(("log", "cfg", msg)),
                )
                self.integration_time = ti_opt
                if hasattr(self, "_acq_state"):
                    self._acq_state["ti"] = ti_opt
                    self._acq_state["optimized"] = True
                self.task_queue.put(("log", "cfg", f"TI optimo para calibracion: {ti_opt * 1000:.2f} ms"))
                self.task_queue.put(("plot", "cfg",
                                     lambda: self._set_calibration_status(
                                         "Tiempo de integracion optimizado",
                                         foreground="green")))
            except Exception as e:
                self.task_queue.put(("error", "cfg", str(e)))

        self._run_in_thread(worker)

    def _cfg_measure_ref(self, ref_type):
        if not self._check_hw("cfg"):
            return
        if ref_type not in {"r0", "r1"}:
            return

        prefix = self._resolve_active_prefix()
        label = "R_0 (oscuridad)" if ref_type == "r0" else "R_1 (referencia blanca)"
        self._append_log("cfg", f"Midiendo {label}...")
        self._update_progress("cfg", 0, 100)

        def worker():
            try:
                ti = self.integration_time or self._get_float(f"{prefix}_ti")
                data = self.spectrometer.take_averaged_series(
                    ti,
                    self._get_int(f"{prefix}_nmed"),
                    self._get_float(f"{prefix}_tespera"),
                    progress_cb=lambda c, t: self.task_queue.put(("progress", "cfg", c, t)),
                    log_cb=lambda msg: self.task_queue.put(("log", "cfg", msg)),
                )
                self.integration_time = ti
                self.calibration_wavelengths = np.array(self.spectrometer.wavelengths, copy=True)

                if ref_type == "r0":
                    self.R_0 = data
                    status_msg = "R_0 medido"
                else:
                    self.R_1 = data
                    status_msg = "R_1 medido"
                if self.R_0 is not None and self.R_1 is not None:
                    self._save_calibration(ti, prefix, log_tab="cfg")

                self.calibration_loaded = self.R_0 is not None and self.R_1 is not None
                self.task_queue.put(("log", "cfg",
                    f"{label} completado. Rango: {data.min():.6f} - {data.max():.6f}"))
                self.task_queue.put(("plot", "cfg", lambda: (
                    self._set_calibration_status(
                        "Calibracion lista" if self.calibration_loaded else status_msg,
                        foreground="green" if self.calibration_loaded else "orange"
                    ),
                    self._refresh_calibration_plot()
                )))
            except Exception as e:
                self.task_queue.put(("error", "cfg", str(e)))

        self._run_in_thread(worker)

    # ================================================================
    # TAB 2: ADQUISICION (Individual + Temporal)
    # ================================================================

    def _build_tab_acquisition(self):
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="2. Adquisicion")

        # Selector de modo arriba
        top_bar = ttk.Frame(tab)
        top_bar.pack(fill="x", padx=8, pady=(6, 2))

        ttk.Label(top_bar, text="Modo:", font=("Segoe UI", 10, "bold")).pack(side="left", padx=(0, 8))
        self.acq_mode_var = tk.StringVar(value="individual")
        ttk.Radiobutton(top_bar, text="Espectro Individual",
                        variable=self.acq_mode_var, value="individual",
                        command=self._acq_mode_changed).pack(side="left", padx=4)
        ttk.Radiobutton(top_bar, text="Serie Temporal",
                        variable=self.acq_mode_var, value="temporal",
                        command=self._acq_mode_changed).pack(side="left", padx=4)

        ttk.Separator(tab, orient="horizontal").pack(fill="x", padx=5, pady=2)

        body = ttk.Frame(tab)
        body.pack(fill="both", expand=True)

        # Panel izquierdo
        left = ttk.Frame(body, width=320)
        left.pack(side="left", fill="y", padx=(8, 4), pady=4)
        left.pack_propagate(False)

        self.acq_instructions = ttk.Label(
            left, text="", wraplength=300, justify="left",
            font=("Segoe UI", 9))
        self.acq_instructions.pack(anchor="w", padx=4, pady=(4, 8))

        lf_cal_info = ttk.LabelFrame(left, text="  Calibracion Activa  ")
        lf_cal_info.pack(fill="x", padx=4, pady=2)
        self.acq_calibration_label = ttk.Label(
            lf_cal_info,
            text="Gestiona R_0 y R_1 en la pestana Configuracion.",
            wraplength=285, justify="left", foreground="gray",
            font=("Segoe UI", 8, "italic")
        )
        self.acq_calibration_label.pack(anchor="w", padx=8, pady=6)

        lf_dest = ttk.LabelFrame(left, text="  Destino del Guardado  ")
        lf_dest.pack(fill="x", padx=4, pady=2)
        self.acq_savepath_label = ttk.Label(
            lf_dest, text="", wraplength=285, justify="left",
            font=("Consolas", 8)
        )
        self.acq_savepath_label.pack(anchor="w", padx=8, pady=6)

        lf_flow = ttk.LabelFrame(left, text="  Flujo de Trabajo  ")
        lf_flow.pack(fill="x", padx=4, pady=2)

        self.acq_btns = {}
        steps = [
            ("rm",  "1. Medir R_M / Serie Temporal"),
            ("calc", "2. Calcular M_R y Guardar"),
        ]
        for key, text in steps:
            btn = ttk.Button(lf_flow, text=text, command=lambda k=key: self._acq_step(k))
            btn.pack(fill="x", padx=6, pady=3)
            self.acq_btns[key] = btn

        self.acq_step_status = ttk.Label(left, text="", justify="left",
                                         font=("Consolas", 9))
        self.acq_step_status.pack(anchor="w", padx=8, pady=(8, 2))

        # Panel derecho
        right = ttk.Frame(body)
        right.pack(side="left", fill="both", expand=True, padx=(4, 8), pady=4)

        self._figs["acq"], self._axes["acq"], self._canvases["acq"] = self._make_canvas(right)
        self._pbars["acq"] = self._make_progress(right)
        self._logs["acq"] = self._make_log(right)

        self._acq_state = {
            "optimized": False,
            "r0_done": False,
            "r1_done": False,
            "rm_done": False,
            "R_0": None,
            "R_1": None,
            "R_M": None,
            "ti": None,
            "intensity_opt": None,
            "temporal_data": None,
        }

        self._acq_mode_changed()

    def _acq_mode_changed(self):
        mode = self.acq_mode_var.get()
        if hasattr(self, "_acq_state"):
            self._acq_state["rm_done"] = False
            self._acq_state["R_M"] = None
            self._acq_state["temporal_data"] = None
        if mode == "individual":
            self.acq_instructions.config(
                text="Adquisicion de un espectro individual M_R.\n\n"
                     "1) Verifica en Configuracion que la calibracion este lista\n"
                     "2) Mide R_M con la muestra\n"
                     "3) Calcula M_R y guarda los datos")
            self.acq_btns["rm"].config(text="1. Medir R_M (muestra)")
            self.acq_btns["calc"].config(text="2. Calcular M_R y Guardar")
        else:
            self.acq_instructions.config(
                text="Adquisicion de serie temporal de espectros.\n\n"
                     "1) Verifica en Configuracion que la calibracion este lista\n"
                     "2) Adquiere la serie temporal de R_M\n"
                     "3) Procesa y guarda la serie")
            self.acq_btns["rm"].config(text="1. Adquirir Serie Temporal R_M")
            self.acq_btns["calc"].config(text="2. Procesar y Guardar Serie")
        self._update_save_preview()
        self._refresh_calibration_plot()
        self._update_acq_status()

    def _update_acq_status(self):
        s = self._acq_state
        if self.calibration_loaded and not s["r0_done"]:
            s["R_0"] = self.R_0
            s["R_1"] = self.R_1
            s["r0_done"] = True
            s["r1_done"] = True
            if self.integration_time:
                s["ti"] = self.integration_time
                s["optimized"] = True

        rm_label = "Serie temporal adquirida" if self.acq_mode_var.get() == "temporal" else "R_M medido"
        checks = [
            ("TI de calibracion listo" if s["optimized"] else "TI de calibracion pendiente", s["optimized"]),
            ("Calibracion cargada" if (s["r0_done"] and s["r1_done"]) else "Calibracion pendiente",
             s["r0_done"] and s["r1_done"]),
            (rm_label if s["rm_done"] else ("Serie temporal pendiente" if self.acq_mode_var.get() == "temporal"
                                            else "R_M pendiente"), s["rm_done"]),
        ]
        lines = []
        for text, done in checks:
            marker = "[OK]" if done else "[  ]"
            lines.append(f"  {marker} {text}")
        self.acq_step_status.config(text="\n".join(lines))
        self._update_save_preview()

    def _acq_step(self, step):
        if not self._check_hw("acq"):
            return

        s = self._acq_state
        mode = self.acq_mode_var.get()
        prefix = "s" if mode == "individual" else "t"

        if step == "opt":
            self._acq_optimize(prefix)
        elif step == "r0":
            self._acq_measure_ref("r0", prefix)
        elif step == "r1":
            if not s["r0_done"]:
                messagebox.showwarning("Flujo", "Primero mide R_0.")
                return
            self._acq_measure_ref("r1", prefix)
        elif step == "rm":
            if not (s["r0_done"] and s["r1_done"]):
                messagebox.showwarning("Flujo", "Primero carga o mide la calibracion en la pestana Configuracion.")
                return
            if mode == "individual":
                self._acq_measure_ref("rm", prefix)
            else:
                self._acq_temporal_series(prefix)
        elif step == "calc":
            if mode == "individual":
                self._acq_calculate_mr_single(prefix)
            else:
                self._acq_calculate_mr_temporal(prefix)

    def _acq_optimize(self, prefix):
        self._append_log("acq", "Optimizando tiempo de integracion...")

        def worker():
            try:
                ti_init = self._get_float(f"{prefix}_ti")
                umbral = 1.0 * self._get_float(f"{prefix}_sat")
                ti_opt, int_opt = self.spectrometer.optimize_integration_time(
                    ti_init, umbral,
                    progress_cb=lambda c, t: self.task_queue.put(("progress", "acq", c, t)),
                    log_cb=lambda msg: self.task_queue.put(("log", "acq", msg)),
                )
                self._acq_state["ti"] = ti_opt
                self._acq_state["intensity_opt"] = int_opt
                self._acq_state["optimized"] = True
                self.integration_time = ti_opt
                self.task_queue.put(("log", "acq", f"TI optimo: {ti_opt*1000:.2f} ms"))

                wl = self.spectrometer.wavelengths
                lmin = self._get_float(f"{prefix}_lmin")
                lmax = self._get_float(f"{prefix}_lmax")
                mask = (wl >= lmin) & (wl <= lmax)

                def do_plot():
                    VisualizationEngine.plot_single_spectrum(
                        self._figs["acq"], self._axes["acq"],
                        wl[mask], int_opt[mask],
                        ylabel="Intensidad",
                        title=f"Espectro optimizado (TI={ti_opt*1000:.2f} ms)",
                        color="dodgerblue"
                    )
                    self._canvases["acq"].draw()
                    self._update_acq_status()

                self.task_queue.put(("plot", "acq", do_plot))
            except Exception as e:
                self.task_queue.put(("error", "acq", str(e)))

        self._run_in_thread(worker)

    def _acq_measure_ref(self, ref_type, prefix):
        labels = {"r0": "R_0 (oscuridad)", "r1": "R_1 (referencia)", "rm": "R_M (muestra)"}
        self._append_log("acq", f"Midiendo {labels[ref_type]}...")

        def worker():
            try:
                ti = self._acq_state.get("ti") or self._get_float(f"{prefix}_ti")
                data = self.spectrometer.take_averaged_series(
                    ti, self._get_int(f"{prefix}_nmed"),
                    self._get_float(f"{prefix}_tespera"),
                    progress_cb=lambda c, t: self.task_queue.put(("progress", "acq", c, t)),
                    log_cb=lambda msg: self.task_queue.put(("log", "acq", msg)),
                )
                key_map = {"r0": "R_0", "r1": "R_1", "rm": "R_M"}
                self._acq_state[key_map[ref_type]] = data
                self._acq_state[f"{ref_type}_done"] = True

                if ref_type == "r0":
                    self.R_0 = data
                    self.calibration_wavelengths = np.array(self.spectrometer.wavelengths, copy=True)
                elif ref_type == "r1":
                    self.R_1 = data
                    self.calibration_wavelengths = np.array(self.spectrometer.wavelengths, copy=True)
                    self._save_calibration(ti, prefix)

                self.task_queue.put(("log", "acq",
                    f"{labels[ref_type]} completado. Rango: {data.min():.6f} - {data.max():.6f}"))

                wl = self.spectrometer.wavelengths
                lmin = self._get_float(f"{prefix}_lmin")
                lmax = self._get_float(f"{prefix}_lmax")
                mask = (wl >= lmin) & (wl <= lmax)

                def do_plot():
                    colors = {"r0": "gray", "r1": "orange", "rm": "green"}
                    VisualizationEngine.plot_single_spectrum(
                        self._figs["acq"], self._axes["acq"],
                        wl[mask], data[mask],
                        ylabel="Intensidad", title=f"Espectro {labels[ref_type]}",
                        color=colors[ref_type]
                    )
                    self._canvases["acq"].draw()
                    self._update_acq_status()

                self.task_queue.put(("plot", "acq", do_plot))
            except Exception as e:
                self.task_queue.put(("error", "acq", str(e)))

        self._run_in_thread(worker)

    def _save_calibration(self, ti, prefix, log_tab="acq"):
        try:
            cal_dir = CALIBRACIONES_DIR
            cal_dir.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(cal_dir / "calibracion_actual.npz",
                                wavelengths=self.calibration_wavelengths
                                if self.calibration_wavelengths is not None
                                else self.spectrometer.wavelengths,
                                R_0=self.R_0, R_1=self.R_1)
            meta = {
                "fecha_calibracion": datetime.now().isoformat(timespec="seconds"),
                "tiempo_integracion_s": float(ti),
                "num_mediciones_promedio": self._get_int(f"{prefix}_nmed"),
            }
            with open(cal_dir / "calibracion_actual.json", "w", encoding="utf-8") as f:
                json.dump(meta, f, indent=2, ensure_ascii=False)
            self.calibration_loaded = True
            self.task_queue.put(("log", log_tab, "Calibracion guardada automaticamente."))
        except Exception as e:
            self.task_queue.put(("log", log_tab, f"Aviso: no se pudo guardar calibracion: {e}"))

    def _acq_temporal_series(self, prefix):
        self._append_log("acq", "Adquiriendo serie temporal...")

        def worker():
            try:
                ti = self._acq_state.get("ti") or self._get_float(f"{prefix}_ti")
                n_esp = self._get_int("t_nesp")
                tiempos, espectros_raw, t_total = self.spectrometer.acquire_temporal_series(
                    ti, n_esp,
                    progress_cb=lambda c, t: self.task_queue.put(("progress", "acq", c, t)),
                    log_cb=lambda msg: self.task_queue.put(("log", "acq", msg)),
                )
                self._acq_state["temporal_data"] = (tiempos, espectros_raw, t_total)
                self._acq_state["rm_done"] = True
                self.task_queue.put(("log", "acq",
                    f"Serie completada: {len(espectros_raw)} espectros en {t_total:.3f} s"))

                def do_plot():
                    wl = self.spectrometer.wavelengths
                    lmin = self._get_float(f"{prefix}_lmin")
                    lmax = self._get_float(f"{prefix}_lmax")
                    mask = (wl >= lmin) & (wl <= lmax)
                    VisualizationEngine.plot_single_spectrum(
                        self._figs["acq"], self._axes["acq"],
                        wl[mask], espectros_raw[-1][mask],
                        ylabel="Intensidad", title="Ultimo espectro de la serie",
                        color="green"
                    )
                    self._canvases["acq"].draw()
                    self._update_acq_status()

                self.task_queue.put(("plot", "acq", do_plot))
            except Exception as e:
                self.task_queue.put(("error", "acq", str(e)))

        self._run_in_thread(worker)

    def _create_subject_folder(self, temporal=False):
        ruta_sujeto, ruta_series, sig_id = self._preview_subject_paths(temporal=temporal)
        ruta_sujeto.parent.mkdir(parents=True, exist_ok=True)
        ruta_series.mkdir(parents=True)

        datos = [
            f"sujeto_id: {sig_id:03d}",
            f"fecha_registro: {datetime.now().strftime('%d/%m/%Y %H:%M:%S')}",
            f"nombre: {self.param_vars.get('suj_nombre', tk.StringVar()).get() or 'No especificado'}",
            f"edad: {self.param_vars.get('suj_edad', tk.StringVar()).get() or 'No especificado'}",
            f"municipio: {self.param_vars.get('suj_municipio', tk.StringVar()).get() or 'No especificado'}",
            f"fitzpatrick: {self.param_vars.get('suj_fitz', tk.StringVar()).get() or 'No especificado'}",
            f"afecciones: {self.param_vars.get('suj_afecciones', tk.StringVar()).get() or 'No especificado'}",
        ]
        (ruta_sujeto / "datos_sujeto.txt").write_text("\n".join(datos), encoding="utf-8")
        self._update_save_preview()
        return ruta_sujeto, ruta_series, sig_id

    def _acq_calculate_mr_single(self, prefix):
        s = self._acq_state
        if not (s["r0_done"] and s["r1_done"] and s["rm_done"]):
            messagebox.showwarning("Flujo", "Completa la calibracion y mide R_M antes de guardar.")
            return

        wl = self.spectrometer.wavelengths
        r0, r1, rm = s["R_0"], s["R_1"], s["R_M"]
        lmin = self._get_float(f"{prefix}_lmin")
        lmax = self._get_float(f"{prefix}_lmax")
        fc = self._get_float(f"{prefix}_bfc")
        orden = self._get_int(f"{prefix}_border")
        mobj = self._get_int(f"{prefix}_mobj")
        r_std = self._get_float(f"{prefix}_rstd")

        sp = self.sp
        wl_t, r0_t = sp.truncate_to_range(wl, r0, lmin, lmax)
        _, r1_t = sp.truncate_to_range(wl, r1, lmin, lmax)
        _, rm_t = sp.truncate_to_range(wl, rm, lmin, lmax)
        r0_b = sp.apply_butterworth(r0_t, fc, orden)
        r1_b = sp.apply_butterworth(r1_t, fc, orden)
        rm_b = sp.apply_butterworth(rm_t, fc, orden)
        mr = sp.compute_reflectance(rm_b, r0_b, r1_b, r_std)

        factor = max(1, len(wl_t) // mobj)
        wl_dec = wl_t[::factor]
        r0_dec = r0_b[::factor]
        r1_dec = r1_b[::factor]
        rm_dec = rm_b[::factor]
        mr_dec = mr[::factor]

        ruta_sujeto, ruta_series, sig_id = self._create_subject_folder(temporal=False)

        pd.DataFrame({"wavelength_nm": wl_dec, "intensity": r0_dec}).to_csv(
            ruta_series / "R_0_data.csv", index=False)
        pd.DataFrame({"wavelength_nm": wl_dec, "intensity": r1_dec}).to_csv(
            ruta_series / "R_1_data.csv", index=False)
        pd.DataFrame({"wavelength_nm": wl_dec, "intensity": rm_dec}).to_csv(
            ruta_series / "R_M_data.csv", index=False)
        pd.DataFrame({"wavelength_nm": wl_dec, "reflectance": mr_dec}).to_csv(
            ruta_series / "M_R_data.csv", index=False)

        fig = self._figs["acq"]
        ax = self._axes["acq"]
        VisualizationEngine.plot_single_spectrum(
            fig, ax, wl_dec, mr_dec,
            ylabel="Reflectancia", title=f"M_R ({len(wl_dec)} puntos)", color="green"
        )
        self._canvases["acq"].draw()

        self._append_log("acq", f"Datos guardados en: {ruta_sujeto}")
        messagebox.showinfo("Listo", f"M_R calculado y guardado en:\n{ruta_series}")

    def _acq_calculate_mr_temporal(self, prefix):
        s = self._acq_state
        if not (s["r0_done"] and s["r1_done"] and s["rm_done"]):
            messagebox.showwarning("Flujo", "Completa la calibracion y la serie temporal antes de guardar.")
            return
        if s["temporal_data"] is None:
            messagebox.showwarning("Flujo", "No hay datos de serie temporal.")
            return

        self._append_log("acq", "Procesando serie temporal...")

        def worker():
            try:
                tiempos, espectros_raw, t_total = s["temporal_data"]
                wl = self.spectrometer.wavelengths
                r0, r1 = s["R_0"], s["R_1"]
                lmin = self._get_float(f"{prefix}_lmin")
                lmax = self._get_float(f"{prefix}_lmax")
                fc = self._get_float(f"{prefix}_bfc")
                orden = self._get_int(f"{prefix}_border")
                mobj = self._get_int(f"{prefix}_mobj")
                r_std = self._get_float(f"{prefix}_rstd")

                refs = SignalProcessing.prepare_processed_references(
                    r0, r1, wl, lmin, lmax, fc, orden, mobj, r_std)

                ruta_sujeto, ruta_series, sig_id = self._create_subject_folder(temporal=True)

                pd.DataFrame({"wavelength_nm": refs["wavelengths_diezmado"],
                               "intensity": refs["r0_diezmado"]}).to_csv(
                    ruta_series / "R_0_data.csv", index=False)
                pd.DataFrame({"wavelength_nm": refs["wavelengths_diezmado"],
                               "intensity": refs["r1_diezmado"]}).to_csv(
                    ruta_series / "R_1_data.csv", index=False)

                rows_raw = []
                for med, (t_rel, esp_raw) in enumerate(zip(tiempos, espectros_raw), 1):
                    n_wl = len(wl)
                    rows_raw.append(np.column_stack([
                        np.full(n_wl, med), np.full(n_wl, t_rel), wl, esp_raw
                    ]))
                df_raw = pd.DataFrame(np.vstack(rows_raw),
                                      columns=["medicion", "tiempo", "lambda", "intensidad"])
                df_raw["medicion"] = df_raw["medicion"].astype(int)
                df_raw.to_csv(ruta_series / "R_M_trueraw_tiempo_data.csv", index=False)

                sp = SignalProcessing
                wl_dec = refs["wavelengths_diezmado"]
                r0_b = refs["r0_butter"]
                denom = refs["denominador"]
                factor_d = refs["factor_diezmado"]

                rows_mr = []
                for med, (t_rel, esp_raw) in enumerate(zip(tiempos, espectros_raw), 1):
                    _, rm_t = sp.truncate_to_range(wl, esp_raw, lmin, lmax)
                    rm_b = sp.apply_butterworth(rm_t, fc, orden)
                    num = rm_b - r0_b
                    mr_vals = r_std * (num / denom)
                    mr_dec = mr_vals[::factor_d]
                    n_pts = len(wl_dec)
                    rows_mr.append(np.column_stack([
                        np.full(n_pts, med), np.full(n_pts, t_rel), wl_dec, mr_dec
                    ]))

                df_mr = pd.DataFrame(np.vstack(rows_mr),
                                     columns=["medicion", "tiempo", "lambda", "reflectancia"])
                df_mr["medicion"] = df_mr["medicion"].astype(int)
                ruta_mr = ruta_series / "M_R_tiempo_data.csv"
                df_mr.to_csv(ruta_mr, index=False)

                self.task_queue.put(("log", "acq",
                    f"Serie guardada: {len(espectros_raw)} espectros en {ruta_sujeto}"))
                self.task_queue.put(("done", "acq", ("temporal_done", str(ruta_mr))))
            except Exception as e:
                self.task_queue.put(("error", "acq", str(e)))

        self._run_in_thread(worker)

    def _on_acq_task_done(self, result):
        tag, data = result
        if tag == "temporal_done":
            messagebox.showinfo("Listo", f"Serie temporal guardada en:\n{data}")
        self._update_acq_status()

    # ================================================================
    # TAB 3: VISUALIZACION M_R SENCILLO
    # ================================================================

    def _build_tab_viz_mr_single(self):
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="3. Viz M_R")

        top = ttk.Frame(tab)
        top.pack(fill="x", padx=8, pady=(6, 2))

        ttk.Label(top, text="Archivo CSV de reflectancia:").pack(side="left", padx=(0, 4))
        self.param_vars["vmrs_csv"] = tk.StringVar(value=str(MR_SINGLE_CSV_DEFAULT))
        ttk.Entry(top, textvariable=self.param_vars["vmrs_csv"], width=55).pack(side="left", padx=2)
        ttk.Button(top, text="...",
                   command=lambda: self._browse_file("vmrs_csv", [("CSV", "*.csv")])).pack(side="left", padx=2)

        cfg = ttk.LabelFrame(tab, text="  Configuracion de Visualizacion  ")
        cfg.pack(fill="x", padx=8, pady=4)
        dv = DEFAULTS["viz"]
        self._make_param_entry(cfg, "Lambda min (nm):", "vmrs_lmin", dv["lambda_vis_min"], 0)
        self._make_param_entry(cfg, "Lambda max (nm):", "vmrs_lmax", dv["lambda_vis_max"], 0, col=1)
        ttk.Button(cfg, text="Graficar", command=self._viz_mr_single_plot).grid(
            row=1, column=0, columnspan=4, sticky="ew", padx=6, pady=(4, 2))
        self._add_confirm_button(cfg)

        self._figs["vmrs"], self._axes["vmrs"], self._canvases["vmrs"] = self._make_canvas(tab, figsize=(10, 5))

    def _viz_mr_single_plot(self):
        csv_path = self.param_vars["vmrs_csv"].get()
        if not csv_path or not Path(csv_path).exists():
            messagebox.showwarning("Aviso", "Selecciona un archivo CSV valido.")
            return
        try:
            df = pd.read_csv(csv_path)
            fig = self._figs["vmrs"]
            ax = self._axes["vmrs"]
            if "wavelength_nm" in df.columns and "reflectance" in df.columns:
                lambdas = df["wavelength_nm"].values
                values = df["reflectance"].values
                mask = (lambdas >= self._get_float("vmrs_lmin")) & (lambdas <= self._get_float("vmrs_lmax"))
                if not np.any(mask):
                    mask = np.ones_like(lambdas, dtype=bool)
                VisualizationEngine.plot_single_spectrum(
                    fig, ax, lambdas[mask], values[mask],
                    ylabel="Reflectancia", title="M_R - Reflectancia Medida", color="green"
                )
            elif "lambda" in df.columns and "reflectancia" in df.columns:
                first = df[df["medicion"] == df["medicion"].min()]
                lambdas = first["lambda"].values
                values = first["reflectancia"].values
                mask = (lambdas >= self._get_float("vmrs_lmin")) & (lambdas <= self._get_float("vmrs_lmax"))
                if not np.any(mask):
                    mask = np.ones_like(lambdas, dtype=bool)
                VisualizationEngine.plot_single_spectrum(
                    fig, ax, lambdas[mask], values[mask],
                    ylabel="Reflectancia", title="M_R - Primer espectro", color="green"
                )
            else:
                messagebox.showwarning("Aviso", f"Columnas no reconocidas: {list(df.columns)}")
                return
            self._canvases["vmrs"].draw()
        except Exception as e:
            messagebox.showerror("Error", str(e))

    # ================================================================
    # TAB 4: VISUALIZACION M_R TEMPORAL
    # ================================================================

    def _build_tab_viz_mr_temporal(self):
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="4. Viz M_R Temporal")

        top = ttk.Frame(tab)
        top.pack(fill="x", padx=8, pady=(6, 2))

        ttk.Label(top, text="CSV temporal:").pack(side="left", padx=(0, 4))
        self.param_vars["vmrt_csv"] = tk.StringVar(value=str(MR_TEMPORAL_CSV_DEFAULT))
        ttk.Entry(top, textvariable=self.param_vars["vmrt_csv"], width=55).pack(side="left", padx=2)
        ttk.Button(top, text="...",
                   command=lambda: self._browse_file("vmrt_csv", [("CSV", "*.csv")])).pack(side="left", padx=2)

        cfg = ttk.LabelFrame(tab, text="  Configuracion de Visualizacion  ")
        cfg.pack(fill="x", padx=8, pady=4)
        dv = DEFAULTS["viz"]
        self._make_param_entry(cfg, "Lambda min (nm):", "vmrt_lmin", dv["lambda_vis_min"], 0)
        self._make_param_entry(cfg, "Lambda max (nm):", "vmrt_lmax", dv["lambda_vis_max"], 0, col=1)
        self._make_param_entry(cfg, "Espectro inicio:", "vmrt_eini", dv["espectro_inicio"], 1)
        self._make_param_entry(cfg, "Espectro fin:", "vmrt_efin", dv["espectro_fin"], 1, col=1)
        self._make_param_entry(cfg, "Promedio waterfall:", "vmrt_avg", dv["promedio_espectros"], 2)
        self._add_confirm_button(cfg)

        btn_bar = ttk.Frame(tab)
        btn_bar.pack(fill="x", padx=8, pady=2)
        ttk.Button(btn_bar, text="Waterfall 3D",
                   command=lambda: self._temporal_viz("vmrt", "red", "Reflectancia",
                                                      "M_R Waterfall 3D")).pack(side="left", padx=4)
        ttk.Button(btn_bar, text="Animacion 2D",
                   command=lambda: self._temporal_anim_setup("vmrt", "red",
                                                              "Reflectancia")).pack(side="left", padx=4)

        self._figs["vmrt"], self._axes["vmrt"], self._canvases["vmrt"] = self._make_canvas(tab, figsize=(10, 5))

        ctrl = ttk.Frame(tab)
        ctrl.pack(fill="x", padx=8, pady=(2, 6))

        ttk.Button(ctrl, text="\u25b6 Play",
                   command=lambda: self._anim_play("vmrt")).pack(side="left", padx=4)
        ttk.Button(ctrl, text="\u23f8 Pausa",
                   command=lambda: self._anim_pause("vmrt")).pack(side="left", padx=4)

        self._anim_states["vmrt"] = {
            "playing": False, "frame": 0, "after_id": None, "data": None,
            "interval": 150,
        }

        slider_var = tk.IntVar(value=1)
        self._anim_states["vmrt"]["slider_var"] = slider_var
        slider_label = ttk.Label(ctrl, text="Medicion: 1", width=20)
        slider_label.pack(side="left", padx=(12, 4))
        self._anim_states["vmrt"]["slider_label"] = slider_label

        slider = ttk.Scale(ctrl, from_=1, to=100, variable=slider_var,
                           command=lambda val: self._anim_on_slider("vmrt", val))
        slider.pack(side="left", fill="x", expand=True, padx=4)
        self._anim_states["vmrt"]["slider"] = slider

        ttk.Label(ctrl, text="Vel:").pack(side="left", padx=(8, 2))
        speed_var = tk.StringVar(value="150")
        speed_combo = ttk.Combobox(ctrl, textvariable=speed_var, width=6, state="readonly",
                                   values=["50", "100", "150", "200", "300", "500"])
        speed_combo.pack(side="left", padx=2)
        speed_combo.bind("<<ComboboxSelected>>",
                         lambda e: self._anim_set_speed("vmrt", int(speed_var.get())))
        ttk.Label(ctrl, text="ms").pack(side="left")

    # ================================================================
    # TAB 5: PROCESAMIENTO IAD
    # ================================================================

    def _build_tab_iad(self):
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="5. Procesamiento IAD")

        top = ttk.Frame(tab)
        top.pack(fill="x", padx=8, pady=(6, 2))

        ttk.Label(top, text="Modo:", font=("Segoe UI", 10, "bold")).pack(side="left", padx=(0, 6))
        self.iad_mode_var = tk.StringVar(value="single")
        ttk.Radiobutton(top, text="Espectro unico", variable=self.iad_mode_var,
                        value="single", command=self._iad_mode_changed).pack(side="left", padx=4)
        ttk.Radiobutton(top, text="Serie temporal", variable=self.iad_mode_var,
                        value="temporal", command=self._iad_mode_changed).pack(side="left", padx=4)

        file_bar = ttk.Frame(tab)
        file_bar.pack(fill="x", padx=8, pady=2)
        ttk.Label(file_bar, text="CSV entrada:").pack(side="left", padx=(0, 4))
        self.param_vars["iad_csv"] = tk.StringVar(value=str(MR_SINGLE_CSV_DEFAULT))
        ttk.Entry(file_bar, textvariable=self.param_vars["iad_csv"], width=60).pack(side="left", padx=2)
        ttk.Button(file_bar, text="...",
                   command=lambda: self._browse_file("iad_csv", [("CSV", "*.csv")])).pack(side="left", padx=2)

        action_bar = ttk.Frame(tab)
        action_bar.pack(fill="x", padx=8, pady=4)
        ttk.Button(action_bar, text="EJECUTAR IAD", command=self._run_iad).pack(side="left", padx=4)
        ttk.Label(action_bar, text="(Parametros IAD configurados en pestana Configuracion)",
                  font=("Segoe UI", 8, "italic"), foreground="gray").pack(side="left", padx=12)

        self._figs["iad"], self._axes["iad"], self._canvases["iad"] = self._make_canvas(tab)
        self._pbars["iad"] = self._make_progress(tab)
        self._logs["iad"] = self._make_log(tab)

    def _iad_mode_changed(self):
        mode = self.iad_mode_var.get()
        if mode == "single":
            self.param_vars["iad_csv"].set(str(MR_SINGLE_CSV_DEFAULT))
        else:
            self.param_vars["iad_csv"].set(str(MR_TEMPORAL_CSV_DEFAULT))

    def _run_iad(self):
        csv_path = self.param_vars["iad_csv"].get()
        rxt_path = self.param_vars["iad_rxt"].get()
        exe_path = self.param_vars["iad_exe"].get()
        out_dir = self.param_vars["iad_outdir"].get()
        mode = self.iad_mode_var.get()

        if not Path(csv_path).exists():
            messagebox.showerror("Error", f"CSV no encontrado:\n{csv_path}")
            return

        proc = IADProcessor(exe_path, rxt_path)
        ok, msg = proc.validate_paths()
        if not ok:
            messagebox.showerror("Error", msg)
            return

        self._append_log("iad", f"Iniciando IAD ({mode})...")
        self._update_progress("iad", 0, 100)

        def worker():
            try:
                use_g_fijo = self.iad_gfijo_var.get()
                g_fijo = self._get_float("iad_g")
                fast = self.iad_fast_var.get()

                if mode == "single":
                    df_res, ruta = proc.run_single_spectrum(
                        csv_path, out_dir, use_g_fijo, g_fijo, fast,
                        progress_cb=lambda c, t: self.task_queue.put(("progress", "iad", c, t)),
                    )
                else:
                    workers = self._get_int("iad_workers")
                    max_med = self._get_int("iad_maxmed")
                    paso = self._get_int("iad_paso")
                    if max_med <= 0:
                        max_med = None
                    df_res, ruta = proc.run_temporal(
                        csv_path, out_dir, use_g_fijo, g_fijo, fast,
                        workers, max_med, paso,
                        progress_cb=lambda c, t: self.task_queue.put(("progress", "iad", c, t)),
                    )

                def do_plot():
                    fig = self._figs["iad"]
                    ax = self._axes["iad"]
                    if "mu_a_mm-1" in df_res.columns:
                        VisualizationEngine.plot_iad_results(fig, ax, df_res)
                    self._canvases["iad"].draw()

                self.task_queue.put(("plot", "iad", do_plot))
                self.task_queue.put(("log", "iad", f"Resultados guardados en: {ruta}"))
                self.task_queue.put(("done", "iad", ("iad_done", str(ruta))))
            except Exception as e:
                self.task_queue.put(("error", "iad", str(e)))

        self._run_in_thread(worker)

    def _on_iad_task_done(self, result):
        tag, data = result
        messagebox.showinfo("IAD Completado", f"Resultados en:\n{data}")

    # ================================================================
    # TAB 6: VISUALIZACION IAD SENCILLO
    # ================================================================

    def _build_tab_viz_iad_single(self):
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="6. Viz IAD")

        top = ttk.Frame(tab)
        top.pack(fill="x", padx=8, pady=(6, 2))

        ttk.Label(top, text="CSV resultados IAD:").pack(side="left", padx=(0, 4))
        self.param_vars["viads_csv"] = tk.StringVar(value=str(IAD_SINGLE_CSV_DEFAULT))
        ttk.Entry(top, textvariable=self.param_vars["viads_csv"], width=55).pack(side="left", padx=2)
        ttk.Button(top, text="...",
                   command=lambda: self._browse_file("viads_csv", [("CSV", "*.csv")])).pack(side="left", padx=2)

        cfg = ttk.LabelFrame(tab, text="  Configuracion de Visualizacion  ")
        cfg.pack(fill="x", padx=8, pady=4)
        dv = DEFAULTS["viz"]
        self._make_param_entry(cfg, "Lambda min (nm):", "viads_lmin", dv["lambda_vis_min"], 0)
        self._make_param_entry(cfg, "Lambda max (nm):", "viads_lmax", dv["lambda_vis_max"], 0, col=1)
        ttk.Label(cfg, text="Parametro:", anchor="e").grid(row=1, column=0, sticky="e", padx=(6, 2), pady=2)
        self.param_vars["viads_param"] = tk.StringVar(value=dv["param_iad"])
        ttk.Combobox(cfg, textvariable=self.param_vars["viads_param"], width=18, state="readonly",
                     values=["mu_a_mm-1", "mu_s_prime_mm-1", "g"]).grid(
                         row=1, column=1, sticky="w", padx=(2, 6), pady=2)
        ttk.Button(cfg, text="Graficar", command=self._viz_iad_single_plot).grid(
            row=1, column=2, columnspan=2, sticky="ew", padx=6, pady=2)
        self._add_confirm_button(cfg)

        self._figs["viads"], self._axes["viads"], self._canvases["viads"] = self._make_canvas(tab, figsize=(10, 5))

    def _viz_iad_single_plot(self):
        csv_path = self.param_vars["viads_csv"].get()
        if not csv_path or not Path(csv_path).exists():
            messagebox.showwarning("Aviso", "Selecciona un archivo CSV valido.")
            return
        try:
            df = pd.read_csv(csv_path)
            fig = self._figs["viads"]
            ax = self._axes["viads"]
            param = self.param_vars["viads_param"].get()
            labels = {"mu_a_mm-1": "mu_a (1/mm)", "mu_s_prime_mm-1": "mu_s' (1/mm)", "g": "g"}
            ylabel = labels.get(param, param)
            if param in df.columns and "lambda_nm" in df.columns:
                lmin = self._get_float("viads_lmin")
                lmax = self._get_float("viads_lmax")
                df_plot = df[(df["lambda_nm"] >= lmin) & (df["lambda_nm"] <= lmax)].copy()
                if df_plot.empty:
                    df_plot = df.copy()
                VisualizationEngine.plot_iad_results(fig, ax, df_plot, param=param, ylabel=ylabel)
                self._canvases["viads"].draw()
            else:
                messagebox.showwarning("Aviso", f"Columna '{param}' no encontrada en el CSV.")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    # ================================================================
    # TAB 7: VISUALIZACION IAD TEMPORAL
    # ================================================================

    def _build_tab_viz_iad_temporal(self):
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="7. Viz IAD Temporal")

        top = ttk.Frame(tab)
        top.pack(fill="x", padx=8, pady=(6, 2))

        ttk.Label(top, text="CSV IAD temporal:").pack(side="left", padx=(0, 4))
        self.param_vars["viadt_csv"] = tk.StringVar(value=str(IAD_TEMPORAL_CSV_DEFAULT))
        ttk.Entry(top, textvariable=self.param_vars["viadt_csv"], width=55).pack(side="left", padx=2)
        ttk.Button(top, text="...",
                   command=lambda: self._browse_file("viadt_csv", [("CSV", "*.csv")])).pack(side="left", padx=2)

        cfg = ttk.LabelFrame(tab, text="  Configuracion de Visualizacion  ")
        cfg.pack(fill="x", padx=8, pady=4)
        dv = DEFAULTS["viz"]
        self._make_param_entry(cfg, "Lambda min (nm):", "viadt_lmin", dv["lambda_vis_min"], 0)
        self._make_param_entry(cfg, "Lambda max (nm):", "viadt_lmax", dv["lambda_vis_max"], 0, col=1)
        self._make_param_entry(cfg, "Espectro inicio:", "viadt_eini", dv["espectro_inicio"], 1)
        self._make_param_entry(cfg, "Espectro fin:", "viadt_efin", dv["espectro_fin"], 1, col=1)
        self._make_param_entry(cfg, "Promedio waterfall:", "viadt_avg", dv["promedio_espectros"], 2)
        ttk.Label(cfg, text="Parametro:", anchor="e").grid(row=2, column=2, sticky="e", padx=(6, 2), pady=2)
        self.param_vars["viadt_param"] = tk.StringVar(value=dv["param_iad"])
        ttk.Combobox(cfg, textvariable=self.param_vars["viadt_param"], width=18, state="readonly",
                     values=["mu_a_mm-1", "mu_s_prime_mm-1", "g"]).grid(
                         row=2, column=3, sticky="w", padx=(2, 6), pady=2)
        self._add_confirm_button(cfg)

        btn_bar = ttk.Frame(tab)
        btn_bar.pack(fill="x", padx=8, pady=2)

        param_labels = {"mu_a_mm-1": "mu_a (1/mm)", "mu_s_prime_mm-1": "mu_s' (1/mm)", "g": "g"}
        ttk.Button(btn_bar, text="Waterfall 3D",
                   command=lambda: self._temporal_viz("viadt", "blue",
                                                      param_labels.get(self.param_vars["viadt_param"].get(),
                                                                       self.param_vars["viadt_param"].get()),
                                                      "IAD Waterfall 3D",
                                                      param_col=self.param_vars["viadt_param"].get())
                   ).pack(side="left", padx=4)
        ttk.Button(btn_bar, text="Animacion 2D",
                   command=lambda: self._temporal_anim_setup("viadt", "blue",
                                                              param_labels.get(self.param_vars["viadt_param"].get(),
                                                                               self.param_vars["viadt_param"].get()),
                                                              param_col=self.param_vars["viadt_param"].get())
                   ).pack(side="left", padx=4)

        self._figs["viadt"], self._axes["viadt"], self._canvases["viadt"] = self._make_canvas(tab, figsize=(10, 5))

        ctrl = ttk.Frame(tab)
        ctrl.pack(fill="x", padx=8, pady=(2, 6))

        ttk.Button(ctrl, text="\u25b6 Play",
                   command=lambda: self._anim_play("viadt")).pack(side="left", padx=4)
        ttk.Button(ctrl, text="\u23f8 Pausa",
                   command=lambda: self._anim_pause("viadt")).pack(side="left", padx=4)

        self._anim_states["viadt"] = {
            "playing": False, "frame": 0, "after_id": None, "data": None,
            "interval": 150,
        }

        slider_var = tk.IntVar(value=1)
        self._anim_states["viadt"]["slider_var"] = slider_var
        slider_label = ttk.Label(ctrl, text="Medicion: 1", width=20)
        slider_label.pack(side="left", padx=(12, 4))
        self._anim_states["viadt"]["slider_label"] = slider_label

        slider = ttk.Scale(ctrl, from_=1, to=100, variable=slider_var,
                           command=lambda val: self._anim_on_slider("viadt", val))
        slider.pack(side="left", fill="x", expand=True, padx=4)
        self._anim_states["viadt"]["slider"] = slider

        ttk.Label(ctrl, text="Vel:").pack(side="left", padx=(8, 2))
        speed_var = tk.StringVar(value="150")
        speed_combo = ttk.Combobox(ctrl, textvariable=speed_var, width=6, state="readonly",
                                   values=["50", "100", "150", "200", "300", "500"])
        speed_combo.pack(side="left", padx=2)
        speed_combo.bind("<<ComboboxSelected>>",
                         lambda e: self._anim_set_speed("viadt", int(speed_var.get())))
        ttk.Label(ctrl, text="ms").pack(side="left")

    # ================================================================
    # METODOS COMPARTIDOS: CARGA TEMPORAL + WATERFALL + ANIMACION
    # ================================================================

    def _load_temporal_data(self, tab_key, param_col=None):
        csv_var = f"{tab_key}_csv"
        csv_path = self.param_vars[csv_var].get()
        if not csv_path or not Path(csv_path).exists():
            messagebox.showwarning("Aviso", "Selecciona un archivo CSV valido.")
            return None

        df = pd.read_csv(csv_path)

        if tab_key == "vmrt":
            pivot = df.pivot_table(index=["medicion", "tiempo"],
                                   columns="lambda", values="reflectancia")
            lmin_key, lmax_key = "vmrt_lmin", "vmrt_lmax"
            eini_key, efin_key = "vmrt_eini", "vmrt_efin"
        else:
            col = param_col or self.param_vars["viadt_param"].get()
            pivot = df.pivot_table(index=["medicion", "tiempo"],
                                   columns="lambda_nm", values=col)
            lmin_key, lmax_key = "viadt_lmin", "viadt_lmax"
            eini_key, efin_key = "viadt_eini", "viadt_efin"
        pivot = pivot.sort_index()

        tiempos = np.array([t for _, t in pivot.index])
        lambdas = pivot.columns.values.astype(float)
        espectros = pivot.values

        lmin = self._get_float(lmin_key)
        lmax = self._get_float(lmax_key)
        mask_l = (lambdas >= lmin) & (lambdas <= lmax)
        lambdas = lambdas[mask_l]
        espectros = espectros[:, mask_l]

        idx_ini = self._get_int(eini_key) - 1
        idx_fin = self._get_int(efin_key)
        if idx_fin <= 0:
            idx_fin = len(tiempos)
        tiempos = tiempos[idx_ini:idx_fin]
        espectros = espectros[idx_ini:idx_fin]

        if len(tiempos) == 0 or len(lambdas) == 0:
            messagebox.showwarning("Aviso", "No hay datos en el rango seleccionado.")
            return None

        return lambdas, tiempos, espectros

    def _temporal_viz(self, tab_key, color_base, ylabel, title, param_col=None):
        data = self._load_temporal_data(tab_key, param_col)
        if data is None:
            return
        lambdas, tiempos, espectros = data

        fig = self._figs[tab_key]
        fig.clf()
        ax = fig.add_subplot(111, projection="3d")
        self._axes[tab_key] = ax

        avg = self._get_int("vmrt_avg" if tab_key == "vmrt" else "viadt_avg")
        VisualizationEngine.create_waterfall_3d(
            fig, ax, lambdas, tiempos, espectros,
            ylabel=ylabel, title=title,
            color_base=color_base, avg_count=avg
        )
        self._canvases[tab_key].draw()

    def _temporal_anim_setup(self, tab_key, color_base, ylabel, param_col=None):
        self._anim_pause(tab_key)

        data = self._load_temporal_data(tab_key, param_col)
        if data is None:
            return
        lambdas, tiempos, espectros = data

        fig = self._figs[tab_key]
        fig.clf()
        ax = fig.add_subplot(111)
        self._axes[tab_key] = ax

        anim_data = VisualizationEngine.setup_animation_2d(
            fig, ax, lambdas, tiempos, espectros,
            ylabel=ylabel, color_base=color_base
        )

        state = self._anim_states[tab_key]
        state["data"] = anim_data
        state["frame"] = 0
        state["playing"] = False
        state["slider_var"].set(1)
        state["slider"].configure(to=anim_data["n"])
        state["slider_label"].config(text=f"Medicion: 1/{anim_data['n']}")

        self._canvases[tab_key].draw()

    def _anim_play(self, tab_key):
        state = self._anim_states[tab_key]
        if state["data"] is None:
            messagebox.showinfo("Info", "Primero carga los datos y presiona 'Animacion 2D'.")
            return
        if state["playing"]:
            return
        state["playing"] = True
        self._anim_tick(tab_key)

    def _anim_pause(self, tab_key):
        state = self._anim_states[tab_key]
        state["playing"] = False
        if state.get("after_id"):
            self.root.after_cancel(state["after_id"])
            state["after_id"] = None

    def _anim_tick(self, tab_key):
        state = self._anim_states[tab_key]
        if not state["playing"] or state["data"] is None:
            return

        data = state["data"]
        idx = state["frame"]

        data["line"].set_ydata(data["espectros"][idx])
        data["line"].set_color(data["colores"][idx])
        data["title"].set_text(
            f"Medicion {idx+1}/{data['n']}  -  t = {data['tiempos'][idx]:.3f} s")
        state["slider_var"].set(idx + 1)
        state["slider_label"].config(text=f"Medicion: {idx+1}/{data['n']}")
        self._canvases[tab_key].draw_idle()

        state["frame"] = (idx + 1) % data["n"]
        state["after_id"] = self.root.after(state["interval"],
                                            lambda: self._anim_tick(tab_key))

    def _anim_on_slider(self, tab_key, val):
        state = self._anim_states[tab_key]
        if state["data"] is None:
            return
        data = state["data"]
        idx = int(float(val)) - 1
        idx = max(0, min(idx, data["n"] - 1))
        state["frame"] = idx

        data["line"].set_ydata(data["espectros"][idx])
        data["line"].set_color(data["colores"][idx])
        data["title"].set_text(
            f"Medicion {idx+1}/{data['n']}  -  t = {data['tiempos'][idx]:.3f} s")
        state["slider_label"].config(text=f"Medicion: {idx+1}/{data['n']}")
        self._canvases[tab_key].draw_idle()

    def _anim_set_speed(self, tab_key, interval_ms):
        self._anim_states[tab_key]["interval"] = interval_ms

    # ── Cleanup ──

    def _on_close(self):
        for key in self._anim_states:
            self._anim_pause(key)
        self.spectrometer.disconnect()
        plt.close("all")
        self.root.destroy()


# ============================================================
# SECCION 7: MAIN
# ============================================================

if __name__ == "__main__":
    app = IADApp()
    app.run()
