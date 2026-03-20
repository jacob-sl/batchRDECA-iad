import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def langmuir_isotherm(Ce, q_max, K_L):
    return (q_max * K_L * Ce) / (1 + K_L * Ce)

Ce = np.array([10, 25, 45, 70, 100, 140, 190])
qe = np.array([5, 12, 20, 30, 38, 43, 47])

initial_guess = (50, 0.1)

popt, pcov = curve_fit(langmuir_isotherm, Ce, qe, p0=initial_guess)

q_max_fit, K_L_fit = popt

print(f"q_max = {q_max_fit:.3f} mg/g")
print(f"K_L = {K_L_fit:.4f} L/mg")

plt.scatter(Ce, qe, color='blue', label='Datos originales')
plt.plot(Ce, langmuir_isotherm(Ce, q_max_fit, K_L_fit), color='red', label='Ajuste Langmuir')

plt.xlabel('Concentración equilibrada (mg/L)')
plt.ylabel('Capacidad adsórfenta (mg/g)')
plt.title('Ajuste del isoterma de Langmuir')
plt.legend()
plt.show()