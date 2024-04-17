# Programa que realiza los calculos iniciales de impedancia

import numpy as np
import math
import matplotlib.pyplot as plt


# frecuencia en gigahertz
frecuencia = 2304 * 10**6
S11_modulo = 0.656
S11_fase = 146.7
S12_modulo = 0.122
S12_fase = 46.1
S21_modulo = 2.3
S21_fase = 44.7
S22_modulo = 0.172
S22_fase = -117.1

# Pasamos a rectangulares
S11 = S11_modulo * np.exp(1j * S11_fase)
S12 = S12_modulo * np.exp(1j * S12_fase)
S21 = S21_modulo * np.exp(1j * S21_fase)
S22 = S22_modulo * np.exp(1j * S22_fase)

# Calculo de estabilidad con Rollet 
delta = S11 * S22 - S12 * S21
delta_modulo = np.abs(delta)
delta_fase = np.angle(delta)
k = (1 - np.abs(S11)**2 - np.abs(S22)**2 + np.abs(delta)**2) / (2 * np.abs(S12 * S21))



# Si k es mayor a 1 es incondicionalmente estable

if k > 1:
    print("El circuito es incondicionalmente estable")
else:
    print("El circuito es condicionalmente estable")
    
# Calculo de coeficientes 
B1 = 1 + np.abs(S11)**2 - np.abs(S22)**2 - np.abs(delta)**2
B2 = 1 + np.abs(S22)**2 - np.abs(S11)**2 - np.abs(delta)**2

C1 = S11 - delta * np.conj(S22)
C2 = S22 - delta * np.conj(S11)

B1_modulo = np.abs(B1)
B2_modulo = np.abs(B2)
C1_modulo = np.abs(C1)
C2_modulo = np.abs(C2)
B1_fase = np.angle(B1)
B2_fase = np.angle(B2)
C1_fase = np.angle(C1)
C2_fase = np.angle(C2)  

# Calculo del modulo del coeficiente de reflexion de entrada Gamma_rm
if B1_modulo > 0:
    Gamma_ms = np.conj(C1) * ((B1 - np.sqrt(B1**2 - 4 * C1_modulo**2)) / (2 * C1_modulo**2))
else:
    Gamma_ms = np.conj(C1) * ((B1 + np.sqrt(B1**2 - 4 * C1_modulo**2)) / (2 * C1_modulo**2))

Gamma_ms_modulo = np.abs(Gamma_ms)
Gamma_ms_fase = np.angle(Gamma_ms)


# Calculo del modulo del coeficiente de reflexion de salida Gamma_ml

if B2_modulo > 0:
    Gamma_ml = np.conj(C2) * ((B2 - np.sqrt(B2**2 - 4 * C2_modulo**2)) / (2 * C2_modulo**2))
else:
    Gamma_ml = np.conj(C2) * ((B2 + np.sqrt(B2**2 - 4 * C2_modulo**2)) / (2 * C2_modulo**2))
    
Gamma_ml_modulo = np.abs(Gamma_ml)
Gamma_ml_fase = np.angle(Gamma_ml)

Gamma_s = Gamma_ms
Gamma_l = Gamma_ml
Gamma_in = np.conj(Gamma_s)
Gamma_in_modulo = np.abs(Gamma_in)
Gamma_in_fase = np.angle(Gamma_in)
Gamma_out = np.conj(Gamma_l)
Gamma_out_modulo = np.abs(Gamma_out)
Gamma_out_fase = np.angle(Gamma_out)

# Calculo de ganancia de transconductancia Gt maxima unilateral
Gtu_max = (S11_modulo/S21_modulo) * (k - np.sqrt(k**2 - 1))

# Caalculo de impedancias

Z0 = 50
Zin = Z0 * (1 + Gamma_in) / (1 - Gamma_in)
Zout = Z0 * (1 + Gamma_out) / (1 - Gamma_out)
Zs = Z0 * (1 + Gamma_s) / (1 - Gamma_s)
Zl = Z0 * (1 + Gamma_l) / (1 - Gamma_l)

# Verificacion de los calculos

if Zs == (np.conj(Zin)) and Zl == (np.conj(Zout)):
    print("Los calculos son correctos")
else:
    print("Los calculos son incorrectos")

# Impresion de resultados. Rectangular y polar
print("Delta: ", delta)
print("Delta_modulo: ", delta_modulo)
print("Delta_fase: ", delta_fase)
print("k: ", k)

print("B1: ", B1)
print("B1_modulo: ", B1_modulo)
print("B1_fase: ", B1_fase)
print("B2: ", B2) 
print("B2_modulo: ", B2_modulo)
print("B2_fase: ", B2_fase)
print("C1: ", C1)
print("C1_modulo: ", C1_modulo)
print("C1_fase: ", C1_fase)
print("C2: ", C2)
print("C2_modulo: ", C2_modulo)
print("C2_fase: ", C2_fase)

print("Gamma_in: ", Gamma_in)
print("Gamma_in_modulo: ", Gamma_in_modulo)
print("Gamma_in_fase: ", Gamma_in_fase)
print("Gamma_out: ", Gamma_out)
print("Gamma_out_modulo: ", Gamma_out_modulo)
print("Gamma_out_fase: ", Gamma_out_fase)

print("Zin: ", Zin)
print("Zout: ", Zout)
print("Zs: ", Zs)
print("Zl: ", Zl)



