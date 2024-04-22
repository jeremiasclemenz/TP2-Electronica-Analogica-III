import cmath
import math
import numpy as np
import skrf as rf
from skrf.media import DistributedCircuit

# Definir las variables
S11 = cmath.rect(0.656, math.radians(146.7))  # Ejemplo de valor complejo
S12 = cmath.rect(0.122, math.radians(46.1))
S21 = cmath.rect(2.3, math.radians(44.7))
S22 = cmath.rect(0.172, math.radians(-117.1))
S11 = (-0.6887813841901239+0.27968511721780415j)
S12 = (0.0329006355170892+0.05124800662046915j)
S21 = (2.0992062681698505+3.981498969443093j)
S22 = (-0.5099988830350469+0.10747362142872341j)
f = 1.6e9
# Definir la función para calcular K, Δ, Γin, Γout
def calcular_parametros(S11, S12, S21, S22):
    delta = S11 * S22 - S12 * S21
    K = (1 + abs(delta)**2 - abs(S11)**2 - abs(S22)**2) / (2 * abs(S12 * S21))
    B1 = 1 + abs(S11)**2 - abs(S22)**2 - abs(delta)**2
    B2 = 1 + abs(S22)**2 - abs(S11)**2 - abs(delta)**2
    C1 = S11 - (delta*S22.conjugate())
    C2 = S22 - (delta*S11.conjugate())
    Gin_magnitude = abs((B1 - math.sqrt(B1**2 - 4*abs(C1)**2)) / (2*abs(C1))) #if B1 > 0 else abs(B1 + math.sqrt(B1**2 - 4*abs(C1)) / (2*abs(C1)))
    Gin = cmath.rect(Gin_magnitude, cmath.phase(C1))
    Gout_magnitude = abs((B2 - math.sqrt(B2**2 - 4*abs(C2)**2)) / (2*abs(C2))) #if B2 > 0 else abs(B2 + math.sqrt(B2**2 - 4*abs(C2)) / (2*abs(C2))))
    Gout = cmath.rect(Gout_magnitude, cmath.phase(C2))
    return K, delta, Gin, Gout, B1, B2, C1, C2

#En formato rectangular
# Calcular los parámetros
K, delta, Gin, Gout, B1, B2, C1, C2 = calcular_parametros(S11, S12, S21, S22)

# Imprimir los resultados
print("K:", K)
print("Δ:", delta)
print("B1", B1)
print("B2", B2)
print("C1", C1)
print("C2", C2)
print("Γin:", Gin)
print("Γout:", Gout)

#En formato polar
# Calcular los parámetros
K, delta, Gin, Gout, B1, B2, C1, C2 = calcular_parametros(S11, S12, S21, S22)

# Convertir a formato polar
delta_polar = cmath.polar(delta)
Gin_polar = cmath.polar(Gin)
Gout_polar = cmath.polar(Gout)

# Convertir los ángulos a grados
delta_polar = (delta_polar[0], math.degrees(delta_polar[1]))
Gin_polar = (Gin_polar[0], math.degrees(Gin_polar[1]))
Gout_polar = (Gout_polar[0], math.degrees(Gout_polar[1]))
print("\n---------- Cálculo de K, Δ, Γin, Γout ----------")
# Imprimir los resultados
print("K:", K)
print("Δ (polar):", delta_polar)
print("Γin (polar):", Gin_polar)
print("Γout (polar):", Gout_polar)

print("\n---------- Cálculo de Zin y Zout ----------")
# Definir la variable Zo
Zo = 50  # Impedancia característica

# Calcular Zin y Zout
Zin = Zo * (1 + Gin) / (1 - Gin)
Zout = Zo * (1 + Gout) / (1 - Gout)

# Imprimir los resultados
print("Zin:", Zin)
print("Zout:", Zout)

print("\n---------- Cálculo de Zs y ZL ----------")
# Calcular Zs y ZL
Zs = Zin.conjugate()
ZL = Zout.conjugate()

# Imprimir los resultados
print("Zs:", Zs)
print("ZL:", ZL)

print("\n---------- Pasaje del modelo serie a paralelo ----------")

# Calcular Rp y Xp
Rp = (Zin.real**2 + Zin.imag**2) / Zin.real #Zin.real*(1+(Zin.imag/Zin.real)**2 )
Xp = (Zin.real**2 + Zin.imag**2) / Zin.imag #Zin.imag*(1+(Zin.real/Zin.imag)**2 )

# Imprimir los resultados
print("Rp:", Rp)
print("Xp:", Xp)

print("\n---------- Cálculo del capacitor ----------")

# Calcular C
C = 1 / (2 * np.pi * f * Xp)

# Imprimir el resultado
print("C:", C)

print("\n---------- Cálculo del conversor λ/4 ----------")
Zgen = 50 #Impedancia de la fuente
# Calcular el conversor λ/4
lambda_4 = np.sqrt(Rp* Zgen)

# Imprimir el resultado
print("Conversor λ/4:", lambda_4)



print("\n---------- Cálculo de Zo ----------")

# Calcular Zo
Zo = np.sqrt(Zout.real* 50)

# Imprimir el resultado
print("Zo:", Zo)

print("\n---------- Cálculo del capacitor para cancelar capacitor de Zout ----------")

print("\n---------- Cálculo de Xc ----------")

# Calcular Xc
Xc = Zo**2 / -Zout.imag

# Imprimir el resultado
print("Xc:", Xc)

print("\n---------- Cálculo del capacitor para la frecuencia utilizada ----------")

# Calcular C
C = 1 / (2 * np.pi * f * Xc)

# Imprimir el resultado
print("C:", C)


