import numpy as np
import skrf as rf
import pandas as pd
import matplotlib.pyplot as plt
import re
import os

# Datos
f = 1.6e9

IC_sel = 110/1000
Vce_sel = 2.5

# Parametros S
dato = r'C:\Users\Jeremias\Documents\GitHub\TP2-Electronica-Analogica-III\Hoja de datos\BFP 450\Infineon-RFTransistor\SPAR\BFP450/BFP450_w_noise_'
str_model = dato+'VCE_'+str(Vce_sel)+'V_IC_'+str(IC_sel)+'A.s2p'

model = rf.Network(str_model)
S11 = model['1.6ghz'].s[0][0][0]
S12 = model['1.6ghz'].s[0][0][1]
S21 = model['1.6ghz'].s[0][1][0]
S22 = model['1.6ghz'].s[0][1][1]

S11mod = np.abs(S11)
S11fase = np.angle(S11) * 180 / np.pi
S12mod = np.abs(S12)
S12fase = np.angle(S12) * 180 / np.pi
S21mod = np.abs(S21)
S21fase = np.angle(S21) * 180 / np.pi
S22mod = np.abs(S22)
S22fase = np.angle(S22) * 180 / np.pi





# Cálculo de Delta
Delta = S11 * S22 - S12 * S21
Deltamod = np.abs(Delta)
Deltafase = np.angle(Delta) * 180 / np.pi

# Cálculo de K
k = (1 - S11mod**2 - S22mod**2 + Deltamod**2) / (2 * np.abs(S12 * S21))

# Estabilidad del sistema
if k > 1 and Deltamod < 1:
    print('SISTEMA INCONDICIONALMENTE ESTABLE\n')
else:
    print('SISTEMA CONDICIONALMENTE ESTABLE\n')

# Cálculo de B1 y B2
B1 = 1 + S11mod**2 - S22mod**2 - Deltamod**2
B2 = 1 + S22mod**2 - S11mod**2 - Deltamod**2

# Cálculo de C1 y C2
C1 = S11 - Delta * np.conj(S22)
C2 = S22 - Delta * np.conj(S11)
C1mod = np.abs(C1)
C1fase = np.angle(C1) * 180 / np.pi
C2mod = np.abs(C2)
C2fase = np.angle(C2) * 180 / np.pi

# Cálculo de MAG
MAG = 10 * np.log10(S21mod / S12mod)

# Cálculo de MSG
if B1 > 0:
    MSG = MAG + 10 * np.log10(k - np.sqrt(k**2 - 1))
    Gammainmod = np.abs((B1 - np.sqrt(B1**2 - 4 * C1mod**2)) / (2 * C1mod))
else:
    MSG = MAG + 10 * np.log10(k + np.sqrt(k**2 - 1))
    Gammainmod = np.abs((B1 + np.sqrt(B1**2 - 4 * C1mod**2)) / (2 * C1mod))

Gammainfase = C1fase
Gammain = Gammainmod * np.exp(1j * np.deg2rad(Gammainfase))

# Cálculo de Gammaout
if B2 > 0:
    Gammaoutmod = np.abs((B2 - np.sqrt(B2**2 - 4 * C2mod**2)) / (2 * C2mod))
else:
    Gammaoutmod = np.abs((B2 + np.sqrt(B2**2 - 4 * C2mod**2)) / (2 * C2mod))

Gammaoutfase = C2fase
Gammaout = Gammaoutmod * np.exp(1j * np.deg2rad(Gammaoutfase))

# Impresión de resultados
print('∆ = {:.2f} + j {:.2f} = {:.2f} ∠ {:.2f}°\n'.format(np.real(Delta), np.imag(Delta), Deltamod, Deltafase))
print('K = {:.2f}\n'.format(k))
print('B1 = {:.2f}\n'.format(B1))
print('B2 = {:.2f}\n'.format(B2))
print('C1 = {:.2f} + j {:.2f} = {:.2f} ∠ {:.2f}°\n'.format(np.real(C1), np.imag(C1), C1mod, C1fase))
print('C2 = {:.2f} + j {:.2f} = {:.2f} ∠ {:.2f}°\n'.format(np.real(C2), np.imag(C2), C2mod, C2fase))
print('MAG = {:.2f} dB\n'.format(MAG))
print('MSG = {:.2f} dB\n\n'.format(MSG))
print('Γin = {:.2f} + j {:.2f} = {:.2f} ∠ {:.2f}°\n'.format(np.real(Gammain), np.imag(Gammain), Gammainmod, Gammainfase))
print('Γout = {:.2f} + j {:.2f} = {:.2f} ∠ {:.2f}°\n'.format(np.real(Gammaout), np.imag(Gammaout), Gammaoutmod, Gammaoutfase))

# Impedancias
Zo = 50
Zin = Zo * (1 + Gammain) / (1 - Gammain)
Zout = Zo * (1 + Gammaout) / (1 - Gammaout)
Zs = np.conj(Zin)
ZL = np.conj(Zout)

# Cálculo de parámetros adicionales
Rsin = np.real(Zin)
Xsin = np.imag(Zin)
Rpin = Rsin * (1 + (Xsin / Rsin)**2)
Xpin = Rsin * Rpin / Xsin
Cin = (1 / (2 * np.pi * Xpin * f)) * 1e12
Zoin = np.sqrt(Zo * Rpin)

Rsout = np.real(Zout)
Xsout = np.imag(Zout)
Zoout = np.sqrt(Zo * Rsout)
Xcout = (Zoout**2 / Xsout) * (-1)
Cout = (1 / (2 * np.pi * f * Xcout)) * 1e12

# Impresión de impedancias y parámetros adicionales
print('Zo = {:.2f} Ω\n'.format(Zo))
print('Zin = {:.2f} + j {:.2f}\n'.format(np.real(Zin), np.imag(Zin)))
print('Zout = {:.2f} + j {:.2f}\n'.format(np.real(Zout), np.imag(Zout)))
print('Zs = {:.2f} + j {:.2f}\n'.format(np.real(Zs), np.imag(Zs)))
print('ZL = {:.2f} + j {:.2f}\n\n'.format(np.real(ZL), np.imag(ZL)))

print('Rsin = {:.2f} Ω\n'.format(Rsin))
print('Xsin = {:.2f} Ω\n'.format(Xsin))
print('Rpin = {:.2f} Ω\n'.format(Rpin))
print('Xpin = {:.2f} Ω\n'.format(Xpin))
print('Cin = {:.2f} pF\n'.format(Cin))
print('Zoin = {:.2f} Ω\n\n'.format(Zoin))

print('Rsout = {:.2f} Ω\n'.format(Rsout))
print('Xsout = {:.2f} Ω\n'.format(Xsout))
print('Zoout = {:.2f} Ω\n'.format(Zoout))
print('Xcout = {:.2f} Ω\n'.format(Xcout))
print('Cout = {:.2f} pF\n\n'.format(Cout))
