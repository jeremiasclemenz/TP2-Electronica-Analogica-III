
import numpy as np
import skrf as rf
import matplotlib.pyplot as plt

# Cargar el archivo de datos
archivo = 'C:\\Users\\Jeremias\\Documents\\GitHub\\TP2-Electronica-Analogica-III\\Hoja de datos\\BFP 450\\Infineon-RFTransistor-SPAR.zip-SM-v02_20-EN\\SPAR\\BFP450\\BFP450_w_noise_'

bfp450_1V_IC_0_10A = rf.Network(archivo + 'VCE_1.0V_IC_0.10A.s2p')

# Graficar los parámetros S, con grid y marca para 1.5 GHz
plt.figure()
bfp450_1V_IC_0_10A.plot_s_db()
plt.grid()
plt.axvline(x=1.5e9, color='r', linestyle='--')
plt.axhline(y=-0.5, color='r', linestyle='--')
plt.show()
