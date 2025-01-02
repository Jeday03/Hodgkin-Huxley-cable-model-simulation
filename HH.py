import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Constantes do modelo
Cm = 1.0        # Capacitância da membrana (uF/cm^2)
gNa = 120.0     # Condutância máxima do Na+ (mS/cm^2)
gK = 36.0       # Condutância máxima do K+ (mS/cm^2)
gL = 0.3        # Condutância de vazamento (mS/cm^2)
ENa = 50.0      # Potencial de reversão do Na+ (mV)
EK = -77.0      # Potencial de reversão do K+ (mV)
EL = -54.387    # Potencial de reversão do vazamento (mV)
Ra = 10.0       # Resistência axial (ohm.cm)

# Parâmetros do cabo
length = 1.0    # Comprimento do cabo (cm)
dx = 0.0001       # Incremento espacial (cm)
Nx = int(length / dx)  # Número de compartimentos
dt = 0.01       # Incremento temporal (ms)
time = 50.0     # Tempo total (ms)

# Condições iniciais
V_rest = -65.0 # Potencial de repouso (mV)
V = np.ones(Nx) * V_rest # Potencial de membrana (mV)
m = np.zeros(Nx) # Variável de gating m
h = np.ones(Nx) # Variável de gating h
n = np.zeros(Nx) # Variável de gating n


#m: Probabilidade de ativação dos canais de sódio.
#h: Probabilidade de inativação dos canais de sódio.
#n: Probabilidade de ativação dos canais de potássio.

# Funções de ativação/inativação
def alpha_m(V):return 0.1 * (V + 40) / (1 - np.exp(-(V + 40) / 10)) # Função de ativação sódio
def beta_m(V): return 4.0 * np.exp(-(V + 65) / 18) # Função de inativação sódio
def alpha_h(V): return 0.07 * np.exp(-(V + 65) / 20) # Função de ativação sódio
def beta_h(V): return 1 / (1 + np.exp(-(V + 35) / 10)) # Função de inativação sódio
def alpha_n(V): return 0.01 * (V + 55) / (1 - np.exp(-(V + 55) / 10)) # Função de ativação potássio
def beta_n(V): return 0.125 * np.exp(-(V + 65) / 80) # Função de inativação potássio

# Correntes iônicas
def I_Na(V, m, h): return gNa * m**3 * h * (V - ENa) # Corrente de sódio
def I_K(V, n): return gK * n**4 * (V - EK) # Corrente de potássio
def I_L(V): return gL * (V - EL) # Corrente de vazamento

# Estímulo externo
I_ext = np.zeros(Nx) # Corrente externa (uA/cm^2)
I_ext[Nx // 2] = 10.0  # Estímulo aplicado no meio do cabo

# Simulação
time_steps = int(time / dt) # Número de passos de tempo
Vs = np.zeros((time_steps, Nx))  # Potencial ao longo do cabo
for t in range(time_steps):
    # Atualizar variáveis de gating
    m += dt * (alpha_m(V) * (1 - m) - beta_m(V) * m)
    h += dt * (alpha_h(V) * (1 - h) - beta_h(V) * h)
    n += dt * (alpha_n(V) * (1 - n) - beta_n(V) * n)

    #print("m", m)
    #print("h", h)
    #print("n", n)

    
    # Correntes iônicas
    Iion = I_Na(V, m, h) + I_K(V, n) + I_L(V)
    

    # Atualizar potencial de membrana
    d2Vdx2 = (np.roll(V, -1) - 2 * V + np.roll(V, 1)) / dx**2
    V += dt / Cm * (-Iion + (1 / (Ra * Cm)) * d2Vdx2 + I_ext)
    #print(Cm * (-Iion + (1 / (Ra * Cm)) * d2Vdx2 + I_ext))
    
    # Salvar o potencial de membrana
    Vs[t, :] = V

#print(V)

print("Max V:", np.max(V))
print("Min V:", np.min(V))

# Selecionar algumas posições ao longo do cabo para plotar
posicoes = [int(Nx * 0.25), int(Nx * 0.5), int(Nx * 0.75)]  # 25%, 50%, 75% do cabo

# Criar a figura
plt.figure(figsize=(10, 8))

for i, pos in enumerate(posicoes, 1):
    plt.subplot(len(posicoes), 1, i)  # Criar subplots
    plt.plot(np.linspace(0, time, time_steps), Vs[:, pos], label=f'Posição {pos*dx:.2f} cm')
    plt.xlabel('Tempo (ms)')
    plt.ylabel('Potencial (mV)')
    plt.legend()
    plt.grid(True)

    # Adicionar título em cada subplot
    plt.title(f'Propagação do Potencial na Posição {i}')

plt.tight_layout()
plt.show()
