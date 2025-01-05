import numpy as np
import matplotlib.pyplot as plt

# Parâmetros da equação do cabo e do modelo Hodgkin-Huxley
C_m = 1.0        # Capacitância da membrana (µF/cm²)
g_Na = 120.0     # Condutância máxima de Na+ (mS/cm²)
g_K = 36.0       # Condutância máxima de K+ (mS/cm²)
g_L = 0.3        # Condutância de vazamento (mS/cm²)
E_Na = 50.0      # Potencial de reversão de Na+ (mV)
E_K = -77.0      # Potencial de reversão de K+ (mV)
E_L = -54.387    # Potencial de reversão de vazamento (mV)

a = 15       # raio do axônio (µm)
R = 5000        # Resistência longitudinal (ohm·cm)

dx = 0.01        # Espaçamento espacial (cm)
dt = 0.01        # Passo de tempo (ms)
L = 3.0          # Comprimento total da fibra (cm)
T = 50.0         # Tempo total de simulação (ms)

# Função para corrente de estímulo aplicada
def I_ap(t, x):
    if  x == dx:
        return 20.0  # µA/cm²
    return 0.0

def safe_exp(x):
    return np.exp(np.clip(x, -100, 100))  # Trunca valores extremos para evitar overflow

def alpha_n(V):
    return 0.01 * (V + 55) / (1 - safe_exp(-(V + 55) / 10))

def beta_n(V):
    return 0.125 * safe_exp(-(V + 65) / 80)

def alpha_m(V):
    return 0.1 * (V + 40) / (1 - safe_exp(-(V + 40) / 10))

def beta_m(V):
    return 4.0 * safe_exp(-(V + 65) / 18)

def alpha_h(V):
    return 0.07 * safe_exp(-(V + 65) / 20)

def beta_h(V):
    return 1 / (1 + safe_exp(-(V + 35) / 10))

# Simulação do modelo Hodgkin-Huxley com equação do cabo
def hodgkin_huxley_1D():
    n_x = int(L / dx)  # Número de pontos espaciais
    n_t = int(T / dt)  # Número de passos temporais

    # Variáveis de estado
    V = np.ones(n_x) * -65.0  # Potencial inicial (mV)
    n = np.zeros(n_x) + 0.3177
    m = np.zeros(n_x) + 0.0529
    h = np.zeros(n_x) + 0.5961

    # Armazenamento para visualização
    V_time = np.zeros((n_t, n_x))

    # Constante difusiva D
    D = (a / (2 * R)) * (dt / dx**2)

    # Iteração temporal
    for t_idx in range(n_t):
        V_new = V.copy()

        V_new[0] = V_new[1]
        V_new[0] = V_new[1]

        n[0] = n[1]
        m[0] = m[1]
        h[0] = h[1]

        n[-1] = n[-2]
        m[-1] = m[-2]
        h[-1] = h[-2]

        for x_idx in range(1, n_x - 1):
            # Correntes iônicas
            I_Na = g_Na * m[x_idx]**3 * h[x_idx] * (V[x_idx] - E_Na)
            I_K = g_K * n[x_idx]**4 * (V[x_idx] - E_K)
            I_L = g_L * (V[x_idx] - E_L)
            I_stim = I_ap(t_idx * dt, x_idx * dx)

            # Atualização do potencial usando a equação do cabo
            dV_dt = (D * (V[x_idx+1] - 2 * V[x_idx] + V[x_idx-1]) - (I_Na + I_K + I_L - I_stim)*dt) / C_m
            V_new[x_idx] = V[x_idx] + dV_dt
            
            # Atualização dos gates
            dn = (alpha_n(V[x_idx]) * (1 - n[x_idx]) - beta_n(V[x_idx]) * n[x_idx]) * dt
            dm = (alpha_m(V[x_idx]) * (1 - m[x_idx]) - beta_m(V[x_idx]) * m[x_idx]) * dt
            dh = (alpha_h(V[x_idx]) * (1 - h[x_idx]) - beta_h(V[x_idx]) * h[x_idx]) * dt

            n[x_idx] += dn
            m[x_idx] += dm
            h[x_idx] += dh

        # Atualizar V e armazenar o tempo
        V = V_new
        V_time[t_idx, :] = V

    return V_time

# Executar simulação
V_time = hodgkin_huxley_1D()

# Plotar o resultado em subplots
fig, axes = plt.subplots(len([0, int(10 / dt), int(20 / dt), int(30 / dt), int(40 / dt)]), 1, figsize=(10, 12), constrained_layout=True)
time_points = [0, int(10 / dt), int(20 / dt), int(30 / dt), int(40 / dt)]  # Tempos para visualizar
x = np.linspace(0, L, int(L/dx))

for ax, t in zip(axes, time_points):
    ax.plot(x, V_time[t, :])
    ax.set_title(f"t = {t*dt:.1f} ms")
    ax.set_ylabel("Potencial (mV)")
    ax.grid()

axes[-1].set_xlabel("Posição x (cm)")
plt.show()
