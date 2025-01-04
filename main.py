import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

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
    if x == dx:
        return 20.0  # µA/cm²
    return 0.0

# Funções auxiliares para as variáveis do modelo
def safe_exp(x):
    return np.exp(np.clip(x, -100, 100))

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
    n_x = int(L / dx)
    n_t = int(T / dt)

    V = np.ones(n_x) * -65.0
    n = np.zeros(n_x) + 0.3177
    m = np.zeros(n_x) + 0.0529
    h = np.zeros(n_x) + 0.5961

    V_time = np.zeros((n_t, n_x))
    D = (a / (2 * R)) * (dt / dx**2)

    for t_idx in range(n_t):
        V_new = V.copy()

        for x_idx in range(1, n_x - 1):
            I_Na = g_Na * m[x_idx]**3 * h[x_idx] * (V[x_idx] - E_Na)
            I_K = g_K * n[x_idx]**4 * (V[x_idx] - E_K)
            I_L = g_L * (V[x_idx] - E_L)
            I_stim = I_ap(t_idx * dt, x_idx * dx)

            dV_dt = (D * (V[x_idx+1] - 2 * V[x_idx] + V[x_idx-1]) -
                    (I_Na + I_K + I_L - I_stim)) / C_m
            V_new[x_idx] += dV_dt * dt

            dn = (alpha_n(V[x_idx]) * (1 - n[x_idx]) - beta_n(V[x_idx]) * n[x_idx]) * dt
            dm = (alpha_m(V[x_idx]) * (1 - m[x_idx]) - beta_m(V[x_idx]) * m[x_idx]) * dt
            dh = (alpha_h(V[x_idx]) * (1 - h[x_idx]) - beta_h(V[x_idx]) * h[x_idx]) * dt

            n[x_idx] += dn
            m[x_idx] += dm
            h[x_idx] += dh

        V = V_new
        V_time[t_idx, :] = V

    return V_time

# Executar simulação
V_time = hodgkin_huxley_1D()

# Criar animação
x = np.linspace(0, L, int(L/dx))

fig, ax = plt.subplots(figsize=(10, 6))
line, = ax.plot([], [], lw=2)
ax.set_xlim(0, L)
ax.set_ylim(-80, 50)
ax.set_xlabel("Posição x (cm)")
ax.set_ylabel("Potencial de membrana V (mV)")
ax.set_title("Propagação do Potencial de Ação")

def init():
    line.set_data([], [])
    return line,

def update(frame):
    line.set_data(x, V_time[frame, :])
    return line,

n_frames = V_time.shape[0]
ani = FuncAnimation(fig, update, frames=n_frames, init_func=init, blit=True)

# Salvar o GIF
ani.save("propagacao_potencial.gif", writer="imagemagick", fps=140)
