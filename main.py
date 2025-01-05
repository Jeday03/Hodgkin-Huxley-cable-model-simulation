import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd
import sys
import configparser

# Função para carregar configurações do arquivo config.txt
def carregar_configuracoes(config_file):
    config = configparser.ConfigParser()
    config.read(config_file)
    
    params = {
        'C_m': float(config['PARAMETROS']['C_m']),
        'g_Na': float(config['PARAMETROS']['g_Na']),
        'g_K': float(config['PARAMETROS']['g_K']),
        'g_L': float(config['PARAMETROS']['g_L']),
        'E_Na': float(config['PARAMETROS']['E_Na']),
        'E_K': float(config['PARAMETROS']['E_K']),
        'E_L': float(config['PARAMETROS']['E_L']),
        'a': float(config['CABO']['a']),
        'R': float(config['CABO']['R']),
        'dx': float(config['SIMULACAO']['dx']),
        'dt': float(config['SIMULACAO']['dt']),
        'L': float(config['SIMULACAO']['L']),
        'T': float(config['SIMULACAO']['T']),
    }
    return params

# Função para corrente de estímulo aplicada
def I_ap(t, x, dx):
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
def hodgkin_huxley_1D(params):
    C_m, g_Na, g_K, g_L = params['C_m'], params['g_Na'], params['g_K'], params['g_L']
    E_Na, E_K, E_L = params['E_Na'], params['E_K'], params['E_L']
    a, R, dx, dt, L, T = params['a'], params['R'], params['dx'], params['dt'], params['L'], params['T']

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
            I_stim = I_ap(t_idx * dt, x_idx * dx, dx)

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

# Função principal
def main():
    if len(sys.argv) != 2:
        print("Uso: python main.py config.txt")
        sys.exit(1)

    config_file = sys.argv[1]
    params = carregar_configuracoes(config_file)

    # Executar simulação
    V_time = hodgkin_huxley_1D(params)

    # Salvar os dados da simulação em um arquivo CSV
    x = np.linspace(0, params['L'], int(params['L']/params['dx']))
    t = np.linspace(0, params['T'], int(params['T']/params['dt']))
    df = pd.DataFrame(V_time, columns=[f"x={xi:.2f}" for xi in x])
    df.insert(0, "Tempo (ms)", t)
    df.to_csv("resultado_simulacao_HH.csv", index=False)

    # Criar animação
    fig, ax = plt.subplots(figsize=(10, 6))
    line, = ax.plot([], [], lw=2)
    ax.set_xlim(0, params['L'])
    ax.set_ylim(-80, 50)
    ax.set_xlabel("Posição x (cm)")
    ax.set_ylabel("Potencial de membrana V (mV)")
    ax.set_title("Propagação do Potencial de Ação Modelo HH em um cabo")

    def init():
        line.set_data([], [])
        return line,

    def update(frame):
        line.set_data(x, V_time[frame, :])
        return line,

    n_frames = V_time.shape[0]
    frames_reduzidos = np.arange(0, n_frames, 10)
    ani = FuncAnimation(fig, update, frames=frames_reduzidos, init_func=init, blit=True)

    # Salvar o GIF
    ani.save("propagacao_potencial_HH_cabo.gif", writer="imagemagick", fps=60)

if __name__ == "__main__":
    main()
