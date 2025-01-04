import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import json
import sys

def safe_exp(x):
    """Evita overflow em cálculos exponenciais."""
    return np.exp(np.clip(x, -500, 500))

def safe_power(base, exp):
    """Evita overflow em cálculos de potência."""
    return np.power(np.clip(base, -1e3, 1e3), exp)

def hodgkin_huxley_cable_model(config_file):
    with open(config_file, 'r') as file:
        params = json.load(file)

    # Parâmetros do modelo
    cm = params.get("cm", 1.0)
    a = params.get("a", 1.0)
    rl = params.get("rl", 1.0)
    gna = params.get("gna", 120.0)
    gk = params.get("gk", 36.0)
    gl = params.get("gl", 0.3)
    ena = params.get("ena", 50.0)
    ek = params.get("ek", -77.0)
    el = params.get("el", -54.4)
    T_max = params.get("T_max", 10.0)
    L_max = params.get("L_max", 1.0)
    dt = params.get("dt", 0.0001)
    dx = params.get("dx", 0.001)
    vm0 = params.get("vm0", -65.0)
    m0 = params.get("m0", 0.05)
    h0 = params.get("h0", 0.6)
    n0 = params.get("n0", 0.32)
    J = np.array(params["J"])
    Mie = np.array(params["Mie"])

    Nx = int(L_max / dx) + 1
    Nt = int(T_max / dt) + 1

    Vm = np.full((Nt, Nx), vm0)
    m = np.full(Nx, m0)
    h = np.full(Nx, h0)
    n = np.full(Nx, n0)

    Vm[0, Nx // 2] = 20.0

    def alpha_m(V): return 0.1 * (V + 40) / (1 - safe_exp(-(V + 40) / 10))
    def beta_m(V): return 4 * safe_exp(-(V + 65) / 18)
    def alpha_h(V): return 0.07 * safe_exp(-(V + 65) / 20)
    def beta_h(V): return 1 / (1 + safe_exp(-(V + 35) / 10))
    def alpha_n(V): return 0.01 * (V + 55) / (1 - safe_exp(-(V + 55) / 10))
    def beta_n(V): return 0.125 * safe_exp(-(V + 65) / 80)

    def calculate_currents(V, m, h, n):
        Ina = gna * safe_power(m, 3) * h * (V - ena)
        Ik = gk * safe_power(n, 4) * (V - ek)
        Il = gl * (V - el)
        return Ina, Ik, Il

    for t in range(1, Nt):
        for x in range(1, Nx - 1):
            d2Vm_dx2 = (Vm[t-1, x+1] - 2 * Vm[t-1, x] + Vm[t-1, x-1]) / dx**2
            if Mie[x] == 1:
                d2Vm_dx2 /= 100

            Ina, Ik, Il = calculate_currents(Vm[t-1, x], m[x], h[x], n[x])
            Iion = Ina + Ik + Il

            Vm[t, x] = Vm[t-1, x] + dt / (cm * a) * (-Iion + rl * d2Vm_dx2 + J[t, x])
            Vm[t, x] = np.clip(Vm[t, x], -100, 50)  # Limitar valores extremos

            m[x] += dt * (alpha_m(Vm[t-1, x]) * (1 - m[x]) - beta_m(Vm[t-1, x]) * m[x])
            h[x] += dt * (alpha_h(Vm[t-1, x]) * (1 - h[x]) - beta_h(Vm[t-1, x]) * h[x])
            n[x] += dt * (alpha_n(Vm[t-1, x]) * (1 - n[x]) - beta_n(Vm[t-1, x]) * n[x])

            # Corrigir valores inválidos
            if np.isnan(m[x]) or np.isinf(m[x]):
                m[x] = 0.05
            if np.isnan(h[x]) or np.isinf(h[x]):
                h[x] = 0.6
            if np.isnan(n[x]) or np.isinf(n[x]):
                n[x] = 0.32

    np.savetxt("results_table.csv", Vm, delimiter=",")

    fig, ax = plt.subplots()
    line, = ax.plot(np.linspace(0, L_max, Nx), Vm[0])
    ax.set_xlim(0, L_max)
    ax.set_ylim(-80, 60)
    ax.set_xlabel("Posição (cm)")
    ax.set_ylabel("Voltagem (mV)")

    def update(frame):
        line.set_ydata(Vm[frame])
        return line,

    ani = FuncAnimation(fig, update, frames=range(0, Nt, 10), blit=True)
    ani.save("voltage_position.gif", fps=60)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python main.py config.txt")
    else:
        config_path = sys.argv[1]
        hodgkin_huxley_cable_model(config_path)
