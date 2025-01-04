import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
import matplotlib.animation as animation

# Configuração do código (mesma lógica anterior)
np.random.seed(1000)
tempo_min = 0.0
tempo_max = 50.0
condutancia_K = 36.0
condutancia_Na = 120.0
condutancia_L = 0.3
capacitancia_membrana = 1.0
potencial_K = -12.0
potencial_Na = 115.0
potencial_L = 10.613
tempo = np.linspace(tempo_min, tempo_max, 1000)

def alfa_n(potencial_membrana):
    return (0.01 * (10.0 - potencial_membrana)) / (np.exp(1.0 - (0.1 * potencial_membrana)) - 1.0)

def beta_n(potencial_membrana):
    return 0.125 * np.exp(-potencial_membrana / 80.0)

def alfa_m(potencial_membrana):
    return (0.1 * (25.0 - potencial_membrana)) / (np.exp(2.5 - (0.1 * potencial_membrana)) - 1.0)

def beta_m(potencial_membrana):
    return 4.0 * np.exp(-potencial_membrana / 18.0)

def alfa_h(potencial_membrana):
    return 0.07 * np.exp(-potencial_membrana / 20.0)

def beta_h(potencial_membrana):
    return 1.0 / (np.exp(3.0 - (0.1 * potencial_membrana)) + 1.0)

def n_infinito(potencial_membrana=0.0):
    return alfa_n(potencial_membrana) / (alfa_n(potencial_membrana) + beta_n(potencial_membrana))

def m_infinito(potencial_membrana=0.0):
    return alfa_m(potencial_membrana) / (alfa_m(potencial_membrana) + beta_m(potencial_membrana))

def h_infinito(potencial_membrana=0.0):
    return alfa_h(potencial_membrana) / (alfa_h(potencial_membrana) + beta_h(potencial_membrana))

def estimulo_corrente(t):
    estimulo = np.zeros_like(t)
    estimulo[(t > 0.0) & (t < 1.0)] = 150.0
    estimulo[(t > 10.0) & (t < 11.0)] = 50.0
    return estimulo

def calcular_derivadas(estados, t0):
    derivadas = np.zeros((4,))
    potencial_membrana = estados[0]
    n = estados[1]
    m = estados[2]
    h = estados[3]

    condutancia_K_total = (condutancia_K / capacitancia_membrana) * np.power(n, 4.0)
    condutancia_Na_total = (condutancia_Na / capacitancia_membrana) * np.power(m, 3.0) * h
    condutancia_L_total = condutancia_L / capacitancia_membrana

    derivadas[0] = (estimulo_corrente(t0) / capacitancia_membrana) - (condutancia_K_total * (potencial_membrana - potencial_K)) - (condutancia_Na_total * (potencial_membrana - potencial_Na)) - (condutancia_L_total * (potencial_membrana - potencial_L))
    derivadas[1] = (alfa_n(potencial_membrana) * (1.0 - n)) - (beta_n(potencial_membrana) * n)
    derivadas[2] = (alfa_m(potencial_membrana) * (1.0 - m)) - (beta_m(potencial_membrana) * m)
    derivadas[3] = (alfa_h(potencial_membrana) * (1.0 - h)) - (beta_h(potencial_membrana) * h)

    return derivadas

estado_inicial = np.array([0.0, n_infinito(), m_infinito(), h_infinito()])
resultados = odeint(calcular_derivadas, estado_inicial, tempo)

# Criar animação do potencial do neurônio
figura, eixo = plt.subplots(figsize=(8, 6))
eixo.set_xlim(tempo_min, tempo_max)
eixo.set_ylim(np.min(resultados[:, 0]) - 10, np.max(resultados[:, 0]) + 10)
eixo.set_xlabel('Tempo (ms)')
eixo.set_ylabel('Vm (mV)')
eixo.set_title('Evolução do Potencial do Neurônio')
linha, = eixo.plot([], [], lw=2)

def atualizar(frame):
    linha.set_data(tempo[:frame], resultados[:frame, 0])
    return linha,

ani = animation.FuncAnimation(figura, atualizar, frames=len(tempo), interval=10, blit=True)
ani.save('potencial_neuronio.gif', writer='pillow')

print("GIF salvo como 'potencial_neuronio.gif'.")
