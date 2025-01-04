import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import json
import sys

def modelo_cabo_hodgkin_huxley(arquivo_configuracao):
    with open(arquivo_configuracao, 'r') as arquivo:
        parametros = json.load(arquivo)

    # Parâmetros
    cm = parametros.get("cm", 1.0)  # Capacitância da membrana (uF/cm^2)
    a = parametros.get("a", 1.0)    # Raio do axônio (cm)
    rl = parametros.get("rl", 1.0)  # Resistência longitudinal (Ohm*cm)
    gna = parametros.get("gna", 120.0)  # Condutância de sódio (mS/cm^2)
    gk = parametros.get("gk", 36.0)  # Condutância de potássio (mS/cm^2)
    gl = parametros.get("gl", 0.3)  # Condutância de vazamento (mS/cm^2)
    ena = parametros.get("ena", 50.0)  # Potencial de reversão do sódio (mV)
    ek = parametros.get("ek", -77.0)  # Potencial de reversão do potássio (mV)
    el = parametros.get("el", -54.4)  # Potencial de reversão do vazamento (mV)
    tempo_total = parametros.get("T_max", 10.0)  # Tempo total de simulação (ms)
    comprimento_max = parametros.get("L_max", 1.0)  # Comprimento do axônio (cm)
    passo_tempo = parametros.get("dt", 0.01)  # Passo de tempo (ms)
    passo_espaco = parametros.get("dx", 0.01)  # Passo de espaço (cm)
    vm_inicial = parametros.get("vm0", -65.0)  # Potencial de membrana inicial (mV)
    m_inicial = parametros.get("m0", 0.05)
    h_inicial = parametros.get("h0", 0.6)
    n_inicial = parametros.get("n0", 0.32)
    corrente_externa = np.array(parametros["J"])
    segmentos_mielinicos = np.array(parametros["Mie"])

    numero_posicoes = int(comprimento_max / passo_espaco) + 1
    numero_tempos = int(tempo_total / passo_tempo) + 1

    potencial_membrana = np.full((numero_tempos, numero_posicoes), vm_inicial)  # Matriz de potencial de membrana
    m = np.full(numero_posicoes, m_inicial)
    h = np.full(numero_posicoes, h_inicial)
    n = np.full(numero_posicoes, n_inicial)

    # Funções auxiliares
    def alfa_m(V):
        return 0.1 * (V + 40) / (1 - np.exp(-(V + 40) / 10))

    def beta_m(V):
        return 4 * np.exp(-(V + 65) / 18)

    def alfa_h(V):
        return 0.07 * np.exp(-(V + 65) / 20)

    def beta_h(V):
        return 1 / (1 + np.exp(-(V + 35) / 10))

    def alfa_n(V):
        return 0.01 * (V + 55) / (1 - np.exp(-(V + 55) / 10))

    def beta_n(V):
        return 0.125 * np.exp(-(V + 65) / 80)

    def calcular_correntes_ionicas(V, m, h, n):
        corrente_na = gna * m**3 * h * (V - ena)
        corrente_k = gk * n**4 * (V - ek)
        corrente_vazamento = gl * (V - el)
        return corrente_na, corrente_k, corrente_vazamento

    # Simulação
    for t in range(1, numero_tempos):
        for x in range(1, numero_posicoes - 1):
            d2Vm_dx2 = (potencial_membrana[t-1, x+1] - 2 * potencial_membrana[t-1, x] + potencial_membrana[t-1, x-1]) / passo_espaco**2
            if segmentos_mielinicos[x] == 1:  # Segmento mielinizado
                d2Vm_dx2 /= 100
            
            corrente_na, corrente_k, corrente_vazamento = calcular_correntes_ionicas(potencial_membrana[t-1, x], m[x], h[x], n[x])
            corrente_ionica = corrente_na + corrente_k + corrente_vazamento

            potencial_membrana[t, x] = potencial_membrana[t-1, x] + passo_tempo / (cm * a) * (-corrente_ionica + rl * d2Vm_dx2 + corrente_externa[t, x])

            m[x] += passo_tempo * (alfa_m(potencial_membrana[t-1, x]) * (1 - m[x]) - beta_m(potencial_membrana[t-1, x]) * m[x])
            h[x] += passo_tempo * (alfa_h(potencial_membrana[t-1, x]) * (1 - h[x]) - beta_h(potencial_membrana[t-1, x]) * h[x])
            n[x] += passo_tempo * (alfa_n(potencial_membrana[t-1, x]) * (1 - n[x]) - beta_n(potencial_membrana[t-1, x]) * n[x])

    # Salvar resultados
    np.savetxt("tabela_resultados.csv", potencial_membrana, delimiter=",")

    # Criar GIF
    figura, eixo = plt.subplots()
    linha, = eixo.plot(np.linspace(0, comprimento_max, numero_posicoes), potencial_membrana[0])
    eixo.set_xlim(0, comprimento_max)
    eixo.set_ylim(-80, 60)
    eixo.set_xlabel("Posição (cm)")
    eixo.set_ylabel("Voltagem (mV)")

    def atualizar(frame):
        linha.set_ydata(potencial_membrana[frame])
        return linha,

    animacao = FuncAnimation(figura, atualizar, frames=numero_tempos, blit=True)
    animacao.save("voltagem_posicao.gif", fps=30)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python main.py config.txt")
    else:
        caminho_config = sys.argv[1]
        modelo_cabo_hodgkin_huxley(caminho_config)
