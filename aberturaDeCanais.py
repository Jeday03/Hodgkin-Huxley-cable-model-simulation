import numpy as np
import matplotlib.pyplot as plt

# Parâmetros do modelo
capacitancia_membrana = 1.0  # Capacitância da membrana (µF/cm²)
condutancia_Na = 120.0  # Condutância máxima de Na+ (mS/cm²)
condutancia_K = 36.0  # Condutância máxima de K+ (mS/cm²)
condutancia_vazamento = 0.3  # Condutância de vazamento (mS/cm²)
potencial_Na = 50.0  # Potencial de reversão de Na+ (mV)
potencial_K = -77.0  # Potencial de reversão de K+ (mV)
potencial_vazamento = -54.387  # Potencial de reversão do vazamento (mV)

# Função para o estímulo
def estimulo_corrente(tempo):
    return 10.0 if 10 <= tempo <= 20 else 0.0  # Estímulo em µA/cm²

# Funções auxiliares para os gates
def taxa_alpha_n(volt):
    return 0.01 * (volt + 55) / (1 - np.exp(-(volt + 55) / 10))

def taxa_beta_n(volt):
    return 0.125 * np.exp(-(volt + 65) / 80)

def taxa_alpha_m(volt):
    return 0.1 * (volt + 40) / (1 - np.exp(-(volt + 40) / 10))

def taxa_beta_m(volt):
    return 4.0 * np.exp(-(volt + 65) / 18)

def taxa_alpha_h(volt):
    return 0.07 * np.exp(-(volt + 65) / 20)

def taxa_beta_h(volt):
    return 1 / (1 + np.exp(-(volt + 35) / 10))

# Simulação do modelo
def simular_hh(duracao, passo_tempo):
    num_passos = int(duracao / passo_tempo)
    tempos = np.linspace(0, duracao, num_passos)

    # Condições iniciais
    potencial = -65.0  # Potencial inicial da membrana (mV)
    gate_n = 0.3177
    gate_m = 0.0529
    gate_h = 0.5961

    # Armazenamento dos resultados
    registro_potencial = np.zeros(num_passos)
    registro_I_Na = np.zeros(num_passos)
    registro_I_K = np.zeros(num_passos)
    registro_I_vazamento = np.zeros(num_passos)
    registro_I_estimulo = np.zeros(num_passos)

    # Iteração temporal
    for passo in range(num_passos):
        corrente_Na = condutancia_Na * gate_m**3 * gate_h * (potencial - potencial_Na)
        corrente_K = condutancia_K * gate_n**4 * (potencial - potencial_K)
        corrente_vazamento = condutancia_vazamento * (potencial - potencial_vazamento)
        corrente_estimulo = estimulo_corrente(tempos[passo])

        # Atualização do potencial de membrana
        delta_potencial = (corrente_estimulo - corrente_Na - corrente_K - corrente_vazamento) / capacitancia_membrana
        potencial += delta_potencial * passo_tempo

        # Atualização dos gates
        gate_n += (taxa_alpha_n(potencial) * (1 - gate_n) - taxa_beta_n(potencial) * gate_n) * passo_tempo
        gate_m += (taxa_alpha_m(potencial) * (1 - gate_m) - taxa_beta_m(potencial) * gate_m) * passo_tempo
        gate_h += (taxa_alpha_h(potencial) * (1 - gate_h) - taxa_beta_h(potencial) * gate_h) * passo_tempo

        # Armazenamento dos valores
        registro_potencial[passo] = potencial
        registro_I_Na[passo] = corrente_Na
        registro_I_K[passo] = corrente_K
        registro_I_vazamento[passo] = corrente_vazamento
        registro_I_estimulo[passo] = corrente_estimulo

    return tempos, registro_potencial, registro_I_Na, registro_I_K, registro_I_vazamento, registro_I_estimulo

# Configuração da simulação
duracao_simulacao = 100.0  # Duração total (ms)
passo_simulacao = 0.001  # Passo de tempo (ms)

# Executa a simulação
resultados = simular_hh(duracao_simulacao, passo_simulacao)

# Plotagem dos resultados
plt.figure(figsize=(12, 8))

# Potencial de membrana
plt.subplot(2, 1, 1)
plt.plot(resultados[0], resultados[1], label="Potencial (mV)")
plt.title("Modelo de Hodgkin-Huxley & Gráficos de Abertura dos Canais Ionicos e Relação Entre Eles")
plt.ylabel("Potencial de Membrana (mV)")
plt.legend()

# Correntes iônicas
plt.subplot(2, 1, 2)
plt.plot(resultados[0], resultados[2], label="Corrente Na+ (µA/cm²)")
plt.plot(resultados[0], resultados[3], label="Corrente K+ (µA/cm²)")
plt.plot(resultados[0], resultados[4], label="Corrente de Vazamento (µA/cm²)")
plt.plot(resultados[0], resultados[5], label="Corrente de Estímulo (µA/cm²)", linestyle='--')
plt.xlabel("Tempo (ms)")
plt.ylabel("Corrente (µA/cm²)")
plt.legend()

plt.tight_layout()
plt.show()
