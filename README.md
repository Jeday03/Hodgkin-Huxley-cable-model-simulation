# Modelo Hodgkin-Huxley em Cabo

## Visão Geral
Este projeto implementa o modelo Hodgkin-Huxley (HH) para a propagação do potencial de ação em um cabo unidimensional. O objetivo é estudar o comportamento dos potenciais de ação em fibras nervosas utilizando simulações numéricas. O sistema gera resultados em formatos como gráficos estáticos, animações e arquivos CSV.

![propagacao_potencial_HH_cabo](https://github.com/user-attachments/assets/c299a7f6-f4d2-4579-ac81-7424745de471)

---

## Estrutura do Código

### 1. **main.py**
Simula o modelo Hodgkin-Huxley ao longo de um cabo unidimensional.

#### Recursos Principais:
- **Simulação do Modelo HH:** Calcula o potencial de membrana ao longo do tempo e da posição.
- **Geração de Arquivo CSV:** Salva os potenciais de membrana para cada ponto do cabo e instante de tempo.
- **Criação de Animação (GIF):** Mostra a propagação do potencial de ação.

#### Execução:
```bash
python main.py config.txt
```
- **Arquivo de Configuração:** Define parâmetros como tamanho do cabo, passo temporal, estímulo aplicado, etc.

---

### 2. **aberturaDeCanais.py**
Explora a dinâmica dos canais iônicos envolvidos no modelo HH, analisando variáveis como condutâncias de sódio (Na⁺) e potássio (K⁺).

#### Recursos Principais:
- **Gráficos das Condutâncias:** Mostra como os canais reagem ao potencial de membrana.
- **Simulação em Pontos Focais:** Foca na dinâmica temporal em um único ponto do cabo.
- 
![Figure_3](https://github.com/user-attachments/assets/87f576ce-913b-4c70-b485-fa84177f5ef5)

---

### 3. **grafico_hh_instantes.py**
Gera gráficos estáticos do potencial de membrana para diferentes instantes ao longo do cabo.

![Figure_1](https://github.com/user-attachments/assets/dc3251a5-daab-41f7-b44e-e82cf25e7b27)

#### Recursos Principais:
- **Gráficos Comparativos:** Apresenta o estado da membrana em quatro momentos distintos durante a simulação.
- **Visualização Detalhada:** Ajuda a compreender como o potencial de ação se propaga ao longo do cabo.
![Figure_2](https://github.com/user-attachments/assets/f5723b1d-70bd-4efe-a79c-68e95c8c7aca)

---

## Parâmetros do Modelo

### Biofísicos
- **Capacitância da membrana (C_m):** 1.0 µF/cm²
- **Condutâncias máximas:**
  - Sódio (g_Na): 120.0 mS/cm²
  - Potássio (g_K): 36.0 mS/cm²
  - Vazamento (g_L): 0.3 mS/cm²
- **Potenciais de reversão:**
  - Sódio (E_Na): 50.0 mV
  - Potássio (E_K): -77.0 mV
  - Vazamento (E_L): -54.387 mV

### Espaciais e Temporais
- **Raio do axônio (a):** 0.1 cm
- **Resistência longitudinal (R):** 100.0 ohm·cm
- **Espaçamento espacial (dx):** 0.01 cm
- **Passo de tempo (dt):** 0.01 ms
- **Comprimento total (L):** 1.0 cm
- **Tempo total de simulação (T):** 50.0 ms

---

## Dependências
Certifique-se de que as seguintes bibliotecas estão instaladas:
- **numpy**
- **matplotlib**
- **pandas**
- **configparser**
- **sys**




---

## Como Executar
1. Clone o repositório:
   ```bash
   git clone <url-do-repositorio>
   ```
2. Configure o arquivo `config.txt` com os parâmetros desejados.
3. Execute o arquivo principal:
   ```bash
   python main.py config.txt
   ```
4. Para análises complementares, utilize os scripts `aberturaDeCanais.py` e `grafico_hh_instantes.py`.

---

## Resultados

### Saídas Geradas:
- **CSV:** Contém os potenciais de membrana para cada posição ao longo do cabo em cada instante de tempo.
- **Gráficos:**
  - Propagação do potencial ao longo do cabo.
  - Dinâmica dos canais iônicos (Na⁺, K⁺).
- **GIF:** Demonstra a propagação do potencial de ação ao longo do tempo.

---
