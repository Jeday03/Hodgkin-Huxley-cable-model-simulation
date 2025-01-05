# Simulação do Modelo Hodgkin-Huxley em um Cabo

## 1. Introdução
Este código implementa a simulação do modelo Hodgkin-Huxley (HH) em um cabo utilizando as equações de condução de potenciais de ação. O modelo é baseado em um sistema unidimensional (“1D”), considerando as propriedades biofísicas da membrana neuronal.

## 2. Dependências
O código utiliza as seguintes bibliotecas do Python:
- **numpy**: Para cálculos matriciais e operações matemáticas.
- **matplotlib**: Para gerar gráficos e animações.
- **pandas**: Para manipulação e exportação de dados.
- **sys**: Para manipulação de argumentos de linha de comando.
- **configparser**: Para leitura de arquivos de configuração.

## 3. Estrutura do Código
O código está organizado nas seguintes funções principais:

### 3.1. `carregar_configuracoes(config_file)`
Carrega os parâmetros da simulação a partir de um arquivo de configuração no formato INI.
- **Parâmetro**:
  - `config_file` (str): Caminho para o arquivo de configuração.
- **Retorno**:
  - Dicionário contendo os parâmetros da simulação.

### 3.2. `I_ap(t, x, dx)`
Define a corrente de estímulo aplicada na membrana.
- **Parâmetros**:
  - `t` (float): Tempo.
  - `x` (float): Posição.
  - `dx` (float): Intervalo espacial.
- **Retorno**:
  - Corrente aplicada (float).

### 3.3. Funções auxiliares
Incluem funções para os valores de transição das variáveis do modelo HH:
- `alpha_n(V), beta_n(V)`
- `alpha_m(V), beta_m(V)`
- `alpha_h(V), beta_h(V)`

### 3.4. `hodgkin_huxley_1D(params)`
Executa a simulação do modelo HH no cabo, atualizando as variáveis ao longo do tempo e espaço.
- **Parâmetro**:
  - `params` (dict): Parâmetros da simulação.
- **Retorno**:
  - Matriz contendo o potencial de membrana em cada instante de tempo e posição.

### 3.5. `main()`
Ponto de entrada do programa.
- Lê o arquivo de configuração.
- Executa a simulação.
- Salva os resultados em um arquivo CSV.
- Gera uma animação da propagação do potencial de ação.

## 4. Arquivo de Configuração (`config.txt`)
O arquivo deve conter as seguintes seções e parâmetros:

```ini
[PARAMETROS]
C_m = 1.0
g_Na = 120.0
g_K = 36.0
g_L = 0.3
E_Na = 50.0
E_K = -77.0
E_L = -54.387

[CABO]
a = 0.1
R = 100.0

[SIMULACAO]
dx = 0.01
dt = 0.01
L = 1.0
T = 50.0
