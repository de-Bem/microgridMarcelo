import numpy as np

# Configurações iniciais
ts = 36                # Tempo de amostragem em segundos
SOC1_INI = 62          # Armazenamento da Bateria 1 (%)
SOC2_INI = 55          # Armazenamento da Bateria 2 (%)

Ndisbat1 = 2.8087 / 2  # Fator discarga bateria 1
Nchbat1 = 2.8087       # Fator carga bateria 1
Ndisbat2 = 1.2         # Fator discarga bateria 2
Nchbat2 = 2            # Fator discarga bateria 2

normaliza = 30 * 60  # O modelo do bordon vem com todas variaveis dividindo por 
#60 para passar de min para segundos. O modelo que tirei tava ajustado com 
#um tempo de amostragem de 30 segundos. Desta forma, é divido o modelo por
#esta constante 'normaliza' para normalizar, e depois multiplo pelo 'ts'
#para colocar no tempo de amostragem utilizado. Você poderia simplemente
#apagar 'normaliza' e reajustar os valores de fator de discarga e carga.

# Matrizes do modelo
microgrid_A = np.eye(2)
microgrid_B = np.array([[Nchbat1, Nchbat1, (Ndisbat1 - 2.8087), (Ndisbat1 - Nchbat1)],
                        [-Ndisbat2, 0, (-Ndisbat2 + Nchbat2), 0]]) * ts / normaliza
microgrid_C = np.eye(2)
microgrid_D = np.zeros((4, 4))
microgrid_E = np.array([[2.8087, (2.8087 / 2)],
                        [0, 0]]) * ts / normaliza

# Parâmetros do controle
Np = 8  # Horizonte de predição
Nu = 2  # Horizonte de controle

# Fatores de ponderação
delta = np.array([1e-7, 1e-7])   # Fator de ponderação do seguimento de referência
lambda_ = np.array([1e-4, 1e-4]) # Fator de ponderação do esforço do incremento de controle
alpha = np.array([1e-4, 4e-4])   # Fator de ponderação do esforço de controle

# Restrições da saída
Ymax = np.array([75, 80])
Ymin = np.array([40, 30])

# Restrições do controle
Umax = np.array([4, 6])
Umin = np.array([-4, 0])

# Restrições do incremento de controle
DeltaUmax = np.array([1, 1])
DeltaUmin = np.array([-1, -1])
Pbatmax = np.array([0.6])
Pbatmin = np.array([-0.6])
DeltaPbatmax = np.array([0.2])
DeltaPbatmin = np.array([-0.2])

# Referência
ref = np.array([80, 60])

# Parâmetros de Simulação
tf = 86400                           # Tempo de simulação  [s]
n = int(tf / ts)                     # Número de amostras  [-]
t = np.arange(0, (tf -ts), ts)       # Tempo               [s]
