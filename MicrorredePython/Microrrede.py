import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm  # Para barra de progresso
from classe import*
import scipy.io
import casadi as ca
from Casadi_MPC import*

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

# Botando tudo na classe microgrid
microgrid = Microgrid(A=microgrid_A, 
                      B=microgrid_B,
                      C=microgrid_C, 
                      D=microgrid_D,
                      E=microgrid_E)

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

#Botando tudo na classe controller
controller = Controller(Np = Np,
                        Nu = Nu,
                        delta = delta,
                        lambda_ = lambda_,
                        alpha = alpha,
                        Ymax = Ymax,
                        Ymin = Ymin,
                        Umax = Umax,
                        Umin = Umin,
                        DeltaUmax = DeltaUmax,
                        DeltaUmin = DeltaUmin,
                        Pbatmax = Pbatmax,
                        Pbatmin = Pbatmin,
                        DeltaPbatmax = DeltaPbatmax,
                        DeltaPbatmin = DeltaPbatmin,
                        ref = ref)

# Matrizes Aumentadas
A = np.block([[microgrid.A, microgrid.E],
              [np.zeros((1, 2)), np.eye(1), np.zeros((1, 1))],
              [np.zeros((1, 3)), np.ones((1, 1))]])

B = np.block([[microgrid.B],
              [np.zeros((2, 4))]])

C = np.block([[microgrid.C, np.zeros((2, 1)), np.zeros((2, 1))]])

#Tamanho do processo
nx = A.shape[1]  # Número de estados
ny = C.shape[0]  # Número de saídas
ne = B.shape[1]  # Número de entradas de controle

# Matrizes Aumentadas para colocar o integrador no sistema (Vide Livro do Bordon)
M = np.block([[A, B],
              [np.zeros((ne, nx)), np.eye(ne)]])
N = np.block([[B],
              [np.eye(ne)]])
Q = np.block([[C, np.zeros((ny, ne))]])

# Matrizes de Ponderação
Qu = np.diag(alpha)
Qx = np.diag(delta)
Qdu = np.diag(lambda_)

#################################################

#controlador=Casadi_MPC(xam,ne,ny,nx,controller.Np,controller.Nu, M,N,Q,Qu,Qx,Qdu,controller.Umin,controller.Umax, controller.DeltaUmin,controller.DeltaUmax,controller.Ymin,controller.Ymax,controller.ref)

#################################################


# Parâmetros de Simulação
tf = 86400                           # Tempo de simulação  [s]
n = int(tf / ts)                     # Número de amostras  [-]
t = np.arange(0, (tf -ts), ts)       # Tempo               [s]

###############################################################

# Simulação 

# Carregar o arquivo .mat
mat_data = scipy.io.loadmat('MicrorredePython/DemandaArray.mat')
Demanda = mat_data['Demanda']

mat_data = scipy.io.loadmat('MicrorredePython/Pot_solar_Array.mat')
Pot_solar = mat_data['Pot_solar']

xaux = np.zeros((4, n))  # Estados aumentados
x = np.zeros((2, n))     # Estados
y = np.zeros((2, n))     # Saídas
u = np.zeros((4, n))     # Controle
Pbat = np.zeros((n))
print(Pbat.shape)
# Condições iniciais
x[:, 0] = [SOC1_INI, SOC2_INI]
y[:, 0] = [SOC1_INI, SOC2_INI]

# Simulação
for k in tqdm(range(n-1), desc="Simulando..."): 
    Pnet_act = Pot_solar[k] - Demanda[k] #Calcula a perturbação (se tivesse potencia do vento seria Psolar+Pvento-Demanda)
    if Pnet_act >= 0:
        zpert = 0
    else:
        zpert = Pnet_act
        Pnet_act = 0
    
    xaux[:,k] = np.hstack((x[:, k], Pnet_act, zpert))

    if (k==0):
       xam= np.hstack((xaux[:,k], 0, 0, 0, 0)) 
    else:
       xam= np.hstack((xaux[:,k], u[:,k-1]))
  
#############################################################

    sol = Casadi_MPC(xam,ne,ny,nx,controller.Np,controller.Nu, M,N,Q,Qu,Qx,Qdu,controller.Umin,controller.Umax, controller.DeltaUmin,controller.DeltaUmax,controller.Ymin,controller.Ymax,controller.ref)

############################################################
    du =  sol['x'][0:4]  # Variação de controle

    if k == 0:
        u[:, k] = du
    else:
        u[:, k] = u[:, k-1] + du

    Pbat[k] = -u[0, k] - u[1, k] - Pnet_act - zpert  # Bateria 1 serve de compensador da rede (vide livro do Bordon)
    x[:, k+1] = microgrid.A @ x[:, k] + microgrid.B @ u[:, k] + microgrid.E @ np.hstack((np.array(Pnet_act),np.array(zpert)))
    y[:, k] = microgrid.C @ x[:, k]

# Plotagem dos resultados
plt.figure()
plt.plot(t, Pbat[:len(t)], label='Bateria 1')
plt.plot(t, Demanda[:len(t)], label='Demanda')
plt.plot(t, Pot_solar[:len(t)], label='Potência Solar')
plt.plot(t, u[0, :-1], label='Bateria 2')
plt.plot(t, u[1, :-1], label='Rede')
plt.legend()
plt.show()

plt.figure()
plt.plot(t, y[0, :-1], label='SOC1')
plt.plot(t, y[1, :-1], label='SOC2')
plt.legend()
plt.show()