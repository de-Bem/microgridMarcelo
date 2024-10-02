import casadi as ca
import numpy as np
import matplotlib.pyplot as plt
from classes import*
from Casadi_MPC import*
from setup import*

# Botando tudo na classe microgrid
microgrid = Microgrid(A=microgrid_A, 
                      B=microgrid_B,
                      C=microgrid_C, 
                      D=microgrid_D,
                      E=microgrid_E)

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

controlador=Casadi_MPC(ne,ny,nx,controller.Np,controller.Nu, M,N,Q,Qu,Qx,Qdu,controller.Umin,controller.Umax, controller.DeltaUmin,controller.DeltaUmax,controller.Ymin,controller.Ymax,controller.ref);

#################################################

# Simulação 

# Carregar o arquivo .mat
mat_data = scipy.io.loadmat('DemandaArray.mat')
Demanda = mat_data['Demanda']

mat_data = scipy.io.loadmat('Pot_solar_Array.mat')
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

    sol = solver(x0=xam)

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