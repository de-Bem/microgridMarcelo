import casadi as ca
import numpy as np


def Casadi_MPC(xam,ne, ns, nx, N, Nc, A, B, C, Qu, Qx, Qdu, u_min, u_max, du_min, du_max, y_min, y_max, x_ref):

    # Inicializa variáveis
    u = ca.MX.sym('u', ne, N+2)        # Ação de controle (inputs)
    du = ca.MX.sym('du', ne, N+1)      # Variação de controle
    x = ca.MX.sym('x', nx+ne, N+2)     # Estados do sistema 
    y = ca.MX.sym('y', ns, N+1)        # Saídas
    Pbat = ca.MX.sym('Pbat', 1, N+1)   # Potência da bateria
    dsoc2 = ca.MX.sym('dsoc2', 1, N+1) # Variável binária para controle do estado de carga
    dgrid = ca.MX.sym('dgrid', 1, N+1) # Variável binária para controle da rede elétrica

    P = ca.MX.sym('u', ne+ne+2, 1) # Parametros que serao passados para o solver [xam,x_ref]

    e = 1e-10                          # Um pequeno valor epsilon para evitar divisão por zero

    # Define a função objetivo e as restrições 
    objective = 0
    rest_1 = ca.MX.sym('rest_1',nx+ne,N)
    rest_2 = ca.MX.sym('rest_2',2,N)
    rest_3 = ca.MX.sym('rest_3',1,N)
    rest_4 = ca.MX.sym('rest_4',1,N)
    x[:,0] = P[0:nx+ne]

    # Monta QP
    for k in range(N):
        # Função objetivo
        if k > Nc + 1:  # Após o horizonte de controle Nc
            objective += (x[0:nx-2, k] - P[nx+ne:]).T @ Qx @ (x[0:nx-2, k] - P[nx+ne:]) + \
                         x[nx:nx+2, k].T @ Qu @ x[nx:nx+2, k]
        else:  # Dentro do horizonte de controle
            objective += (x[0:nx-2, k] - P[nx+ne:]).T @ Qx @ (x[0:nx-2, k] - P[nx+ne:]) + \
                         x[nx:nx+2, k+1].T @ Qu @ x[nx:nx+2, k+1] + \
                         du[0:2, k].T @ Qdu @ du[0:2, k]
            
        rest_1[:,k] = A @ x[:,k] + B @ du[:,k] - x[:,k+1]
        rest_2[:,k] = C @ x[:,k] - y[:,k]
        rest_3[:,k] = -x[nx,k]-x[nx-1,k]-x[nx+1,k+1]-x[nx+2,k+1] - Pbat[0,k]
        rest_4[:,k] = Pbat[0,k]

    g = ca.vertcat(
        rest_1,
        rest_2,
        rest_3,
        rest_4
    )

    OPT_variables = ca.vertcat(
    x.reshape((-1, 1)),
    du.reshape((-1, 1))
    )

    nlp_prob = {
    'f': objective,
    'x': OPT_variables,
    'g': g,
    'p': P
    }

    opts = {
    'ipopt': {
        'max_iter': 2000,
        'print_level': 0,
        'acceptable_tol': 1e-8,
        'acceptable_obj_change_tol': 1e-6
    },
    'print_time': 0
    }
    # Função de solução otimizada
    return