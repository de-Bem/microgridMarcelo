import casadi as ca
import numpy as np

def Casadi_MPC(ne, ns, nx, N, Nc, A, B, C, Qu, Qx, Qdu, u_min, u_max, du_min, du_max, y_min, y_max, x_ref):


    # Inicializa variáveis
    u = ca.MX.sym('u', ne, N+2)        # Ação de controle (inputs)
    du = ca.MX.sym('du', ne, N+1)      # Variação de controle
    x = ca.MX.sym('x', nx+ne, N+2)     # Estados do sistema 
    y = ca.MX.sym('y', ns, N+1)        # Saídas
    Pbat = ca.MX.sym('Pbat', 1, N+1)   # Potência da bateria
    dsoc2 = ca.MX.sym('dsoc2', 1, N+1) # Variável binária para controle do estado de carga
    dgrid = ca.MX.sym('dgrid', 1, N+1) # Variável binária para controle da rede elétrica
    X0 = ca.SX.sym('X0', nx+ne)        # Estado inicial

    e = 1e-10                          # Um pequeno valor epsilon para evitar divisão por zero

    # Define a função objetivo e as restrições 
    objective = 0
    constraints = []

    constraints.append(x[:,1] == X0)
    # Monta QP
    for k in range(1, N+1):
        # Função objetivo
        if k > Nc:  # Após o horizonte de controle Nc
            objective += ca.mtimes([(x[0:nx-2, k] - x_ref).T, Qx, (x[0:nx-2, k] - x_ref)]) + \
                         ca.mtimes([x[nx:nx+2, k].T, Qu, x[nx:nx+2, k]])
        else:  # Dentro do horizonte de controle
            objective += ca.mtimes([(x[0:nx-2, k] - x_ref).T, Qx, (x[0:nx-2, k] - x_ref)]) + \
                         ca.mtimes([x[nx:nx+2, k+1].T, Qu, x[nx:nx+2, k+1]]) + \
                         ca.mtimes([du[0:2, k].T, Qdu, du[0:2, k]])

        # Dinâmica do sistema
        constraints.append(x[:, k+1] == ca.mtimes(A, x[:, k]) + ca.mtimes(B, du[:, k]))
        constraints.append(y[:, k] == ca.mtimes(C, x[:, k]))
        
        ### Ver certinho esse fmax e fmin

        # Restrições de estados 
        constraints.append(ca.fmin(y_min, x[0:nx-2, k]) == y_min)  
        constraints.append(ca.fmax(y_max, x[0:nx-2, k]) == y_max)

        # Restrições para variações de controle
        constraints.append(ca.fmax(np.hstack((du_min,du_min)), du[:, k]) == du[:, k])
        constraints.append(ca.fmin(np.hstack((du_max,du_max)), du[:, k]) == du[:, k])

        ### Ver se não é nx+1

        # Restrições no controle 
        constraints.append(ca.fmax(u_min, x[nx:nx+2, k+1]) == x[nx:nx+2, k+1])
        constraints.append(ca.fmin(u_max, x[nx:nx+2, k+1]) == x[nx:nx+2, k+1])

        # Dinamica da bateria 1
        constraints.append(Pbat[0, k] == -x[nx-1, k] - x[nx-2, k] - x[nx, k+1] - x[nx+1, k+1])

        # Restrições na potência da bateria
        constraints.append(Pbat[0, k] <= 0.4)
        constraints.append(Pbat[0, k] >= -0.4)

        # Definição de Delta SOC2 (estado de carga)
        constraints.append(x[nx, k+1] <= u_max[0] * (1 - dsoc2[0, k]) + e)
        constraints.append(u_min[0] * dsoc2[0, k] + e <= x[nx, k+1])

        # Definição de zSoc2
        constraints.append(x[nx+2, k+1] <= u_max[0] * dsoc2[0, k])
        constraints.append(u_min[0] * dsoc2[0, k] <= x[nx+2, k+1])
        constraints.append(x[nx+2, k+1] <= x[nx, k+1] - u_min[0] * (1 - dsoc2[0, k]))
        constraints.append(x[nx, k+1] - u_max[0] * (1 - dsoc2[0, k]) <= x[nx+2, k+1])

        # Definição de Delta grid
        constraints.append(x[nx+1, k+1] <= u_max[1] * (1 - dgrid[0, k]) + e)
        constraints.append(u_min[1] * dgrid[0, k] + e <= x[nx+1, k+1])

        # Definição de zgrid (essa variavel nem era necessaria, pq restringi para o
        # grid ser sempre maior que 0 (no Brasil não se pode vender pra rede),
        # mas deixo aqui como exemplo se quiserem implementar algo no futuro

        constraints.append(x[nx+3, k+1] <= u_max[1] * dgrid[0, k])
        constraints.append(u_min[1] * dgrid[0, k] <= x[nx+3, k+1])
        constraints.append(x[nx+3, k+1] <= x[nx+1, k+1] - u_min[1] * (1 - dgrid[0, k]))
        constraints.append(x[nx+1, k+1] - u_max[1] * (1 - dgrid[0, k]) <= x[nx+3, k+1])
    oo = ca.vertcat(ca.reshape(du, -1, 1),  ca.reshape(u, -1, 1), 
                                            ca.reshape(x, -1, 1), ca.reshape(Pbat, -1, 1), ca.reshape(y, -1, 1),
                                            ca.reshape(dsoc2, -1, 1) , ca.reshape(dgrid, -1, 1)),
    
    print(oo)
    # Montagem da função de otimização
    prob = {'f': objective, 'x': ca.vertcat(ca.reshape(du, -1, 1),  ca.reshape(u, -1, 1), 
                                            ca.reshape(x, -1, 1), ca.reshape(Pbat, -1, 1), ca.reshape(y, -1, 1),
                                            ca.reshape(dsoc2, -1, 1) , ca.reshape(dgrid, -1, 1)), 
            'x0': ca.reshape(X0, -1, 1), 'g': ca.vertcat(*constraints)}

    # Definir as opções para o solver Bonmin
    opts = {'bonmin': {'print_level': 5}}
    
    # Solver Bonmin
    solver = ca.nlpsol('solver', 'bonmin', prob, opts)

    # Função de solução otimizada
    return solver