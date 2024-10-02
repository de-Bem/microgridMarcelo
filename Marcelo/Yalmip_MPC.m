function [yalmip] = Yalmip_MPC(ne,ns,nx,N,Nc,A,B,C,Qu,Qx,Qdu,u_min,u_max,du_min,du_max,y_min,y_max,x_ref)

% Inicializa Variaveis 
u = sdpvar(ne,N+2);
du = sdpvar(ne,N+1);
x = sdpvar(nx+ne,N+2);
y = sdpvar(ns,N+1);
Pbat=sdpvar(1,N+1);
dsoc2 = binvar(1,N+1);
dgrid = binvar(1,N+1);
dpert = binvar(1,N+1);
e=1e-10;

%Monta QP
constraints = [];    % Constraints
objective = 0;       % Objective function

for k = 2:N+1  
    if k>Nc+1
        objective = objective + (x(1:(nx-2),k) - x_ref)'*Qx*(x(1:(nx-2),k) - x_ref)+x((nx+1):end-2,k)'*Qu*x((nx+1):end-2,k);
    else
        objective = objective +  x((nx+1):end-2,k+1)'*Qu*x((nx+1):end-2,k+1) + (x(1:(nx-2),k) -x_ref)'*Qx*(x(1:(nx-2),k) - x_ref) + (du(1:2,k))'*Qdu*(du(1:2,k));
    end
    constraints = [constraints, x(:,k+1) == A*x(:,k) + B*du(:,k)];
    constraints = [constraints, y(:,k) == C*x(:,k)];
    constraints = [constraints, y_min <= x(1:(nx-2),k) <=  y_max];
    constraints = [constraints, [du_min;du_min] <= du(:,k) <= [du_max;du_max]];
    constraints = [constraints, u_min <= x(nx+1:nx+2,k+1) <= u_max];
    constraints =[constraints, Pbat(1,k)==-x(nx,k)-x(nx-1,k)-x(nx+1,k+1)-x(nx+2,k+1)];
    %constraints = [constraints, -0.4 <= Pbat(1,k)-Pbat(1,k-1) <= 0.4];   
    constraints = [constraints, -0.4 <= Pbat(1,k) <= 0.4];
    
    %Define Delta Soc2
    constraints = [constraints, x(nx+1,k+1)<=u_max(1,1)*(1-dsoc2(1,k))+e];
    constraints = [constraints, u_min(1,1)*dsoc2(1,k)+e<=x(nx+1,k+1)];
    
    %Define zSoc2
    constraints = [constraints, x(nx+3,k+1)<=u_max(1,1)*dsoc2(1,k)];
    constraints = [constraints, u_min(1,1)*dsoc2(1,k)<=x(nx+3,k+1)];
    constraints = [constraints, x(nx+3,k+1)<=x(nx+1,k+1)-u_min(1,1)*(1-dsoc2(1,k))];
    constraints = [constraints, x(nx+1,k+1)-u_max(1,1)*(1-dsoc2(1,k))<=x(nx+3,k+1)];
    
    %Define Delta grid
    constraints = [constraints, x(nx+2,k+1)<=u_max(2,1)*(1-dgrid(1,k))+e];
    constraints = [constraints, u_min(2,1)*dgrid(1,k)+e<=x(nx+2,k+1)];
    
    %Define zgrid (essa variavel nem era necessaria, pq restringi para o
    %grid ser sempre maior que 0 (no Brasil não se pode vender pra rede),
    %mas deixo aqui como exemplo se quiserem implementar algo no futuro
    constraints = [constraints, x(nx+4,k+1)<=u_max(2,1)*dgrid(1,k)];
    constraints = [constraints, u_min(2,1)*dgrid(1,k)<=x(nx+4,k+1)];
    constraints = [constraints, x(nx+4,k+1)<=x(nx+2,k+1)-u_min(2,1)*(1-dgrid(1,k))];
    constraints = [constraints, x(nx+2,k+1)-u_max(2,1)*(1-dgrid(1,k))<=x(nx+4,k+1)];  
end

% Monta Controlador
options = sdpsettings('verbose',1,'solver','gurobi'); 
Parameters = [x(:,2)];
solutions_out = {[objective],[du(:,2)],[Pbat(1,2)]};
yalmip = optimizer(constraints, objective,options,Parameters,solutions_out);



