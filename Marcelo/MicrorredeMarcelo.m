%% Microgrid 

addpath("C:\gurobi1103\win64\matlab\")
addpath(genpath("C:\Users\Pedro\Documents\MATLAB\YALMIP-master\"))

clear all;
close all hidden;
close all;
clc;

bar=waitbar(0, 'Preparando simulação...', 'Name', 'Microgrid');


ts = 36;%Tempo de amostragem em segunds (Foi escolhido assim para bater com a amostragem da coleta de dados das energias renovaveis)

SOC1_INI = 62; %Armazenamento da Bateria 1 (%)
SOC2_INI = 55; %Armazenamento da Bateria 2 (%)

Ndisbat1=2.8087/2; %Fator discarga bateria 1
Nchbat1=2.8087; %Fator carga bateria 1
Ndisbat2=1.2; %Fator discarga bateria 2
Nchbat2=2; %Fator discarga bateria 2

normaliza=30*60; %O modelo do bordon vem com todas variaveis dividindo por 
%60 para passar de min para segundos. O modelo que tirei tava ajustado com 
%um tempo de amostragem de 30 segundos. Desta forma, é divido o modelo por
%esta constante 'normaliza' para normalizar, e depois multiplo pelo 'ts'
%para colocar no tempo de amostragem utilizado. Você poderia simplemente
%apagar 'normaliza' e reajustar os valores de fator de discarga e carga.

%Matrizes do modelo
microgrid.A = eye(2);
microgrid.B = [Nchbat1 Nchbat1 Ndisbat1-2.8087 Ndisbat1-Nchbat1;
               -Ndisbat2 0 -Ndisbat2+Nchbat2 0]*ts/normaliza;                
microgrid.C = eye(2);
microgrid.D = zeros(4);
microgrid.E = [2.8087 2.8087/2;
               0 0]*ts/normaliza;
           

%% Parametros do controle
% Horizontes de predição Np e de controle Nu
controller.Np = 8;
controller.Nu = 2;

%Fatores e ponderação
controller.delta = [1e-7 1e-7];         % - Fator de ponderação do seguimento de referência
controller.lambda = [1e-4 1e-4];    % - Fator de ponderação do esforço do incremento de controle
controller.alpha = [1e-4 4e-4];       % - Fator de ponderação do esforço de controle

% Restrições da saída
controller.Ymax = [75; 80];
controller.Ymin = [40; 30];

% Restrições do controle
controller.Umax = [4; 6];
controller.Umin = [-4; 0];

% Restrições do incremento de controle
controller.DeltaUmax = [1; 1];
controller.DeltaUmin = [-1; -1];
controller.Pbatmax=[0.6];
controller.Pbatmin=[-0.6];
controller.DeltaPbatmax=[0.2];
controller.DeltaPbatmin=[-0.2];

% Referência
controller.ref = [80; 60];        % [SOc1-%, SOC2-%]

%Matrizes Aumentadas para colocar a perturbação nos estados (Vide Livro do Bordon)

A = [microgrid.A microgrid.E; zeros(1,2) eye(1) 0;zeros(1,3) ones(1,1)]; 
B = [microgrid.B; zeros(2,4)];
C = [microgrid.C zeros(2,1) zeros(2,1)];


% Tamanho do processo
nx = size(A,2);         % Numero de estados
ny = size(C,1);         % Numero de saidas
ne = size(B,2);        % Numero de entradas de controle

%Matrizes Aumentadas para colocar o integrador no sistema (Vide Livro do Bordon)
M = [A B; zeros(ne,nx) eye(ne)];
N = [B; eye(ne)];
Q = [C zeros(ny,ne)];

%% MONTANDO AS MATRIZES de Ponderação
Qu=diag(controller.alpha);
Qx=diag(controller.delta);
Qdu=diag(controller.lambda);   

%% Montando o controlador com o Yalmip

controlador=Yalmip_MPC(ne,ny,nx,controller.Np,controller.Nu, M,N,Q,Qu,Qx,Qdu,controller.Umin,controller.Umax, controller.DeltaUmin,controller.DeltaUmax,controller.Ymin,controller.Ymax,controller.ref);

%% Parâmetros de Simulação
% Simulação
tf = 86400;                % Tempo de simulação  [s]
n  = tf/ts;              % Número de amostras  [ ]
t  = 0:ts:(tf-ts);       % Tempo               [s]    


%***********************************************************************    

%% Simulação 
load Demanda %Carregando a demanda (arquivo externo)
%Colocando o tamanho da demanda para bater com o tamanho das amostras da energia renovavel
j=1;
Demanda=[];
for i=1:27
    Demanda=[Demanda;ts_power.Data(i)/1000*ones(89,1)];
end

load P_sunny %Carregando Potencia renovavel (aqui é solar, mas poderia ser qualquer outra)
Pot_solar = ts_power.Data(:)*0.350/1000; %Potencia solar em kW

xaux     = zeros(4, n);   % Estados aumentados
x     = zeros(2, n);      % Estados
y     = zeros(2, n);      % Saídas
u     = zeros(4, n);      % Controle


x(:,1) = [SOC1_INI; SOC2_INI];
y(:,1) = [SOC1_INI; SOC2_INI];


microgrid.A = eye(2);
microgrid.B = [2.8087 2.8087 2.8087/2-2.8087 2.8087/2-2.8087;
               -1.2 0 -2+1.2 0]*ts/(30*60);                
microgrid.C = eye(2);
microgrid.D = zeros(4);
microgrid.E = [2.8087 2.8087/2;
               0 0]*ts/(30*60);
           
for k = 1:n
    Pnet_act=Pot_solar(k)-Demanda(k); %Calcula a perturbação (se tivesse potencia do vento seria Psolar+Pvento-Demanda)
    if Pnet_act>=0
        zpert=0;
    else
        zpert=Pnet_act;
        Pnet_act=0;
    end
    xaux(:,k)=[x(:,k);Pnet_act;zpert];
    if k==1
       xam=[xaux(:,k);0;0;0;0]; 
    else
        xam=[xaux(:,k);u(:,k-1)];
    end
    sol=controlador{xam};
    du=sol{2};
    if k==1
        u(:,k)=sol{2};
    else
        u(:,k) = u(:,k-1)+sol{2};
    end
    disp(u(:,k));
    Pbat(1,k)=-u(1,k)-u(2,k)-Pnet_act-zpert;  %A bateria 1 serve de compensador da rede (vide livro do Bordon) 
    x(:, k+1) = microgrid.A*x(:,k) + microgrid.B*u(:,k)+microgrid.E*[Pnet_act;zpert];
    y(:, k)   = microgrid.C*x(:,k);
    waitbar(k/n,bar, 'Simulando...');
end
close(bar);
close all hidden;

figure
hold on
plot(t,Pbat);
plot(t,Demanda(1:2400));
plot(t,Pot_solar(1:2400));
plot(t,u(1,:));
plot(t,u(2,:));
legend('Bateria 1','Demanda', 'Potência Solar', 'Bateria 2', 'Rede');

figure
hold on
plot(t,y(1,:));
plot(t,y(2,:));
legend('SOC1','SOC2');