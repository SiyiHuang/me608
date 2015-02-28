% This script is written to solve HW3P1. 
% The problem to solve is the non dimensionalized equation. The parameters are inputted at beginning, then non dimensionalized.
% The original problem will be reformed after the solution. 

%k_1 = 20; % Conduction Coef of Material 1
%k_2 = 80; % Conduction Coef of Material 2
%rho_1 = 3e3; % Density of M 1
%rho_2 = 2e3; % Density of M 2
%Cp_1 = 400; % Specific Heat of Material 1
%Cp_2 = 400; % Specific Heat of Material 2
%D = 1; % Dimension of the Geometry
%S = 10; % Heat Generation Rate
%Cell = 11; % Do not make change of this! The mesh configuration manually generated in this case. 

alpha12 = 1; % Ratio of thermal diffusivity
k12 = 0.25; % Ratio of conductivity

% Mesh Generation : 11*11 control volumes only. 
% 5*5 for region 2, else for region 1.

delta_t = 1e-2; % Define Time Step; Note that the estimate convergence time is 1.

Loc_x = zeros(11);
for i = 1:3
Loc_x(:,i) = 1/24*(2*i-1);
end
for i = 4:8
Loc_x(:,i) = 0.05*(2*(i-3)-1) + 0.25;
end
for i = 9:11
Loc_x(:,i) = 1/24*(2*(i-8)-1) + .75;
end
Loc_y = Loc_x';

AE = ones(11);
AE(:,11) = 0;
AE([1:3,9:11],[3,8]) = 10/11;
AE([1:3,9:11],4:7) = 10/12;
AE(4:8,[3,8]) = 12/11 * 2/(1+alpha12);
AE(4:8,4:7) = 2/(1+alpha12);
AE(4:8,[1,2,9,10]) = 1.2;
AW = ones(11);
AW(:,1) = 0;
AW([1:3,9:11],[4,9]) = 10/11;
AW([1:3,9:11],5:8) = 10/12;
AW(4:8,[4,9]) = 12/11 * 2/(1+alpha12);
AW(4:8,5:8) = 2/(1+alpha12);
AW(4:8,[2,3,10,11]) = 1.2;
AS = AW';
AN = AE';
AP0 = ones(11) /120;
AP0([1:3,9:11],[1:3,9:11]) = 1/144;
AP0(4:8,4:8) = 0.01;
AP0 = AP0 / delta_t;
AP = AE+AW+AN+AS+AP0;
b = zeros(11);
b(4:8,4:8) = k12/alpha12;

iter = 1;
ITER = iter;
Iter = 1;
T = ones(11);
T_c = T(6,6);
T0 = T;
T0(1,1) = 2*T(1,1); %a random T0 to walk in the loop
while abs(max(max(T-T0)) - min(min(T-T0)))> 1e-3
    T0 = T;
    T_ = 2*T; %a random T0 to walk in the loop
while max(max(abs(T-T_))) > 1e-3
    T_ = T;
    T = lblTDMA(AP,AW,AE,AS,AN,T,b+AP0.*T0);
    iter = iter + 1;
end
T_c = [T_c, T(6,6)];
ITER = [ITER,iter];
Iter = Iter + 1;
end

figure(1);
mesh(Loc_x,Loc_y,T);
xlabel('x');
ylabel('y');
zlabel('T*');
figure(2);
plot(delta_t*(0:Iter-1),T_c);
xlabel('\tau');
ylabel('T*_c');

