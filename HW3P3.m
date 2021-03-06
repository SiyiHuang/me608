% This script is to solve a 1D transient diffusion with Dirichlet boundary
% condition. 
% The schemes will be generalized represented in this case, i.e. a factor f
% will be involved. and f = 0,1,0.5 will represent explicit, implicit and
% Crank-Nicholson respectly. 

k = 150; % Thermal Conductivity
C_p = 700; % Specific Heat
rho = 2300; % Density
alpha = k/C_p/rho; % Thermal Diffusivity
L = 1; % Size of 1D region. 
T_i = 300; % Initial Temperature
T_w = 500; % Boundary Condition: Wall Temperature

cellNum = 100; % Cell Number

d_x = L / cellNum; % dist btwn Ctrs / size of cell
d_t = L^2 / alpha / 100; % Time Step
f = 1;

AE = ones(cellNum,1) * alpha / d_x;
AW = AE;
AE(cellNum) = 0;
AW(1) = 0;
AP0 = d_x/d_t;
AP = f * (AE + AW) + AP0;
AP([1,cellNum]) = AP([1,cellNum]) + alpha / d_x / 2;
b = zeros(cellNum,1);
b([1,cellNum]) = alpha / d_x / 2 * T_w;

iter = 1;
ITER = iter;
Iter = 1;
T = T_i*ones(cellNum,1);
Theta = 1;
Theta0 = 0; % Random walkin
THETA = Theta;
Nu = cellNum * 2;
while abs(Theta - Theta0) > 1e-5;
   Theta0 = Theta;
   TE = [T(2:cellNum);0];
   TW = [0;T(1:cellNum-1)];
   b_ = b + (1-f)*(AE.*TE + AW.*TW) + (AP0 - (1-f) * (AE+AW)).*T;
   T = TDMA(AP,AW*f,AE*f,T,b_);
   T_avg = sum(T)/cellNum;
   Theta = (T(fix(cellNum/2)) - T_w) / (T_avg - T_w);
   THETA = [THETA;Theta];
   Nu = [Nu; cellNum * 2 * (T_w - T(1)) / (T_w - T_avg)];
end

t = d_t * (1:size(THETA));
tau = alpha * t/ L^2;
figure(1);
plot(tau,THETA);
xlabel('\tau=\alphat/L^2')
ylabel('\Theta')
figure(2);
plot(t,Nu);
xlabel('t')
ylabel('Nu')
