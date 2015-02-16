% This is not a independent piece of script. 
% It must be run after HW2P2.m
% Or else the result data must be imported with corresponding .dat file.

kappa = 1e-2; % Thermal conductivity
rho = 10; % Density of fluid
C_p = 100; % Thermal Capacity of the fluid
q_0 = 10; % Heat generation rate

U_ = U;
for i = 1:Cell2
    U_(i,:) = U_(i,:) * Delta_x_2;
    U_(:,i) = U_(:,i) * Delta_x_2;
end
for i = (Cell2+1):Cell
    U_(i,:) = U_(i,:) * Delta_x_1;
    U_(:,i) = U_(:,i) * Delta_x_1;
end
U_int = sum(sum(U_));

q = q_0 * (ones(Cell) + sin(pi * Location_x / L) + sin(pi * Location_y / L));

q_ = q;
for i = 1:Cell2
    q_(i,:) = q_(i,:) * Delta_x_2;
    q_(:,i) = q_(:,i) * Delta_x_2;
end
for i = (Cell2+1):Cell
    q_(i,:) = q_(i,:) * Delta_x_1;
    q_(:,i) = q_(:,i) * Delta_x_1;
end
q_int = sum(sum(q_));

dTdx = q_int / U_int / rho / C_p;

% Governing Equation
AE = zeros(Cell);
AW = zeros(Cell);
AS = zeros(Cell);
AN = zeros(Cell);

for i = 2:Cell2
    AS(i,:) = kappa * Delta_x_2/ delta_x_2;
    AW(:,i) = kappa * Delta_x_2/ delta_x_2;
end
AS(Cell2+1,:) = kappa * Delta_x_1/ delta_x_c;
AW(:,Cell2+1) = kappa * Delta_x_1/ delta_x_c;
for i = (Cell2+2):Cell
    AS(i,:) = kappa * Delta_x_1/ delta_x_1;
    AW(:,i) = kappa * Delta_x_1/ delta_x_1;
end

for i = 1:(Cell2-1)
    AN(i,:) = kappa * Delta_x_2/ delta_x_2;
    AE(:,i) = kappa * Delta_x_2/ delta_x_2;
end
AN(Cell2,:) = kappa * Delta_x_2/ delta_x_c;
AE(:,Cell2) = kappa * Delta_x_2/ delta_x_c;
for i = (Cell2+1):(Cell-1)
    AN(i,:) = kappa * Delta_x_1/ delta_x_1;
    AE(:,i) = kappa * Delta_x_1/ delta_x_1;
end

B = q_ - U_ * rho * C_p * dTdx;

T = 100*ones(Cell);

AP = AE + AW + AS + AN;
AP(1:Cell2,Cell2+1) = AP(1:Cell2,Cell2+1) - kappa * Delta_x_1/ delta_x_c;
AP(Cell2+1,1:Cell2) = AP(Cell2+1,1:Cell2) - kappa * Delta_x_1/ delta_x_c;
AP(1:Cell2,1:Cell2) = 1;
AE(1:Cell2,1:Cell2) = 0;
AW(1:Cell2,1:Cell2+1) = 0;
AS(1:Cell2+1,1:Cell2) = 0;
AN(1:Cell2,1:Cell2) = 0;
B(1:Cell2,1:Cell2) = T(1:Cell2,1:Cell2);

% AP(1:Cell2,Cell2+1) = AP(1:Cell2,Cell2+1); 
% AP(Cell2+1,1:Cell2) = AP(Cell2+1,1:Cell2);
% AP(1:Cell2,1:Cell2) = 1;
% AE(1:Cell2,1:Cell2) = 0;
% AW(1:Cell2,1:Cell2) = 0;
% AS(1:Cell2,1:Cell2) = 0;
% AN(1:Cell2,1:Cell2) = 0;
% B(1:Cell2,1:Cell2) = 0;

% AP(Cell2+5,Cell2+5) = 1;
% AW(Cell2+5,Cell2+5) = 0;
% AE(Cell2+5,Cell2+5) = 0;
% AN(Cell2+5,Cell2+5) = 0;
% AS(Cell2+5,Cell2+5) = 0;
% B(Cell2+5,Cell2+5) = T(Cell2+5,Cell2+5);


iter = 1;
T0 = T;
alpha = 0.3; 
T = lblTDMA(AP,AW,AE,AS,AN,T,B);
% Under-Relaxation
while max(max(abs(T-T0))) > 0.0001
    mod = ones(Cell);
    mod(1:Cell2,1:Cell2) = 0;
    T = T - (T(Cell,Cell) - 100)*mod;
    T0 = T;
    T = lblTDMA(AP / alpha,AW,AE,AS,AN,T,B+(1-alpha)/alpha*AP.*T0);
    iter = iter+1;
    T(Cell2+5,Cell2+5)  
end

T_b_int1 = U.*T;
for i = 1:Cell2
    T_b_int1(i,:) = T_b_int1(i,:) * Delta_x_2;
    T_b_int1(:,i) = T_b_int1(:,i) * Delta_x_2;
end
for i = (Cell2+1):Cell
    T_b_int1(i,:) = T_b_int1(i,:) * Delta_x_1;
    T_b_int1(:,i) = T_b_int1(:,i) * Delta_x_1;
end
uT_int = sum(sum(T_b_int1));
T_b = uT_int / U_int;
T_n = (T - T_b) / (q_0 * D_h^2 / kappa);