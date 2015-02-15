%This piece of script is intended to solve Problem 2 in Homework 2.
%The reduced size model is taken from the north east part of the annulus
%due to symmetry.

%Set Parameters
L = 0.1; % Outer Dimension of duct
l = 0.02; % Annulus thickness
C = L / 2 - l; % half size of center
mu = 1e-3; % Dynamic viscosity
Cell1 = 20;% Cell number across the annulus
Cell2 = 30;% Cell number across half size of center
Cell = Cell1+Cell2; % Cell number across half L
Delta_x_1 = l / Cell1; % Cell size across the annulus
Delta_x_2 = C / Cell2; % Cell size along center block
delta_x_1 = l / Cell1; % Centroid distance across the annulus
delta_x_2 = C / Cell2; % Centroid distance along center block
delta_x_b_1 = l / Cell1 / 2; % C Disntance at annulus boundary
delta_x_b_2 = C / Cell2 / 2; % C D at non annulus block
delta_x_c = (delta_x_1 + delta_x_2) / 2; % C dist at interface
dpdx = 100 ; %Pressure Gradient in x dir
% Mesh generation
Location_x = zeros(Cell);
Location_y = zeros(Cell);
for i = 1:(Cell)
    for j = 1:(Cell)
        if i <= Cell2
            Location_x(i,j) = i * Delta_x_2 - delta_x_b_2;
        else 
            Location_x(i,j) = C + (i-Cell2) * Delta_x_1 - delta_x_b_1;
        end
        if j <= Cell2
            Location_y(i,j) = j * Delta_x_2 - delta_x_b_2;
        else 
            Location_y(i,j) = C + (j-Cell2) * Delta_x_1 - delta_x_b_1;
        end
    end
end

U = ones(Cell);

% Governing Equation
AE = zeros(Cell);
AW = zeros(Cell);
AS = zeros(Cell);
AN = zeros(Cell);

for i = 2:Cell2
    AS(i,:) = mu * Delta_x_2/ delta_x_2;
    AW(:,i) = mu * Delta_x_2/ delta_x_2;
end
AS(Cell2+1,:) = mu * Delta_x_1/ delta_x_c;
AW(:,Cell2+1) = mu * Delta_x_1/ delta_x_c;
for i = (Cell2+2):Cell
    AS(i,:) = mu * Delta_x_1/ delta_x_1;
    AW(:,i) = mu * Delta_x_1/ delta_x_1;
end

for i = 1:(Cell2-1)
    AN(i,:) = mu * Delta_x_2/ delta_x_2;
    AE(:,i) = mu * Delta_x_2/ delta_x_2;
end
AN(Cell2,:) = mu * Delta_x_2/ delta_x_c;
AE(:,Cell2) = mu * Delta_x_2/ delta_x_c;
for i = (Cell2+1):(Cell-1)
    AN(i,:) = mu * Delta_x_1/ delta_x_1;
    AE(:,i) = mu * Delta_x_1/ delta_x_1;
end
AB2 = mu * Delta_x_1 / delta_x_b_2;
AB1 = mu * Delta_x_1 / delta_x_b_1;
B = dpdx * ones(Cell);
for i = 1:Cell2
    B(i,:) = B(i,:) * Delta_x_2;
    B(:,i) = B(:,i) * Delta_x_2;
end
for i = (Cell2+1):Cell
    B(i,:) = B(i,:) * Delta_x_1;
    B(:,i) = B(:,i) * Delta_x_1;
end

AP = AE + AW + AS + AN;

% Segragate Blocks and BC reset.
AP(1:Cell2,Cell2+1) = AP(1:Cell2,Cell2+1) + AB2;
AP(1:Cell2,Cell) = AP(1:Cell2,Cell) + AB2;
AP(Cell2+1,1:Cell2) = AP(Cell2+1,1:Cell2) + AB2;
AP(Cell,1:Cell2) = AP(Cell,1:Cell2) + AB2;
AP(Cell,(Cell2+1):Cell) = AP(Cell,(Cell2+1):Cell) + AB1;
AP((Cell2+1):Cell,Cell) = AP((Cell2+1):Cell,Cell) + AB1;

AP(1:Cell2,1:Cell2) = 1;
AE(1:Cell2,1:Cell2) = 0;
AW(1:Cell2,1:Cell2) = 0;
AS(1:Cell2,1:Cell2) = 0;
AN(1:Cell2,1:Cell2) = 0;
B(1:Cell2,1:Cell2) = 0;

iter = 1;
U0 = U;
U = lblTDMA(AP,AW,AE,AS,AN,U,B);
while max(max(abs(U-U0))) > 1e-5
    U0 = U;
    U = lblTDMA(AP,AW,AE,AS,AN,U,B);
    iter = iter+1;
end

% non-dimensionalize:
D_h = 4 * (L^2 - (2*C)^2) / (4*L + 8*C);
U_n = U * mu / (D_h^2 * dpdx);
