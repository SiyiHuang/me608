% THis is the script to solve Problem4 of ME608 HW4.

% Properties: 
rho = 2; % Density
D = 3; % Diffusion Constant
L = 3; % Geometry
phi_s = 0.5; % scalar generated from source. 
% Mesh Generation:
n_cell = 50; % Generating n_cell * n_cell uniform size
delta = L / n_cell; % Size of cell
x_ctr = zeros(n_cell);
y_ctr = zeros(n_cell);
x_w = zeros(n_cell);
y_s = zeros(n_cell);
x_e = zeros(n_cell);
y_n = zeros(n_cell);
for i = 1:n_cell
x_ctr(:,i) = delta * (2*i-1) / 2 ;
y_ctr(i,:) = delta * (2*i-1) / 2 ;
x_w(:,i) = delta * (i-1) ;
y_s(i,:) = delta * (i-1) ;
x_e(:,i) = delta * i;
y_n(i,:) = delta * i;
end
% XX Fields of mesh
u = x_ctr.*x_ctr + 1;
v = y_ctr.*y_ctr + 2;
u_w = x_w.^2 + 1;
u_e = x_e.^2 + 1;
v_s = y_s.^2 + 2;
v_n = y_n.^2 + 2;
S = 4 * (x_ctr + y_ctr);
S_phi = S * phi_s;
% Coef Matrix and source Matrix
F_e = rho*u_e*delta;
AE = D*ones(n_cell);
AE(:,n_cell) = 0;
F_w = rho*u_w*delta;
AW = D*ones(n_cell) + F_w;
AW(:,1) = 0;
F_n = rho*v_n*delta;
AN = D*ones(n_cell);
AN(n_cell,:) = 0;
F_s = rho*v_s*delta;
AS = D*ones(n_cell) + F_s;
AS(1,:) = 0;
AP = AE+AW+AN+AS+S*delta*delta;
AP(:,1) = AP(:,1) + D*2 + F_w(:,1); % Dirichlet Boundary Condition
AP(1,:) = AP(1,:) + D*2 + F_s(1,:);
B = S_phi * delta * delta;
B(:,1) = B(:,1)+0*ones(n_cell,1).*(F_s(:,1)+2*D); % Dirichlet Boundary Condition
B(1,:) = B(1,:)+1*ones(1,n_cell).*(F_s(1,:)+2*D);
% Initial guess of phi field:
phi = 0.5*ones(n_cell);

% Initialize tolerance
tol = 1;

% Start Solving
while tol > 1e-5
phi_old = phi;
% QUICK can be treated as a modification of UDS. 
% In implementation, an additional term with previous phi is placed in
% source term. 
% Velocity field is given. For e face, W,P,E are involved, while for w face, 
% WW,W,P are involved.
% Similarly for s and n. 
% Thus, first two, and last row/collumn are not applicable for corresponding
% QUICK scheme. Upwind is used instead.
B_s = B; 
B_s(3:(n_cell),:) = B(3:(n_cell),:) + F_w(3:(n_cell),:)...
.*( (phi(3:(n_cell),:) + phi(2:(n_cell-1),:))/2 - ...
(phi(3:(n_cell),:) + phi(1:(n_cell-2),:) - 2 * phi(2:(n_cell-1),:)) / 8 ...
- phi(2:(n_cell-1),:) );
B_s(2:(n_cell-1),:) = B(2:(n_cell-1),:) + F_e(2:(n_cell-1),:)...
.*( (phi(3:(n_cell),:) + phi(2:(n_cell-1),:))/2 - ...
(phi(3:(n_cell),:) + phi(1:(n_cell-2),:) - 2 * phi(2:(n_cell-1),:)) / 8 ...
- phi(2:(n_cell-1),:) );
B_s(:,3:(n_cell)) = B(:,3:(n_cell)) + F_s(:,3:(n_cell))...
.*( (phi(:,3:(n_cell)) + phi(:,2:(n_cell-1)))/2 - ...
(phi(:,3:(n_cell)) + phi(:,1:(n_cell-2)) - 2 * phi(:,2:(n_cell-1))) / 8 ...
- phi(:,2:(n_cell-1)) );
B_s(:,2:(n_cell-1)) = B(:,2:(n_cell-1)) + F_n(:,2:(n_cell-1))...
.*( (phi(:,3:(n_cell)) + phi(:,2:(n_cell-1)))/2 - ...
(phi(:,3:(n_cell)) + phi(:,1:(n_cell-2)) - 2 * phi(:,2:(n_cell-1))) / 8 ...
- phi(:,2:(n_cell-1)) );

phi = lblTDMA(AP, AW, AE, AS, AN, phi, B_s);
tol = max(max(abs(phi-phi_old)));
% mesh(phi_old);  % This two lines are used for debugging
% pause;
end
ctr_lx = delta:delta:L;
ctr_l1 = (phi(n_cell/2,:)+phi(n_cell/2+1,:))/2;
ctr_l2 = (phi(:,n_cell/2)+phi(:,n_cell/2+1))/2;
plot(ctr_lx,ctr_l1,ctr_lx,ctr_l2);
xlabel('L');
ylabel('\phi (50x50)');
legend('hrz ctr line','vtc ctr line');
figure;
c = contour(phi);
clabel(c);
xlabel('x');
ylabel('y');