% This script is to solve the rhomboidal duct velocity profile. 

L = 5e-2; % The length of the equilateral triangle. 
D_H = sqrt(3) / 2 * L; %Hydrolic Diameter.
mu = 1e-3; % Dynamic Viscosity.
dPdz = 10; % Pressure drop

% For 8 cell schematics, size:
A_f = L/2; 
D_xi = sqrt(3)*L/6;
D_xi_b = sqrt(3)*L/12;
V_0 = sqrt(3) / 16 * L*2;

%Symmetric
b = ones(4,1) * dPdz * V_0;
b = b + [0;0;0;0]; % Boundary condition. Though 0 in this case. 
A = [1 -1 0 0; -1 3 -1 -1; 0 -1 1 0; 0 -1 0 1] * mu * A_f / D_xi;
A = A + diag([2 0 1 1]) * mu * A_f / D_xi_b;

w = zeros(4,1);
w0 = ones(4,1);
while max(max(abs(w-w0))) > 1e-5
    for i = 1:4
        w0 = w;
        w(i) = (1/A(i,i))*(b(i)-A(i,:)*w + A(i,i)*w(i));
    end
end

w1 = pinv(A)*b; % for small cell number, use pinv function to check G-S correctness. 

% Symmetric;

w = [w; zeros(4,1)];
w(8:-1:5) = w(1:4);

w_m = sum(w) / 8;
fRe = 2 * dPdz * D_H^2 / w_m / mu;

w_dl = w / w_m;

x_diag1 = [-3:-1,1:3]*D_xi;
w_diag1 = [0,w_dl(8),w_dl(7),w_dl(2),w_dl(1),0];
x_diag2 = [-2,-1,1,2]*L/4;
w_diag2 = [0, (w_dl(3)+w_dl(5))/2, (w_dl(6)+w_dl(4))/2, 0];
figure(1);
plot(x_diag1,w_diag1,x_diag2,w_diag2);
legend('Long Diag', 'Short Diag');
xlabel('Dist from CTR');
ylabel('T*');


% Fitting
Loc_x = [0,1.25,1,0.75,1,0.5,0.75,0.5,0.25,1.5,0.5,1];
Loc_y = [0 5 4 5 2 4 1 2 1 6,6,0]*sqrt(3)/12;
w_dl = [0;w_dl;0;0;0];
f = fit([Loc_x',Loc_y'],w_dl,'poly23');
figure(2);
plot(f,[Loc_x',Loc_y'],w_dl);
zlabel('T*');


x = (0:24)/16;
y1 = sqrt(3)/3*x;
x_ = 0.5+(0:20)/40;
y2 = sqrt(3)*(1-x_);
p00 = -0.05882;       %The following coefs are from printing f=fit(). 
p10 =       2.965;  
p01 =       1.712; 
p20 =      -2.918;  
p11 =        3.26;  
p02 =        -4.8;  
p21 =   9.375e-15;  
p12 =  -2.172e-15;  
p03 =   -1.54e-14;
w_diag1 = p00 + p10*x + p01*y1 + p20*x.^2 + p11*x.*y1 + p02*y1.^2 + p21*x.^2.*y1 + p12*x.*y1.^2 + p03*y1.^3;
w_diag2 = p00 + p10*x_ + p01*y2 + p20*x_.^2 + p11*x_.*y2 + p02*y2.^2 + p21*x_.^2.*y2 + p12*x_.*y2.^2 + p03*y2.^3;
figure(3);
plot((-12:12)*sqrt(3)*L/24,w_diag1,(-10:10)/20*L,w_diag2);
legend('Fit long diag','Fit short diag');
xlabel('Dist from CTR');
ylabel('T*');
