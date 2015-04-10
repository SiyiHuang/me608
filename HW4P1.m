% This is a matlab script to solve Problem 1 in ME608 HW4. 

% Properties and Constants

rho = 1; % Density
p_1 = 28; % Pressure at cell centroid 1
p_3 = 0;  % Pressure at cell centroid 3
A = [3,1];  % Cross sectional area
u_1 = 0;  

% The momentum equation is discretized using UDS. 

% Initial guess 
m_dot = 5; % m_dot = rho * u * A;
u = m_dot / rho ./ A; % Guess of velicity
p_2 = 25; % I'm curious how Patankar came up with this? 
up = zeros(1,2);

% Start SIMPLE:
tol = 1; % Initialize tolerance. 
iter = 0;
while tol > 1e-2
%%#1 Guess of p. 
p_2s = p_2;
%%#2 Solve u_Ps in staggered grid by guessed pressure source. 

a_W = [0, rho*m_dot]; % UDS
a_E = [0, 0];
a_P = a_W + a_E + [rho*m_dot, 0];
b_u = [m_dot*u_1 + A(1)*(p_1 - p_2s), A(2)*(p_2s-p_3)];
us = TDMA(a_P, a_W, a_E, u, b_u);

%%#3 Find F* and form pressure correction sourse term. 
F_s = rho.*us.*A;
b = F_s(1) - F_s(2);
%%#4 Solve pressure correction term p'
aw = A(1)/u(1);
ae = A(2)/u(2);
p_2p = b / (aw+ae);
%%#5 Correct Pressure and velocity 
p_2 = p_2s + p_2p;
up = TDMA(a_P, a_W, a_E, up, A.*[-p_2p, p_2p]);
alpha = 0.24;
tol = max(abs(us+alpha*up-u));
u = us+alpha*up;
m_dot = u(1)*A(1);
%u
%pause
iter = iter+1;
end

u
p_2
iter
