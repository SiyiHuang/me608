Cell = 50;
T1 = 100;
T2 = 200;

AW = ones(Cell);
AE = ones(Cell);
AS = ones(Cell);
AN = ones(Cell);
AW(:,1) = 0;
AE(:,Cell) = 0;
AS(1,:) = 0;
AN(Cell,:) = 0;

AP = AW + AE + AS + AN;
AP(:,1) = AP(:,1) + 2;
AP(:,Cell) = AP(:,Cell) + 2;
AP(1,:) = AP(1,:) + 2;
AP(Cell,:) = AP(Cell,:) + 2;

B = zeros(Cell);
B(1,:) = B(1,:) + 2 * T1;
B(:,1) = B(:,1) + 2 * T1;
B(:,Cell) = B(:,Cell) + 2 * T1;
B(Cell,:) = B(Cell,:) + 2 * T2;

T = T1 * ones(Cell);
T0 = T;
iter = 0;
T = lblTDMA(AP,AW,AE,AS,AN,T,B);
while max(max(abs(T-T0))) > 1e-5
    T0 = T;
    T = lblTDMA(AP,AW,AE,AS,AN,T,B);
    iter = iter + 1;
end
