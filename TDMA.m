function output = TDMA(AP, AW, AE,X, B)
%This function is a TDMA solver of tridiagnal matrix.
%The input are supposed to be column vectors.
N = length(B);

for i = 2:N
    r = -AW(i)/AP(i-1);
    AP(i) = AP(i) + r * AE(i-1);
    B(i) = B(i) - r * B(i-1);
end
X(N) = B(N) / AP(N);
for i = (N-1):(-1):1
    X(i) = (B(i) + AE(i) * X(i+1))/ AP(i);
end
output = X;
end