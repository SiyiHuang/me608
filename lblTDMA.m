function output = lblTDMA(AP, AW, AE, AS, AN, X, B)
%This is the line by line TDMA solver for 2D problem.
%To achieve bias balance, A same process is conducted on
%two ways.
%This iterative method requires a input of initial guess of X
[NJ,NI] = size(X);
for j = 1:NJ
    AP1D = AP(j,:);
    AE1D = AE(j,:);
    AW1D = AW(j,:);
    B1D = B(j,:);
    X1D = X(j,:);
    for i = 1:NI
        if j > 1
            B1D(i) = B1D(i) + AS(j,i) * X(j-1,i);
        end
        if j < NJ
            B1D(i) = B1D(i) + AN(j,i) * X(j+1,i);
        end
    end
    X1D = TDMA(AP1D,AW1D,AE1D,X1D,B1D);
    X(j,:) = X1D;    
end

% for j = NJ:(-1):1
%     AP1D = AP(:,j);
%     AE1D = AE(:,j);
%     AW1D = AW(:,j);
%     B1D = B(:,j);
%     X1D = X(:,j);
%     for i = 1:NI
%         if j > 1
%             B1D(i) = B1D(i) - AS(i,j) * X(i,j-1);
%         end
%         if j < NJ
%             B1D(i) = B1D(i) - AN(i,j) * X(i,j+1);
%         end
%     end
%     X1D = TDMA(AP1D,AW1D,AE1D,X1D,B1D);
%     X(:,j) = X1D;    
% end
output = X;
end