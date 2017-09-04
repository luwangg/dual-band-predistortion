%2D-DPD PDmodel
function [U1,U2] = PDmodel2(x_low,x_high,K,Q)

% Q = 2;%记忆深度
% K = 5;%多项式阶数
N = length(x_high)-Q;
U1=zeros(N,(Q+1)*sum(1:K)); U2=zeros(N,(Q+1)*sum(1:K));

for q = 0:Q
    for k = 0:K-1
        for j = 0:k
            U1(:,q*sum(1:K)+sum(0:k)+j+1) = x_low(q+1:N+q).*(abs(x_low(q+1:N+q)).^(k-j)).*(abs(x_high(q+1:N+q)).^j);
        end
    end
end
for q = 0:Q
    for k = 0:K-1
        for j = 0:k
            U2(:,q*sum(1:K)+sum(0:k)+j+1) = x_high(q+1:N+q).*(abs(x_low(q+1:N+q)).^(k-j)).*(abs(x_high(q+1:N+q)).^j);
        end
    end
end
end