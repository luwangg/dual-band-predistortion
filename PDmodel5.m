%2D-DPDS2 PDmodel
function [U1,U2] = PDmodel5(x_low,x_high,K,Q)

% Q = 2;%记忆深度
% K = 5;%多项式阶数
N = length(x_high)-Q;
%U1=zeros(N,(Q+1)*9); U2=zeros(N,(Q+1)*9);
U1=zeros(N,(Q+1)*6); U2=zeros(N,(Q+1)*6);

for q = 0:Q
    U1(:,q+1) = x_low(q+1:N+q);
    U1(:,(Q+1)+q+1) = x_low(q+1:N+q).*abs(x_low(q+1:N+q).^2);
    U1(:,(Q+1)*2+q+1) = x_low(q+1:N+q).*abs(x_high(q+1:N+q).^2);
    U1(:,(Q+1)*3+q+1) = x_low(q+1:N+q).*abs(x_low(q+1:N+q).^4);
    U1(:,(Q+1)*4+q+1) = x_low(q+1:N+q).*abs(x_low(q+1:N+q).^2).*abs(x_high(q+1:N+q).^2);
    U1(:,(Q+1)*5+q+1) = x_low(q+1:N+q).*abs(x_high(q+1:N+q).^4);
%     U1(:,(Q+1)*6+q+1) = (x_low(q+1:N+q).^2).*conj(x_high(q+1:N+q));
%     U1(:,(Q+1)*7+q+1) = (x_low(q+1:N+q).^2).*abs(x_low(q+1:N+q).^2).*conj(x_high(q+1:N+q));
%     U1(:,(Q+1)*8+q+1) = (x_low(q+1:N+q).^2).*abs(x_high(q+1:N+q).^2).*conj(x_high(q+1:N+q));
end

for q = 0:Q
    U2(:,q+1) = x_high(q+1:N+q);
    U2(:,(Q+1)+q+1) = x_high(q+1:N+q).*abs(x_low(q+1:N+q).^2);
    U2(:,(Q+1)*2+q+1) = x_high(q+1:N+q).*abs(x_high(q+1:N+q).^2);
    U2(:,(Q+1)*3+q+1) = x_high(q+1:N+q).*abs(x_low(q+1:N+q).^4);
    U2(:,(Q+1)*4+q+1) = x_high(q+1:N+q).*abs(x_low(q+1:N+q).^2).*abs(x_high(q+1:N+q).^2);
    U2(:,(Q+1)*5+q+1) = x_high(q+1:N+q).*abs(x_high(q+1:N+q).^4);
%     U2(:,(Q+1)*6+q+1) = (x_high(q+1:N+q).^2).*conj(x_low(q+1:N+q));
%     U2(:,(Q+1)*7+q+1) = (x_high(q+1:N+q).^2).*abs(x_high(q+1:N+q).^2).*conj(x_low(q+1:N+q));
%     U2(:,(Q+1)*8+q+1) = (x_high(q+1:N+q).^2).*abs(x_low(q+1:N+q).^2).*conj(x_low(q+1:N+q));
end
end