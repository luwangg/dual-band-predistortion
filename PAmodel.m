function [y]=PAmodel(x)



M=3;N = length(x);
y=zeros(1,N);v1=zeros(1,N);v2=zeros(1,N);v3=zeros(1,N);y=zeros(1,N);

n=1:N;%ÐÅºÅ³¤¶È
y(1:M-1)=x(1:M-1);
%hammerstein
for i=M:N
    %v1(i) = saleh(x(i));
    v1(i) = 0.2*( x(i) + (x(i))^2  + (x(i))^3 + (x(i))^4 + (x(i))^5);%+(x(i))^6+(x(i))^7) ;
    y(i) = v1(i) + 0.1*v1(i-1)+0.1*y(i-2);

end
%wiener
% for i=M:N
% 
%     v2(i) = x(i) + 0.1*x(i-2) + 0.1*v2(i-1);%IIR
%     %y(i) = saleh(v2(i));
%     y(i) = v2(i) + 0.1*v2(i)^2 + 0.1*v2(i)^3;
% end
%y = y./1.17;
%wiener-hammerstein
% for i=M:N
%     v3(i) = x(i) + 0.1*x(i-2) + 0.1*v3(i-1);
%     v4(i) = saleh(v3(i));
%     y(i) = v4(i) + 0.1*v4(i-2) + 0.1*y(i-1);
% end

% py=10*log10(abs(fft(y)));
% L=length(py);
% p1y=circshift(py',L/2)';
% figure,plot(n,p1y);

%figure,plot(10*log10(pwelch(hamm)));

end