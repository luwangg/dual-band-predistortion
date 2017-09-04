%ls algorithm 2D-DPD
close all;clear;clc;
pic = 1;
if pic==1
    load('input5-10.mat');
    fsig_o = fft(sig_o);
    xhigh_fft = [fsig_o(1:1229);zeros(6963,1)];%input16
    xlow_fft = [zeros(6963,1);fsig_o(6964:8192)];
end
if pic==2
    load('input10-15.mat');
    fsig_o = fft(sig_o);
    xhigh_fft=[zeros(409,1);fsig_o(410:1638);zeros(6554,1)];%input10-15
    xlow_fft=[zeros(6554,1);fsig_o(6555:7783);zeros(409,1)];
end
if pic==3
    load('input15-20.mat');
    fsig_o = fft(sig_o);
    xhigh_fft = [zeros(819,1);fsig_o(820:2048);zeros(6144,1)]; %input20
    xlow_fft = [zeros(6144,1);fsig_o(6145:7373);zeros(819,1)];
end
if pic==4
    load('input20-25.mat');
    fsig_o = fft(sig_o);
    xhigh_fft = [zeros(1228,1);fsig_o(1229:2458);zeros(5734,1)]; %input20
    xlow_fft = [zeros(5734,1);fsig_o(5735:6964);zeros(1228,1)];
end
if pic==5
    load('input25-30.mat');
    fsig_o = fft(sig_o);
    xhigh_fft = [zeros(1637,1);fsig_o(1638:2867);zeros(5325,1)]; %input20
    xlow_fft = [zeros(5325,1);fsig_o(5326:6555);zeros(1637,1)];
end
if pic==6
    load('input30-35.mat');
    fsig_o = fft(sig_o);
    xhigh_fft=[zeros(1638,1);fsig_o(1639:3686);zeros(4506,1)];
    xlow_fft=[zeros(4506,1);fsig_o(4507:6554);zeros(1638,1)];
end
if pic==7
    load('input35-40.mat');
    fsig_o = fft(sig_o);
    xhigh_fft=[zeros(2048,1);fsig_o(2049:4096);zeros(4096,1)];
    xlow_fft=[zeros(4096,1);fsig_o(4097:6144);zeros(2048,1)];
    nfft =length(sig_o);
    f = -50:100/nfft:49.99;
    plot(f,circshift(abs(fft(sig_o)),nfft/2),'k');
 %   axis([-30 30 0 120]);
end
if pic==8
    load('input30-40.mat');
    fsig_o = fft(sig_o);
    xhigh_fft=[zeros(1638,1);fsig_o(1639:4096);zeros(4096,1)];
    xlow_fft=[zeros(4096,1);fsig_o(4097:6554);zeros(1638,1)];
end
if pic==9
    load('input325-375.mat');
    fsig_o = fft(sig_o);
    xhigh_fft=[zeros(2252,1);fsig_o(2253:3482);zeros(4710,1)];
    xlow_fft=[zeros(4710,1);fsig_o(4711:5940);zeros(2252,1)];
end
if pic==10
    load('input14-16.mat');
    fsig_o = fft(sig_o);
    xhigh_fft=[zeros(983,1);fsig_o(984:1475);zeros(6717,1)];
    xlow_fft=[zeros(6717,1);fsig_o(6718:7209);zeros(983,1)];
end

x_low = ifft(xlow_fft);
x_high = ifft(xhigh_fft);
Q = 2;%memory depth
K = 5;%polynomial order
num = length(x_low) - Q;
%%
%无预失真输出
input = sig_o;%dual-band input
input = input./max(abs(input));
%output without DPD
output = PAmodel2(input);
out_fft=fft(output);
if pic==1
    outhigh_fft = [out_fft(1:1229),zeros(1,6963)];
    outlow_fft = [zeros(1,6963),out_fft(6964:8192)];
end
if pic==7
%     outhigh_fft=[zeros(1,2866),out_fft(2867:3277),zeros(1,4915)];
%     outlow_fft=[zeros(1,4915),out_fft(4916:5326),zeros(1,2866)];
    
    outhigh_fft=[zeros(1,2048),out_fft(2049:4096),zeros(1,4096)];
    outlow_fft=[zeros(1,4096),out_fft(4097:6144),zeros(1,2048)];
end
if pic==10
    outhigh_fft=[zeros(1,819),out_fft(820:1638),zeros(1,6554)];
    outlow_fft=[zeros(1,6554),out_fft(6555:7373),zeros(1,819)];
end
out_low = ifft(outlow_fft);out_low = out_low./max(abs(out_low));
out_high = ifft(outhigh_fft);out_high = out_high./max(abs(out_high));
y =output;
%%
%参数设置
x1 = max(abs(x_low));
x2 = max(abs(x_high));
temp = 0.5*(x1+x2);
x_low = x_low./x1;%normalization
x_high = x_high./x2;

% f_sigo = circshift(fsig_o,length(sig_o)/2);
% figure,plot(fo(1638:6554),abs(f_sigo(1638:6554)),'k');
% axis([-60 60 0 120]);

%初始预失真器参数为c1=c2=[1,0,0,...]
z1 = x_low(1+Q:num+Q);
z2 = x_high(1+Q:num+Q);
z = input;
rate = 1;choose=1;

ya = output;ya = ya./max(abs(ya));
za = input(1+Q:num+Q);

for q=0 : Q
    for k=1 : K
        Ua(:,(Q-q)*K+k)=ya(q+1:num+q).*abs(ya(q+1:num+q)).^(k-1);

    end
end
for q=0 : Q
    for k=1 : K
        Uxa(:,(Q-q)*K+k)=input(q+1:num+q).*abs(input(q+1:num+q)).^(k-1);

    end
end
    ca=inv(Ua'*Ua)*(Ua')*za;
    za = Uxa*ca;
    za=za./max(abs(za));
    
    ya(1+Q:num+Q)=PAmodel2(za);

%    y=ya; z(1+Q:end)=za;
for iter = 1:5
%%
%数据预处理
%     dd=1;
%     y_00 = y(dd:length(input)+dd-1);%从y中取与原始信号等长的一段数据
%     y_00=y_00./max(abs(y_00));%归一化
%     y_00=y_00.';
%     
%     ck=xcorr(input,y_00);%用相关操作进行预对齐
%     [cm,C_indx]=max(ck);
%     %C_indx = C_indx mod 4000;
%     y_00=circshift(y_00,[C_indx,0]);
% 
% %     y_0 = y_00;
%     y_0 = align_signals(y_00.',input.');
%          align_signals(y_00.',sig_o.');
    
    y_fft = fft(y);
    if pic==1
        yhigh_fft = [y_fft(1:1229),zeros(1,6963)];
        ylow_fft = [zeros(1,6963),y_fft(6964:8192)];
    end
    if pic==2
        yhigh_fft=[y_fft(1:2048),zeros(1,6144)];%input10-15
        ylow_fft=[zeros(1,6144),y_fft(6145:8192)];
    end
    if pic==3
        yhigh_fft = [zeros(1,819),y_fft(820:2048),zeros(1,6144)]; %input20
        ylow_fft = [zeros(1,6144),y_fft(6145:7373),zeros(1,819)];
    end
    if pic==4
        yhigh_fft = [zeros(1,1228),y_fft(1229:2458),zeros(1,5734)]; %input20
        ylow_fft = [zeros(1,5734),y_fft(5735:6964),zeros(1,1228)];
    end
    if pic==5
%         yhigh_fft = [zeros(1,1637),y_fft(1638:2867),zeros(1,5325)]; 
%         ylow_fft = [zeros(1,5325),y_fft(5326:6555),zeros(1,1637)];
        yhigh_fft = [y_fft(1:4096),zeros(1,4096)]; 
        ylow_fft = [zeros(1,4096),y_fft(4097:8192)];
    end
    if pic==6
        yhigh_fft=[zeros(1,1638),y_fft(1639:3686),zeros(1,4506)];
        ylow_fft=[zeros(1,4506),y_fft(4507:6554),zeros(1,1638)];
    end
    if pic==7
    yhigh_fft=[zeros(1,2048),y_fft(2049:4096),zeros(1,4096)];
    ylow_fft=[zeros(1,4096),y_fft(4097:6144),zeros(1,2048)];
    end
    if pic==8
    yhigh_fft=[zeros(1,1638),y_fft(1639:4096),zeros(1,4096)];
    ylow_fft=[zeros(1,4096),y_fft(4097:6554),zeros(1,1638)];
    end
    if pic==9
    yhigh_fft=[zeros(1,2252),y_fft(2253:3482),zeros(1,4710)];
    ylow_fft=[zeros(1,4710),y_fft(4711:5940),zeros(1,2252)];
    end
    if pic==10
%     yhigh_fft=[zeros(1,983),y_fft(984:1475),zeros(1,6717)];
%     ylow_fft=[zeros(1,6717),y_fft(6718:7209),zeros(1,983)];
    yhigh_fft=[zeros(1,819),y_fft(820:1638),zeros(1,6554)];
    ylow_fft=[zeros(1,6554),y_fft(6555:7373),zeros(1,819)];
    end

    y_low = ifft(ylow_fft);  y_low = y_low./max(abs(y_low)); 
    y_high = ifft(yhigh_fft);  y_high = y_high./max(abs(y_high)); 
% if choose == 1
%     y_high = y_high.*sqrt(rate);
% end
% if choose == 0
%     y_low = y_low.*sqrt(1/rate);
% end
%%
%预失真器建模
tic
    [Ux1,Ux2]=PDmodel2(x_low,x_high,K,Q);
    [Uy1,Uy2]=PDmodel2(y_low,y_high,K,Q);

%%
%LS算法求解预失真器参数 
    c1=inv(Uy1'*Uy1)*(Uy1')*z1;
    c2=inv(Uy2'*Uy2)*(Uy2')*z2;
toc    
%预失真器输出相加进入PA
    z1 = Ux1*c1;%z1 = z1/rate;
    z2 = Ux2*c2;%z2 = z2/1.01;
    z(1+Q:end) = z1+z2;
    z=z./max(abs(z));
    
    y=PAmodel2(z);
    if pic ==1
        [power_hf_inband,power_lf_inband,rate,choose] =powertest1(y);
    end
    if pic==2
        [power_hf_inband,power_lf_inband,rate,choose] =powertest2(y);
    end
    if pic==7
        [power_hf_inband,power_lf_inband,rate,choose] =powertest7(y);
    end
    if pic==10
        [power_hf_inband,power_lf_inband,rate,choose] =powertest10(y);
    end

end

%     yhigh_fft=[zeros(1,2866),y_fft(2867:3277),zeros(1,4915)];
%     ylow_fft=[zeros(1,4915),y_fft(4916:5326),zeros(1,2866)];
 y_fft=fft(ya);
        yhigh_fft = [y_fft(1:1229),zeros(1,6963)];%filter_input14/16_正好3阶带宽
        ylow_fft = [zeros(1,6963),y_fft(6964:8192)];
    y_low = ifft(ylow_fft);  y_low = y_low./max(abs(y_low)); 
    y_high = ifft(yhigh_fft);  y_high = y_high./max(abs(y_high));
%  y_fft = fft(ya);
%  yhigh_fft = [y_fft(1:1229),zeros(1,6963)];%filter_input14/16_正好3阶带宽
%  ylow_fft = [zeros(1,6963),y_fft(6964:8192)];
%     y_low = ifft(ylow_fft);  y_low = y_low./max(abs(y_low)); 
%     y_high = ifft(yhigh_fft);  y_high = y_high./max(abs(y_high));

%NMSE result
% e_low = x_low - y_low.';  e_high = x_high - y_high.';
% NMSE_low=db(sum(e_low.^2)/sum(abs(x_low).^2));disp('low band NMSE=');disp(NMSE_low);
% NMSE_high=db(sum(e_high.^2)/sum(abs(x_high).^2));disp('high band NMSE=');disp(NMSE_high);

%PSD curve
%load('pxx_a_16.mat');

nfft = 1024;
pxx_out = pwelch(output,1024); p_out = 10*log10(circshift(pxx_out,nfft/2));
pxx_y = pwelch(y,1024); p_y = 10*log10(circshift(pxx_y,nfft/2));
pxx_in = pwelch(input,1024); p_in = 10*log10(circshift(pxx_in,nfft/2));
pxx_a = pwelch(ya,1024); p_a = 10*log10(circshift(pxx_a,nfft/2));
f = -50:100/nfft:49.99;%input12
f=f+1950;
figure, plot(f,p_out,'k');
hold on; plot(f,p_y,'ro-');
hold on; plot(f,p_in,'--g');
%hold on;plot(f,10*log10(circshift(pxx_a,nfft/2)),'y*-');
xlabel('频率(MHz)');ylabel('功率密度(dB/Hz)');
%axis([-15 15 -60 0]);
legend('无预失真','2D-DPD模型','输入信号');
grid on;
axis([1950-15 1950+15 -60 0])
pxx_out = pwelch(output);
pxx_y = pwelch(y);
pxx_in = pwelch(input);
nfft = length(pxx_in);
f = -50:100/nfft:49.99;%input12
%acpr
% if pic==1
% [acpr_mpdpd_lf_left,acpr_mpdpd_lf_right,acpr_mpdpd_hf_left,acpr_mpdpd_hf_right] = acpr14(f,pxx_a);
% [acpr_2ddpd_lf_left,acpr_2ddpd_lf_right,acpr_2ddpd_hf_left,acpr_2ddpd_hf_right] = acpr14(f,pxx_y);
% end
if pic==7
nfft=length(pxx_y);
pxx_y = 10*log10(pxx_y);
power_hf_inband = sum(pxx_y(758:778))/20;
power_lf_inband = sum(pxx_y(1270:1290))/20;
power_hf_outbandLeft = sum(pxx_y(655:675))/20;
power_hf_outbandRight= sum(pxx_y(860:880))/20;
acpr_hf_left = power_hf_outbandLeft - power_hf_inband;
acpr_hf_right = power_hf_outbandRight - power_hf_inband;

power_lf_outbandLeft = sum(pxx_y(1167:1187))/20;
power_lf_outbandRight=sum(pxx_y(1372:1392))/20;
acpr_lf_left = power_lf_outbandLeft - power_lf_inband;
acpr_lf_right = power_lf_outbandRight - power_lf_inband;
%[acpr_withoutdpd_lf_left,acpr_withoutdpd_lf_right,acpr_withoutdpd_hf_left,acpr_withoutdpd_hf_right] = acpr7(f,pxx_out);
%  disp(sprintf('\r\n2D-DPD Low Band ACPR_left %2.2f',acpr_2ddpd_lf_left));
%  disp(sprintf('\r\n2D-DPD Low Band ACPR_right %2.2f',acpr_2ddpd_lf_right));
%  disp(sprintf('\r\n2D-DPD High Band ACPR_left %2.2f',acpr_2ddpd_hf_left));
%  disp(sprintf('\r\n2D-DPD High Band ACPR_right %2.2f',acpr_2ddpd_hf_right));
end
if pic==10
nfft=length(pxx_y);
pxx_y = 10*log10(pxx_y);
power_hf_inband = sum(pxx_y(297:317))/20;
power_lf_inband = sum(pxx_y(1731:1751))/20;
power_hf_outbandLeft = sum(pxx_y(256:276))/20;
power_hf_outbandRight= sum(pxx_y(338:358))/20;
acpr_hf_left = power_hf_outbandLeft - power_hf_inband;
acpr_hf_right = power_hf_outbandRight - power_hf_inband;

power_lf_outbandLeft = sum(pxx_y(1690:1710))/20;
power_lf_outbandRight=sum(pxx_y(1772:1792))/20;
acpr_lf_left = power_lf_outbandLeft - power_lf_inband;
acpr_lf_right = power_lf_outbandRight - power_lf_inband;
end
%AMAM & AMPM
figure,
subplot(121);plot(abs(x_low(15:length(x_high)-15)),abs(out_low(15:length(x_high)-15)),'b.');hold on;plot(abs(x_low(15:length(x_high)-15)),abs(y_low(15:length(x_high)-15)),'r.');xlabel('低频段归一化输入信号幅度');ylabel('低频段信号输出信号幅度');legend('无预失真输出','预失真后输出');axis([0 1 0 1]);grid on;
subplot(122);plot(abs(x_low(15:length(x_high)-15)),angle(out_low(15:length(x_high)-15).'./x_low(15:length(x_high)-15)),'b.');hold on;plot(abs(x_low(15:length(x_high)-15)),angle(y_low(15:length(x_high)-15).'./x_low(15:length(x_high)-15)),'r.');xlabel('低频段归一化输入信号幅度');ylabel('低频段信号输入输出相位差');legend('无预失真输出','预失真后输出');axis([0 1 -1 1]);grid on;
% subplot(223);plot(abs(x_high(15:length(x_high)-10)),abs(y_high(15:length(x_high)-10)),'.');hold on;plot(abs(x_high(15:length(x_high)-10)),abs(out_high(15:length(x_high)-10)),'b.');xlabel('高频段归一化输入信号幅度');ylabel('高频段信号输出信号幅度');legend('双频段预失真后输出','无预失真输出');axis([0 1 -1 1]);
% subplot(224);plot(abs(x_high(11:length(x_high)-10)),angle(out_high(11:length(x_high)-10).'./x_high(11:length(x_high)-10)),'y.');hold on;plot(abs(x_high(15:length(x_high)-10)),angle(out_high(15:length(x_high)-10).'./x_low(15:length(x_high)-10)),'b.');xlabel('高频段归一化输入信号幅度');ylabel('高频段信号输入输出相位差');legend('双频段预失真后输出','无预失真输出');axis([0 1 -1 1]);
box on;
figure,
subplot(121);plot(abs(x_high(15:length(x_high)-10)),abs(out_high(15:length(x_high)-10)),'b.');hold on;plot(abs(x_high(15:length(x_high)-10)),abs(y_high(15:length(x_high)-10)),'r.');xlabel('高频段归一化输入信号幅度');ylabel('高频段信号输出信号幅度');legend('无预失真输出','预失真后输出');axis([0 1 0 1]);grid on;
subplot(122);plot(abs(x_high(15:length(x_high)-10)),angle(out_high(15:length(x_high)-10).'./x_high(15:length(x_high)-10)),'b.');hold on;plot(abs(x_high(15:length(x_high)-10)),angle(y_high(15:length(x_high)-10).'./x_high(15:length(x_high)-10)),'r.');xlabel('高频段归一化输入信号幅度');ylabel('高频段信号输入输出相位差');legend('无预失真输出','预失真后输出');axis([0 1 -1 1]);grid on;
box on;
% figure,
% subplot(121);plot(abs(x_low(15:length(x_low)-15)),abs(y_low(15:length(x_low)-15)),'.');hold on;plot(abs(x_low(15:length(x_low)-15)),angle(y_low(15:length(x_low)-15).'./x_low(15:length(x_low)-15)),'y.');xlabel('低频段归一化输入信号');ylabel('低频段信号输出输入相位差');axis([0 1 -1 1]);
% subplot(122);hold on;plot(abs(x_high(15:length(x_low)-15)),abs(y_high(15:length(x_low)-15)),'.');hold on;plot(abs(x_high(15:length(x_high)-15)),angle(y_high(15:length(x_high)-15).'./x_high(15:length(x_high)-15)),'y.');xlabel('高频段归一化输入信号');ylabel('高频段信号输出输入相位差');axis([0 1 -1 1]);
