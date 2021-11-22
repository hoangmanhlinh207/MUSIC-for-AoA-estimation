%%
% Code name: 2D MUSIC algorithm 
% Author: Hoang Manh Linh
% hoangmanhlinhbachkhoa@gmail.com
% Reference: 
% Introduction to Direction-of-Arrival Estimation-Zhizhang Chen, Gopal Gokeda, Yiqiang Yu
% DOA estimation based on MUSIC algorithm-Honghao Tang
clc
clear 
close all
format long
N=200;fs=2e11;
doa=[40 60]/180*pi;
w=[pi/4 pi/4]'*95e9;
M=10;
Msub=3;
P=length(w);
c=3e8;
lambda=c*2*pi/w(1);
deltad=lambda/2;
% deltad=lambda/1.5;
snr=10;
D=zeros(P,M);
for k=1:P
D(k,:)=exp(-1i*2*[0:M-1]*pi*deltad*sin(doa(k))/lambda);
end
s=2*exp(1i*(w*[1:N]));
x=D'*s;
x=x+awgn(x,snr);
figure,
%% Without spatial smoothing
R_old=x*x';
J=fliplr(eye(M));
R_old=R_old+J*conj(R_old)*J;
[N,~]=eig(R_old);
NN=N(:,1:M-P);
theta=-90:0.5:90;
for ii=1:length(theta)
SS=zeros(1,length(M));
for jj=0:M-1
SS(1+jj)=exp(-1i*2*jj*pi*deltad*sin(theta(ii)/180*pi)/lambda);
end
PP= SS*NN*NN'*SS';
Pmusic_im(ii)=abs(1/PP);
end
Pmusic_im=10*log10(Pmusic_im/max(Pmusic_im));
plot(theta,Pmusic_im,'c');
hold on
%% When we use only one sub sub-array
R_sub_ma=[]; 
for t=1:M-Msub+1
x_sub=x(t:t+Msub-1,:);
R_sub=x_sub*x_sub';
R_sub_ma(t,:,:)=R_sub;
end
%% When we use spatial smoothing
R_fin=reshape(mean(R_sub_ma,1),[Msub,Msub]);
%J=fliplr(eye(Msub));
%R=R+J*conj(R)*J;       %Forward backward averaging, apply if interested
[N,V]=eig(R_fin);
NN=N(:,1:Msub-P);
theta=-90:0.5:90;
for ii=1:length(theta)
    SS=zeros(1,length(M));
    for jj=0:Msub-1
        SS(1+jj)=exp(-1i*2*jj*pi*deltad*sin(theta(ii)/180*pi)/lambda);
    end
    PP= SS*NN*NN'*SS';
    Pmusic_im(ii)=abs(1/PP);
end
Pmusic_im=10*log10(Pmusic_im/max(Pmusic_im));
plot(theta,Pmusic_im,'-k')
xlabel('angle \theta/degree')
ylabel('spectrum function P(\theta) /dB')
title('DOA estimation based on improved MUSIC algorithm')
grid on
% axis([-91 91 min(Pmusic_im)-2 0]);
legend('no spatial smoothing','with spatial smoothing');
legend boxoff
% % End of code