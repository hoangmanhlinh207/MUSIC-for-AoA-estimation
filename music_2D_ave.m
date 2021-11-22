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
%% Actual DOA
phi=[20 30]/180*pi;    % Azimuth
theta=[10 40]/180*pi;  % Elevation
P=length(theta);               % 
%% Signal frequency
f=1/8*95e9;                  % frequency of incoming signals
snapshot=500;
fs=2e11;                     % sampling frequency
%% Array characteristic
M=4;    % 4x4 array
N=4;
c=3e8;
lambda=c/f;
deltad=lambda/2;
%% Signal parameters
S=200;  % Number of signal samples
snr=10; % SNR
QAM=4;
A=2;    % Ampliture
s1_mod=qammod(randi(QAM,snapshot,1)'-1,QAM);  % minus 1 because randi(M) generate integer from 1 to M, whereas qamod(x,M)
s1=s1_mod.*exp(1i*2*pi*f*[1:snapshot]/fs);
s2_mod=qammod(randi(QAM,snapshot,1)'-1,QAM);  % minus 1 because randi(M) generate integer from 1 to M, whereas qamod(x,M)
                                     % takes value of x from 0 to M)
s2=s2_mod.*exp(1i*2*pi*f*[1:snapshot]/fs);
s=[s1;s1];                           % Totally coherent sources
D_vec_save=[];
for t=1:P
ami=(exp(-1i*2*pi*deltad*sin(theta(t))*cos(phi(t))/lambda*[0:M-1]))';
vi=exp(-1i*2*pi*deltad*sin(theta(t))*sin(phi(t))/lambda*[0:N-1]);
D=ami*vi;                                                               % matrix in 2D
D_vec=D(:);                                                             % Change 2D into 1D
D_vec_save=[D_vec_save D_vec];                                          % Save the matrix
end
%% Received signal generator
x=D_vec_save*s;              % Received signal matrix without noise
x=x+awgn(x,snr);             % received signal matrix with noise
%% Perform MUSIC
R=x*x';
J=fliplr(eye(M*N));
R=0.5*(R+J*conj(R)*J);  % Forward/Backward (FB) averaging
[vec,~]=eig(R);
NN=vec(:,1:M*N-P);
%% Create the search angles
theta_search=(-90:0.5:90)/180*pi;  % Elevation search vector
phi_search=(-90:0.5:90)/180*pi;
tic
term1=-1i*2*pi*deltad/lambda*[0:M-1];    % save for less re-calculation
term2=-1i*2*pi*deltad/lambda*[0:N-1];    % save for less re-calculation
for kk=1:length(phi_search)
    for tt=1:length(theta_search)
        ami_temp=(exp(term1*sin(theta_search(tt))*cos(phi_search(kk))))';
        vi_temp=exp(term2*sin(theta_search(tt))*sin(phi_search(kk)));
        D_temp=ami_temp*vi_temp;
        D_vec_temp=D_temp(:);
        PP(tt)= D_vec_temp'*NN*NN'*D_vec_temp;
        Pmusic_im(tt)=abs(1/PP(tt));
    end
    Pmusic_im_all(:,kk)=Pmusic_im;
end
toc
%% Plot the result
Pmusic_im_all=10*log10(Pmusic_im_all/max(max(Pmusic_im_all)));
figure,
fts=25;
contour(theta_search*180/pi,phi_search*180/pi,Pmusic_im_all)
ylabel('Elevation angle \theta [degree]','Fontsize',fts)
xlabel('Azimuth angle \phi [degree]','Fontsize',fts)
title('With Forward/Backward (FB) averaging','fontsize',fts)
set(gca,'fontsize',fts);
grid on

%% Perform MUSIC
R=x*x';
[vec,~]=eig(R);
NN=vec(:,1:M*N-P);
%% Create the search angles
theta_search=(-90:0.5:90)/180*pi;  % Elevation search vector
phi_search=(-90:0.5:90)/180*pi;
tic
term1=-1i*2*pi*deltad/lambda*[0:M-1];    % save for less re-calculation
term2=-1i*2*pi*deltad/lambda*[0:N-1];    % save for less re-calculation
for kk=1:length(phi_search)
    for tt=1:length(theta_search)
        ami_temp=(exp(term1*sin(theta_search(tt))*cos(phi_search(kk))))';
        vi_temp=exp(term2*sin(theta_search(tt))*sin(phi_search(kk)));
        D_temp=ami_temp*vi_temp;
        D_vec_temp=D_temp(:);
        PP(tt)= D_vec_temp'*NN*NN'*D_vec_temp;
        Pmusic_im(tt)=abs(1/PP(tt));
    end
    Pmusic_im_all(:,kk)=Pmusic_im;
end
toc
%% Plot the result
Pmusic_im_all=10*log10(Pmusic_im_all/max(max(Pmusic_im_all)));
figure,
fts=25;
contour(theta_search*180/pi,phi_search*180/pi,Pmusic_im_all)
ylabel('Elevation angle \theta [degree]','Fontsize',fts)
xlabel('Azimuth angle \phi [degree]','Fontsize',fts)
title('Without Forward/Backward (FB) averaging','fontsize',fts)
set(gca,'fontsize',fts);
grid on

% End of code