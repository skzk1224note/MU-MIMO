%% Riceで固有モード伝送とBD法の伝送レートを計算
% Programed by K. Nishimori 2016
clear; close all;
clc;

%% SNRと伝送レートの関係 (IEEE802.11ac, 20MHz)

TRT(1) = 7.2;
TRT(2) = 14.4;
TRT(3) = 21.7;
TRT(4) = 28.9;
TRT(5) = 43.3;
TRT(6) = 57.8;
TRT(7) = 65.0;
TRT(8) = 72.2;
TRT(9) = 86.7;

SNRT(1) = 6.637744;
SNRT(2) = 10.24946;
SNRT(3) = 12.07158;
SNRT(4) = 15.32538;
SNRT(5) = 18.35141;
SNRT(6) = 22.51627;
SNRT(7) = 23.88286;
SNRT(8) = 25.31453;
SNRT(9) = 29.25163;

% TRT(1) = 15; %IEEE802.11ac 40MHz伝送レート
% TRT(2) = 30;
% TRT(3) = 45;
% TRT(4) = 60;
% TRT(5) = 90;
% TRT(6) = 120;
% TRT(7) = 135;
% TRT(8) = 150;
% TRT(9) = 180;clear all;close all;
%%

% testfile1 = 'TR2x32x4u_BD_CDFSNR20dB.csv';
testfile2 = 'EIG_1x32x16u_BD_CDFSNR20dB_K10.csv';

% パラメータ
% パラメータ条件 NT >= NR*NU
SN_tar  = 20;        % CDF表示のためのターゲットSNR [dB]
SN_max = 40;         % 最大SNR[dB]
Ntri   = 1000;      % 伝搬チャネル行列の発生回数
NT     = 16;          % 送信素子数
NR     = 2;          % 受信素子数
NU     = 8;          % ユーザ数
I      = eye(NT,NT); % NTxNTの単位行列
K_dB   = 10;    % RicianのKファクタ
K      = 10^(K_dB/10);
d_t      = 0.5; % 送信アンテナ間隔（in wavelength)
d_r      = 0.5; % 受信アンテナ間隔（in wavelength)

derad = pi/180;

%%
% 伝搬チャネル行列
% Rayleigh
%H = (randn(Ntri,NR*NU,NT)+1i*randn(Ntri,NR*NU,NT))/sqrt(2);

E_BD   = zeros(Ntri, NR);              % BDの固有値
E_TDMA = zeros(Ntri, NR);              % TDMAの固有値
E_UB   = zeros(Ntri, NR*NU);           % Upper boundのチャネル容量

%%
%for SN_tar = 10:5:30

%for SN_tar = 0:30

H_los = zeros(NU*NR,NT);

for itr = 1:Ntri              % 試行回数のループ
    
    % iid Rayleigh チャネル行列
    H_iid = (randn(NR*NU,NT)+1j*randn(NR*NU,NT))/sqrt(2);
    
    
    H_los = zeros(NU*NR,NT);
   % LOS チャネル
    Theta_t = (rand-0.5)*360;   % ユーザ毎の送信角 (-180deg - 180deg)
    Theta_r = (rand(1,NU)-0.5)*360; % ユーザ毎の受信角 (-180deg - 180deg)
    a_t = exp(-1j*2*pi*d_t*(0:NT-1).'*sin(Theta_t*derad));  % 送信モードベクトル
    for n = 1 : NU
        a_r = exp(-1j*2*pi*d_r*(0:NR-1).'*sin(Theta_r(1,n)*derad)); % ユーザ毎の受信モードベクトル
        H_los((n-1)*NR+1:(n-1)*NR+NR,:) = a_r*a_t.';    % ユーザ毎のLOSチャネル行列
    end
 
    
    H0 = sqrt(K/(K+1))*H_los + sqrt(1/(K+1))*H_iid;
    
    % BD algorithm
    [W_BD,U_BD,S_BD] = bd(NT,NR,NU,H0); % function, bd.m を使用
    
    % SVD
    [U,S,V] = svd(H0);
    
    % 固有値分布
    Pn = 1/(10^(SN_tar/10));
    for nuser = 1 : NU
        if NR==1
            E_BD(itr,:,nuser) = 10*log10(S_BD(1,1,nuser).^2/(NT*Pn));
        else
            E_BD(itr,:,nuser) = 10*log10(diag(S_BD(:,:,nuser)).^2/(NT*Pn));
        end
    end
    E_BD_max(itr) = max(max(E_BD(itr,:,:)));
    E_BD_min(itr) = min(min(E_BD(itr,:,:)));
    
    
    %% Transmission rate
    TR_MU = zeros(NR,NU);
    for n = 1: length(TRT)
        TR_MU(squeeze(E_BD(itr,:,:)).' >= SNRT(n)) = TRT(n);
    end
    sumTR_MU(itr,1)=sum(sum(TR_MU,1));
    TR_MU_1st(itr,1) = TR_MU(1,1);
end

P = zeros(Ntri,3);
P(:,1) = [1/Ntri:1/Ntri:1].'*100;
P(:,2) = sort(E_BD_max);
P(:,3) = sort(E_BD_min);
csvwrite(testfile2,P);

fprintf('%d %f %f\n', NU, median(E_BD_max),median(E_BD_min));


sumTR_MUs = sort(sumTR_MU(:,1));



