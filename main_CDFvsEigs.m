% cap_bd_muser.m
% NU-ユーザでのBD法の計算
% TDMA, Upper boundと比較

clear;%close;

% 出力ファイル名
% testfile1 = 'Capacity2x8x4u_BD.csv';
% testfile2 = 'Capacity2x8x4u_BD_CDFSNR20dB.csv';
% testfile3 = 'Eig2x8x4u_BD_CDFSNR20dB.csv';

% パラメータ
% パラメータ条件 NT >= NR*NU
SN_tar  = 30;        % CDF表示のためのターゲットSNR [dB]
%SN_max = 40;         % 最大SNR[dB]
Ntri   = 1000;       % 伝搬チャネル行列の発生回数
NT     = 16;          % 送信素子数
NR     = 2;          % 受信素子数(=2に固定)
NU     = 8;          % ユーザ数
I      = eye(NT,NT); % NTxNTの単位行列

% 擬似雑音
a = NT/(10^(SN_tar/10));
%a2 = NR*NT/(10^(SN_tar/10));    % 擬似雑音 for BMSN2
        
% 所望のチャネル行列 for BMSN
for nuser = 1:NU
    T(:,:,nuser) = eye(NR,NR);
end

if SN_tar < 10
    target_snr=strcat('SNR= ',num2str(SN_tar,'%01d'),'dB');
else
    target_snr=strcat('SNR=',num2str(SN_tar,'%02d'),'dB');
end

% 伝搬チャネル行列
% i.i.d. Rayleigh
H = (randn(Ntri,NR*NU,NT)+1j*randn(Ntri,NR*NU,NT))/sqrt(2);

E_ZF = zeros(Ntri, NR, NU);             % ZF-CIの固有値
E_BD = zeros(Ntri, NR, NU);             % BDの固有値
E_MMSE = zeros(Ntri, NR, NU);           % MMSE-CIの固有値
E_BMSN_BF = zeros(Ntri, NR, NU);        % BMSN-BFの固有値
E_BMSN_GE = zeros(Ntri, NR, NU);        % BMSN-GEの固有値

for k = 1:Ntri              % 試行回数のループ

    H0 = squeeze(H(k,:,:)); % k番目の試行回数での伝搬チャネル行列
    
    % ZF-CI algorithm
    [W_ZF,U_ZF,S_ZF,~,~] = zf(NT,NR,NU,H0); % function zfd.m を使用
    
    % BD algorithm
    [W_BD,U_BD,S_BD,~,~] = bd(NT,NR,NU,H0); % function bd.m を使用
    
    % MMSE-CI algorithm
    [W_MMSE,U_MMSE,S_MMSE,RIPM,~] = mmse(NT,NR,NU,H0,a); % function mmse.m を使用
    
    % BMSN-BF algorithm
    [W_BMSN_BF,U_BMSN_BF,S_BMSN_BF,RIPBF,~] = bmsn_bf(NT,NR,NU,H0,a,T); % function bmsn_bf.m を使用
      
    % BMSN-GE algorithm
    [W_BMSN_GE,U_BMSN_GE,S_BMSN_GE,RIPGE,~] = bmsn_gev(NT,NR,NU,H0,a); % function bmsn_gev.m を使用
    
    
    % ユーザ毎の固有値分布
    snt = 1/(10^(SN_tar/10));
    for nuser=1:NU
        if NR==1
            E_ZF(k,:,nuser) = 10*log10(S_ZF(1,1,nuser).^2/(NT*snt));
            E_BD(k,:,nuser) = 10*log10(S_BD(1,1,nuser).^2/(NT*snt));
            E_MMSE(k,:,nuser) = 10*log10((S_MMSE(1,1,nuser).^2)./(RIPM(1,nuser)+NT*snt));
            E_BMSN_BF(k,:,nuser) = 10*log10((S_BMSN_BF(1,1,nuser).^2)./(RIPBF(1,nuser)+NT*snt));
            E_BMSN_GE(k,:,nuser) = 10*log10((S_BMSN_GE(1,1,nuser).^2)./(RIPGE(1,nuser)+NT*snt));
        else 
            E_ZF(k,:,nuser) = 10*log10(diag(S_ZF(:,:,nuser)).^2/(NT*snt));   
            E_BD(k,:,nuser) = 10*log10(diag(S_BD(:,:,nuser)).^2/(NT*snt));
            E_MMSE(k,:,nuser) = 10*log10((diag(S_MMSE(:,:,nuser)).^2)./(RIPM(:,nuser)+NT*snt));
            E_BMSN_BF(k,:,nuser) = 10*log10((diag(S_BMSN_BF(:,:,nuser)).^2)./(RIPBF(:,nuser)+NT*snt));
            E_BMSN_GE(k,:,nuser) = 10*log10((diag(S_BMSN_GE(:,:,nuser)).^2)./(RIPGE(:,nuser)+NT*snt));
        end
    end

    
end

% CDF of Eigenvalue at Target SNR
rr2(:,1) = (1/Ntri:1/Ntri:1).'*100; 
Y = rr2;
rr3 = zeros(Ntri,NR*5);
for nuser=1:NU
    E_ZF(:,:,nuser) = sort(E_ZF(:,:,nuser),1);
    E_BD(:,:,nuser) = sort(E_BD(:,:,nuser),1);
    E_MMSE(:,:,nuser) = sort(E_MMSE(:,:,nuser),1);
    E_BMSN_BF(:,:,nuser) = sort(E_BMSN_BF(:,:,nuser),1);
    E_BMSN_GE(:,:,nuser) = sort(E_BMSN_GE(:,:,nuser),1);
end
% ユーザ平均
E_ZFm = mean(E_ZF,3);
E_BDm = mean(E_BD,3);
E_MMSEm = mean(E_MMSE,3);
E_BMSN_BFm = mean(E_BMSN_BF,3);
E_BMSN_GEm = mean(E_BMSN_GE,3);

for nn=1:NR
    rr3(:,nn) = E_BMSN_BFm(:,nn);
    rr3(:,nn+NR) = E_BMSN_GEm(:,nn);
    rr3(:,nn+2*NR) = E_MMSEm(:,nn);
    rr3(:,nn+3*NR) = E_BDm(:,nn);
    rr3(:,nn+4*NR) = E_ZFm(:,nn);
end


%% CDF of EGV at Target SNR 
figure;
mycol = [1 0 1;1 0 1;0 0 0;0 0 0;
      1 0 0;1 0 0;0 0 1;0 0 1;
      1 0 0;1 0 0;
      0 0 0;0 0 0];
set(groot,'defaultAxesColorOrder',mycol)
plot(rr3(:,1),Y,'r-d','MarkerIndices',10:100:length(Y),'Linewidth',2);
axis([-20 30 0 100]);
grid on;
hold on;
set(gca,'XTick',-20:5:30,'Fontsize',14,'Fontname','Arial')
xlabel('SINR of eigenvalue [dB]','Fontsize',16,'Fontname','Arial');
ylabel('CDF [%]','Fontsize',16,'Fontname','Arial');
plot(rr3(:,2),Y,'r--d','MarkerIndices',50:100:length(Y),'Linewidth',2);
plot(rr3(:,3),Y,'b-s','MarkerIndices',10:100:length(Y),'Linewidth',2);
plot(rr3(:,4),Y,'b--s','MarkerIndices',50:100:length(Y),'Linewidth',2);
plot(rr3(:,5),Y,'g-o','MarkerIndices',10:100:length(Y),'Linewidth',2);
plot(rr3(:,6),Y,'g--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
plot(rr3(:,7),Y,'c-x','MarkerIndices',10:100:length(Y),'Linewidth',2);
plot(rr3(:,8),Y,'c--x','MarkerIndices',50:100:length(Y),'Linewidth',2);
plot(rr3(:,9),Y,'m-^','MarkerIndices',10:100:length(Y),'Linewidth',2);
plot(rr3(:,10),Y,'m--^','MarkerIndices',50:100:length(Y),'Linewidth',2);
%plot(rr3(:,6),rr3(:,4),'m');
legend('BMSN-BF \lambda_1','BMSN-BF \lambda_2','BMSN-GE \lambda_1','BMSN-GE \lambda_2','MMSE-CI \lambda_1','MMSE-CI \lambda_2',...
    'BD \lambda_1','BD \lambda_2','ZF-CI \lambda_1','ZF-CI \lambda_2','Location','southeast');
title(target_snr);



