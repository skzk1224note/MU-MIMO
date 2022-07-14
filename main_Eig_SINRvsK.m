% K_lambda1,lambda2

clear;

% K[dB]の範囲設定
K_min = -20;        % 最小K [dB]
K_max = 20;         % 最大K [dB]

CDF = 10;   % 着目CDF

Nt = 16;    % 送信素子数
Nr = 2;     % 受信素子数 (アンテナ選択の場合1, そうでない場合は2)
Nu = 8;     % ユーザ数

SIMU = 1000; %試行回数（通常 1000）

I = eye(Nt,Nt); % Nt*Ntの単位行列
Nru = Nr*Nu;

SNR_tar = 10;   % ターゲットSNR[dB]

% SNR_max = 30; %最大SNR[dB]
d_t      = 0.5;      % 送信アンテナ間隔（in wavelength)
d_r      = 0.5;      % 受信アンテナ間隔（in wavelength)
derad = pi/180;      % degree -> rad

%T 所望のチャネル行列
T = zeros(Nr,Nr);
for inu = 1:Nu
    T(:,:,inu) = eye(Nr,Nr);
end

% 雑音と擬似雑音
snt = 1/(10^(SNR_tar/10));
a = Nt*snt;

if CDF < 10
    target_CDF=strcat('CDF= ',num2str(CDF,'%01d'),'%');
else
    target_CDF=strcat('CDF= ',num2str(CDF,'%02d'),'%');
end

if SNR_tar < 10
    target_SNR=strcat('SNR= ',num2str(SNR_tar,'%01d'),'dB');
else
    target_SNR=strcat('SNR= ',num2str(SNR_tar,'%02d'),'dB');
end

K_box=(K_min:5:K_max).'; % figureの横軸のためのKの箱
LK=length(K_box);        % Kの箱の大きさ

sigma2 = 1/(10^(SNR_tar/10)); % noise power
a = sigma2*Nt;

% 出力ファイル名 with alpha and SNR in dB
%evfile1 = strcat(folder,fn1,num2str(snr(isnr),'%02d'),'dB_1000itr.csv');

E_ZF = zeros(SIMU, Nr, Nu);             % ZF-CIの固有値のSINR
E_BD = zeros(SIMU, Nr, Nu);             % BDの固有値のSINR
E_MMSE = zeros(SIMU, Nr, Nu);           % MMSE-CIの固有値のSINR
E_GMI1 = zeros(SIMU, Nr, Nu);           % GMI1の固有値のSINR
E_GMI2 = zeros(SIMU, Nr, Nu);           % GMI2の固有値のSINR
E_BMSN_BF = zeros(SIMU, Nr, Nu);        % BMSN-BFの固有値のSINR
E_BMSN_GE = zeros(SIMU, Nr, Nu);        % BMSN-GEの固有値のSINR

Eigs_ZF = zeros(SIMU, Nr, Nu);             % ZF-CIの固有値
Eigs_BD = zeros(SIMU, Nr, Nu);             % BDの固有値
Eigs_MMSE = zeros(SIMU, Nr, Nu);           % MMSE-CIの固有値
Eigs_GMI1 = zeros(SIMU, Nr, Nu);           % GMI1の固有値のSINR
Eigs_GMI2 = zeros(SIMU, Nr, Nu);           % GMI2の固有値のSINR
Eigs_BMSN_BF = zeros(SIMU, Nr, Nu);        % BMSN-BFの固有値
Eigs_BMSN_GE = zeros(SIMU, Nr, Nu);        % BMSN-GEの固有値

    
for ik = 1:LK
    
    K_tar = K_box(ik);

    for isimu = 1:SIMU
    
    % 伝搬チャネル行列のマルチパス波成分(NLoSチャネル)
    H_iid = (randn(Nu*Nr,Nt)+1j*randn(Nu*Nr,Nt))/sqrt(2);
    %伝搬チャネル行列の直接波成分(LoSチャネル)
    H_los = zeros(Nu*Nr,Nt);
    
    Theta_t = (rand-0.5)*360;   % ユーザ毎の送信角 (-180deg - 180deg)
    Theta_r = (rand(1,Nu)-0.5)*360; % ユーザ毎の受信角 (-180deg - 180deg)
    a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t*derad));  % 送信モードベクトル
    
    for n = 1 : Nu
        a_r = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad)); % ユーザ毎の受信モードベクトル
        H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t.';                % ユーザ毎のLoSチャネル行列
    end
      
    K = 10^(K_tar/10);
         
    % 伝搬チャネル行列 H=[sqrt(K/(K+1))*(LOS チャネル)]...
    %                   .+[sqrt(1/(K+1))*(NLOS チャネル)]
    H = sqrt(K/(K+1))*H_los + sqrt(1/(K+1))*H_iid;
    
    % ZF-CI algorithm
    [W_ZF,U_ZF,S_ZF,~,~] = zf(Nt,Nr,Nu,H); % function zf.m を使用
    
    % BD algorithm
    [W_BD,U_BD,S_BD,~,~] = bd(Nt,Nr,Nu,H); % function bd.m を使用
    
    % MMSE-CI algorithm
    [W_MMSE,U_MMSE,S_MMSE,RIPM,~] = mmse(Nt,Nr,Nu,H,a); % function mmse.m を使用
    
    % GMI1 algorithm
    [W_GMI1,U_GMI1,S_GMI1,RIPGM1,~] = gmmse_m1(Nt,Nr,Nu,H,a); % function gmmse.m1 を使用
    
    % GMI2 algorithm
    [W_GMI2,U_GMI2,S_GMI2,RIPGM2,~] = gmmse_m2(Nt,Nr,Nu,H,a); % function gmmse.m2 を使用
    
    % BMSN-BF algorithm
    [W_BMSN_BF,U_BMSN_BF,S_BMSN_BF,RIPBF,~] = bmsn_bf(Nt,Nr,Nu,H,a,T); % function bmsn_bf.m を使用
      
    % BMSN-GE algorithm
    [W_BMSN_GE,U_BMSN_GE,S_BMSN_GE,RIPGE,~] = bmsn_gev(Nt,Nr,Nu,H,a); % function bmsn_gev.m を使用
    
    
    % ユーザ毎の固有値分布
     for nuser=1:Nu
        if Nr==1
            E_ZF(isimu,:,nuser) = 10*log10(S_ZF(1,1,nuser).^2/(Nt*snt));
            E_BD(isimu,:,nuser) = 10*log10(S_BD(1,1,nuser).^2/(Nt*snt));
            E_MMSE(isimu,:,nuser) = 10*log10((S_MMSE(1,1,nuser).^2)./(RIPM(1,nuser)+Nt*snt));
            E_GMI1(isimu,:,nuser) = 10*log10((S_GMI1(1,1,nuser).^2)./(RIPGM1(1,nuser)+Nt*snt));
            E_GMI2(isimu,:,nuser) = 10*log10((S_GMI2(1,1,nuser).^2)./(RIPGM2(1,nuser)+Nt*snt));
            E_BMSN_BF(isimu,:,nuser) = 10*log10((S_BMSN_BF(1,1,nuser).^2)./(RIPBF(1,nuser)+Nt*snt));
            E_BMSN_GE(isimu,:,nuser) = 10*log10((S_BMSN_GE(1,1,nuser).^2)./(RIPGE(1,nuser)+Nt*snt));
        else 
            E_ZF(isimu,:,nuser) = 10*log10(diag(S_ZF(:,:,nuser)).^2/(Nt*snt));   
            E_BD(isimu,:,nuser) = 10*log10(diag(S_BD(:,:,nuser)).^2/(Nt*snt));
            E_MMSE(isimu,:,nuser) = 10*log10((diag(S_MMSE(:,:,nuser)).^2)./(RIPM(:,nuser)+Nt*snt));
            E_GMI1(isimu,:,nuser) = 10*log10((diag(S_GMI1(:,:,nuser)).^2)./(RIPGM1(:,nuser)+Nt*snt));
            E_GMI2(isimu,:,nuser) = 10*log10((diag(S_GMI2(:,:,nuser)).^2)./(RIPGM2(:,nuser)+Nt*snt));
            E_BMSN_BF(isimu,:,nuser) = 10*log10((diag(S_BMSN_BF(:,:,nuser)).^2)./(RIPBF(:,nuser)+Nt*snt));
            E_BMSN_GE(isimu,:,nuser) = 10*log10((diag(S_BMSN_GE(:,:,nuser)).^2)./(RIPGE(:,nuser)+Nt*snt));
        end % Nr end
    end % nuser end  
    end % isimu end
    
    % ソーティング（昇順）
    for nuser=1:Nu
        E_ZF(:,:,nuser) = sort(E_ZF(:,:,nuser),1);
        E_BD(:,:,nuser) = sort(E_BD(:,:,nuser),1);
        E_MMSE(:,:,nuser) = sort(E_MMSE(:,:,nuser),1);
        E_GMI1(:,:,nuser) = sort(E_GMI1(:,:,nuser),1);
        E_GMI2(:,:,nuser) = sort(E_GMI2(:,:,nuser),1);
        E_BMSN_BF(:,:,nuser) = sort(E_BMSN_BF(:,:,nuser),1);
        E_BMSN_GE(:,:,nuser) = sort(E_BMSN_GE(:,:,nuser),1);
    end
    % ユーザ平均
    E_ZFm = mean(E_ZF,3);
    E_BDm = mean(E_BD,3);
    E_MMSEm = mean(E_MMSE,3);
    E_GMI1m = mean(E_GMI1,3);
    E_GMI2m = mean(E_GMI2,3);
    E_BMSN_BFm = mean(E_BMSN_BF,3);
    E_BMSN_GEm = mean(E_BMSN_GE,3);
    
    % ターゲットCDF値の固有値を抽出
    for nn=1:Nr
        rr3(ik,nn) = E_BMSN_BFm(round(CDF*SIMU/100),nn);     % BMSN_BFの2つの固有値が入った箱
        rr3(ik,nn+Nr) = E_BMSN_GEm(round(CDF*SIMU/100),nn);  % BMSN_GEの2つの固有値が入った箱
        rr3(ik,nn+2*Nr) = E_MMSEm(round(CDF*SIMU/100),nn);   % MMSEの2つの固有値が入った箱
        rr3(ik,nn+3*Nr) = E_GMI1m(round(CDF*SIMU/100),nn);
        rr3(ik,nn+4*Nr) = E_GMI2m(round(CDF*SIMU/100),nn);
        rr3(ik,nn+5*Nr) = E_BDm(round(CDF*SIMU/100),nn);     % BDの2つの固有値が入った箱
        rr3(ik,nn+6*Nr) = E_ZFm(round(CDF*SIMU/100),nn);     % ZFの2つの固有値が入った箱
    end
    fprintf('K = %d dB\n',K_box(ik));  % 解析中の計算過程を見るためのK = ??dBの表示
end% k_dB end

%% 
figure;
mycol = [0 0 1;
         0 0 1;
         1 0 0;
         1 0 0;
         0 0.7 0;
         0 0.7 0;
         1 0 1;
         1 0 1;
         0 0 0;
         0 0 0]; % グラフの色
set(groot,'defaultAxesColorOrder',mycol) % figureの書式を設定
axis([K_min K_max SNR_tar-60 SNR_tar]);              % figureの軸のとり方を設定
grid on;
hold on;
set(gca,'XTick',K_min:5:K_max,'Fontsize',14,'Fontname','Arial') % figureの書式を設定
plot(K_box,rr3(:,1),'r-d','LineWidth',2);
plot(K_box,rr3(:,2),'r--d','LineWidth',2);
plot(K_box,rr3(:,3),'b-s','LineWidth',2);
plot(K_box,rr3(:,4),'b--s','LineWidth',2);
plot(K_box,rr3(:,5),'g-o','LineWidth',2);
plot(K_box,rr3(:,6),'g--o','LineWidth',2);
plot(K_box,rr3(:,7),'y-+','LineWidth',2);
plot(K_box,rr3(:,8),'y--+','LineWidth',2);
plot(K_box,rr3(:,9),'k-*','LineWidth',2);
plot(K_box,rr3(:,10),'k--*','LineWidth',2);
plot(K_box,rr3(:,11),'c-x','LineWidth',2);
plot(K_box,rr3(:,12),'c--x','LineWidth',2);
plot(K_box,rr3(:,13),'m-^','LineWidth',2);
plot(K_box,rr3(:,14),'m--^','LineWidth',2);

% plot(K_box,rr3(:,1),'r-d','MarkerIndices',5:5:length(rr3(:,1)),'LineWidth',2);
% plot(K_box,rr3(:,2),'r--d','MarkerIndices',10:5:length(rr3(:,2)),'LineWidth',2);
% plot(K_box,rr3(:,3),'b-s','MarkerIndices',5:5:length(rr3(:,3)),'LineWidth',2);
% plot(K_box,rr3(:,4),'b--s','MarkerIndices',10:5:length(rr3(:,4)),'LineWidth',2);
% plot(K_box,rr3(:,5),'g-o','MarkerIndices',5:5:length(rr3(:,5)),'LineWidth',2);
% plot(K_box,rr3(:,6),'g--o','MarkerIndices',10:5:length(rr3(:,6)),'LineWidth',2);
% plot(K_box,rr3(:,7),'c-x','MarkerIndices',5:5:length(rr3(:,7)),'LineWidth',2);
% plot(K_box,rr3(:,8),'c--x','MarkerIndices',10:5:length(rr3(:,8)),'LineWidth',2);
% plot(K_box,rr3(:,9),'m-^','MarkerIndices',5:5:length(rr3(:,9)),'LineWidth',2);
% plot(K_box,rr3(:,10),'m--^','MarkerIndices',10:5:length(rr3(:,10)),'LineWidth',2);

legend('BMSN-BF \lambda_1','BMSN-BF \lambda_2','BMSN-GE \lambda_1','BMSN-GE \lambda_2','MMSE-CI \lambda_1','MMSE-CI \lambda_2',...
    'GMI1 \lambda_1','GMI1 \lambda_2','GMI2 \lambda_1','GMI2 \lambda_2','BD \lambda_1','BD \lambda_2','ZF-CI \lambda_1','ZF-CI \lambda_2','Location','southeast');
xlabel('K [dB]','Fontsize',16,'Fontname','Arial');
ylabel('SINR of eigenvalue [dB]','Fontsize',16,'Fontname','Arial');
set(gca,'Fontsize',16,'Fontname','Arial');
title(strcat(target_CDF,", ",target_SNR));
grid on;
hold on;

%End