% cap_bd_muser.m
% NU-ユーザでのBD法の計算
% TDMA, Upper boundと比較

clear;%close;
overBD_switch = 0;% 1:on,0:off overBD figureの有無

% パラメータ条件 NT >= NR*NU
K_min = -20;      % 最小K[dB]
K_max = 20;       % 最大K[dB]

CDF = 10;         % ABRのターゲットCDF
SNR_tar = 10;     % ターゲットSNR[dB]
Nt = 16;          % 送信素子数
Nr = 2;           % 受信素子数 (アンテナ選択の場合1, そうでない場合は2)
Nu = 8;           % ユーザ数
SIMU = 10;      % 試行回数（通常 1000)

An = 2; % 指向性関数の係数(cosパターンの場合は2)

I = eye(Nt,Nt);
Nru = Nr*Nu;

d_t = 0.5;   % 送信アンテナ間隔（in wavelength)
d_r = 0.5;        % 受信アンテナ間隔（in wavelength)
derad = pi/180;   % degree -> rad
K_box=(K_min:5:K_max).'; % figureの横軸のためのKの箱
LK=length(K_box);        % Kの箱の大きさ

%T 所望のチャネル行列
T = zeros(Nr,Nr);
for inu = 1:Nu
    T(:,:,inu) = eye(Nr,Nr);
end
Algorithms = ["BMSN-BF","BMSN-GE","MMSE-CI","GMI1","GMI2","BD","ZF-CI","BMSN-GE3"];
MSt = zeros(SIMU, Nr*Nu, numel(Algorithms)-1);
if Nr > 1
    MSt_stream1 = zeros(SIMU, (Nr-1)*Nu);
else
    MSt_stream1 = zeros(SIMU, Nr*Nu);
end
Q = zeros(SIMU, Nu, numel(Algorithms));
QmC = zeros(LK, numel(Algorithms));
QmCcos = zeros(LK, numel(Algorithms));
% 入力ファイル基幹名
% folder='CSV3/';
% fn1 = 'Eig16x2x8u_BD_SNR';
% fn2 = 'Eig16x2x8u_BMSN1(GEV)_SNR';
% fn3 = 'Eig16x2x8u_BMSN2a_SNR';  % alpha = 1e-2
% fn4 = 'Eig16x2x8u_BMSN2b_SNR';  % alpha = 1e-6
% fn5 = 'Eig16x2x8u_BMSN3_SNR';   % alpha = sigma2Nt equal to MMSE
% fn6 = 'Eig16x2x8u_ZF_SNR';
% fn7 = 'Eig16x2x8u_MMSE_SNR';

% 入力ファイル基幹名
% folder='CSV4/';
% fn1 = 'Eig16x2x8u_ZF_SNR';
% fn2 = 'Eig16x2x8u_BD_SNR';
% %fn3 = 'Eig16x2x8u_BMSN2_SNR'; % 2: alpha = 1e-2, 3: alpha=sigma2xNT
% %fn4 = 'Eig16x2x8u_BDAS_SNR';  % antenna selection
% %fn5 = 'Eig16x2x8u_BMSNAS2_SNR';  % alpha = 1e-2 and antenna selection
% fn3 = 'Eig16x2x8u_BMSN3_SNR';   % alpha = sigma2Nt equal to MMSE
% fn4 = 'Eig16x2x8u_BMSN3(GEV)_SNR';   % alpha = sigma2Nt equal to MMSE
% fn5 = 'Eig16x2x8u_MMSE_SNR';

% 出力ファイル名
%folder= 'CSV/';
%cdfn1 = 'ABRCDFvsSNR_16x2x8u_BD_MMSE_BMSNs_USERave_CDF';
%cdfn2 = 'Eig_CDF_16x2x8u_BMSN_SNR';
%%
if CDF < 10
   target_CDF=strcat('CDF=',num2str(CDF,'%01d'),'%');
else
   target_CDF=strcat('CDF=',num2str(CDF,'%02d'),'%');
end

if SNR_tar < 10
   target_SNR=strcat('SNR=',num2str(SNR_tar,'%01d'),'dB');
else
   target_SNR=strcat('SNR=',num2str(SNR_tar,'%02d'),'dB');
end

if d_t < 10
   target_d_t=strcat('d_t=',num2str(d_t,'%.3g'),'\lambda');
else
   target_d_t=strcat('d_t=',num2str(d_t,'%02d'),'\lambda');
end

if d_r < 10
   target_d_r=strcat('d_r=',num2str(d_r,'%.3g'),'\lambda');
else
   target_d_r=strcat('d_r=',num2str(d_r,'%02d'),'\lambda');
end
% 出力ファイル名 with SNR in dB
%cdffile1 = strcat(folder,cdfn1,num2str(CDF,'%02d'),'_1000itr.csv');
%cdffile2 = strcat(folder,cdfn2,num2str(SN_tar,'%02d'),'dB_1000itr.csv');
%%
H_los = zeros(Nu*Nr,Nt); % 伝搬チャネル行列の直接波成分(LOS チャネル)
sigma2 = 1/(10^(SNR_tar/10)); % noise power
a = sigma2*Nt; % 擬似雑音

St = zeros(Nr,Nu);
St_GE3 = zeros(Nr-1,Nu);
% MSt_BMSN_GE3 = zeros(SIMU, Nr*Nu);
% MSt_BMSN_BF = zeros(SIMU, Nr*Nu);
% MSt_BMSN_GE = zeros(SIMU, Nr*Nu);
% MSt_MMSE = zeros(SIMU, Nr*Nu);
% MSt_GMI1 = zeros(SIMU, Nr*Nu);
% MSt_GMI2 = zeros(SIMU, Nr*Nu);
% MSt_BD = zeros(SIMU, Nr*Nu);
% MSt_ZF = zeros(SIMU, Nr*Nu);
%%
tic;
for Directivity_switch = 0:1 % 1:on,0:off 送受信素子の指向性考慮の有無
    for ik = 1:LK
    
    K_tar = K_box(ik);

% 出力ファイル名 with alpha and SNR in dB
%evfile1 = strcat(folder,fn1,num2str(snr(isnr),'%02d'),'dB_1000itr.csv');
    
    for isimu = 1:SIMU
        
        % 伝搬チャネル行列のマルチパス成分 (i.i.d.Rayleigh)
        H_iid = (randn(Nr*Nu,Nt)+1j*randn(Nr*Nu,Nt))/sqrt(2);
        % LOS チャネル
        
        if Directivity_switch == 1
            Theta_t = (rand(1,Nu)-0.5)*180; % ユーザ毎の送信角 指向性:(-90deg - 90deg)
            Theta_r = (rand(1,Nu)-0.5)*180; % ユーザ毎の受信角 指向性:(-90deg - 90deg)
            for n = 1 : Nu
                a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad))*An*cos(Theta_t(1,n)*derad); % ユーザ毎の送信モードベクトル
                a_r = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad))*An*cos(Theta_r(1,n)*derad); % ユーザ毎の受信モードベクトル
                H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t';                 % ユーザ毎のLOSチャネル行列
            end
        else
            Theta_t = (rand(1,Nu)-0.5)*360; % ユーザ毎の送信角 等方性:(-180deg - 180deg) 
            Theta_r = (rand(1,Nu)-0.5)*360; % ユーザ毎の受信角 等方性:(-180deg - 180deg) 
            for n = 1 : Nu
            a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad)); % ユーザ毎の送信モードベクトル
            a_r = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad)); % ユーザ毎の受信モードベクトル
            H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t';                 % ユーザ毎のLOSチャネル行列
            end
        end
    % Kを真値にする
    K = 10^(K_tar/10);
      
    % H=[sqrt(K/(K+1))*(LOS チャネル)]+[sqrt(1/(K+1))*(NLOS チャネル)]
    H = sqrt(K/(K+1))*H_los + sqrt(1/(K+1))*H_iid;
    
    %% BMSN-BF
    [~,~,STT,RIP,~] = bmsn_bf(Nt,Nr,Nu,H,a,T);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    % MSt_BMSN_BF(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    MSt(isimu,:,strcmp("BMSN-BF",Algorithms))=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
 
    %% BMSN-GE
    [~,~,STT,RIP,~] = bmsn_gev(Nt,Nr,Nu,H,a);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    % MSt_BMSN_GE(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    MSt(isimu,:,strcmp("BMSN-GE",Algorithms))=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
% Stge = St
%     RIP
%     reshape(St,[Nru,1])
%     reshape(RIP,[Nru,1])
%     ((reshape(St,[Nru,1]).').^2)
% (reshape(RIP,[Nru,1]).'+Nt*sigma2)
    %% MMSE-CI
    [~,~,STT,RIP,~] = mmse(Nt,Nr,Nu,H,a);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    % MSt_MMSE(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    MSt(isimu,:,strcmp("MMSE-CI",Algorithms))=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    %% GMI1
    [~,~,STT,RIP,~] = gmmse_m1(Nt,Nr,Nu,H,a);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    % MSt_GMI1(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    MSt(isimu,:,strcmp("GMI1",Algorithms))=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    %% GMI2
    [~,~,STT,RIP,~] = gmmse_m2(Nt,Nr,Nu,H,a);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    % MSt_GMI2(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    MSt(isimu,:,strcmp("GMI2",Algorithms))=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    %% BD
    [~,~,STT,RIP,~] = bd(Nt,Nr,Nu,H);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    % MSt_BD(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    MSt(isimu,:,strcmp("BD",Algorithms))=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);  
    %% ZF-CI
    [~,~,STT,RIP,~] = zf(Nt,Nr,Nu,H);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end    
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    % MSt_ZF(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    MSt(isimu,:,strcmp("ZF-CI",Algorithms))=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    %% BMSN-GE3
    [~,~,STT,RIP,~] = bmsn_ge3(Nt,Nr,Nu,H,a);
    for inu=1:Nu

% St_3 = St
%     RIP
%     ((St_GE3.').^2)
%     RIP(1,:).'
% RIP(1,:)+Nt*sigma2
    % fprintf('Iteration = %d / %d\n',isimu, SIMU);    
    end % isimu
    
    %% 各アルゴリズムのChannel Capacity
    for nuser=1:Nu
        ns = Nr*(nuser-1)+1:Nr*nuser;
        ns_GE3 = nuser; % BMSN-GE3用のアンテナ数
        
        for n_alg = 1:numel(Algorithms)-1
            Q(:,nuser,n_alg)=sort(sum(log2(1 + MSt(:,ns,n_alg)),2)); % 各アルゴリズムのChannel capacityを計算
            % Q(:,nuser,strcmp(" ",Algorithms))=sort(sum(log2(1 + MSt(:,ns,strcmp(" ",Algorithms))),2)); % 各アルゴリズムのChannel capacityを計算
        end
        Q(:,nuser,strcmp("BMSN-GE3",Algorithms))=sort(sum(log2(1 + MSt_stream1(:,ns_GE3)),2)); % BMSN-GE3のChannel capacityを計算
    end % nuser end

    Qm = mean(Q,2); % Qのユーザ回数平均
    if Directivity_switch == 0
        QmC(ik,:)=Qm(round(CDF*SIMU/100),:); % Channel capacity
    else
        QmCcos(ik,:)=Qm(round(CDF*SIMU/100),:); % cos Channel capacity
    end
    
    QmC_overBD=QmC(:,:)./QmC(:,strcmp("BD",Algorithms));% C / C_BD

    fprintf('K = %d dB\n',K_box(ik));
    end % ik end
    
    %csvwrite(cdffile1,[K,QmC]);
end % Directivity end
toc;
%% グラフ表示 figure ABRvsK
figure;
mycol = [1 0 0;0 0 1;0 1 0;0 1 1;1 0 1;0 0 0;1 1 0;0.85 0.325 0.098
         1 0 0;0 0 1;0 1 0;0 1 1;1 0 1;0 0 0;1 1 0;0.85 0.325 0.098]; % 色
set(groot,'defaultAxesColorOrder',mycol)
for n_alg = 1:numel(Algorithms)
    plot(K_box,QmC(:,n_alg),'--','Linewidth',2);
    grid on; hold on;
end
for n_alg = 1:numel(Algorithms)
    plot(K_box,QmCcos(:,n_alg),'-','Linewidth',2);
    grid on; hold on;
end
%plot(K_box,QmC(:,1),'r-',K_box,QmC(:,2),'b-',K_box,QmC(:,3),'g-',K_box,QmC(:,4),'y-',K_box,QmC(:,5),'k-',K_box,QmC(:,6),'c-',K_box,QmC(:,7),'m-','Linewidth',2);

axis([K_min,K_max,0,max(max([QmC QmCcos]))+2])
set(gca,'Fontsize',18,'Fontname','Times New Roman');
lgd = legend;
lgd.NumColumns = numel(Algorithms)/4; % 凡例の列数を指定
legend([Algorithms strcat(Algorithms,"cos")],'Location','Northwest');

xlabel('{\it{K}} [dB]','Fontsize',6,'Fontname','Times New Roman');
ylabel('Channel capacity [bits/s/Hz]','Fontsize',6,'Fontname','Times New Roman');
set(gca,'Fontsize',18,'Fontname','Times New Roman');
title(strcat(target_CDF,',',target_SNR,',',target_d_t,',',target_d_r));
grid on; hold on;
%% グラフ表示 figure Channel Capacity overBD vs K
if overBD_switch == 1
    figure;
    mycol = [1 0 1;0 1 0;1 0 0;0 0 1;0 1 1;0 0 0]; % 色
    set(groot,'defaultAxesColorOrder',mycol)
    %plot(K,QmC(:,1),K,QmC(:,2),K,QmC(:,3),K,QmC(:,4),'Linewidth',2);
    
    plot(K_box,QmC_overBD(:,1),'r-',K_box,QmC_overBD(:,2),'b-',K_box,QmC_overBD(:,3),'g-',K_box,QmC_overBD(:,4),'y-',K_box,QmC_overBD(:,5),'k-',K_box,QmC_overBD(:,7),'m-','Linewidth',2);
    
    axis([K_min,K_max,0,5])
    set(gca,'Fontsize',18,'Fontname','Times New Roman');
    lgd = legend;
    lgd.NumColumns = 3; % 凡例の列数を指定
    legend(Algorithms,'Location','Northwest');
    xlabel('{\it{K}} [dB]','Fontsize',6,'Fontname','Times New Roman');
    ylabel('C / C_{BD}','Fontsize',6,'Fontname','Times New Roman');
    set(gca,'Fontsize',18,'Fontname','Times New Roman');
    title(strcat(target_CDF,',',target_SNR,',',target_d_t,',',target_d_r));
    grid on; hold on;
end
% End