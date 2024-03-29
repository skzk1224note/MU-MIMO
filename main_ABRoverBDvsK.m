% cap_bd_muser.m
% NU-ユーザでのBD法の計算
% TDMA, Upper boundと比較

clear;%close;

% パラメータ条件 NT >= NR*NU
K_min = -20;      % 最小K[dB]
K_max = 20;       % 最大K[dB]

CDF = 10;         % ABRのターゲットCDF
SNR_tar = 20;     % ターゲットSNR[dB]
Nt = 16;          % 送信素子数
Nr = 2;           % 受信素子数 (アンテナ選択の場合1, そうでない場合は2)
Nu = 8;
SIMU = 100;      % 試行回数（通常 5000）

I = eye(Nt,Nt);
Nru = Nr*Nu;

d_t = 0.5;        % 送信アンテナ間隔（in wavelength)
d_r = 0.5;        % 受信アンテナ間隔（in wavelength)
derad = pi/180;   % degree -> rad


%T 所望のチャネル行列
T = zeros(Nr,Nr);
for inu = 1:Nu
    T(:,:,inu) = eye(Nr,Nr);
end

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
folder= 'CSV/';
cdfn1 = 'ABRCDFvsSNR_16x2x8u_BD_MMSE_BMSNs_USERave_CDF';
%cdfn2 = 'Eig_CDF_16x2x8u_BMSN_SNR';

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

K_box=(K_min:5:K_max).'; % figureの横軸のためのKの箱
LK=length(K_box);        % Kの箱の大きさ

% 出力ファイル名 with SNR in dB
cdffile1 = strcat(folder,cdfn1,num2str(CDF,'%02d'),'_1000itr.csv');
%cdffile2 = strcat(folder,cdfn2,num2str(SN_tar,'%02d'),'dB_1000itr.csv');

% 伝搬チャネル行列のマルチパス成分 (i.i.d.Rayleigh)
H_iid = (randn(Nr*Nu,Nt)+1j*randn(Nr*Nu,Nt))/sqrt(2);
% 伝搬チャネル行列の直接波成分(LOS チャネル)
H_los = zeros(Nu*Nr,Nt);
    
sigma2 = 1/(10^(SNR_tar/10)); % noise power
a = sigma2*Nt;

for ik = 1:LK
    
    K_tar = K_box(ik);

% 出力ファイル名 with alpha and SNR in dB
%evfile1 = strcat(folder,fn1,num2str(snr(isnr),'%02d'),'dB_1000itr.csv');
    
    for isimu = 1:SIMU
        Theta_t = (rand(1,Nu)-0.5)*360;   % ユーザ毎の送信角 (-180deg - 180deg)
        % Theta_t = (rand-0.5)*360;   % ユーザ毎の送信角 (-180deg - 180deg)
        Theta_r = (rand(1,Nu)-0.5)*360; % ユーザ毎の受信角 (-180deg - 180deg)
        % a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t*derad));  % 送信モードベクトル
        for n = 1 : Nu
            a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad));  % 送信モードベクトル
            a_r = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad)); % ユーザ毎の受信モードベクトル
            H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t.';                % ユーザ毎のLOSチャネル行列
        end
        
    % Kを真値に  
    K = 10^(K_tar/10);
      
    % H=[sqrt(K/(K+1))*(LOS チャネル)]+[sqrt(1/(K+1))*(NLOS チャネル)]
    H = sqrt(K/(K+1))*H_los + sqrt(1/(K+1))*H_iid;
    
    Hu = peruser(H,Nu);
        
    % Hから1ユーザのチャネル行列を除いた行列
    He = zeros((Nu-1)*Nr,Nt,Nu);    

    % BMSN-BF
    [~,~,STT,RIP,~] = bmsn_bf(Nt,Nr,Nu,H,a,T);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    MSt_BMSN1(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
 
    % BMSN-GE
    [~,~,STT,RIP,~] = bmsn_gev(Nt,Nr,Nu,H,a);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    MSt_BMSN2(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
 
    % MMSE-CI
    [~,~,STT,RIP,~] = mmse(Nt,Nr,Nu,H,a);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    MSt_MMSE(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    
    %GMI1
    [~,~,STT,RIP,~] = gmmse_m1(Nt,Nr,Nu,H,a);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    MSt_GMI1(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    
    %GMI2
    [~,~,STT,RIP,~] = gmmse_m2(Nt,Nr,Nu,H,a);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    MSt_GMI2(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);

    % BD
    [~,~,STT,RIP,~] = bd(Nt,Nr,Nu,H);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    MSt_BD(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
        
    % ZF-CI
    [~,~,STT,RIP,~] = zf(Nt,Nr,Nu,H);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    MSt_ZF(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    
    % fprintf('Iteration = %d / %d\n',isimu, SIMU);
    
    end % isimu

    % 各アルゴリズムのABR
    for nuser=1:Nu
        if Nr==2
            Q(:,1,nuser)=sort(sum(log2(1 + MSt_BMSN1(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
            Q(:,2,nuser)=sort(sum(log2(1 + MSt_BMSN2(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
            Q(:,3,nuser)=sort(sum(log2(1 + MSt_MMSE(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
            Q(:,4,nuser)=sort(sum(log2(1 + MSt_GMI1(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
            Q(:,5,nuser)=sort(sum(log2(1 + MSt_GMI2(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
            Q(:,6,nuser)=sort(sum(log2(1 + MSt_BD(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
            Q(:,7,nuser)=sort(sum(log2(1 + MSt_ZF(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
        end
        
        if Nr==1
            Q(:,1,nuser)=sort(log2(1 + MSt_BMSN1(:,nuser*Nr)));
            Q(:,2,nuser)=sort(log2(1 + MSt_BMSN2(:,nuser*Nr)));
            Q(:,3,nuser)=sort(log2(1 + MSt_MMSE(:,nuser*Nr)));
            Q(:,4,nuser)=sort(log2(1 + MSt_GMI1(:,nuser*Nr)));
            Q(:,5,nuser)=sort(log2(1 + MSt_GMI2(:,nuser*Nr)));
            Q(:,6,nuser)=sort(log2(1 + MSt_BD(:,nuser*Nr)));
            Q(:,7,nuser)=sort(log2(1 + MSt_ZF(:,nuser*Nr)));
        end
    end % nuser end

    Qm = mean(Q,3);
    
    % ABR/ABR_BD
    QmC(ik,:)=Qm(round(CDF*SIMU/100),:);
    
    QmC_overBD=QmC(:,:)./QmC(:,6);

    fprintf('K = %d dB\n',K_box(ik));
end % ik end

%csvwrite(cdffile1,[K,QmC]);

%% CDF of EGV at Target SNR

% グラフ表示 figure ABRvsK
figure;
mycol = [1 0 1;
         0 1 0;
         1 0 0;
         0 0 1;
         0 1 1;
         0 0 0]; % 色
set(groot,'defaultAxesColorOrder',mycol)
%plot(K,QmC(:,1),K,QmC(:,2),K,QmC(:,3),K,QmC(:,4),'Linewidth',2);

plot(K_box,QmC(:,1),'r-d',K_box,QmC(:,2),'b-s',K_box,QmC(:,3),'g-o',K_box,QmC(:,4),'y-d',K_box,QmC(:,5),'y--d',K_box,QmC(:,6),'c-x',K_box,QmC(:,7),'m-^','Linewidth',2);

axis([K_min,K_max,0,15])
set(gca,'Fontsize',14,'Fontname','Times New Roman');
legend('BMSN-BF','BMSN-GE','MMSE-CI','GMI1','GMI2','BD','ZF-CI','Location','Northwest');
xlabel('{\it{K}} [dB]','Fontsize',16,'Fontname','Times New Roman');
ylabel('Achievable bit rate [bits/s/Hz]','Fontsize',16,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
title(strcat(target_CDF,',',target_SNR));
grid on;
hold on;
 
% グラフ表示 figure ABRoverBDvsK
figure;
mycol = [1 0 1;
         0 1 0;
         1 0 0;
         0 0 1;
         0 1 1;
         0 0 0]; % 色
set(groot,'defaultAxesColorOrder',mycol)
%plot(K,QmC(:,1),K,QmC(:,2),K,QmC(:,3),K,QmC(:,4),'Linewidth',2);

plot(K_box,QmC_overBD(:,1),'r-d',K_box,QmC_overBD(:,2),'b-s',K_box,QmC_overBD(:,3),'g-x',K_box,QmC_overBD(:,5),'m-^','Linewidth',2);

axis([K_min,K_max,0,3])
set(gca,'Fontsize',14,'Fontname','Times New Roman');
legend('BMSN-BF','BMSN-GE','MMSE-CI','ZF-CI','Location','Northwest');
xlabel('{\it{K}} [dB]','Fontsize',16,'Fontname','Times New Roman');
ylabel('ABR / ABR_B_D','Fontsize',16,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
title(strcat(target_CDF,',',target_SNR));
grid on;
hold on;

% End