% K_lambda1,lambda2
% 固有モード伝送後の特異値を固有値としている

clear;

% K[dB]の範囲設定
K_min = -20;        % 最小K [dB]
K_max = 20;         % 最大K [dB]
weight=0;   % weight figure on=0,off=1

Nt = 16;    % 送信素子数
Nr = 1;     % 受信素子数 (アンテナ選択の場合1, そうでない場合は2)
Nu = 16;     % ユーザ数

SIMU = 1000; % 試行回数（通常 1000）

I = eye(Nt,Nt); % Nt*Ntの単位行列

SNR_tar = 10;   % ターゲットSNR[dB]

d_t      = 0.5;      % 送信アンテナ間隔（in wavelength)
d_r      = 0.5;      % 受信アンテナ間隔（in wavelength)
derad = pi/180;      % degree -> rad

%T 所望のチャネル行列
T = zeros(Nr,Nr);
for inu = 1:Nu
    T(:,:,inu) = eye(Nr,Nr);
end
%%
if SNR_tar < 10
    target_SNR=strcat('SNR= ',num2str(SNR_tar,'%01d'),'dB');
else
    target_SNR=strcat('SNR= ',num2str(SNR_tar,'%02d'),'dB');
end
H_isimu=zeros(Nu*Nr,Nt,SIMU);
H_iid_isimu=zeros(Nu*Nr,Nt,SIMU);
H_los_isimu=zeros(Nu*Nr,Nt,SIMU);

Norm_H=zeros(Nu*Nr,1,SIMU);      %ユーザチャネルのノルム
Norm_H_los=zeros(Nu*Nr,1,SIMU);
Norm_H_iid=zeros(Nu*Nr,1,SIMU);

HW_BMSN_BF=zeros(Nr*Nu,Nr*Nu,SIMU);%HWの生成
HW_BMSN_GE=zeros(Nr*Nu,Nr*Nu,SIMU);
HW_MMSE=zeros(Nr*Nu,Nr*Nu,SIMU);
HW_BD=zeros(Nr*Nu,Nr*Nu,SIMU);
HW_ZF=zeros(Nr*Nu,Nr*Nu,SIMU);

HW_BMSN_BFn=zeros(1,SIMU);%HWノルムの生成
HW_BMSN_GEn=zeros(1,SIMU);
HW_MMSEn=zeros(1,SIMU);
HW_BDn=zeros(1,SIMU);
HW_ZFn=zeros(1,SIMU);

Norm_W_ZF_whole=zeros(1,SIMU);%全体Wノルムの生成
Norm_W_BD_whole=zeros(1,SIMU);
Norm_W_MMSE_whole=zeros(1,SIMU);
Norm_W_BMSN_BF_whole=zeros(1,SIMU);
Norm_W_BMSN_GE_whole=zeros(1,SIMU);

Norm_W_ZF_user=zeros(Nu,SIMU);%ユーザWノルムの生成
Norm_W_BD_user=zeros(Nu,SIMU);
Norm_W_MMSE_user=zeros(Nu,SIMU);
Norm_W_BMSN_BF_user=zeros(Nu,SIMU);
Norm_W_BMSN_GE_user=zeros(Nu,SIMU);

Norm_W_ZF_user_avsimu=zeros(1,Nu);%ユーザWノルム試行回数平均の生成
Norm_W_BD_user_avsimu=zeros(1,Nu);
Norm_W_MMSE_user_avsimu=zeros(1,Nu);
Norm_W_BMSN_BF_user_avsimu=zeros(1,Nu);
Norm_W_BMSN_GE_user_avsimu=zeros(1,Nu);
%%
K_box=(K_min:5:K_max).'; % figureの横軸のためのKの箱
LK=length(K_box);        % Kの箱の大きさ
sigma2 = 1/(10^(SNR_tar/10)); % noise power
a = sigma2*Nt;
 dd3=zeros(LK,5); hh3=zeros(LK,5); ii3=zeros(LK,5); 

% 出力ファイル名 with alpha and SNR in dB
%evfile1 = strcat(folder,fn1,num2str(snr(isnr),'%02d'),'dB_1000itr.csv');

W_ZF_isimu = zeros(Nt,Nr,Nu,SIMU);     %Wの作成        
W_BD_isimu = zeros(Nt,Nr,Nu,SIMU);             
W_MMSE_isimu = zeros(Nt,Nr,Nu,SIMU);           
W_BMSN_BF_isimu = zeros(Nt,Nr,Nu,SIMU);       
W_BMSN_GE_isimu = zeros(Nt,Nr,Nu,SIMU);     

W_ZF_whole=zeros(Nt,Nr*Nu,SIMU);            %W全体の作成
W_BD_whole=zeros(Nt,Nr*Nu,SIMU);
W_MMSE_whole=zeros(Nt,Nr*Nu,SIMU);
W_BMSN_BF_whole=zeros(Nt,Nr*Nu,SIMU);
W_BMSN_GE_whole=zeros(Nt,Nr*Nu,SIMU);
%%
for ik = 1:LK
    K_tar = K_box(ik);
    
    for isimu = 1:SIMU
    
        H_iid = (randn(Nu*Nr,Nt)+1j*randn(Nu*Nr,Nt))/sqrt(2); % 伝搬チャネル行列のマルチパス波成分(NLoSチャネル)
        H_los = zeros(Nu*Nr,Nt);                              % 伝搬チャネル行列の直接波成分(LoSチャネル)
    
        Theta_t = (rand(1,Nu)-0.5)*360;   % ユーザ毎の送信角 (-180deg - 180deg)
        Theta_r = (rand(1,Nu)-0.5)*360;   % ユーザ毎の受信角 (-180deg - 180deg)  % a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t*derad));  % 送信モードベクトル
    
    
        for n = 1 : Nu
            a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad)); % ユーザ毎の送信モードベクトル
            a_r = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad)); % ユーザ毎の受信モードベクトル
            H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t';                 % ユーザ毎のLOSチャネル行列
        end
      
    K = 10^(K_tar/10);
         
    H = sqrt(K/(K+1))*H_los + sqrt(1/(K+1))*H_iid; % 伝搬チャネル行列 H=[sqrt(K/(K+1))*(LOS チャネル)]+[sqrt(1/(K+1))*(NLOS チャネル)]
   
    % ZF-CI algorithm
    [W_ZF,U_ZF,S_ZF,~,~] = zf(Nt,Nr,Nu,H); % function zf.m を使用
    
    % BD algorithm
    [W_BD,U_BD,S_BD,~,~] = bd(Nt,Nr,Nu,H); % function bd.m を使用
    
    % MMSE-CI algorithm
    [W_MMSE,U_MMSE,S_MMSE,RIPM,~] = mmse(Nt,Nr,Nu,H,a); % function mmse.m を使用
    
    % BMSN-BF algorithm
    [W_BMSN_BF,U_BMSN_BF,S_BMSN_BF,RIPBF,~] = bmsn_bf(Nt,Nr,Nu,H,a,T); % function bmsn_bf.m を使用
      
    % BMSN-GE algorithm
    [W_BMSN_GE,U_BMSN_GE,S_BMSN_GE,RIPGE,~] = bmsn_gev(Nt,Nr,Nu,H,a); % function bmsn_gev.m を使用
    
    W_ZF_isimu(:,:,:,isimu)=W_ZF;          
    W_BD_isimu(:,:,:,isimu) = W_BD;             
    W_MMSE_isimu(:,:,:,isimu) = W_MMSE;           
    W_BMSN_BF_isimu(:,:,:,isimu) =W_BMSN_BF;       
    W_BMSN_GE_isimu(:,:,:,isimu) =W_BMSN_GE;
    
    H_isimu(:,:,isimu)=H;
    H_iid_isimu(:,:,isimu)=H_iid;
    H_los_isimu(:,:,isimu)=H_los;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nuser=1:Nu
    ns = Nr*(nuser-1)+1:Nr*nuser;
    HW_BMSN_BF(ns,ns,isimu)=H_isimu(ns,:,isimu)*W_BMSN_BF_isimu(:,:,nuser,isimu);
    HW_BMSN_GE(ns,ns,isimu)=H_isimu(ns,:,isimu)*W_BMSN_GE_isimu(:,:,nuser,isimu);
    HW_MMSE(ns,ns,isimu)=H_isimu(ns,:,isimu)*W_MMSE_isimu(:,:,nuser,isimu);
    HW_BD(ns,ns,isimu)=H_isimu(ns,:,isimu)*W_BD_isimu(:,:,nuser,isimu);
    HW_ZF(ns,ns,isimu)=H_isimu(ns,:,isimu)*W_ZF_isimu(:,:,nuser,isimu);
end

for n = 1 : Nu*Nr
        Norm_H(n,:,isimu)=norm(H_isimu(n,:,isimu));   % HWの作成
        Norm_H_iid(n,:,isimu)=norm(H_iid_isimu(n,:,isimu));
        Norm_H_los(n,:,isimu)=norm(H_los_isimu(n,:,isimu));
end
    HW_BMSN_BFn(1,isimu)=norm(HW_BMSN_BF(:,:,isimu));  % HWのノルム
    HW_BMSN_GEn(1,isimu)=norm(HW_BMSN_GE(:,:,isimu));
    HW_MMSEn(1,isimu)=norm(HW_MMSE(:,:,isimu));
    HW_BDn(1,isimu)=norm(HW_BD(:,:,isimu));
    HW_ZFn(1,isimu)=norm(HW_ZF(:,:,isimu));
    
for nuser=1:Nu
    ns = Nr*(nuser-1)+1:Nr*nuser;    
    W_ZF_whole(:,ns,isimu)=W_ZF_isimu(:,:,nuser,isimu);
    W_BD_whole(:,ns,isimu)=W_BD_isimu(:,:,nuser,isimu);
    W_MMSE_whole(:,ns,isimu)=W_MMSE_isimu(:,:,nuser,isimu);
    W_BMSN_BF_whole(:,ns,isimu)=W_BMSN_BF_isimu(:,:,nuser,isimu);
    W_BMSN_GE_whole(:,ns,isimu)=W_BMSN_GE_isimu(:,:,nuser,isimu);
end
    Norm_W_ZF_whole(1,isimu)=norm(W_ZF_whole(:,:,isimu));% W全体のノルム
    Norm_W_BD_whole(1,isimu)=norm(W_BD_whole(:,:,isimu));
    Norm_W_MMSE_whole(1,isimu)=norm(W_MMSE_whole(:,:,isimu));
    Norm_W_BMSN_BF_whole(1,isimu)=norm(W_BMSN_BF_whole(:,:,isimu));
    Norm_W_BMSN_GE_whole(1,isimu)=norm(W_BMSN_GE_whole(:,:,isimu));
    
for nuser=1:Nu    
    Norm_W_ZF_user(nuser,isimu)=norm(W_ZF_isimu(:,:,nuser,isimu));% ユーザごとのWノルム
    Norm_W_BD_user(nuser,isimu)=norm(W_BD_isimu(:,:,nuser,isimu));
    Norm_W_MMSE_user(nuser,isimu)=norm(W_MMSE_isimu(:,:,nuser,isimu));
    Norm_W_BMSN_BF_user(nuser,isimu)=norm(W_BMSN_BF_isimu(:,:,nuser,isimu));
    Norm_W_BMSN_GE_user(nuser,isimu)=norm(W_BMSN_GE_isimu(:,:,nuser,isimu));
end
    end   % isimu end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    W_ZF_isimu_av=mean(W_ZF_isimu,4);% 試行回数の平均
    W_BD_isimu_av=mean(W_BD_isimu,4);
    W_MMSE_isimu_av=mean(W_MMSE_isimu,4);
    W_BMSN_BF_isimu_av=mean(W_BMSN_BF_isimu,4);
    W_BMSN_GE_isimu_av=mean(W_BMSN_GE_isimu,4);
    
    W_ZF_user_av=mean(W_ZF_isimu_av,3);% ユーザの平均
    W_BD_user_av=mean(W_BD_isimu_av,3);
    W_MMSE_user_av=mean(W_MMSE_isimu_av,3);
    W_BMSN_BF_user_av=mean(W_BMSN_BF_isimu_av,3);
    W_BMSN_GE_user_av=mean(W_BMSN_GE_isimu_av,3);
    
    Norm_W_ZF=mean(Norm_W_ZF_whole,2);% W全体のノルムの試行回数平均
    Norm_W_BD=mean(Norm_W_BD_whole,2);
    Norm_W_MMSE=mean(Norm_W_MMSE_whole,2);
    Norm_W_BMSN_BF=mean(Norm_W_BMSN_BF_whole,2);
    Norm_W_BMSN_GE=mean(Norm_W_BMSN_GE_whole,2);
      
    Norm_W_ZF_user_avsimu=mean(Norm_W_ZF_user,2);% ユーザごとのWノルムの試行回数平均
    Norm_W_BD_user_avsimu=mean(Norm_W_BD_user,2);
    Norm_W_MMSE_user_avsimu=mean(Norm_W_MMSE_user,2);
    Norm_W_BMSN_BF_user_avsimu=mean(Norm_W_BMSN_BF_user,2);
    Norm_W_BMSN_GE_user_avsimu=mean(Norm_W_BMSN_GE_user,2);
    
    Norm_W_ZF_avsimu_user=mean(Norm_W_ZF_user_avsimu);% ユーザごとのWノルムの試行回数平均のあとさらにユーザ平均
    Norm_W_BD_avsimu_user=mean(Norm_W_BD_user_avsimu);
    Norm_W_MMSE_avsimu_user=mean(Norm_W_MMSE_user_avsimu);
    Norm_W_BMSN_BF_avsimu_user=mean(Norm_W_BMSN_BF_user_avsimu);
    Norm_W_BMSN_GE_avsimu_user=mean(Norm_W_BMSN_GE_user_avsimu);        
 
    Norm_HW_BMSN_BFn=mean(HW_BMSN_BFn,2);
    Norm_HW_BMSN_GEn=mean(HW_BMSN_GEn,2);
    Norm_HW_MMSEn=mean(HW_MMSEn,2);
    Norm_HW_BDn=mean(HW_BDn,2);
    Norm_HW_ZFn=mean(HW_ZFn,2);
    
    Norm_H_av=mean(Norm_H,3);
    Norm_H_iid_av=mean(Norm_H_iid,3);
    Norm_H_los_av=mean(Norm_H_los,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
        dd3(ik,1) = Norm_W_BMSN_GE_avsimu_user;     % 送信ウエイトBMSN_BFの平均
        dd3(ik,2) = Norm_W_BMSN_BF_avsimu_user;  % 送信ウエイトBMSN_GEの平均
        dd3(ik,3) = Norm_W_MMSE_avsimu_user;   % 送信ウエイトMMSEの平均
        dd3(ik,4) = Norm_W_BD_avsimu_user;     % 送信ウエイトBDの平均
        dd3(ik,5) = Norm_W_ZF_avsimu_user;     % 送信ウエイトZFの平均            
       
        hh3(ik,1) = Norm_HW_BMSN_BFn;     % BMSN_BFのHW
        hh3(ik,2) = Norm_HW_BMSN_GEn;  % BMSN_GEのHW
        hh3(ik,3) = Norm_HW_MMSEn;   % MMSEのHW
        hh3(ik,4) = Norm_HW_BDn;     % BDのHW
        hh3(ik,5) = Norm_HW_ZFn;     % ZFのHW
        
        ii3(ik,1) = Norm_W_BMSN_BF;     % BMSN_BFのW
        ii3(ik,2) = Norm_W_BMSN_GE;  % BMSN_GEのW
        ii3(ik,3) = Norm_W_MMSE;   % MMSEのW
        ii3(ik,4) = Norm_W_BD;     % BDのW
        ii3(ik,5) = Norm_W_ZF;     % ZFのW
    
    fprintf('K = %d dB\n',K_box(ik));  % 解析中の計算過程を見るためのK = ??dBの表示
 end% ik end

if weight==0
%% グラフ表示　送信ウエイト(ユーザあたり)
figure;
mycol = [0 0 1;1 0 0;0 0.7 0;1 0 1;0 0 0]; % グラフの色
set(groot,'defaultAxesColorOrder',mycol) % figureの書式を設定
axis([K_min K_max 0 2]);              % figureの軸のとり方を設定
grid on;
hold on;
set(gca,'XTick',K_min:5:K_max,'Fontsize',14,'Fontname','Times New Roman') % figureの書式を設定
plot(K_box,dd3(:,1),'r-d','LineWidth',4);
plot(K_box,dd3(:,2),'b-s','LineWidth',4);
plot(K_box,dd3(:,3),'g-o','LineWidth',4);
plot(K_box,dd3(:,4),'c-x','LineWidth',4);
plot(K_box,dd3(:,5),'m-^','LineWidth',4);

legend('BMSN-BF','BMSN-GE','MMSE-CI','BD','ZF-CI','Location','southeast');
xlabel('{\it{K}} [dB]','Fontsize',36,'Fontname','Times New Roman');
ylabel('norm {\bf{W}}_{user}','Fontsize',36,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
title(target_SNR);
grid on;
hold on;
%% グラフ表示　送信ウエイト
figure;
mycol = [0 0 1;1 0 0;0 0.7 0;1 0 1;0 0 0]; % グラフの色
set(groot,'defaultAxesColorOrder',mycol) % figureの書式を設定
axis([K_min K_max 0 3]);              % figureの軸のとり方を設定
grid on;
hold on;
set(gca,'XTick',K_min:5:K_max,'Fontsize',14,'Fontname','Times New Roman') % figureの書式を設定
plot(K_box,ii3(:,1),'r-d','LineWidth',4);
plot(K_box,ii3(:,2),'b-s','LineWidth',4);
plot(K_box,ii3(:,3),'g-o','LineWidth',4);
plot(K_box,ii3(:,4),'c-x','LineWidth',4);
plot(K_box,ii3(:,5),'m-^','LineWidth',4);

legend('BMSN-BF','BMSN-GE','MMSE-CI','BD','ZF-CI','Location','southeast');
xlabel('{\it{K}} [dB]','Fontsize',36,'Fontname','Times New Roman');
ylabel('norm \bf{W}','Fontsize',36,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
title(target_SNR);
grid on;
hold on;

end
%% グラフ表示
figure;
mycol = [0 0 1;1 0 0;0 0.7 0;1 0 1;0 0 0]; % グラフの色
set(groot,'defaultAxesColorOrder',mycol) % figureの書式を設定
axis([K_min K_max 0 6]);              % figureの軸のとり方を設定
grid on;
hold on;
set(gca,'XTick',K_min:5:K_max,'Fontsize',14,'Fontname','Times New Roman') % figureの書式を設定
plot(K_box,hh3(:,1),'r-d','LineWidth',4);
plot(K_box,hh3(:,2),'b-s','LineWidth',4);
plot(K_box,hh3(:,3),'g-o','LineWidth',4);
plot(K_box,hh3(:,4),'c-x','LineWidth',4);
plot(K_box,hh3(:,5),'m-^','LineWidth',4);

legend('BMSN-BF','BMSN-GE','MMSE-CI','BD','ZF-CI','Location','southeast');
xlabel('{\it{K}} [dB]','Fontsize',36,'Fontname','Times New Roman');
ylabel('norm \bf{HW}','Fontsize',36,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
title(target_SNR);
grid on;
hold on;
% End