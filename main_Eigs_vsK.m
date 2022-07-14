% K_lambda1,lambda2
% 固有モード伝送後の特異値を固有値としている

clear;

% K[dB]の範囲設定
K_min = -20;        % 最小K [dB]
K_max = 20;         % 最大K [dB]

CDF = 10;   % 注目CDF
wantCDF = 1;  % 0ならoff,1ならon(figure用)

Nt = 16;    % 送信素子数
Nr = 2;     % 受信素子数 (アンテナ選択の場合1, そうでない場合は2)
Nu = 8;     % ユーザ数

SIMU = 10; % 試行回数（通常 1000）
Directivity_switch = 1;% 1:on,0:off 送受信素子の指向性考慮の有無
An = 2; % 指向性関数の係数

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
%%
K_box=(K_min:5:K_max).'; % figureの横軸のためのKの箱
LK=length(K_box);        % Kの箱の大きさ
sigma2 = 1/(10^(SNR_tar/10)); % noise power
a = sigma2*Nt;
aa3=zeros(LK,10); bb3=zeros(LK,10); cc3=zeros(LK,5); dd3=zeros(LK,10); 
ee3=zeros(LK,10); ff3=zeros(LK,5);

% 出力ファイル名 with alpha and SNR in dB
%evfile1 = strcat(folder,fn1,num2str(snr(isnr),'%02d'),'dB_1000itr.csv');

E_ZF = zeros(SIMU, Nr, Nu);             % ZF-CIの固有値のSINR
E_BD = zeros(SIMU, Nr, Nu);             % BDの固有値のSINR
E_MMSE = zeros(SIMU, Nr, Nu);           % MMSE-CIの固有値のSINR
E_BMSN_BF = zeros(SIMU, Nr, Nu);        % BMSN-BFの固有値のSINR
E_BMSN_GE = zeros(SIMU, Nr, Nu);        % BMSN-GEの固有値のSINR

Eigs_ZF = zeros(SIMU, Nr, Nu);             % ZF-CIの固有値
Eigs_BD = zeros(SIMU, Nr, Nu);             % BDの固有値
Eigs_MMSE = zeros(SIMU, Nr, Nu);           % MMSE-CIの固有値
Eigs_BMSN_BF = zeros(SIMU, Nr, Nu);        % BMSN-BFの固有値
Eigs_BMSN_GE = zeros(SIMU, Nr, Nu);        % BMSN-GEの固有値

for ik = 1:LK
    K_tar = K_box(ik);
    
    for isimu = 1:SIMU
    
        H_iid = (randn(Nu*Nr,Nt)+1j*randn(Nu*Nr,Nt))/sqrt(2); % 伝搬チャネル行列のマルチパス波成分(NLoSチャネル)
        H_los = zeros(Nu*Nr,Nt);                              % 伝搬チャネル行列の直接波成分(LoSチャネル)
    
        Theta_t = (rand(1,Nu)-0.5)*360;   % ユーザ毎の送信角 (-180deg - 180deg)
        Theta_r = (rand(1,Nu)-0.5)*360;   % ユーザ毎の受信角 (-180deg - 180deg)  % a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t*derad));  % 送信モードベクトル
    
    
        if Directivity_switch == 1
            for n = 1 : Nu
                a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad))*An*cos(Theta_t(1,n)*derad); % ユーザ毎の送信モードベクトル
                a_r = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad))*An*cos(Theta_r(1,n)*derad); % ユーザ毎の受信モードベクトル
                H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t';                 % ユーザ毎のLOSチャネル行列
            end
        else
            for n = 1 : Nu
            a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad)); % ユーザ毎の送信モードベクトル
            a_r = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad)); % ユーザ毎の受信モードベクトル
            H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t';                 % ユーザ毎のLOSチャネル行列
            end
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
        
    for nuser=1:Nu         % ユーザ毎の固有値分布
        if Nr==1
            E_ZF(isimu,:,nuser) = 10*log10(S_ZF(1,1,nuser).^2/(Nt*sigma2));
            E_BD(isimu,:,nuser) = 10*log10(S_BD(1,1,nuser).^2/(Nt*sigma2));
            E_MMSE(isimu,:,nuser) = 10*log10((S_MMSE(1,1,nuser).^2)./(Nt*sigma2));
            E_BMSN_BF(isimu,:,nuser) = 10*log10((S_BMSN_BF(1,1,nuser).^2)./(Nt*sigma2));
            E_BMSN_GE(isimu,:,nuser) = 10*log10((S_BMSN_GE(1,1,nuser).^2)./(Nt*sigma2));
        else 
            E_ZF(isimu,:,nuser) = 10*log10(diag(S_ZF(:,:,nuser)).^2/(Nt*sigma2));   
            E_BD(isimu,:,nuser) = 10*log10(diag(S_BD(:,:,nuser)).^2/(Nt*sigma2));
            E_MMSE(isimu,:,nuser) = 10*log10((diag(S_MMSE(:,:,nuser)).^2)./(Nt*sigma2));
            E_BMSN_BF(isimu,:,nuser) = 10*log10((diag(S_BMSN_BF(:,:,nuser)).^2)./(Nt*sigma2));
            E_BMSN_GE(isimu,:,nuser) = 10*log10((diag(S_BMSN_GE(:,:,nuser)).^2)./(Nt*sigma2));
        end % Nr end
    end % nuser end
    
  
    for nuser=1:Nu          %　ここから固有値そのものの値を解析
        if Nr==1
            Eigs_ZF(isimu,:,nuser) = S_ZF(1,1,nuser);
            Eigs_BD(isimu,:,nuser) = S_BD(1,1,nuser);
            Eigs_MMSE(isimu,:,nuser) = S_MMSE(1,1,nuser);
            Eigs_BMSN_BF(isimu,:,nuser) = S_BMSN_BF(1,1,nuser);
            Eigs_BMSN_GE(isimu,:,nuser) = S_BMSN_GE(1,1,nuser);
        else 
            Eigs_ZF(isimu,:,nuser) = diag(S_ZF(:,:,nuser));
            Eigs_BD(isimu,:,nuser) = diag(S_BD(:,:,nuser));
            Eigs_MMSE(isimu,:,nuser) = diag(S_MMSE(:,:,nuser));
            Eigs_BMSN_BF(isimu,:,nuser) = diag(S_BMSN_BF(:,:,nuser));
            Eigs_BMSN_GE(isimu,:,nuser) = diag(S_BMSN_GE(:,:,nuser));
        end % Nr end
    end % nuser end    
    end % isimu end
        
    for nuser=1:Nu    % ソーティング（昇順）
        E_ZF(:,:,nuser) = sort(E_ZF(:,:,nuser),1);
        E_BD(:,:,nuser) = sort(E_BD(:,:,nuser),1);
        E_MMSE(:,:,nuser) = sort(E_MMSE(:,:,nuser),1);
        E_BMSN_BF(:,:,nuser) = sort(E_BMSN_BF(:,:,nuser),1);
        E_BMSN_GE(:,:,nuser) = sort(E_BMSN_GE(:,:,nuser),1);% 固有値のSINR
        
        Eigs_ZF(:,:,nuser) = sort( Eigs_ZF(:,:,nuser),1);
        Eigs_BD(:,:,nuser) = sort( Eigs_BD(:,:,nuser),1);
        Eigs_MMSE(:,:,nuser) = sort( Eigs_MMSE(:,:,nuser),1);
        Eigs_BMSN_BF(:,:,nuser) = sort( Eigs_BMSN_BF(:,:,nuser),1);
        Eigs_BMSN_GE(:,:,nuser) = sort( Eigs_BMSN_GE(:,:,nuser),1);% 固有値そのもの
    end
   %%  % ユーザ平均
    
    E_ZFm = mean(E_ZF,3);
    E_BDm = mean(E_BD,3);
    E_MMSEm = mean(E_MMSE,3);
    E_BMSN_BFm = mean(E_BMSN_BF,3);
    E_BMSN_GEm = mean(E_BMSN_GE,3);% 固有値のSINR
    
    % 全試行回数の平均
    E_ZF_all = mean(E_ZFm,1);
    E_BD_all = mean(E_BDm,1);
    E_MMSE_all = mean(E_MMSEm,1);
    E_BMSN_BF_all = mean(E_BMSN_BFm,1);
    E_BMSN_GE_all = mean(E_BMSN_GEm,1);
    
    % 固有値SINRのEigs1とEigs2の平均
    E_ZF_whole = mean(E_ZF_all);
    E_BD_whole = mean(E_BD_all);
    E_MMSE_whole = mean(E_MMSE_all);
    E_BMSN_BF_whole = mean(E_BMSN_BF_all);
    E_BMSN_GE_whole = mean(E_BMSN_GE_all);
   %% ユーザ平均
    Eigs_ZFm = mean(Eigs_ZF,3);
    Eigs_BDm = mean(Eigs_BD,3);
    Eigs_MMSEm = mean(Eigs_MMSE,3);
    Eigs_BMSN_BFm = mean(Eigs_BMSN_BF,3);
    Eigs_BMSN_GEm = mean(Eigs_BMSN_GE,3);% 固有値そのもの
    
    % 全試行回数の平均
    Eigs_ZF_all = mean(Eigs_ZFm,1);
    Eigs_BD_all = mean(Eigs_BDm,1);
    Eigs_MMSE_all = mean(Eigs_MMSEm,1);
    Eigs_BMSN_BF_all = mean(Eigs_BMSN_BFm,1);
    Eigs_BMSN_GE_all = mean(Eigs_BMSN_GEm,1);
    
    % 固有値そのもののEigs1とEigs2の平均
    Eigs_ZF_whole = mean(Eigs_ZF_all);
    Eigs_BD_whole = mean(Eigs_BD_all);
    Eigs_MMSE_whole = mean(Eigs_MMSE_all);
    Eigs_BMSN_BF_whole = mean(Eigs_BMSN_BF_all);
    Eigs_BMSN_GE_whole = mean(Eigs_BMSN_GE_all);
   %%   ターゲットCDF値の固有値を抽出
    
    for nn=1:Nr
        aa3(ik,nn) = Eigs_BMSN_BF_all(:,nn);     % BMSN_BFの2つの固有値の試行回数平均が入った箱
        aa3(ik,nn+Nr) = Eigs_BMSN_GE_all(:,nn);  % BMSN_GEの2つの固有値の試行回数平均が入った箱
        aa3(ik,nn+2*Nr) = Eigs_MMSE_all(:,nn);   % MMSEの2つの固有値の試行回数平均が入った箱
        aa3(ik,nn+3*Nr) = Eigs_BD_all(:,nn);     % BDの2つの固有値の試行回数平均が入った箱
        aa3(ik,nn+4*Nr) = Eigs_ZF_all(:,nn);     % ZFの2つの固有値の試行回数平均が入った箱
                
        bb3(ik,nn) = E_BMSN_BF_all(:,nn);     % BMSN_BFの2つの固有値SINRの試行回数平均が入った箱
        bb3(ik,nn+Nr) = E_BMSN_GE_all(:,nn);  % BMSN_GEの2つの固有値SINRの試行回数平均が入った箱
        bb3(ik,nn+2*Nr) = E_MMSE_all(:,nn);   % MMSEの2つの固有値SINRの試行回数平均が入った箱
        bb3(ik,nn+3*Nr) = E_BD_all(:,nn);     % BDの2つの固有値SINRの試行回数平均が入った箱
        bb3(ik,nn+4*Nr) = E_ZF_all(:,nn);     % ZFの2つの固有値SINRの試行回数平均が入った箱
        
        cc3(ik,1) = E_BMSN_BF_whole;     % BMSN_BFの2つの固有値SINRの平均
        cc3(ik,2) = E_BMSN_GE_whole;  % BMSN_GEの2つの固有値SINRの平均
        cc3(ik,3) = E_MMSE_whole;   % MMSEの2つの固有値SINRの平均
        cc3(ik,4) = E_BD_whole;     % BDの2つの固有値SINRの平均
        cc3(ik,5) = E_ZF_whole;     % ZFの2つの固有値SINRの平均
               
        dd3(ik,nn) = Eigs_BMSN_BFm(round(CDF*SIMU/100),nn);     % BMSN_BFの2つの固有値が入った箱
        dd3(ik,nn+Nr) = Eigs_BMSN_GEm(round(CDF*SIMU/100),nn);  % BMSN_GEの2つの固有値が入った箱
        dd3(ik,nn+2*Nr) = Eigs_MMSEm(round(CDF*SIMU/100),nn);   % MMSEの2つの固有値が入った箱
        dd3(ik,nn+3*Nr) = Eigs_BDm(round(CDF*SIMU/100),nn);     % BDの2つの固有値が入った箱
        dd3(ik,nn+4*Nr) = Eigs_ZFm(round(CDF*SIMU/100),nn);     % ZFの2つの固有値が入った箱
        
        ee3(ik,nn) = E_BMSN_BFm(round(CDF*SIMU/100),nn);     % BMSN_BFの2つの固有値のSINRが入った箱
        ee3(ik,nn+Nr) = E_BMSN_GEm(round(CDF*SIMU/100),nn);  % BMSN_GEの2つの固有値のSINRが入った箱
        ee3(ik,nn+2*Nr) = E_MMSEm(round(CDF*SIMU/100),nn);   % MMSEの2つの固有値のSINRが入った箱
        ee3(ik,nn+3*Nr) = E_BDm(round(CDF*SIMU/100),nn);     % BDの2つの固有値のSINRが入った箱
        ee3(ik,nn+4*Nr) = E_ZFm(round(CDF*SIMU/100),nn);     % ZFの2つの固有値のSINRが入った箱
        
        ff3(ik,1) = Eigs_BMSN_BF_whole;     % BMSN_BFの2つの固有値の平均
        ff3(ik,2) = Eigs_BMSN_GE_whole;  % BMSN_GEの2つの固有値の平均
        ff3(ik,3) = Eigs_MMSE_whole;   % MMSEの2つの固有値の平均
        ff3(ik,4) = Eigs_BD_whole;     % BDの2つの固有値の平均
        ff3(ik,5) = Eigs_ZF_whole;     % ZFの2つの固有値の平均
    end
    fprintf('K = %d dB\n',K_box(ik));  % 解析中の計算過程を見るためのK = ??dBの表示
 end% ik end

%% CDFに注目した場合のグラフ
if wantCDF==1
%% Kに対する固有値[dB]グラフ（CDFに注目）
   figure;
   mycol = [0 0 1;0 0 1;1 0 0;1 0 0;0 0.7 0;0 0.7 0;1 0 1;1 0 1;0 0 0;0 0 0]; % グラフの色
   set(groot,'defaultAxesColorOrder',mycol) % figureの書式を設定
   axis([K_min K_max SNR_tar-50 SNR_tar+10]);              % figureの軸のとり方を設定
   grid on;
   hold on;
   set(gca,'XTick',K_min:5:K_max,'Fontsize',14,'Fontname','Times New Roman') % figureの書式を設定
   if Nr==2
      plot(K_box,ee3(:,1),'r-d','LineWidth',4);
      plot(K_box,ee3(:,2),'r--d','LineWidth',4);
      plot(K_box,ee3(:,3),'b-s','LineWidth',4);
      plot(K_box,ee3(:,4),'b--s','LineWidth',4);
      plot(K_box,ee3(:,5),'g-o','LineWidth',4);
      plot(K_box,ee3(:,6),'g--o','LineWidth',4);
      plot(K_box,ee3(:,7),'c-x','LineWidth',4);
      plot(K_box,ee3(:,8),'c--x','LineWidth',4);
      plot(K_box,ee3(:,9),'m-^','LineWidth',4);
      plot(K_box,ee3(:,10),'m--^','LineWidth',4);
   end
   if Nr==1
      plot(K_box,ee3(:,1),'r-d','LineWidth',4);
      plot(K_box,ee3(:,2),'b-s','LineWidth',4);
      plot(K_box,ee3(:,3),'g-o','LineWidth',4);
      plot(K_box,ee3(:,4),'c-x','LineWidth',4);
      plot(K_box,ee3(:,5),'m-^','LineWidth',4);
   end
   if Nr==1
   legend('BMSN-BF \lambda_1','BMSN-GE \lambda_1','MMSE-CI \lambda_1','BD \lambda_1','ZF-CI \lambda_1','Location','southeast');
   end
   if Nr==2
   legend('BMSN-BF \lambda_1','BMSN-BF \lambda_2','BMSN-GE \lambda_1','BMSN-GE \lambda_2','MMSE-CI \lambda_1','MMSE-CI \lambda_2',...
    'BD \lambda_1','BD \lambda_2','ZF-CI \lambda_1','ZF-CI \lambda_2','Location','southeast');
   end
   xlabel('{\it{K}} [dB]','Fontsize',16,'Fontname','Times New Roman');
   ylabel('Eigenvalues [dB]','Fontsize',16,'Fontname','Times New Roman');
   set(gca,'Fontsize',16,'Fontname','Times New Roman');
   title(strcat(target_CDF,", ",target_SNR));
   grid on;
   hold on;
%% Kに対する固有値そのもののグラフ（CDFに注目）
   figure;
   mycol = [0 0 1;0 0 1;1 0 0;1 0 0;0 0.7 0;0 0.7 0;1 0 1;1 0 1;0 0 0;0 0 0]; % グラフの色
   set(groot,'defaultAxesColorOrder',mycol) % figureの書式を設定
   axis([K_min K_max 0 7]);              % figureの軸のとり方を設定
   grid on;
   hold on;
   set(gca,'XTick',K_min:5:K_max,'Fontsize',14,'Fontname','Times New Roman') % figureの書式を設定
   if Nr==2
      plot(K_box,dd3(:,1),'r-d','LineWidth',4);
      plot(K_box,dd3(:,2),'r--d','LineWidth',4);
      plot(K_box,dd3(:,3),'b-s','LineWidth',4);
      plot(K_box,dd3(:,4),'b--s','LineWidth',4);
      plot(K_box,dd3(:,5),'g-o','LineWidth',4);
      plot(K_box,dd3(:,6),'g--o','LineWidth',4);
      plot(K_box,dd3(:,7),'c-x','LineWidth',4);
      plot(K_box,dd3(:,8),'c--x','LineWidth',4);
      plot(K_box,dd3(:,9),'m-^','LineWidth',4);
      plot(K_box,dd3(:,10),'m--^','LineWidth',4);
      legend('BMSN-BF \lambda_1','BMSN-BF \lambda_2','BMSN-GE \lambda_1','BMSN-GE \lambda_2','MMSE-CI \lambda_1','MMSE-CI \lambda_2',...
       'BD \lambda_1','BD \lambda_2','ZF-CI \lambda_1','ZF-CI \lambda_2','Location','southeast');
   end
   if Nr==1
      plot(K_box,dd3(:,1),'r-d','LineWidth',4);
      plot(K_box,dd3(:,2),'b-s','LineWidth',4);
      plot(K_box,dd3(:,3),'g-o','LineWidth',4);
      plot(K_box,dd3(:,4),'c-x','LineWidth',4);
      plot(K_box,dd3(:,5),'m-^','LineWidth',4);
   end
   if Nr==1
   legend('BMSN-BF \lambda_1','BMSN-GE \lambda_1','MMSE-CI \lambda_1','BD \lambda_1','ZF-CI \lambda_1','Location','southeast');
   end
   xlabel('{\it{K}} [dB]','Fontsize',16,'Fontname','Times New Roman');
   ylabel('Eigenvalues','Fontsize',16,'Fontname','Times New Roman');
   set(gca,'Fontsize',16,'Fontname','Times New Roman');
   title(strcat(target_CDF,", ",target_SNR));
   grid on;
   hold on;
else
%% Kに対するEigs1とEigs2の平均のグラフ（全試行回数平均）
if Nr==2
   figure;
   mycol = [0 0 1;1 0 0;0 0.7 0;1 0 1;0 0 0]; % グラフの色
   set(groot,'defaultAxesColorOrder',mycol) % figureの書式を設定
   axis([K_min K_max 0 6]);              % figureの軸のとり方を設定
   grid on;
   hold on;
   set(gca,'XTick',K_min:5:K_max,'Fontsize',14,'Fontname','Times New Roman') % figureの書式を設定
   plot(K_box,ff3(:,1),'r-d','LineWidth',4);
   plot(K_box,ff3(:,2),'b-s','LineWidth',4);
   plot(K_box,ff3(:,3),'g-o','LineWidth',4);
   plot(K_box,ff3(:,4),'c-x','LineWidth',4);
   plot(K_box,ff3(:,5),'m-^','LineWidth',4);

   legend('BMSN-BF \lambda_{Av}','BMSN-GE \lambda_{Av}','MMSE-CI \lambda_{Av}','BD \lambda_{Av}','ZF-CI \lambda_{Av}','Location','southeast');
   xlabel('{\it{K}} [dB]','Fontsize',36,'Fontname','Times New Roman');
   ylabel('Eigenvalues','Fontsize',36,'Fontname','Times New Roman');
   set(gca,'Fontsize',16,'Fontname','Times New Roman');
   title(target_SNR);
   grid on;
   hold on;
end
%% Kに対するEigs1とEigs2それぞれのグラフ（全試行回数平均）
figure;
mycol = [0 0 1;1 0 0;0 0.7 0;1 0 1;0 0 0]; % グラフの色
set(groot,'defaultAxesColorOrder',mycol) % figureの書式を設定
axis([K_min K_max 0 6]);              % figureの軸のとり方を設定
grid on;
hold on;
set(gca,'XTick',K_min:5:K_max,'Fontsize',14,'Fontname','Times New Roman') % figureの書式を設定
if Nr==1
   plot(K_box,aa3(:,1),'r-d','LineWidth',4);
   plot(K_box,aa3(:,2),'b-s','LineWidth',4);
   plot(K_box,aa3(:,3),'g-o','LineWidth',4);
   plot(K_box,aa3(:,4),'c-x','LineWidth',4);
   plot(K_box,aa3(:,5),'m-^','LineWidth',4);
end
if Nr==2
   plot(K_box,aa3(:,1),'r-d','LineWidth',4);
   plot(K_box,aa3(:,2),'r--d','LineWidth',4);
   plot(K_box,aa3(:,3),'b-s','LineWidth',4);
   plot(K_box,aa3(:,4),'b--s','LineWidth',4);
   plot(K_box,aa3(:,5),'g-o','LineWidth',4);
   plot(K_box,aa3(:,6),'g--o','LineWidth',4);
   plot(K_box,aa3(:,7),'c-x','LineWidth',4);
   plot(K_box,aa3(:,8),'c--x','LineWidth',4);
   plot(K_box,aa3(:,9),'m-^','LineWidth',4);
   plot(K_box,aa3(:,10),'m--^','LineWidth',4);
end

if Nr==1
   legend('BMSN-BF \lambda_1','BMSN-GE \lambda_1','MMSE-CI \lambda_1','BD \lambda_1','ZF-CI \lambda_1','Location','southeast');
end
if Nr>=2
   legend('BMSN-BF \lambda_1','BMSN-BF \lambda_2','BMSN-GE \lambda_1','BMSN-GE \lambda_2','MMSE-CI \lambda_1','MMSE-CI \lambda_2',...
   'BD \lambda_1','BD \lambda_2','ZF-CI \lambda_1','ZF-CI \lambda_2','Location','southeast');
end
xlabel('{\it{K}} [dB]','Fontsize',36,'Fontname','Times New Roman');
ylabel('Eigenvalues','Fontsize',36,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
title(target_SNR);
grid on;
hold on;
%% Kに対するEigs1[dB]とEigs2[dB]それぞれのグラフ（全試行回数平均）
figure;
mycol = [0 0 1;0 0 1;1 0 0;1 0 0;0 0.7 0;0 0.7 0;1 0 1;1 0 1;0 0 0;0 0 0]; % グラフの色
set(groot,'defaultAxesColorOrder',mycol) % figureの書式を設定
axis([K_min K_max SNR_tar-50 SNR_tar+10]);              % figureの軸のとり方を設定
grid on;
hold on;
set(gca,'XTick',K_min:5:K_max,'Fontsize',14,'Fontname','Times New Roman') % figureの書式を設定
if Nr==1
   plot(K_box,bb3(:,1),'r-d','LineWidth',4);
   plot(K_box,bb3(:,2),'b-s','LineWidth',4);
   plot(K_box,bb3(:,3),'g-o','LineWidth',4);
   plot(K_box,bb3(:,4),'c-x','LineWidth',4);
   plot(K_box,bb3(:,5),'m-^','LineWidth',4);
end
if Nr==2
   plot(K_box,bb3(:,1),'r-d','LineWidth',4);
   plot(K_box,bb3(:,2),'r--d','LineWidth',4);
   plot(K_box,bb3(:,3),'b-s','LineWidth',4);
   plot(K_box,bb3(:,4),'b--s','LineWidth',4);
   plot(K_box,bb3(:,5),'g-o','LineWidth',4);
   plot(K_box,bb3(:,6),'g--o','LineWidth',4);
   plot(K_box,bb3(:,7),'c-x','LineWidth',4);
   plot(K_box,bb3(:,8),'c--x','LineWidth',4);
   plot(K_box,bb3(:,9),'m-^','LineWidth',4);
   plot(K_box,bb3(:,10),'m--^','LineWidth',4);
end

if Nr==1
   legend('BMSN-BF \lambda_1','BMSN-GE \lambda_1','MMSE-CI \lambda_1','BD \lambda_1','ZF-CI \lambda_1','Location','southeast');
end
if Nr>=2
   legend('BMSN-BF \lambda_1','BMSN-BF \lambda_2','BMSN-GE \lambda_1','BMSN-GE \lambda_2','MMSE-CI \lambda_1','MMSE-CI \lambda_2',...
   'BD \lambda_1','BD \lambda_2','ZF-CI \lambda_1','ZF-CI \lambda_2','Location','southeast');
end
xlabel('{\it{K}} [dB]','Fontsize',16,'Fontname','Times New Roman');
ylabel('Eigenvalues [dB]','Fontsize',16,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
title(target_SNR);
grid on;
hold on;
%% Kに対するEigs1[dB]とEigs2[dB]との平均のグラフ（全試行回数平均）
figure;
mycol = [0 0 1;0 0 1;1 0 0;1 0 0;0 0.7 0;0 0.7 0;1 0 1;1 0 1;0 0 0;0 0 0]; % グラフの色
set(groot,'defaultAxesColorOrder',mycol) % figureの書式を設定
axis([K_min K_max SNR_tar-50 SNR_tar+10]);              % figureの軸のとり方を設定
grid on;
hold on;
set(gca,'XTick',K_min:5:K_max,'Fontsize',14,'Fontname','Times New Roman') % figureの書式を設定
plot(K_box,cc3(:,1),'r-d','LineWidth',4);
plot(K_box,cc3(:,2),'b-s','LineWidth',4);
plot(K_box,cc3(:,3),'g-o','LineWidth',4);
plot(K_box,cc3(:,4),'c-x','LineWidth',4);
plot(K_box,cc3(:,5),'m-^','LineWidth',4);

legend('BMSN-BF \lambda_{Av}','BMSN-GE \lambda_{Av}','MMSE-CI \lambda_{Av}','BD \lambda_{Av}','ZF-CI \lambda_{Av}','Location','southeast');
xlabel('{\it{K}} [dB]','Fontsize',36,'Fontname','Times New Roman');
ylabel('Eigenvalues[dB]','Fontsize',36,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
title(target_SNR);
grid on;
hold on;
end
%End