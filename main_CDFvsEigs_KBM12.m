% cap_bd_muser.m
% NU-ユーザでのBD法の計算
% TDMA, Upper boundと比較

clear;%close;

% 出力ファイル名
% testfile1 = 'Capacity2x8x4u_BD.csv';
% testfile2 = 'Capacity2x8x4u_BD_CDFSNR20dB.csv';
% testfile3 = 'Eig2x8x4u_BD_CDFSNR20dB.csv';

% パラメータ条件 NT >= NR*NU
SN_tar  = 10;        % CDF表示のためのターゲットSNR [dB]
SIMU   = 1000;       % 伝搬チャネル行列の発生回数
Nt     = 16;         % 送信素子数
Nr     = 2;          % 受信素子数(=2に固定)
Nu     = 8;          % ユーザ数
I      = eye(Nt,Nt); % NTxNTの単位行列

K_dB   = -20;
K = 10^(K_dB/10);    % Kの真値
An = 2; % 指向性関数の係数(cosパターンの場合は2)
 
d_t      = 0.5;      % 送信アンテナ間隔（in wavelength)
d_r      = 0.5;      % 受信アンテナ間隔（in wavelength)
derad = pi/180;      % degree -> rad

% 擬似雑音
a = Nt/(10^(SN_tar/10));
%a2 = NR*NT/(10^(SN_tar/10));    % 擬似雑音 for BMSN2
        
T = zeros(Nr,Nr); % 所望のチャネル行列 for BMSN
for nuser = 1:Nu
    T(:,:,nuser) = eye(Nr,Nr);
end

if SN_tar < 10
    target_snr=strcat('SNR=',num2str(SN_tar,'%01d'),'dB');
else
    target_snr=strcat('SNR=',num2str(SN_tar,'%02d'),'dB');
end

if K_dB < 10
   target_K=strcat('K=',num2str(K_dB,'%01d'),'dB');
else
   target_K=strcat('K=',num2str(K_dB,'%02d'),'dB');
end
%伝搬チャネル行列の直接波成分(LOS チャネル)
H_los = zeros(Nu*Nr,Nt);

E_ZF = zeros(SIMU, Nr, Nu);             % ZF-CIの固有値
E_BD = zeros(SIMU, Nr, Nu);             % BDの固有値
E_MMSE = zeros(SIMU, Nr, Nu);           % MMSE-CIの固有値
E_GMI1 = zeros(SIMU, Nr, Nu);           % GMI1の固有値
E_GMI2 = zeros(SIMU, Nr, Nu);           % GMI2の固有値
E_BMSN_BF = zeros(SIMU, Nr, Nu);        % BMSN-BFの固有値
E_BMSN_GE = zeros(SIMU, Nr, Nu);        % BMSN-GEの固有値
E_BMSN_BF1 = zeros(SIMU, Nr, Nu);        % BMSN-BF1の固有値
E_BMSN_GE1 = zeros(SIMU, Nr, Nu);        % BMSN-GE1の固有値
E_BMSN_BF2 = zeros(SIMU, Nr, Nu);        % BMSN-BF2の固有値
E_BMSN_GE2 = zeros(SIMU, Nr, Nu);        % BMSN-GE2の固有値
E_BMSN_GE3 = zeros(SIMU, Nr, Nu);        % BMSN-GE3の固有値

for Directivity_switch = 0:1
for isimu = 1:SIMU              % 試行回数のループ
      
    % LOS チャネル      
        if Directivity_switch == 1 % 送受信素子の指向性考慮有り
            Theta_t = (rand(1,Nu)-0.5)*180; % ユーザ毎の送信角 指向性:(-90deg - 90deg)
            Theta_r = (rand(1,Nu)-0.5)*180; % ユーザ毎の受信角 指向性:(-90deg - 90deg)
            for n = 1 : Nu
                a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad))*An*cos(Theta_t(1,n)*derad); % ユーザ毎の送信モードベクトル
                a_r = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad))*An*cos(Theta_r(1,n)*derad); % ユーザ毎の受信モードベクトル
                H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t'; % ユーザ毎のLOSチャネル行列
%% LoSチャネル　候補1
%                 a_t_iid = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad)); % ユーザ毎のNLoS送信モードベクトル
%                 a_r_iid = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad)); % ユーザ毎のNLoS受信モードベクトル
%                 Theta_t_iid = (rand(1,Nt)-0.5)*360; % ユーザ毎の送信角 指向性:(-90deg - 90deg)
%                 Theta_r_iid = (rand(1,Nr)-0.5)*360; % ユーザ毎の受信角 指向性:(-90deg - 90deg)
%                 g_theta_t_iid = An*cos(Theta_t_iid*derad);
%                 g_theta_r_iid = An*cos(Theta_r_iid*derad);
%                 a_t_iid = a_t_iid.* g_theta_t_iid';
%                 a_r_iid = a_r_iid.* g_theta_r_iid';
%                 H_iid((n-1)*Nr+1:(n-1)*Nr+Nr,:) = (a_r_iid*a_t_iid'); % ユーザ毎のNLOSチャネル行列
%% LoSチャネル　候補2
%     a_t_iid = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin((rand(1,Nu)-0.5)*360*derad)); % ユーザ毎のNLoS送信モードベクトル
%     a_r_iid = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin((rand(1,Nu)-0.5)*360*derad)); % ユーザ毎のNLoS受信モードベクトル
%     Theta_t_iid = (rand(1,Nt)-0.5)*360; % ユーザ毎の送信角 指向性:(-180deg ～ 180deg)
%     Theta_r_iid = (rand(1,Nr)-0.5)*360; % ユーザ毎の受信角 指向性:(-180deg ～ 180deg)
%     g_theta_t_iid = An*cos(Theta_t_iid*derad);
%     g_theta_r_iid = An*cos(Theta_r_iid*derad);
%     a_t_iid = a_t_iid.* g_theta_t_iid';
%     a_r_iid = a_r_iid.* g_theta_r_iid';
%     H_iid((n-1)*Nr+1:(n-1)*Nr+Nr,:) = (a_r_iid*a_t_iid')/2; % ユーザ毎のNLOSチャネル行列
            end
            Theta_iid = (rand(Nr*Nu,Nt)-0.5)*180; % ユーザ毎の送信角 指向性:(-90deg - 90deg)
            g_theta_iid = An*cos(Theta_iid*derad);
            H_iid = (randn(Nr*Nu,Nt)+1j*randn(Nr*Nu,Nt))/sqrt(2); % 伝搬チャネル行列のマルチパス成分 (i.i.d. Rayleigh , NLOS チャネル)
            H_iid = H_iid * g_theta_iid;
            % Normcos = norm(H_los/sqrt(2),'fro');
        elseif Directivity_switch == 0 % 送受信素子の指向性考慮無し
            H_iid = (randn(Nr*Nu,Nt)+1j*randn(Nr*Nu,Nt))/sqrt(2); % 伝搬チャネル行列のマルチパス成分 (i.i.d. Rayleigh , NLOS チャネル)
            Theta_t = (rand(1,Nu)-0.5)*360; % ユーザ毎の送信角 等方性:(-180deg - 180deg) 
            Theta_r = (rand(1,Nu)-0.5)*360; % ユーザ毎の受信角 等方性:(-180deg - 180deg) 
            for n = 1 : Nu
                a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad)); % ユーザ毎の送信モードベクトル
                a_r = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad)); % ユーザ毎の受信モードベクトル
                H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t';                 % ユーザ毎のLOSチャネル行列
            end
            % Normiso = norm(H_iid,'fro');
        end
    

% 伝搬チャネル行列=[sqrt(K/(K+1))*(LOS チャネル)]...
%                   .+[sqrt(1/(K+1))*(NLOS チャネル)]

    H0 = sqrt(K/(K+1))*H_los + sqrt(1/(K+1))*H_iid;
    
    %H0 = squeeze(H(k,:,:)); % k番目の試行回数での伝搬チャネル行列
    %%
%     % ZF-CI algorithm
%     [W_ZF,U_ZF,S_ZF,~,~] = zf(Nt,Nr,Nu,H0); % function zf.m を使用    
%     % BD algorithm
%     [W_BD,U_BD,S_BD,~,~] = bd(Nt,Nr,Nu,H0); % function bd.m を使用    
%     % MMSE-CI algorithm
%     [W_MMSE,U_MMSE,S_MMSE,RIPM,~] = mmse(Nt,Nr,Nu,H0,a); % function mmse.m を使用    
%     % GMI1 algorithm
%     [W_GMI1,U_GMI1,S_GMI1,RIPGM1,~] = gmmse_m1(Nt,Nr,Nu,H0,a); % function gmmse.m1 を使用    
%     % GMI2 algorithm
%     [W_GMI2,U_GMI2,S_GMI2,RIPGM2,~] = gmmse_m2(Nt,Nr,Nu,H0,a); % function gmmse.m2 を使用
%     
%     % BMSN-BF algorithm
%     [W_BMSN_BF,U_BMSN_BF,S_BMSN_BF,RIPBF,~] = bmsn_bf(Nt,Nr,Nu,H0,a,T); % function bmsn_bf.m を使用      
%     % BMSN-GE algorithm
%     [W_BMSN_GE,U_BMSN_GE,S_BMSN_GE,RIPGE,~] = bmsn_gev(Nt,Nr,Nu,H0,a); % function bmsn_gev.m を使用 
    %%
    % BMSN-BF1 algorithm
    [W_BMSN_BF1,U_BMSN_BF1,S_BMSN_BF1,RIPBF1,~] = bmsn_bf1(Nt,Nr,Nu,H0,a,T); % function bmsn_bf1.m を使用      
    % BMSN-GE1 algorithm
    [W_BMSN_GE1,U_BMSN_GE1,S_BMSN_GE1,RIPGE1,~] = bmsn_ge1(Nt,Nr,Nu,H0,a); % function bmsn_ge1.m を使用    
    % BMSN-BF2 algorithm
    [W_BMSN_BF2,U_BMSN_BF2,S_BMSN_BF2,RIPBF2,~] = bmsn_bf2(Nt,Nr,Nu,H0,a); % function bmsn_bf2.m を使用      
    % BMSN-GE2 algorithm
    [W_BMSN_GE2,U_BMSN_GE2,S_BMSN_GE2,RIPGE2,~] = bmsn_ge2(Nt,Nr,Nu,H0,a); % function bmsn_ge2.m を使用
    % BMSN-GE3 algorithm
    [W_BMSN_GE3,U_BMSN_GE3,S_BMSN_GE3,RIPGE3,~] = bmsn_ge3(Nt,Nr,Nu,H0,a); % function bmsn_ge3.m を使用
    
    % ユーザ毎の固有値分布
    snt = 1/(10^(SN_tar/10));
    for nuser=1:Nu
        if Nr==1
%             E_ZF(isimu,:,nuser) = 10*log10(S_ZF(1,1,nuser).^2/(Nt*snt));
%             E_BD(isimu,:,nuser) = 10*log10(S_BD(1,1,nuser).^2/(Nt*snt));
%             E_MMSE(isimu,:,nuser) = 10*log10((S_MMSE(1,1,nuser).^2)./(RIPM(1,nuser)+Nt*snt));
%             E_GMI1(isimu,:,nuser) = 10*log10((S_GMI1(1,1,nuser).^2)./(RIPGM1(1,nuser)+Nt*snt));
%             E_GMI2(isimu,:,nuser) = 10*log10((S_GMI2(1,1,nuser).^2)./(RIPGM2(1,nuser)+Nt*snt));
%             E_BMSN_BF(isimu,:,nuser) = 10*log10((S_BMSN_BF(1,1,nuser).^2)./(RIPBF(1,nuser)+Nt*snt));
%             E_BMSN_GE(isimu,:,nuser) = 10*log10((S_BMSN_GE(1,1,nuser).^2)./(RIPGE(1,nuser)+Nt*snt));
            E_BMSN_BF1(isimu,:,nuser) = 10*log10((S_BMSN_BF1(1,1,nuser).^2)./(RIPBF1(1,nuser)+Nt*snt));
            E_BMSN_GE1(isimu,:,nuser) = 10*log10((S_BMSN_GE1(1,1,nuser).^2)./(RIPGE1(1,nuser)+Nt*snt));
            E_BMSN_BF2(isimu,:,nuser) = 10*log10((S_BMSN_BF2(1,1,nuser).^2)./(RIPBF2(1,nuser)+Nt*snt));
            E_BMSN_GE2(isimu,:,nuser) = 10*log10((S_BMSN_GE2(1,1,nuser).^2)./(RIPGE2(1,nuser)+Nt*snt));
            E_BMSN_GE3(isimu,:,nuser) = 10*log10((S_BMSN_GE3(1,1,nuser).^2)./(RIPGE3(1,nuser)+Nt*snt));
        else 
%             E_ZF(isimu,:,nuser) = 10*log10(diag(S_ZF(:,:,nuser)).^2/(Nt*snt));   
%             E_BD(isimu,:,nuser) = 10*log10(diag(S_BD(:,:,nuser)).^2/(Nt*snt));
%             E_MMSE(isimu,:,nuser) = 10*log10((diag(S_MMSE(:,:,nuser)).^2)./(RIPM(:,nuser)+Nt*snt));
%             E_GMI1(isimu,:,nuser) = 10*log10((diag(S_GMI1(:,:,nuser)).^2)./(RIPGM1(:,nuser)+Nt*snt));
%             E_GMI2(isimu,:,nuser) = 10*log10((diag(S_GMI2(:,:,nuser)).^2)./(RIPGM2(:,nuser)+Nt*snt));
%             E_BMSN_BF(isimu,:,nuser) = 10*log10((diag(S_BMSN_BF(:,:,nuser)).^2)./(RIPBF(:,nuser)+Nt*snt));
%             E_BMSN_GE(isimu,:,nuser) = 10*log10((diag(S_BMSN_GE(:,:,nuser)).^2)./(RIPGE(:,nuser)+Nt*snt));
            E_BMSN_BF1(isimu,:,nuser) = 10*log10((diag(S_BMSN_BF1(:,:,nuser)).^2)./(RIPBF1(:,nuser)+Nt*snt));
            E_BMSN_GE1(isimu,:,nuser) = 10*log10((diag(S_BMSN_GE1(:,:,nuser)).^2)./(RIPGE1(:,nuser)+Nt*snt));
            E_BMSN_BF2(isimu,:,nuser) = 10*log10((diag(S_BMSN_BF2(:,:,nuser)).^2)./(RIPBF2(:,nuser)+Nt*snt));
            E_BMSN_GE2(isimu,:,nuser) = 10*log10((diag(S_BMSN_GE2(:,:,nuser)).^2)./(RIPGE2(:,nuser)+Nt*snt));
            E_BMSN_GE3(isimu,:,nuser) = 10*log10((diag(S_BMSN_GE3(:,:,nuser)).^2)./(RIPGE3(:,nuser)+Nt*snt));
        end
    end
    
     %　ここから固有値そのものの値を解析
    for nuser=1:Nu
        if Nr==1
%             Eigs_ZF(isimu,:,nuser) = S_ZF(1,1,nuser);
%             Eigs_BD(isimu,:,nuser) = S_BD(1,1,nuser);
%             Eigs_MMSE(isimu,:,nuser) = S_MMSE(1,1,nuser);
%             Eigs_GMI1(isimu,:,nuser) = S_GMI1(1,1,nuser);
%             Eigs_GMI2(isimu,:,nuser) = S_GMI2(1,1,nuser);
%             Eigs_BMSN_BF(isimu,:,nuser) = S_BMSN_BF(1,1,nuser);
%             Eigs_BMSN_GE(isimu,:,nuser) = S_BMSN_GE(1,1,nuser);
            Eigs_BMSN_BF1(isimu,:,nuser) = S_BMSN_BF1(1,1,nuser);
            Eigs_BMSN_GE1(isimu,:,nuser) = S_BMSN_GE1(1,1,nuser);
            Eigs_BMSN_BF2(isimu,:,nuser) = S_BMSN_BF2(1,1,nuser);
            Eigs_BMSN_GE2(isimu,:,nuser) = S_BMSN_GE2(1,1,nuser);
            Eigs_BMSN_GE3(isimu,:,nuser) = S_BMSN_GE3(1,1,nuser);
        else 
%             Eigs_ZF(isimu,:,nuser) = diag(S_ZF(:,:,nuser));
%             Eigs_BD(isimu,:,nuser) = diag(S_BD(:,:,nuser));
%             Eigs_MMSE(isimu,:,nuser) = diag(S_MMSE(:,:,nuser));
%             Eigs_GMI1(isimu,:,nuser) = diag(S_GMI1(:,:,nuser));
%             Eigs_GMI2(isimu,:,nuser) = diag(S_GMI2(:,:,nuser));
%             Eigs_BMSN_BF(isimu,:,nuser) = diag(S_BMSN_BF(:,:,nuser));
%             Eigs_BMSN_GE(isimu,:,nuser) = diag(S_BMSN_GE(:,:,nuser));
            Eigs_BMSN_BF1(isimu,:,nuser) = diag(S_BMSN_BF1(:,:,nuser));
            Eigs_BMSN_GE1(isimu,:,nuser) = diag(S_BMSN_GE1(:,:,nuser));
            Eigs_BMSN_BF2(isimu,:,nuser) = diag(S_BMSN_BF2(:,:,nuser));
            Eigs_BMSN_GE2(isimu,:,nuser) = diag(S_BMSN_GE2(:,:,nuser));
            Eigs_BMSN_GE3(isimu,:,nuser) = diag(S_BMSN_GE3(:,:,nuser));
        end % NR end
    end % nuser end
end

% CDF of Eigenvalue at Target SNR
rr2(:,1) = (1/SIMU:1/SIMU:1).'*100; 
Y = rr2;
% rr3 = zeros(SIMU,Nr*10);
% rr3cos = zeros(SIMU,Nr*9);
for nuser=1:Nu
%     E_ZF(:,:,nuser) = sort(E_ZF(:,:,nuser),1);
%     E_BD(:,:,nuser) = sort(E_BD(:,:,nuser),1);
%     E_MMSE(:,:,nuser) = sort(E_MMSE(:,:,nuser),1);
%     E_GMI1(:,:,nuser) = sort(E_GMI1(:,:,nuser),1);
%     E_GMI2(:,:,nuser) = sort(E_GMI2(:,:,nuser),1);
%     E_BMSN_BF(:,:,nuser) = sort(E_BMSN_BF(:,:,nuser),1);
%     E_BMSN_GE(:,:,nuser) = sort(E_BMSN_GE(:,:,nuser),1);% 固有値のSINR
    E_BMSN_BF1(:,:,nuser) = sort(E_BMSN_BF1(:,:,nuser),1);
    E_BMSN_GE1(:,:,nuser) = sort(E_BMSN_GE1(:,:,nuser),1);
    E_BMSN_BF2(:,:,nuser) = sort(E_BMSN_BF2(:,:,nuser),1);
    E_BMSN_GE2(:,:,nuser) = sort(E_BMSN_GE2(:,:,nuser),1);
    E_BMSN_GE3(:,:,nuser) = sort(E_BMSN_GE3(:,:,nuser),1);
    
end

% ユーザ平均
% E_ZFm = mean(E_ZF,3);
% E_BDm = mean(E_BD,3);
% E_MMSEm = mean(E_MMSE,3);
% E_GMI1m = mean(E_GMI1,3);
% E_GMI2m = mean(E_GMI2,3);
% E_BMSN_BFm = mean(E_BMSN_BF,3);
% E_BMSN_GEm = mean(E_BMSN_GE,3);% 固有値のSINR
E_BMSN_BF1m = mean(E_BMSN_BF1,3);
E_BMSN_GE1m = mean(E_BMSN_GE1,3);
E_BMSN_BF2m = mean(E_BMSN_BF2,3);
E_BMSN_GE2m = mean(E_BMSN_GE2,3);
E_BMSN_GE3m = mean(E_BMSN_GE3,3);

% 試行回数平均
if Directivity_switch == 0
    E_BMSN_BF1meansimu = mean(E_BMSN_BF1m,1);
    E_BMSN_GE1meansimu = mean(E_BMSN_GE1m,1);
    E_BMSN_BF2meansimu = mean(E_BMSN_BF2m,1);
    E_BMSN_GE2meansimu = mean(E_BMSN_GE2m,1);
    E_BMSN_GE3meansimu = mean(E_BMSN_GE3m,1);
elseif Directivity_switch == 1
    E_BMSN_BF1cosmeansimu = mean(E_BMSN_BF1m,1);
    E_BMSN_GE1cosmeansimu = mean(E_BMSN_GE1m,1);
    E_BMSN_BF2cosmeansimu = mean(E_BMSN_BF2m,1);
    E_BMSN_GE2cosmeansimu = mean(E_BMSN_GE2m,1);
    E_BMSN_GE3cosmeansimu = mean(E_BMSN_GE3m,1);
end

   for nn=1:Nr
       if Directivity_switch == 0
           %     rr3(:,nn) = E_BMSN_BFm(:,nn);
           %     rr3(:,nn+Nr) = E_BMSN_GEm(:,nn);
           %     rr3(:,nn+2*Nr) = E_MMSEm(:,nn);
           %     rr3(:,nn+3*Nr) = E_GMI1m(:,nn);
           %     rr3(:,nn+4*Nr) = E_GMI2m(:,nn);
           rr3(:,nn+5*Nr) = E_BMSN_BF1m(:,nn);
           rr3(:,nn+6*Nr) = E_BMSN_GE1m(:,nn);
           rr3(:,nn+7*Nr) = E_BMSN_BF2m(:,nn);
           rr3(:,nn+8*Nr) = E_BMSN_GE2m(:,nn); 
           rr3(:,nn+9*Nr) = E_BMSN_GE3m(:,nn);           
       elseif Directivity_switch == 1
           rr3cos(:,nn+5*Nr) = E_BMSN_BF1m(:,nn);
           rr3cos(:,nn+6*Nr) = E_BMSN_GE1m(:,nn);
           rr3cos(:,nn+7*Nr) = E_BMSN_BF2m(:,nn);
           rr3cos(:,nn+8*Nr) = E_BMSN_GE2m(:,nn);
           rr3cos(:,nn+9*Nr) = E_BMSN_GE3m(:,nn);
       end
    end
end
    
%% CDF of EGV at Target SNR 
figure;
mycol = [1 0 1;1 0 1;0 0 0;0 0 0;1 0 0;1 0 0;0 0 1;0 0 1;1 0 0;1 0 0;0 0 0;0 0 0];
set(groot,'defaultAxesColorOrder',mycol)
if Nr==1
    Holizon_min = round(min(min(rr3cos(:,9))))-5;
else
    Holizon_min = round(min(min(rr3cos(:,18))))-5;
end
Holizon_max = round(max(max(rr3cos)))+5;
axis([Holizon_min Holizon_max 0 100]);
grid on;
hold on;
set(gca,'XTick',Holizon_min:5:Holizon_max,'Fontsize',14,'Fontname','Times New Roman')
xlabel('SINR of eigenvalue [dB]','Fontsize',16,'Fontname','Times New Roman');
ylabel('CDF [%]','Fontsize',16,'Fontname','Times New Roman');
if Nr==1
% plot(rr3(:,1),Y,'r-d','MarkerIndices',10:100:length(Y),'Linewidth',2);
% plot(rr3(:,2),Y,'b-s','MarkerIndices',50:100:length(Y),'Linewidth',2);
plot(rr3(:,6),Y,'g--','MarkerIndices',10:100:length(Y),'Linewidth',2);
plot(rr3(:,7),Y,'y--','MarkerIndices',50:100:length(Y),'Linewidth',2);
plot(rr3(:,8),Y,'k--','MarkerIndices',10:100:length(Y),'Linewidth',2);
plot(rr3(:,9),Y,'c--','MarkerIndices',50:100:length(Y),'Linewidth',2);
plot(rr3(:,10),Y,'r--','MarkerIndices',10:100:length(Y),'Linewidth',2);

% plot(rr3cos(:,1),Y,'r-d','MarkerIndices',10:100:length(Y),'Linewidth',2);
% plot(rr3cos(:,2),Y,'b-s','MarkerIndices',50:100:length(Y),'Linewidth',2);
plot(rr3cos(:,6),Y,'g-o','MarkerIndices',10:100:length(Y),'Linewidth',2);
plot(rr3cos(:,7),Y,'y-o','MarkerIndices',50:100:length(Y),'Linewidth',2);
plot(rr3cos(:,8),Y,'k-o','MarkerIndices',10:100:length(Y),'Linewidth',2);
plot(rr3cos(:,9),Y,'c-o','MarkerIndices',50:100:length(Y),'Linewidth',2);
plot(rr3cos(:,10),Y,'r-o','MarkerIndices',10:100:length(Y),'Linewidth',2);
legend('BMSN-BF1 \lambda_1','BMSN-GE1 \lambda_1','BMSN-BF2 \lambda_1','BMSN-GE2 \lambda_1','BMSN-GE3 \lambda_1',...
    'BMSN-BF1cos \lambda_1','BMSN-GE1cos \lambda_1','BMSN-BF2cos \lambda_1','BMSN-GE2cos \lambda_1','BMSN-GE3cos \lambda_1','Location','southeast');
end

if Nr==2
% plot(rr3(:,1),Y,'r-d','MarkerIndices',10:100:length(Y),'Linewidth',2);
% plot(rr3(:,2),Y,'r--d','MarkerIndices',20:100:length(Y),'Linewidth',2);
% plot(rr3(:,3),Y,'b-s','MarkerIndices',30:100:length(Y),'Linewidth',2);
% plot(rr3(:,4),Y,'b--s','MarkerIndices',40:100:length(Y),'Linewidth',2);
% plot(rr3(:,11),Y,'g-o','MarkerIndices',10:100:length(Y),'Linewidth',2);
% plot(rr3(:,12),Y,'g--o','MarkerIndices',10:100:length(Y),'Linewidth',2);
plot(rr3(:,13),Y,'b--','MarkerIndices',50:100:length(Y),'Linewidth',2);
plot(rr3(:,14),Y,'b:','MarkerIndices',50:100:length(Y),'Linewidth',2);
% plot(rr3(:,15),Y,'k-*','MarkerIndices',10:100:length(Y),'Linewidth',2);
% plot(rr3(:,16),Y,'k--*','MarkerIndices',10:100:length(Y),'Linewidth',2);
plot(rr3(:,17),Y,'k--','MarkerIndices',50:100:length(Y),'Linewidth',2);
plot(rr3(:,18),Y,'k:','MarkerIndices',50:100:length(Y),'Linewidth',2);
plot(rr3(:,19),Y,'r--','MarkerIndices',30:100:length(Y),'Linewidth',2);
% plot(rr3(:,20),Y,'r--o','MarkerIndices',30:100:length(Y),'Linewidth',2);

% plot(rr3cos(:,11),Y,'g-o','MarkerIndices',10:100:length(Y),'Linewidth',2);
% plot(rr3cos(:,12),Y,'g--o','MarkerIndices',10:100:length(Y),'Linewidth',2);
plot(rr3cos(:,13),Y,'b-o','MarkerIndices',50:100:length(Y),'Linewidth',2);
plot(rr3cos(:,14),Y,'b--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
% plot(rr3cos(:,15),Y,'k-*','MarkerIndices',10:100:length(Y),'Linewidth',2);
% plot(rr3cos(:,16),Y,'k--*','MarkerIndices',10:100:length(Y),'Linewidth',2);
plot(rr3cos(:,17),Y,'k-o','MarkerIndices',50:100:length(Y),'Linewidth',2);
plot(rr3cos(:,18),Y,'k--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
plot(rr3cos(:,19),Y,'r-o','MarkerIndices',30:100:length(Y),'Linewidth',2);

%plot(rr3(:,6),rr3(:,4),'m');
legend('BMSN-GE1 \lambda_1','BMSN-GE1 \lambda_2','BMSN-GE2 \lambda_1','BMSN-GE2 \lambda_2','BMSN-GE3 \lambda_1',...
    'BMSN-GE1cos \lambda_1','BMSN-GE1cos \lambda_2','BMSN-GE2cos \lambda_1','BMSN-GE2cos \lambda_2','BMSN-GE3cos \lambda_1','Location','southeast');
end
title(strcat(target_snr,',',target_K));
% End