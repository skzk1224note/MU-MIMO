% cap_bd_muser.m
% NU-ユーザでのBD法の計算
% TDMA, Upper boundと比較

clear;

% 出力ファイル名
% testfile1 = 'Capacity2x8x4u_BD.csv';
% testfile2 = 'Capacity2x8x4u_BD_CDFSNR20dB.csv';
% testfile3 = 'Eig2x8x4u_BD_CDFSNR20dB.csv';

% パラメータ条件 NT >= NR*NU

SN_tar = 10;         % CDF表示のためのターゲットSNR [dB]
SIMU   = 1000;       % 試行回数
Nt     = 16;         % 送信素子数
Nr     = 2;          % 受信素子数(=2に固定)
Nu     = 8;          % ユーザ数
NL     = 8;          % マルチパス波の素波数
I      = eye(Nt,Nt); % NTxNTの単位行列

K_dB   = -20;
K = 10^(K_dB/10);    % Kの真値
 
d_t = 0.5;      % 送信アンテナ間隔（in wavelength)
d_r = 0.5;      % 受信アンテナ間隔（in wavelength)
derad = pi/180;      % degree -> rad


An = 2; % 指向性関数の係数
figure_switch = 0;   % figure表示用のスイッチ 1:on, 0:off（onなら固有値そのもののグラフを表示）

a = Nt/(10^(SN_tar/10)); % 擬似雑音
        
T = zeros(Nr,Nr); % 所望のチャネル行列 for BMSN
for nuser = 1:Nu
    T(:,:,nuser) = eye(Nr,Nr);
end
Algorithms = ["BMSN-BF","BMSN-GE","MMSE","GMI1","GMI2","BD","ZF"]; %,"BMSN-GE3"
S = zeros(Nr, Nr, Nu, numel(Algorithms));
E = zeros(SIMU, Nr, Nu, numel(Algorithms));
Eigs = zeros(SIMU, Nr, Nu, numel(Algorithms));
EigsAll = zeros(SIMU, Nr, Nu, numel(Algorithms));
MSt = zeros(SIMU, Nr, Nu, numel(Algorithms));
C = zeros(SIMU, Nu, numel(Algorithms));
a_t_iid = zeros(Nt, 1, NL);
a_r_iid = zeros(Nr, 1, NL);
X = zeros(Nr, Nt, NL);
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

for Directivity_switch = 0:1 % 送受信素子の指向性考慮の有無 0:無,1:有
    for isimu = 1:SIMU % 試行回数のループ
        if Directivity_switch == 1 % 送受信素子の指向性考慮有り
            Theta_t = (rand(1,Nu)-0.5)*180; % ユーザ毎の送信角 指向性:(-90deg - 90deg)
            Theta_r = (rand(1,Nu)-0.5)*180; % ユーザ毎の受信角 指向性:(-90deg - 90deg)
            % LoS チャネル
            for n = 1 : Nu
                a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad)); % ユーザ毎の送信モードベクトル An*cos(Theta_t(1,n)*derad) *
                a_r = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad)); % ユーザ毎の受信モードベクトル An*cos(Theta_r(1,n)*derad) *
                H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t'; % ユーザ毎のLOSチャネル行列
                % 幾何学NLoSチャネル
                for L = 1:NL
                    a_t_iid(:,:,L) = sqrt(1/Nt) * exp(-1j*2*pi*d_t*(0:Nt-1).' * sin(Theta_t(1,n)*derad)); % ユーザ毎のNLoS送信モードベクトル An*cos(Theta_t(1,n)*derad) *
                    a_r_iid(:,:,L) = sqrt(1/Nr) * exp(-1j*2*pi*d_r*(0:Nr-1).' * sin(Theta_r(1,n)*derad)); % ユーザ毎のNLoS受信モードベクトル An*cos(Theta_r(1,n)*derad) *
                    X(:,:,L) = exp(1j*(rand(Nr,Nt)-0.5)*360*derad) .* (a_r_iid(:,:,L) * a_t_iid(:,:,L)');
                end
                H_iid((n-1)*Nr+1:(n-1)*Nr+Nr,:) = sum(X,3) * sqrt(Nt*Nr/NL);
%                 g_theta_t_iid = 1; % An*cos(Theta_t_iid*derad);
%                 g_theta_r_iid = 1; % An*cos(Theta_r_iid*derad);
%                 a_t_iid = a_t_iid.* g_theta_t_iid';
%                 a_r_iid = a_r_iid.* g_theta_r_iid';
                % H_iid((n-1)*Nr+1:(n-1)*Nr+Nr,:) = (a_r_iid*a_t_iid'); %*sqrt(1/NL); % ユーザ毎のNLOSチャネル行列
            end
        elseif Directivity_switch == 0 % 送受信素子の指向性考慮無し
            H_iid = (randn(Nr*Nu,Nt)+1j*randn(Nr*Nu,Nt))/sqrt(2); % 伝搬チャネル行列のマルチパス成分 (i.i.d. Rayleigh , NLOS チャネル)
            Theta_t = (rand(1,Nu)-0.5)*180; % ユーザ毎の送信角 等方性:(-180deg - 180deg) 
            Theta_r = (rand(1,Nu)-0.5)*180; % ユーザ毎の受信角 等方性:(-180deg - 180deg) 
            Theta_iid = (rand(Nr*Nu,Nt)-0.5)*180; % ユーザ毎の送信角 指向性:(-90deg - 90deg)
            g_theta_iid = 1; %An*cos(Theta_iid*derad);
            H_iid = H_iid * g_theta_iid;
            for n = 1 : Nu
                a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad)); % ユーザ毎の送信モードベクトル
                a_r = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad)); % ユーザ毎の受信モードベクトル
                H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t';                 % ユーザ毎のLOSチャネル行列
            end
            % Normiso = norm(H_iid,'fro');
        end
        
        % 伝搬チャネル行列のマルチパス成分 (i.i.d. Rayleigh, NLOS チャネル)
%         H_iid = (randn(Nr*Nu,Nt)+1j*randn(Nr*Nu,Nt))/sqrt(2);
        % 伝搬チャネル行列=[sqrt(K/(K+1))*(LOS チャネル)]+[sqrt(1/(K+1))*(NLOS チャネル)]
        H0 = sqrt(K/(K+1))*H_los + sqrt(1/(K+1))*H_iid;           
        
        %H0 = squeeze(H(k,:,:)); % k番目の試行回数での伝搬チャネル行列
        
        % BMSN-BF algorithm
        [W_BMSN_BF,U_BMSN_BF,S_BMSN_BF,RIPBF,~] = bmsn_bf(Nt,Nr,Nu,H0,a,T); % function bmsn_bf.m を使用
        
        % BMSN-GE algorithm
        [W_BMSN_GE,U_BMSN_GE,S_BMSN_GE,RIPGE,~] = bmsn_gev(Nt,Nr,Nu,H0,a); % function bmsn_gev.m を使用
        
        % MMSE-CI algorithm
        [W_MMSE,U_MMSE,S_MMSE,RIPM,~] = mmse(Nt,Nr,Nu,H0,a); % function mmse.m を使用
        
        % GMI1 algorithm
        [W_GMI1,U_GMI1,S_GMI1,RIPGM1,~] = gmmse_m1(Nt,Nr,Nu,H0,a); % function mmse.m を使用
        
        % GMI2 algorithm
        [W_GMI2,U_GMI2,S_GMI2,RIPGM2,~] = gmmse_m2(Nt,Nr,Nu,H0,a); % function mmse.m を使用
        
        % BD algorithm
        [W_BD,U_BD,S_BD,RIPBD,~] = bd(Nt,Nr,Nu,H0); % function bd.m を使用
        
        % ZF-CI algorithm
        [W_ZF,U_ZF,S_ZF,RIPZF,~] = zf(Nt,Nr,Nu,H0); % function zf.m を使用
        
        RIP = [RIPBF; RIPGE; RIPM; RIPGM1; RIPGM2; RIPBD; RIPZF];
        S(:,:,:,1) = S_BMSN_BF; S(:,:,:,2) = S_BMSN_GE; S(:,:,:,3) = S_MMSE; S(:,:,:,4) = S_GMI1; 
        S(:,:,:,5) = S_GMI2;S(:,:,:,7) = S_ZF;
        if Nr == 1
            S(:,:,:,6) = S_BD(1,1,:);
        else
            S(:,:,:,6) = S_BD;
        end
        % ユーザ毎の固有値分布
        snt = 1/(10^(SN_tar/10));
        for nalg = 1:numel(Algorithms)
            for nuser=1:Nu
                ns = Nr*(nuser-1)+1:Nr*nuser;
                E(isimu,:,nuser,nalg) = 10*log10((diag(S(:,:,nuser,nalg)).^2)./(RIP(Nr*nalg,nuser)+Nt*snt));
                Eigs(isimu,:,nuser,nalg) = diag(S(:,:,nuser,nalg)); % 固有値そのものの値
                MSt(isimu,:,nuser,nalg)=(diag(S(:,:,nuser,nalg)).^2)./(RIP(Nr*nalg,nuser)+Nt*snt);
                C(isimu,nuser,nalg)=sum(log2(1 + MSt(isimu,:,nuser,nalg)),2); % 各アルゴリズムのChannel capacityを計算
            end            
        end
    end
    for nuser=1:Nu
        E(:,:,nuser,:) = sort(E(:,:,nuser,:),1);
        Eigs(:,:,nuser,:) = sort(Eigs(:,:,nuser,:),1);
        C(:,nuser,:) = sort(C(:,nuser,:),1);
    end
    % ユーザ平均
    E = mean(E,3);
    Eigs = mean(Eigs,3);
    C = mean(C,2);
    %さらにEigs1とEigs2の平均
    EigsAll = mean(Eigs,2);
    
    % CDF of Eigenvalue at Target SNR
    Y = (1/SIMU:1/SIMU:1).'*100;
    if Directivity_switch == 1 % 送受信素子の指向性考慮の有
        rr3_cos = zeros(SIMU,Nr*numel(Algorithms));
        ee3_cos = zeros(SIMU,Nr*numel(Algorithms));        
        for nalg = 1:numel(Algorithms)
            for nn=1:Nr            
                rr3_cos(:,nn+(nalg-1)*Nr) = E(:,nn,nalg);
                ee3_cos(:,nn+(nalg-1)*Nr) = Eigs(:,nn,nalg);
            end
        end
        cc3_cos = C;
        ss3_cos = EigsAll;
        
    else % 送受信素子の指向性考慮の無
        rr3 = zeros(SIMU,Nr*numel(Algorithms));
        ee3 = zeros(SIMU,Nr*numel(Algorithms)); 
        for nalg = 1:numel(Algorithms)
            for nn=1:Nr            
                rr3(:,nn+(nalg-1)*Nr) = E(:,nn,nalg);
                ee3(:,nn+(nalg-1)*Nr) = Eigs(:,nn,nalg);
            end
        end
        cc3 = C;
        ss3 = EigsAll;
    end
end
%% CDF of EGV at Target SNR
figure;
mycol = [1 0 0;0 0 1;0 1 0;1 1 0;0 0 0;0 1 1;1 0 1
        1 0 0;0 0 1;0 1 0;1 1 0;0 0 0;0 1 1;1 0 1]; % 色
set(groot,'defaultAxesColorOrder',mycol)
Holizon_min = round(min(min(rr3)))-5;
Holizon_max = round(max(max(rr3_cos)))+5;
axis([Holizon_min Holizon_max 0 100]);
grid on;
hold on;
set(gca,'XTick',Holizon_min:5:Holizon_max,'Fontsize',14,'Fontname','Times New Roman')
xlabel('SINR of eigenvalue [dB]','Fontsize',16,'Fontname','Times New Roman');
ylabel('CDF [%]','Fontsize',16,'Fontname','Times New Roman');
if Nr==1
    plot(rr3(:,1),Y,'r-','Linewidth',2); plot(rr3(:,2),Y,'b-','Linewidth',2);
    plot(rr3(:,3),Y,'g-','Linewidth',2); plot(rr3(:,4),Y,'y-','Linewidth',2);
    plot(rr3(:,5),Y,'k-','Linewidth',2); plot(rr3(:,6),Y,'c-','Linewidth',2);
    plot(rr3(:,7),Y,'m-','Linewidth',2);
    plot(rr3_cos(:,1),Y,'r-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,2),Y,'b-o','MarkerIndices',50:100:length(Y),'Linewidth',2); 
    plot(rr3_cos(:,3),Y,'g-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,4),Y,'y-o','MarkerIndices',50:100:length(Y),'Linewidth',2); 
    plot(rr3_cos(:,5),Y,'k-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,6),Y,'c-o','MarkerIndices',50:100:length(Y),'Linewidth',2);
    plot(rr3_cos(:,7),Y,'m-o','MarkerIndices',10:100:length(Y),'Linewidth',2);
    legend('BMSN-BF \lambda_1','BMSN-GE \lambda_1','MMSE-CI \lambda_1','GMI1 \lambda_1','GMI2 \lambda_1','BD \lambda_1','ZF-CI \lambda_1',...
        'BMSN-BFcos \lambda_1','BMSN-GEcos \lambda_1','MMSE-CIcos \lambda_1','GMI1cos \lambda_1','GMI2cos \lambda_1','BDcos \lambda_1','ZF-CIcos \lambda_1','Location','southeast');
end

if Nr==2
    plot(rr3(:,1),Y,'r-','Linewidth',2); plot(rr3(:,2),Y,'r--','Linewidth',2);
    plot(rr3(:,3),Y,'b-','Linewidth',2); plot(rr3(:,4),Y,'b--','Linewidth',2);
    plot(rr3(:,5),Y,'g-','Linewidth',2); plot(rr3(:,6),Y,'g--','Linewidth',2);
    plot(rr3(:,7),Y,'y-','Linewidth',2); plot(rr3(:,8),Y,'Y--','Linewidth',2);
    plot(rr3(:,9),Y,'k-','Linewidth',2); plot(rr3(:,10),Y,'k--','Linewidth',2);
    plot(rr3(:,11),Y,'c-','Linewidth',2); plot(rr3(:,12),Y,'c--','Linewidth',2);
    plot(rr3(:,13),Y,'m-','Linewidth',2); plot(rr3(:,14),Y,'m--','Linewidth',2);
    
    plot(rr3_cos(:,1),Y,'r-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,2),Y,'r--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
    plot(rr3_cos(:,3),Y,'b-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,4),Y,'b--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
    plot(rr3_cos(:,5),Y,'g-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,6),Y,'g--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
    plot(rr3_cos(:,7),Y,'y-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,8),Y,'Y--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
    plot(rr3_cos(:,9),Y,'k-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,10),Y,'k--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
    plot(rr3_cos(:,11),Y,'c-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,12),Y,'c--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
    plot(rr3_cos(:,13),Y,'m-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,14),Y,'m--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
    lgd = legend;
    lgd.NumColumns = 2; % 凡例の列数を指定
    legend('BMSN-BF \lambda_1','BMSN-BF \lambda_2','BMSN-GE \lambda_1','BMSN-GE \lambda_2','MMSE-CI \lambda_1','MMSE-CI \lambda_2',...
        'GMI1 \lambda_1','GMI1 \lambda_2','GMI2 \lambda_1','GMI2 \lambda_2','BD \lambda_1','BD \lambda_2','ZF-CI \lambda_1','ZF-CI \lambda_2',...
        'BMSN-BFcos \lambda_1','BMSN-BFcos \lambda_2','BMSN-GEcos \lambda_1','BMSN-GEcos \lambda_2','MMSE-CIcos \lambda_1','MMSE-CIcos \lambda_2',...
        'GMI1cos \lambda_1','GMI1cos \lambda_2','GMI2cos \lambda_1','GMI2cos \lambda_2','BDcos \lambda_1','BDcos \lambda_2',...
        'ZF-CIcos \lambda_1','ZF-CIcos \lambda_2','Location','southeast');
end
title(strcat(target_snr,',',target_K));

 %% Capacity
    figure;
    mycol = [1 0 0;0 0 1;0 1 0;1 1 0;0 0 0;0 1 1;1 0 1]; % 色
        set(groot,'defaultAxesColorOrder',mycol)        
        axis([0 round(max(max(cc3_cos)))+0.5 0 100]);        
        for nalg = 1:numel(Algorithms)
            plot(cc3(:,nalg),Y,'--','Linewidth',2);
            grid on; hold on;
        end
        for nalg = 1:numel(Algorithms)
            plot(cc3_cos(:,nalg),Y,'-o','MarkerIndices',50:100:length(Y),'Linewidth',2);
            grid on; hold on;
        end
        set(gca,'XTick',0:1:round(max(max(C)))+0.5,'Fontsize',14,'Fontname','Times New Roman')
        xlabel('Average Capacity','Fontsize',16,'Fontname','Times New Roman');
        ylabel('CDF [%]','Fontsize',16,'Fontname','Times New Roman');
        lgd = legend;
        lgd.NumColumns = 2; % 凡例の列数を指定
        legend([Algorithms strcat(Algorithms,"cos")],'Location','southeast');
        title(strcat(target_snr,',',target_K));
%% CDF of EGV_nonSINR at Target SNR 
if figure_switch == 1
    figure;
    if Nr==2
        mycol = [1 0 0;0 0 1;0 1 0;1 1 0;0 0 0;0 1 1;1 0 1
            1 0 0;0 0 1;0 1 0;1 1 0;0 0 0;0 1 1;1 0 1]; % 色
        set(groot,'defaultAxesColorOrder',mycol)        
        for nalg = 1:2:numel(Algorithms)*2-1
            plot(ee3(:,nalg),Y,'-','Linewidth',2);
            grid on; hold on;
        end
        for nalg = 2:2:numel(Algorithms)*2
            plot(ee3(:,nalg),Y,'--','Linewidth',2);
            grid on; hold on;
        end
        lgd = legend;
        lgd.NumColumns = 2; % 凡例の列数を指定
        legend([strcat(Algorithms," \lambda_1") strcat(Algorithms," \lambda_2")],'Location','southeast');
    end
    if Nr==1        
        mycol = [1 0 0;0 0 1;0 1 0;1 1 0;0 0 0;0 1 1;1 0 1]; % 色
        set(groot,'defaultAxesColorOrder',mycol)
        for nalg = 1:numel(Algorithms)
            plot(ee3(:,nalg),Y,'-','Linewidth',2);
            grid on; hold on;
        end
        legend(strcat(Algorithms," \lambda_1"),'Location','southeast');
    end
    axis([0 round(max(max(ee3)))+3 0 100]);
    set(gca,'XTick',0:1:round(max(max(ee3)))+3,'Fontsize',14,'Fontname','Times New Roman')
    xlabel('Eigenvalue','Fontsize',16,'Fontname','Times New Roman');
    ylabel('CDF [%]','Fontsize',16,'Fontname','Times New Roman');
    title(strcat(target_snr,',',target_K));   
    %% CDF of EGV at Target SNR 
    if Nr==2
        figure;
        mycol = [1 0 0;0 0 1;0 1 0;1 1 0;0 0 0;0 1 1;1 0 1]; % 色
        set(groot,'defaultAxesColorOrder',mycol)        
        axis([0 round(max(max(EigsAll)))+3 0 100]);
        grid on; hold on;
        for nalg = 1:numel(Algorithms)
            plot(EigsAll(:,nalg),Y,'-','Linewidth',2);
        end
        set(gca,'XTick',0:1:round(max(max(EigsAll)))+3,'Fontsize',14,'Fontname','Times New Roman')
        xlabel('Average Eigenvalue','Fontsize',16,'Fontname','Times New Roman');
        ylabel('CDF [%]','Fontsize',16,'Fontname','Times New Roman');  
        legend(strcat(Algorithms,"avr"),'Location','southeast');
        title(strcat(target_snr,',',target_K));
    end
end
% End