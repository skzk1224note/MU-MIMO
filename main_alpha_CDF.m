clear;
% パラメータ条件 NT >= NR*NU

SN_tar = 30;         % CDF表示のためのターゲットSNR [dB]
SIMU   = 1000;       % 試行回数
Nt     = 16;         % 送信素子数
Nr     = 2;          % 受信素子数(=2に固定)
Nu     = 8;          % ユーザ数
NL     = 20;         % マルチパス波の素波数
I      = eye(Nt,Nt); % NTxNTの単位行列
CDF = 10;
K_dB   = 20;
K = 10^(K_dB/10);    % Kの真値

d_t = 0.5;      % 送信アンテナ間隔（in wavelength)
d_r = 0.5;      % 受信アンテナ間隔（in wavelength)
derad = pi/180;      % degree -> rad

% cosθのn乗（nは整数値）
cos_n = 0;
if rem(cos_n, 2) == 0 % nが偶数のとき
    % 二重階乗
    a = 2:2:cos_n;
    aa = prod(a);
    b = 1:2:cos_n+1;
    bb = prod(b);
    D = bb/aa;
else                  % nが奇数のとき
    % 二重階乗
    a = 1:2:cos_n;
    aa = prod(a);
    b = 2:2:cos_n+1;
    bb = prod(b);
    D = (2*bb) / (pi*aa);
end

An = sqrt(D)*sqrt(1); % 半球のときは1，全球のときは2
%An = 2; % 指向性関数の係数

%a = Nt/(10^(SN_tar/10)); % 擬似雑音Nt * sigma^2 , sigma^2 = 電力　1/(10^(SNR/10))

T = zeros(Nr,Nr); % 所望のチャネル行列 for BMSN
for nuser = 1:Nu
    T(:,:,nuser) = eye(Nr,Nr);
end
Algorithms = ["BMSN-BF","BMSN-GE","BD"]; %,"BMSN-GE3"
S = zeros(Nr, Nr, Nu, numel(Algorithms));
E_cos = zeros(SIMU, Nr, Nu, numel(Algorithms));
a_t_iid_cos = zeros(Nt, 1, NL);a_r_iid_cos = zeros(Nr, 1, NL);
a_t_iid_iso = zeros(Nt, 1, NL);a_r_iid_iso = zeros(Nr, 1, NL);
X_cos = zeros(Nr, Nt, NL);X_iso = zeros(Nr, Nt, NL);
if CDF < 10
   target_CDF=strcat('CDF=',num2str(CDF,'%01d'),'%');
else
   target_CDF=strcat('CDF=',num2str(CDF,'%02d'),'%');
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
a_box = 10^(-4):10^(-4):1;
La = length(a_box);
rr3_cos = zeros(La, numel(Algorithms)*Nr);
tic;
for ia = 1:La
    a = a_box(ia);
for isimu = 1:SIMU % 試行回数のループ
    Theta_t = (rand(1,Nu)-0.5)*180; % ユーザ毎の送信角 指向性:(-90deg - 90deg)
    Theta_r = (rand(1,Nu)-0.5)*180; % ユーザ毎の受信角 指向性:(-90deg - 90deg)
    for n = 1 : Nu
        % 幾何学レイリーチャネル
        for L = 1:NL
            Thetat = (rand(1,1)-0.5)*180; % ユーザ毎の送信角 指向性:(-90deg - 90deg)An*cos(Thetat*derad)
            Thetar = (rand(1,1)-0.5)*180; % ユーザ毎の送信角 指向性:(-90deg - 90deg)An*cos(Thetar*derad)
            gtheta_t = An * ((cos(Thetat*derad))^(cos_n/2));
            gtheta_r = An * ((cos(Thetar*derad))^(cos_n/2));
            a_t_iid_cos(:,:,L) = gtheta_t * sqrt(1/Nt) * exp(1j*Thetat*derad) * exp(-1j*2*pi*d_t*(0:Nt-1).' * sin(Thetat*derad)); % ユーザ毎のNLoS送信モードベクトル An*cos(Theta_t(1,n)*derad) *
            a_r_iid_cos(:,:,L) = gtheta_r * sqrt(1/Nr) * exp(-1j*2*pi*d_r*(0:Nr-1).' * sin(Thetar*derad)); % ユーザ毎のNLoS受信モードベクトル An*cos(Theta_r(1,n)*derad) *
            X_cos(:,:,L) = (a_r_iid_cos(:,:,L) * a_t_iid_cos(:,:,L)');
%             a_t_iid_iso(:,:,L) = sqrt(1/Nt) * exp(1j*Thetat*derad) * exp(-1j*2*pi*d_t*(0:Nt-1).' * sin(Thetat*derad)); % ユーザ毎のNLoS送信モードベクトル An*cos(Theta_t(1,n)*derad) *
%             a_r_iid_iso(:,:,L) = sqrt(1/Nr) * exp(-1j*2*pi*d_r*(0:Nr-1).' * sin(Thetar*derad)); % ユーザ毎のNLoS受信モードベクトル An*cos(Theta_r(1,n)*derad) *
%             X_iso(:,:,L) = (a_r_iid_iso(:,:,L) * a_t_iid_iso(:,:,L)');
        end
        H_iid((n-1)*Nr+1:(n-1)*Nr+Nr,:) = sum(X_cos,3) * sqrt(Nt*Nr/NL);
        a_t = sqrt(1/Nt) * exp(1j*Theta_t(1,n)*derad) * exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad))*An*(cos(Theta_t(1,n)*derad))^(cos_n/2); % ユーザ毎の送信モードベクトル
        a_r = sqrt(1/Nr) * exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad))*An*(cos(Theta_r(1,n)*derad))^(cos_n/2); % ユーザ毎の受信モードベクトル
        H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t' * sqrt(Nt*Nr);                 % ユーザ毎のLOSチャネル行列
    end
    beta = (norm(H_los,'fro'))^2/(norm(H_iid,'fro'))^2; % 直接波とレイリー波を等電力に
    %H_iid = sqrt(beta) * H_iid;
        % H=[sqrt(K/(K+1))*(LOS チャネル)]+[sqrt(1/(K+1))*(NLOS チャネル)]
        H0 = sqrt(K/(K+1))*H_los + sqrt(1/(K+1))*H_iid;
%         
%         (norm(H_los,'fro'))^2
%         (norm(H_iid,'fro'))^2
        
        
        % BMSN-BF algorithm
        [W_BMSN_BF,U_BMSN_BF,S_BMSN_BF,RIPBF,~] = bmsn_bf(Nt,Nr,Nu,H0,a,T); % function bmsn_bf.m を使用
        
        % BMSN-GE algorithm
        [W_BMSN_GE,U_BMSN_GE,S_BMSN_GE,RIPGE,~] = bmsn_gev(Nt,Nr,Nu,H0,a); % function bmsn_gev.m を使用
                
        % BD algorithm
        [W_BD,U_BD,S_BD,~,~] = bd(Nt,Nr,Nu,H0); % function bd.m を使用
        
        
        
        S(:,:,:,1) = S_BMSN_BF; S(:,:,:,2) = S_BMSN_GE;
        if Nr == 1
            S(:,:,:,3) = S_BD(1,1,:);
        else
            S(:,:,:,3) = S_BD;
        end
        % ユーザ毎の固有値分布
        snt = 1/(10^(SN_tar/10));
        
            for nuser=1:Nu
                if Nr==1
                    E_cos(isimu,:,nuser,strcmp("BMSN-BF",Algorithms)) = 10*log10((S_BMSN_BF(1,1,nuser).^2)./(RIPBF(1,nuser)+Nt*snt));
                    E_cos(isimu,:,nuser,strcmp("BMSN-GE",Algorithms)) = 10*log10((S_BMSN_GE(1,1,nuser).^2)./(RIPGE(1,nuser)+Nt*snt));                    
                    E_cos(isimu,:,nuser,strcmp("BD",Algorithms)) = 10*log10(S_BD(1,1,nuser).^2/(Nt*snt));
                else
                    E_cos(isimu,:,nuser,strcmp("BMSN-BF",Algorithms)) = 10*log10((diag(S_BMSN_BF(:,:,nuser)).^2)./(RIPBF(:,nuser)+Nt*snt));
                    E_cos(isimu,:,nuser,strcmp("BMSN-GE",Algorithms)) = 10*log10((diag(S_BMSN_GE(:,:,nuser)).^2)./(RIPGE(:,nuser)+Nt*snt));                    
                    E_cos(isimu,:,nuser,strcmp("BD",Algorithms)) = 10*log10(diag(S_BD(:,:,nuser)).^2/(Nt*snt));
                end
            end
end


for nuser=1:Nu
    E_cos(:,:,nuser,:) = sort(E_cos(:,:,nuser,:),1);
end
% ユーザ平均
E_cos = mean(E_cos,3);

% CDF of Eigenvalue at Target SNR
for nalg = 1:numel(Algorithms)
    for nn=1:Nr            
        rr3_cos(ia,nn+(nalg-1)*Nr) = E_cos(round(CDF*SIMU/100),nn,nalg);
    end
end
end
toc;
%% CDF of EGV at Target SNR
figure;
mycol = [1 0 0;0 0 1;0 1 0;1 1 0;0 0 0;0 1 1;1 0 1
        1 0 0;0 0 1;0 1 0;1 1 0;0 0 0;0 1 1;1 0 1]; % 色
set(groot,'defaultAxesColorOrder',mycol)
Holizon_min = round(min(min(rr3_cos)))-5;
Holizon_max = round(max(max(rr3_cos)))+5;
if Nr==1
    semilogx(a_box,rr3_cos(:,1),'r-',a_box,rr3_cos(:,2),'b-',a_box,rr3_cos(:,3),'k-','Linewidth',2);
    legend('BMSN-BF \lambda_1','BMSN-GE \lambda_1','BD \lambda_1','Location','southeast');
end

if Nr==2
    semilogx(a_box,rr3_cos(:,1),'r-',a_box,rr3_cos(:,2),'r--',a_box,rr3_cos(:,3),'b-',...
        a_box,rr3_cos(:,4),'b--',a_box,rr3_cos(:,5),'k-',a_box,rr3_cos(:,6),'k--','Linewidth',2);
    lgd = legend;
    lgd.NumColumns = 1; % 凡例の列数を指定
    legend('BMSN-BF \lambda_1','BMSN-BF \lambda_2','BMSN-GE \lambda_1','BMSN-GE \lambda_2','BD \lambda_1','BD \lambda_2','Location','southeast');
end
axis([1e-4 1e0 Holizon_min Holizon_max]);
set(gca,'XTick',[1e-4, 1e-3, 1e-2, 1e-1, 1e0],'YTick',Holizon_min:10:Holizon_max,'Fontsize',8,'Fontname','Times New Roman')
xlabel('Pseudo noise','Fontsize',16,'Fontname','Times New Roman');
ylabel('SINR of eigenvalue [dB]','Fontsize',16,'Fontname','Times New Roman');
title(strcat(target_CDF,',',target_snr,',',target_K));
grid on;
hold on;
% End