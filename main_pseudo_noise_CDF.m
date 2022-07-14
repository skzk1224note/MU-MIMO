clear;
% パラメータ条件 NT >= NR*NU

SN_tar = 10;         % CDF表示のためのターゲットSNR [dB]
SIMU   = 1000;       % 試行回数
Nt     = 16;         % 送信素子数
Nr     = 2;          % 受信素子数(=2に固定)
Nu     = 8;          % ユーザ数
NL     = 20;         % マルチパス波の素波数
I      = eye(Nt,Nt); % NTxNTの単位行列

d_t = 0.5;      % 送信アンテナ間隔（in wavelength)
d_r = 0.5;      % 受信アンテナ間隔（in wavelength)
derad = pi/180;      % degree -> rad

% cosθのn乗（nは整数値）
cos_n = 0;
if rem(cos_n, 2) == 0 % nが偶数のとき
    % 二重階乗
    a_box = 2:2:cos_n;
    aa = prod(a_box); 
    b = 1:2:cos_n+1;
    bb = prod(b);
    D = bb/aa;
else                  % nが奇数のとき
    % 二重階乗
    a_box = 1:2:cos_n;
    aa = prod(a_box);
    b = 2:2:cos_n+1;
    bb = prod(b);
    D = (2*bb) / (pi*aa);
end

An = sqrt(D)*sqrt(1); % 半球のときは1，全球のときは2
%An = 2; % 指向性関数の係数

a = Nt/(10^(SN_tar/10)); % 擬似雑音Nt * sigma^2 , sigma^2 = 電力　1/(10^(SNR/10))
T = zeros(Nr,Nr); % 所望のチャネル行列 for BMSN
for nuser = 1:Nu
    T(:,:,nuser) = eye(Nr,Nr);
end
Algorithms = ["BMSN-BF","BMSN-GE","BD"]; %,"BMSN-GE3","MMSE","GMI1","GMI2","ZF"
S = zeros(Nr, Nr, Nu, numel(Algorithms));
E_cos = zeros(SIMU, Nr, Nu, numel(Algorithms));
E_iso = zeros(SIMU, Nr, Nu, numel(Algorithms));
a_t_iid_cos = zeros(Nt, 1, NL);a_r_iid_cos = zeros(Nr, 1, NL);
a_t_iid_iso = zeros(Nt, 1, NL);a_r_iid_iso = zeros(Nr, 1, NL);
X_cos = zeros(Nr, Nt, NL);X_iso = zeros(Nr, Nt, NL);
if SN_tar < 10
        target_snr=strcat('SNR=',num2str(SN_tar,'%01d'),'dB');
else
        target_snr=strcat('SNR=',num2str(SN_tar,'%02d'),'dB');
end

for isimu = 1:SIMU % 試行回数のループ
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
            a_t_iid_iso(:,:,L) = sqrt(1/Nt) * exp(1j*Thetat*derad) * exp(-1j*2*pi*d_t*(0:Nt-1).' * sin(Thetat*derad)); % ユーザ毎のNLoS送信モードベクトル An*cos(Theta_t(1,n)*derad) *
            a_r_iid_iso(:,:,L) = sqrt(1/Nr) * exp(-1j*2*pi*d_r*(0:Nr-1).' * sin(Thetar*derad)); % ユーザ毎のNLoS受信モードベクトル An*cos(Theta_r(1,n)*derad) *
            X_iso(:,:,L) = (a_r_iid_iso(:,:,L) * a_t_iid_iso(:,:,L)');
        end
        H_iid_cos((n-1)*Nr+1:(n-1)*Nr+Nr,:) = sum(X_cos,3) * sqrt(Nt*Nr/NL);
        H_iid_iso((n-1)*Nr+1:(n-1)*Nr+Nr,:) = sum(X_iso,3) * sqrt(Nt*Nr/NL);
    end
    
    for Directivity_switch = 0:1 % 送受信素子の指向性考慮の有無 0:無,1:有
        if Directivity_switch == 1
            H0 = H_iid_cos;
        else
            H0 = H_iid_iso;
        end
        
        % BMSN-BF algorithm 疑似雑音補正
        [W_BMSN_BF,U_BMSN_BF,S_BMSN_BF,RIPBF,~] = bmsn_bf_pseudo_noise(Nt,Nr,Nu,H0,T); % function bmsn_bf.m を使用
        
        % BMSN-BF algorithm 疑似雑音定数
        [W_BMSN_GE,U_BMSN_GE,S_BMSN_GE,RIPGE,~] = bmsn_bf(Nt,Nr,Nu,H0,a,T); % function bmsn_bf.m を使用
        
        % BMSN-GE algorithm
        %[W_BMSN_GE,U_BMSN_GE,S_BMSN_GE,RIPGE,~] = bmsn_gev_pseudo_noise(Nt,Nr,Nu,H0); % function bmsn_gev.m を使用
        
%         % MMSE-CI algorithm
%         [W_MMSE,U_MMSE,S_MMSE,RIPM,~] = mmse(Nt,Nr,Nu,H0,a); % function mmse.m を使用
%         
%         % GMI1 algorithm
%         [W_GMI1,U_GMI1,S_GMI1,RIPGM1,~] = gmmse_m1(Nt,Nr,Nu,H0,a); % function mmse.m を使用
%         
%         % GMI2 algorithm
%         [W_GMI2,U_GMI2,S_GMI2,RIPGM2,~] = gmmse_m2(Nt,Nr,Nu,H0,a); % function mmse.m を使用
        
        % BD algorithm
        [W_BD,U_BD,S_BD,~,~] = bd(Nt,Nr,Nu,H0); % function bd.m を使用
        
%         % ZF-CI algorithm
%         [W_ZF,U_ZF,S_ZF,~,~] = zf(Nt,Nr,Nu,H0); % function zf.m を使用
        
        S(:,:,:,1) = S_BMSN_BF; S(:,:,:,2) = S_BMSN_GE; %S(:,:,:,3) = S_MMSE; S(:,:,:,4) = S_GMI1; 
%         S(:,:,:,5) = S_GMI2;S(:,:,:,7) = S_ZF;
        if Nr == 1
            S(:,:,:,3) = S_BD(1,1,:);
        else
            S(:,:,:,3) = S_BD;
        end
        % ユーザ毎の固有値分布
        snt = 1/(10^(SN_tar/10));
        
        if Directivity_switch == 1
            for nuser=1:Nu
                if Nr==1
                    E_cos(isimu,:,nuser,strcmp("BMSN-BF",Algorithms)) = 10*log10((S_BMSN_BF(1,1,nuser).^2)./(RIPBF(1,nuser)+Nt*snt));
                    E_cos(isimu,:,nuser,strcmp("BMSN-GE",Algorithms)) = 10*log10((S_BMSN_GE(1,1,nuser).^2)./(RIPGE(1,nuser)+Nt*snt));
%                     E_cos(isimu,:,nuser,strcmp("MMSE",Algorithms)) = 10*log10((S_MMSE(1,1,nuser).^2)./(RIPM(1,nuser)+Nt*snt));
%                     E_cos(isimu,:,nuser,strcmp("GMI1",Algorithms)) = 10*log10((S_GMI1(1,1,nuser).^2)./(RIPGM1(1,nuser)+Nt*snt));
%                     E_cos(isimu,:,nuser,strcmp("GMI2",Algorithms)) = 10*log10((S_GMI2(1,1,nuser).^2)./(RIPGM2(1,nuser)+Nt*snt));
                    E_cos(isimu,:,nuser,strcmp("BD",Algorithms)) = 10*log10(S_BD(1,1,nuser).^2/(Nt*snt));
%                     E_cos(isimu,:,nuser,strcmp("ZF",Algorithms)) = 10*log10(S_ZF(1,1,nuser).^2/(Nt*snt));    
                else
                    E_cos(isimu,:,nuser,strcmp("BMSN-BF",Algorithms)) = 10*log10((diag(S_BMSN_BF(:,:,nuser)).^2)./(RIPBF(:,nuser)+Nt*snt));
                    E_cos(isimu,:,nuser,strcmp("BMSN-GE",Algorithms)) = 10*log10((diag(S_BMSN_GE(:,:,nuser)).^2)./(RIPGE(:,nuser)+Nt*snt));
%                     E_cos(isimu,:,nuser,strcmp("MMSE",Algorithms)) = 10*log10((diag(S_MMSE(:,:,nuser)).^2)./(RIPM(:,nuser)+Nt*snt));
%                     E_cos(isimu,:,nuser,strcmp("GMI1",Algorithms)) = 10*log10((diag(S_GMI1(:,:,nuser)).^2)./(RIPGM1(:,nuser)+Nt*snt));
%                     E_cos(isimu,:,nuser,strcmp("GMI2",Algorithms)) = 10*log10((diag(S_GMI2(:,:,nuser)).^2)./(RIPGM2(:,nuser)+Nt*snt));
                    E_cos(isimu,:,nuser,strcmp("BD",Algorithms)) = 10*log10(diag(S_BD(:,:,nuser)).^2/(Nt*snt));
%                     E_cos(isimu,:,nuser,strcmp("ZF",Algorithms)) = 10*log10(diag(S_ZF(:,:,nuser)).^2/(Nt*snt));
                end
            end
        else
            for nuser=1:Nu
                if Nr==1
                    E_iso(isimu,:,nuser,strcmp("BMSN-BF",Algorithms)) = 10*log10((S_BMSN_BF(1,1,nuser).^2)./(RIPBF(1,nuser)+Nt*snt));
                    E_iso(isimu,:,nuser,strcmp("BMSN-GE",Algorithms)) = 10*log10((S_BMSN_GE(1,1,nuser).^2)./(RIPGE(1,nuser)+Nt*snt));
%                     E_iso(isimu,:,nuser,strcmp("MMSE",Algorithms)) = 10*log10((S_MMSE(1,1,nuser).^2)./(RIPM(1,nuser)+Nt*snt));
%                     E_iso(isimu,:,nuser,strcmp("GMI1",Algorithms)) = 10*log10((S_GMI1(1,1,nuser).^2)./(RIPGM1(1,nuser)+Nt*snt));
%                     E_iso(isimu,:,nuser,strcmp("GMI2",Algorithms)) = 10*log10((S_GMI2(1,1,nuser).^2)./(RIPGM2(1,nuser)+Nt*snt));
                    E_iso(isimu,:,nuser,strcmp("BD",Algorithms)) = 10*log10(S_BD(1,1,nuser).^2/(Nt*snt));
%                     E_iso(isimu,:,nuser,strcmp("ZF",Algorithms)) = 10*log10(S_ZF(1,1,nuser).^2/(Nt*snt));    
                else
                    E_iso(isimu,:,nuser,strcmp("BMSN-BF",Algorithms)) = 10*log10((diag(S_BMSN_BF(:,:,nuser)).^2)./(RIPBF(:,nuser)+Nt*snt));
                    E_iso(isimu,:,nuser,strcmp("BMSN-GE",Algorithms)) = 10*log10((diag(S_BMSN_GE(:,:,nuser)).^2)./(RIPGE(:,nuser)+Nt*snt));
%                     E_iso(isimu,:,nuser,strcmp("MMSE",Algorithms)) = 10*log10((diag(S_MMSE(:,:,nuser)).^2)./(RIPM(:,nuser)+Nt*snt));
%                     E_iso(isimu,:,nuser,strcmp("GMI1",Algorithms)) = 10*log10((diag(S_GMI1(:,:,nuser)).^2)./(RIPGM1(:,nuser)+Nt*snt));
%                     E_iso(isimu,:,nuser,strcmp("GMI2",Algorithms)) = 10*log10((diag(S_GMI2(:,:,nuser)).^2)./(RIPGM2(:,nuser)+Nt*snt));
                    E_iso(isimu,:,nuser,strcmp("BD",Algorithms)) = 10*log10(diag(S_BD(:,:,nuser)).^2/(Nt*snt));
%                     E_iso(isimu,:,nuser,strcmp("ZF",Algorithms)) = 10*log10(diag(S_ZF(:,:,nuser)).^2/(Nt*snt));
                end
            end
        end
    end
end

for nuser=1:Nu
    E_cos(:,:,nuser,:) = sort(E_cos(:,:,nuser,:),1);
    E_iso(:,:,nuser,:) = sort(E_iso(:,:,nuser,:),1);
end
% ユーザ平均
E_cos = mean(E_cos,3);
E_iso = mean(E_iso,3);

% CDF of Eigenvalue at Target SNR
Y = (1/SIMU:1/SIMU:1).'*100;
rr3_cos = zeros(SIMU,Nr*numel(Algorithms));
rr3_iso = zeros(SIMU,Nr*numel(Algorithms));

for nalg = 1:numel(Algorithms)
    for nn=1:Nr            
        rr3_cos(:,nn+(nalg-1)*Nr) = E_cos(:,nn,nalg);
        rr3_iso(:,nn+(nalg-1)*Nr) = E_iso(:,nn,nalg);
    end
end

%% CDF of EGV at Target SNR
figure;
mycol = [1 0 0;0 0 1;0 1 0;1 1 0;0 0 0;0 1 1;1 0 1
        1 0 0;0 0 1;0 1 0;1 1 0;0 0 0;0 1 1;1 0 1]; % 色
set(groot,'defaultAxesColorOrder',mycol)
Holizon_min = round(min(min(rr3_cos)))-5;
Holizon_max = round(max(max(rr3_iso)))+5;
axis([Holizon_min Holizon_max 0 100]);
grid on;
hold on;
set(gca,'XTick',Holizon_min-5:5:Holizon_max+5,'Fontsize',8,'Fontname','Times New Roman')
xlabel('SINR of eigenvalue [dB]','Fontsize',16,'Fontname','Times New Roman');
ylabel('CDF [%]','Fontsize',16,'Fontname','Times New Roman');
if Nr==1
    plot(rr3_iso(:,1),Y,'r-','Linewidth',2); plot(rr3_iso(:,2),Y,'b-','Linewidth',2);
%     plot(rr3(:,3),Y,'g-','Linewidth',2); plot(rr3(:,4),Y,'y-','Linewidth',2);
%     plot(rr3(:,5),Y,'k-','Linewidth',2);
    plot(rr3_iso(:,6),Y,'k-','Linewidth',2);
%     plot(rr3(:,7),Y,'m-','Linewidth',2);
    plot(rr3_cos(:,1),Y,'r-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,2),Y,'b-o','MarkerIndices',50:100:length(Y),'Linewidth',2); 
%     plot(rr3_cos(:,3),Y,'g-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,4),Y,'y-o','MarkerIndices',50:100:length(Y),'Linewidth',2); 
%     plot(rr3_cos(:,5),Y,'k-o','MarkerIndices',10:100:length(Y),'Linewidth',2);
    plot(rr3_cos(:,6),Y,'k-o','MarkerIndices',50:100:length(Y),'Linewidth',2);
%     plot(rr3_cos(:,7),Y,'m-o','MarkerIndices',10:100:length(Y),'Linewidth',2);
    legend('BMSN-BF \lambda_1','BMSN-GE \lambda_1','BD \lambda_1',...
        'BMSN-BFcos \lambda_1','BMSN-GEcos \lambda_1','BDcos \lambda_1','Location','southeast');
end

if Nr==2
    plot(rr3_iso(:,1),Y,'r-','Linewidth',2); plot(rr3_iso(:,2),Y,'r--','Linewidth',2);
    plot(rr3_iso(:,3),Y,'b-','Linewidth',2); plot(rr3_iso(:,4),Y,'b--','Linewidth',2);
%     plot(rr3(:,5),Y,'g-','Linewidth',2); plot(rr3(:,6),Y,'g--','Linewidth',2);
%     plot(rr3(:,7),Y,'y-','Linewidth',2); plot(rr3(:,8),Y,'Y--','Linewidth',2);
%     plot(rr3(:,9),Y,'k-','Linewidth',2); plot(rr3(:,10),Y,'k--','Linewidth',2);
    plot(rr3_iso(:,5),Y,'k-','Linewidth',2); plot(rr3_iso(:,6),Y,'k--','Linewidth',2);
%     plot(rr3(:,13),Y,'m-','Linewidth',2); plot(rr3(:,14),Y,'m--','Linewidth',2);
    
    plot(rr3_cos(:,1),Y,'r-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,2),Y,'r--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
    plot(rr3_cos(:,3),Y,'b-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,4),Y,'b--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
%     plot(rr3_cos(:,5),Y,'g-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,6),Y,'g--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
%     plot(rr3_cos(:,7),Y,'y-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,8),Y,'Y--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
%     plot(rr3_cos(:,9),Y,'k-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,10),Y,'k--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
    plot(rr3_cos(:,5),Y,'k-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,6),Y,'k--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
%     plot(rr3_cos(:,13),Y,'m-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,14),Y,'m--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
    lgd = legend;
    lgd.NumColumns = 1; % 凡例の列数を指定
    legend('BMSN-BF \lambda_1','BMSN-BF \lambda_2','BMSN-GE \lambda_1','BMSN-GE \lambda_2','BD \lambda_1','BD \lambda_2',...
        'BMSN-BFcos \lambda_1','BMSN-BFcos \lambda_2','BMSN-GEcos \lambda_1','BMSN-GEcos \lambda_2',...
        'BDcos \lambda_1','BDcos \lambda_2','Location','southeast');
end
title(strcat(target_snr));

% End