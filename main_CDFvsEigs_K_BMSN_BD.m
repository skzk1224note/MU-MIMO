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
cos_n = 4;
D = 2 * (cos_n + 1);

An = sqrt(D);

a = Nt/(10^(SN_tar/10)); % 擬似雑音Nt * sigma^2 , sigma^2 = 電力　1/(10^(SNR/10))
T = zeros(Nr,Nr); % 所望のチャネル行列 for BMSN
for nuser = 1:Nu
    T(:,:,nuser) = eye(Nr,Nr);
end
Algorithms = ["BMSN-BF","BMSN-GE","BD"]; %,"BMSN-GE3","MMSE","GMI1","GMI2","ZF"
S = zeros(Nr, Nr, Nu, numel(Algorithms));
E_cos4 = zeros(SIMU, Nr, Nu, numel(Algorithms));
E_cos0 = zeros(SIMU, Nr, Nu, numel(Algorithms));
E = zeros(SIMU, Nr, Nu, numel(Algorithms));
a_t_iid_cos4 = zeros(Nt, 1, NL);a_r_iid_cos4 = zeros(Nr, 1, NL);
a_t_iid_cos0 = zeros(Nt, 1, NL);a_r_iid_cos0 = zeros(Nr, 1, NL);
a_t_iid = zeros(Nt, 1, NL);a_r_iid = zeros(Nr, 1, NL);
X_cos4 = zeros(Nr, Nt, NL);X_cos0 = zeros(Nr, Nt, NL);X = zeros(Nr, Nt, NL);
if SN_tar < 10
        target_snr=strcat('SNR=',num2str(SN_tar,'%01d'),'dB');
else
        target_snr=strcat('SNR=',num2str(SN_tar,'%02d'),'dB');
end

for isimu = 1:SIMU % 試行回数のループ
    for n = 1 : Nu
        % 幾何学レイリーチャネル
        for L = 1:NL
            Thetat = (rand(1,1)-0.5)*90; % ユーザ毎の送信角 指向性:(-90deg - 90deg)An*cos(Thetat*derad)
            Thetar = (rand(1,1)-0.5)*90; % ユーザ毎の送信角 指向性:(-90deg - 90deg)An*cos(Thetar*derad)
            Thetat_iso = (rand(1,1)-0.5)*180; % ユーザ毎の送信角 指向性:(-90deg - 90deg)An*cos(Thetat*derad)
            Thetar_iso = (rand(1,1)-0.5)*180; % ユーザ毎の送信角 指向性:(-90deg - 90deg)An*cos(Thetar*derad)
            gtheta_t = An * ((cos(Thetat*derad))^(cos_n/2));
            % gtheta_r = An * ((cos(Thetar*derad))^(cos_n/2));
            a_t_iid_cos4(:,:,L) = gtheta_t * sqrt(1/Nt) * exp(1j*Thetat*derad) * exp(-1j*2*pi*d_t*(0:Nt-1).' * sin(Thetat*derad)); % ユーザ毎のNLoS送信モードベクトル An*cos(Theta_t(1,n)*derad) *
            a_r_iid_cos4(:,:,L) = sqrt(1/Nr) * exp(-1j*2*pi*d_r*(0:Nr-1).' * sin(Thetar_iso*derad)); % ユーザ毎のNLoS受信モードベクトル An*cos(Theta_r(1,n)*derad) *
            X_cos4(:,:,L) = (a_r_iid_cos4(:,:,L) * a_t_iid_cos4(:,:,L)');
            
            a_t_iid_cos0(:,:,L) = sqrt(2) * sqrt(1/Nt) * exp(1j*Thetat*derad) * exp(-1j*2*pi*d_t*(0:Nt-1).' * sin(Thetat*derad)); % ユーザ毎のNLoS送信モードベクトル An*cos(Theta_t(1,n)*derad) *
            a_r_iid_cos0(:,:,L) = sqrt(1/Nr) * exp(-1j*2*pi*d_r*(0:Nr-1).' * sin(Thetar_iso*derad)); % ユーザ毎のNLoS受信モードベクトル An*cos(Theta_r(1,n)*derad) *
            X_cos0(:,:,L) = (a_r_iid_cos0(:,:,L) * a_t_iid_cos0(:,:,L)');
            
            a_t_iid(:,:,L) = sqrt(1/Nt) * exp(1j*Thetat_iso*derad) * exp(-1j*2*pi*d_t*(0:Nt-1).' * sin(Thetat_iso*derad)); % ユーザ毎のNLoS送信モードベクトル
            a_r_iid(:,:,L) = sqrt(1/Nr) * exp(-1j*2*pi*d_r*(0:Nr-1).' * sin(Thetar_iso*derad)); % ユーザ毎のNLoS受信モードベクトル
            X(:,:,L) = (a_r_iid(:,:,L) * a_t_iid(:,:,L)');
        end
        H_iid_cos4((n-1)*Nr+1:(n-1)*Nr+Nr,:) = sum(X_cos4,3) * sqrt(Nt*Nr/NL);
        H_iid_cos0((n-1)*Nr+1:(n-1)*Nr+Nr,:) = sum(X_cos0,3) * sqrt(Nt*Nr/NL);
        H_iid((n-1)*Nr+1:(n-1)*Nr+Nr,:) = sum(X,3) * sqrt(Nt*Nr/NL);
    end
    
    for Directivity_switch = 0:2 % 送受信素子の指向性考慮の有無 0:無,1:有
        if Directivity_switch == 1
            H0 = H_iid_cos4;
        elseif Directivity_switch == 0
            H0 = H_iid_cos0;
        else
            H0 = H_iid;
        end
        
        % BMSN-BF algorithm
        [W_BMSN_BF,U_BMSN_BF,S_BMSN_BF,RIPBF,~] = bmsn_bf(Nt,Nr,Nu,H0,a,T); % function bmsn_bf.m を使用
        
        % BMSN-GE algorithm
        [W_BMSN_GE,U_BMSN_GE,S_BMSN_GE,RIPGE,~] = bmsn_gev(Nt,Nr,Nu,H0,a); % function bmsn_gev.m を使用
        
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
        
        S(:,:,:,1) = S_BMSN_BF; S(:,:,:,2) = S_BMSN_GE;
%         S(:,:,:,3) = S_MMSE; S(:,:,:,4) = S_GMI1; 
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
                    E_cos4(isimu,:,nuser,strcmp("BMSN-BF",Algorithms)) = 10*log10((S_BMSN_BF(1,1,nuser).^2)./(RIPBF(1,nuser)+Nt*snt));
                    E_cos4(isimu,:,nuser,strcmp("BMSN-GE",Algorithms)) = 10*log10((S_BMSN_GE(1,1,nuser).^2)./(RIPGE(1,nuser)+Nt*snt));
                    E_cos4(isimu,:,nuser,strcmp("BD",Algorithms)) = 10*log10(S_BD(1,1,nuser).^2/(Nt*snt));
%                     E_cos(isimu,:,nuser,strcmp("MMSE",Algorithms)) = 10*log10((S_MMSE(1,1,nuser).^2)./(RIPM(1,nuser)+Nt*snt));
%                     E_cos(isimu,:,nuser,strcmp("GMI1",Algorithms)) = 10*log10((S_GMI1(1,1,nuser).^2)./(RIPGM1(1,nuser)+Nt*snt));
%                     E_cos(isimu,:,nuser,strcmp("GMI2",Algorithms)) = 10*log10((S_GMI2(1,1,nuser).^2)./(RIPGM2(1,nuser)+Nt*snt));
                    
                    E_cos4(isimu,:,nuser,strcmp("ZF",Algorithms)) = 10*log10(S_ZF(1,1,nuser).^2/(Nt*snt));    
                else
                    E_cos4(isimu,:,nuser,strcmp("BMSN-BF",Algorithms)) = 10*log10((diag(S_BMSN_BF(:,:,nuser)).^2)./(RIPBF(:,nuser)+Nt*snt));
                    E_cos4(isimu,:,nuser,strcmp("BMSN-GE",Algorithms)) = 10*log10((diag(S_BMSN_GE(:,:,nuser)).^2)./(RIPGE(:,nuser)+Nt*snt));
                    E_cos4(isimu,:,nuser,strcmp("BD",Algorithms)) = 10*log10(diag(S_BD(:,:,nuser)).^2/(Nt*snt));
%                     E_cos(isimu,:,nuser,strcmp("MMSE",Algorithms)) = 10*log10((diag(S_MMSE(:,:,nuser)).^2)./(RIPM(:,nuser)+Nt*snt));
%                     E_cos(isimu,:,nuser,strcmp("GMI1",Algorithms)) = 10*log10((diag(S_GMI1(:,:,nuser)).^2)./(RIPGM1(:,nuser)+Nt*snt));
%                     E_cos(isimu,:,nuser,strcmp("GMI2",Algorithms)) = 10*log10((diag(S_GMI2(:,:,nuser)).^2)./(RIPGM2(:,nuser)+Nt*snt));
%                     E_cos(isimu,:,nuser,strcmp("ZF",Algorithms)) = 10*log10(diag(S_ZF(:,:,nuser)).^2/(Nt*snt));
                end
            end
        elseif Directivity_switch == 0
            for nuser=1:Nu
                if Nr==1
                    E_cos0(isimu,:,nuser,strcmp("BMSN-BF",Algorithms)) = 10*log10((S_BMSN_BF(1,1,nuser).^2)./(RIPBF(1,nuser)+Nt*snt));
                    E_cos0(isimu,:,nuser,strcmp("BMSN-GE",Algorithms)) = 10*log10((S_BMSN_GE(1,1,nuser).^2)./(RIPGE(1,nuser)+Nt*snt));
                    E_cos0(isimu,:,nuser,strcmp("BD",Algorithms)) = 10*log10(S_BD(1,1,nuser).^2/(Nt*snt));
%                     E_iso(isimu,:,nuser,strcmp("MMSE",Algorithms)) = 10*log10((S_MMSE(1,1,nuser).^2)./(RIPM(1,nuser)+Nt*snt));
%                     E_iso(isimu,:,nuser,strcmp("GMI1",Algorithms)) = 10*log10((S_GMI1(1,1,nuser).^2)./(RIPGM1(1,nuser)+Nt*snt));
%                     E_iso(isimu,:,nuser,strcmp("GMI2",Algorithms)) = 10*log10((S_GMI2(1,1,nuser).^2)./(RIPGM2(1,nuser)+Nt*snt));
%                     E_iso(isimu,:,nuser,strcmp("ZF",Algorithms)) = 10*log10(S_ZF(1,1,nuser).^2/(Nt*snt));    
                else
                    E_cos0(isimu,:,nuser,strcmp("BMSN-BF",Algorithms)) = 10*log10((diag(S_BMSN_BF(:,:,nuser)).^2)./(RIPBF(:,nuser)+Nt*snt));
                    E_cos0(isimu,:,nuser,strcmp("BMSN-GE",Algorithms)) = 10*log10((diag(S_BMSN_GE(:,:,nuser)).^2)./(RIPGE(:,nuser)+Nt*snt));
                    E_cos0(isimu,:,nuser,strcmp("BD",Algorithms)) = 10*log10(diag(S_BD(:,:,nuser)).^2/(Nt*snt));
%                     E_iso(isimu,:,nuser,strcmp("MMSE",Algorithms)) = 10*log10((diag(S_MMSE(:,:,nuser)).^2)./(RIPM(:,nuser)+Nt*snt));
%                     E_iso(isimu,:,nuser,strcmp("GMI1",Algorithms)) = 10*log10((diag(S_GMI1(:,:,nuser)).^2)./(RIPGM1(:,nuser)+Nt*snt));
%                     E_iso(isimu,:,nuser,strcmp("GMI2",Algorithms)) = 10*log10((diag(S_GMI2(:,:,nuser)).^2)./(RIPGM2(:,nuser)+Nt*snt));
%                     E_iso(isimu,:,nuser,strcmp("ZF",Algorithms)) = 10*log10(diag(S_ZF(:,:,nuser)).^2/(Nt*snt));
                end
            end
        else
            for nuser=1:Nu
                if Nr==1
                    E(isimu,:,nuser,strcmp("BMSN-BF",Algorithms)) = 10*log10((S_BMSN_BF(1,1,nuser).^2)./(RIPBF(1,nuser)+Nt*snt));
                    E(isimu,:,nuser,strcmp("BMSN-GE",Algorithms)) = 10*log10((S_BMSN_GE(1,1,nuser).^2)./(RIPGE(1,nuser)+Nt*snt));
                    E(isimu,:,nuser,strcmp("BD",Algorithms)) = 10*log10(S_BD(1,1,nuser).^2/(Nt*snt));
%                     E(isimu,:,nuser,strcmp("MMSE",Algorithms)) = 10*log10((S_MMSE(1,1,nuser).^2)./(RIPM(1,nuser)+Nt*snt));
%                     E(isimu,:,nuser,strcmp("GMI1",Algorithms)) = 10*log10((S_GMI1(1,1,nuser).^2)./(RIPGM1(1,nuser)+Nt*snt));
%                     E(isimu,:,nuser,strcmp("GMI2",Algorithms)) = 10*log10((S_GMI2(1,1,nuser).^2)./(RIPGM2(1,nuser)+Nt*snt));
%                     E(isimu,:,nuser,strcmp("ZF",Algorithms)) = 10*log10(S_ZF(1,1,nuser).^2/(Nt*snt));    
                else
                    E(isimu,:,nuser,strcmp("BMSN-BF",Algorithms)) = 10*log10((diag(S_BMSN_BF(:,:,nuser)).^2)./(RIPBF(:,nuser)+Nt*snt));
                    E(isimu,:,nuser,strcmp("BMSN-GE",Algorithms)) = 10*log10((diag(S_BMSN_GE(:,:,nuser)).^2)./(RIPGE(:,nuser)+Nt*snt));
                    E(isimu,:,nuser,strcmp("BD",Algorithms)) = 10*log10(diag(S_BD(:,:,nuser)).^2/(Nt*snt));
%                     E(isimu,:,nuser,strcmp("MMSE",Algorithms)) = 10*log10((diag(S_MMSE(:,:,nuser)).^2)./(RIPM(:,nuser)+Nt*snt));
%                     E(isimu,:,nuser,strcmp("GMI1",Algorithms)) = 10*log10((diag(S_GMI1(:,:,nuser)).^2)./(RIPGM1(:,nuser)+Nt*snt));
%                     E(isimu,:,nuser,strcmp("GMI2",Algorithms)) = 10*log10((diag(S_GMI2(:,:,nuser)).^2)./(RIPGM2(:,nuser)+Nt*snt));
%                     E(isimu,:,nuser,strcmp("ZF",Algorithms)) = 10*log10(diag(S_ZF(:,:,nuser)).^2/(Nt*snt));
                end
            end
        end
    end
end

for nuser=1:Nu
    E_cos4(:,:,nuser,:) = sort(E_cos4(:,:,nuser,:),1);
    E_cos0(:,:,nuser,:) = sort(E_cos0(:,:,nuser,:),1);
    E(:,:,nuser,:) = sort(E(:,:,nuser,:),1);
end
% ユーザ平均
E_cos4 = mean(E_cos4,3);
E_cos0 = mean(E_cos0,3);
E = mean(E,3);

% CDF of Eigenvalue at Target SNR
Y = (1/SIMU:1/SIMU:1).'*100;
rr3_cos4 = zeros(SIMU,Nr*numel(Algorithms));
rr3_cos = zeros(SIMU,Nr*numel(Algorithms));
rr3 = zeros(SIMU,Nr*numel(Algorithms));

for nalg = 1:numel(Algorithms)
    for nn=1:Nr            
        rr3_cos4(:,nn+(nalg-1)*Nr) = E_cos4(:,nn,nalg);
        rr3_cos(:,nn+(nalg-1)*Nr) = E_cos0(:,nn,nalg);
        rr3(:,nn+(nalg-1)*Nr) = E(:,nn,nalg);
    end
end

%% CDF of EGV at Target SNR
figure;
mycol = [1 0 0;0 0 1;0 1 0;1 1 0;0 0 0;0 1 1;1 0 1
        1 0 0;0 0 1;0 1 0;1 1 0;0 0 0;0 1 1;1 0 1]; % 色
set(groot,'defaultAxesColorOrder',mycol)
Holizon_min = round(min(min(rr3(:,6))))-5;
Holizon_max = round(max(max(rr3)))+5;
axis([Holizon_min Holizon_max 0 100]);
grid on;
hold on;
set(gca,'XTick',Holizon_min-5:5:Holizon_max+5,'Fontsize',8,'Fontname','Times New Roman')
xlabel('SINR of eigenvalue [dB]','Fontsize',16,'Fontname','Times New Roman');
ylabel('CDF [%]','Fontsize',16,'Fontname','Times New Roman');
if Nr==1
    plot(rr3_cos(:,1),Y,'r-','Linewidth',2); plot(rr3_cos(:,2),Y,'b-','Linewidth',2);
    plot(rr3_cos(:,3),Y,'k-','Linewidth',2);
    plot(rr3_cos4(:,1),Y,'r-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos4(:,2),Y,'b-o','MarkerIndices',50:100:length(Y),'Linewidth',2); 
    plot(rr3_cos4(:,3),Y,'k-o','MarkerIndices',50:100:length(Y),'Linewidth',2);
    legend('BMSN-BF \lambda_1','BMSN-GE \lambda_1','BD \lambda_1',...
        'BMSN-BFcos \lambda_1','BMSN-GEcos \lambda_1','BDcos \lambda_1','Location','southeast');
end

if Nr==2
    plot(rr3(:,1),Y,'r-','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3(:,2),Y,'r--','MarkerIndices',50:100:length(Y),'Linewidth',2);
    plot(rr3(:,3),Y,'b-','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3(:,4),Y,'b--','MarkerIndices',50:100:length(Y),'Linewidth',2);
    plot(rr3(:,5),Y,'k-','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3(:,6),Y,'k--','MarkerIndices',50:100:length(Y),'Linewidth',2);
%     plot(rr3_iso(:,1),Y,'r-','Linewidth',2); plot(rr3_iso(:,2),Y,'r--','Linewidth',2);
%     plot(rr3_iso(:,3),Y,'b-','Linewidth',2); plot(rr3_iso(:,4),Y,'b--','Linewidth',2);
%     plot(rr3_iso(:,5),Y,'b-','Linewidth',2); plot(rr3_iso(:,6),Y,'b--','Linewidth',2);
%     
%     plot(rr3_cos(:,1),Y,'r-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,2),Y,'r--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
%     plot(rr3_cos(:,3),Y,'b-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,4),Y,'b--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
%     plot(rr3_cos(:,5),Y,'r-','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,6),Y,'r--','MarkerIndices',50:100:length(Y),'Linewidth',2);
    
%       plot(rr3(:,1),Y,'r-^','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3(:,2),Y,'r--^','MarkerIndices',50:100:length(Y),'Linewidth',2);
%       plot(rr3(:,3),Y,'b-^','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3(:,4),Y,'b--^','MarkerIndices',50:100:length(Y),'Linewidth',2);
      
    
    lgd = legend;
    lgd.NumColumns = 1; % 凡例の列数を指定
     legend('BMSN-BFiso \lambda_1','BMSN-BFiso \lambda_2',...
        'BMSN-GEiso \lambda_1','BMSN-GEiso \lambda_2',...
        'BDiso \lambda_1','BDiso \lambda_2',...
        'Location','southeast');
%         'BMSN-BFcos^4 \lambda_1','BMSN-BFcos^4 \lambda_2',...
%         'BMSN-BFiso \lambda_1','BMSN-BFiso \lambda_2',...
%         'BMSN-GEiso \lambda_1','BMSN-GEiso \lambda_2',...
%         'BMSN-GEcos^0 \lambda_1','BMSN-GEcos^0 \lambda_2',...
%         'BMSN-GEcos^4 \lambda_1','BMSN-GEcos^4 \lambda_2',...
%         'BDiso \lambda_1','BDiso \lambda_2',...
%         'BDcos^0 \lambda_1','BDcos^0 \lambda_2',...  
%         'BDcos^4 \lambda_1','BDcos^4 \lambda_2',...       
%         'Location','southeast');
        
        
%         'BMSN-BFcos^0 \lambda_1','BMSN-BFcos^0 \lambda_2',...
%         'BMSN-BFcos^4 \lambda_1','BMSN-BFcos^4 \lambda_2',...
%         'BMSN-BFiso \lambda_1','BMSN-BFiso \lambda_2',...
% %         
        % 'BDcos^0 \lambda_1','BDcos^0 \lambda_2',...        
%          'BDcos^4 \lambda_1','BDcos^4 \lambda_2',...        
%          
end
title(strcat(target_snr));

% End