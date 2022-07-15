clear;%close;

% パラメータ
% パラメータ条件 NT >= NR*NU
%SN_tar  = 30;        % CDF表示のためのターゲットSNR [dB]
%SN_max = 40;         % 最大SNR[dB]
SIMU   = 100;       % 伝搬チャネル行列の発生回数
NT     = 24;          % 送信素子数
NR     = 3;          % 受信素子数(=2に固定)
NU     = 8;          % ユーザ数
I      = eye(NT,NT); % NTxNTの単位行列
SN_min = 0;
SN_max = 30;
SNR=SN_min:5:SN_max;
LSNR=length(SNR);

NRs = NR-1; % Antenna Selection
E_BDASms=zeros(LSNR,NRs);
E_BDms = zeros(LSNR,NR);
E_BMSN_GE_TDms = zeros(LSNR,NR);
E_BMSN_BFms = zeros(LSNR,NR);
E_BMSN_GEms = zeros(LSNR,NR);
% 擬似雑音
%a = NT/(10^(SN_tar/10));
%a2 = NR*NT/(10^(SN_tar/10));    % 擬似雑音 for BMSN2
        
% 所望のチャネル行列 for BMSN
for nuser = 1:NU
    T(:,:,nuser) = eye(NR,NR);
end

%if SN_tar < 10
 %   target_snr=strcat('SNR= ',num2str(SN_tar,'%01d'),'dB');
%else
%    target_snr=strcat('SNR=',num2str(SN_tar,'%02d'),'dB');
%end

% 伝搬チャネル行列
% i.i.d. Rayleigh
H = (randn(SIMU,NR*NU,NT)+1j*randn(SIMU,NR*NU,NT))/sqrt(2);

E_BDAS = zeros(SIMU, NRs, NU);           % BD-ASの固有値
E_BD = zeros(SIMU, NR, NU);             % BDの固有値
E_BMSN_GE_TD = zeros(SIMU, NR, NU);   % BMSN-GE-TDの固有値
E_BMSN_BF = zeros(SIMU, NR, NU);        % BMSN-BFの固有値
E_BMSN_GE = zeros(SIMU, NR, NU);        % BMSN-GEの固有値


for isnr=1:LSNR          % 試行回数のループ
    SN = SNR(isnr);
    sigma2 = 1/(10^(SN/10));
    a = sigma2*NT;
    for k = 1:SIMU 
    H0 = squeeze(H(k,:,:)); % k番目の試行回数での伝搬チャネル行列
    
    % BD-AS algorithm
    [W_BDAS,U_BDAS,S_BDAS,~,~] = bd_as(NT,NR,NRs,NU,H0); % function bd_as.m を使用
    
    % BD algorithm
    [W_BD,U_BD,S_BD,~,~] = bd(NT,NR,NU,H0); % function bd.m を使用
    
    % BMSN-GE-TD algorithm
    [W_BMSN_GE_TD,U_BMSN_GE_TD,S_BMSN_GE_TD,RIPGETD,~] = bmsn_getd_3(NT,NR,NU,H0,a); % function bmsn_getd_3.m を使用
    
    % BMSN-BF algorithm
    [W_BMSN_BF,U_BMSN_BF,S_BMSN_BF,RIPBF,~] = bmsn_bf(NT,NR,NU,H0,a,T); % function bmsn_bf.m を使用
      
    % BMSN-GE algorithm
    [W_BMSN_GE,U_BMSN_GE,S_BMSN_GE,RIPGE,~] = bmsn_gev(NT,NR,NU,H0,a); % function bmsn_gev.m を使用
    
    
    % ユーザ毎の固有値分布

    for nuser=1:NU
        if NR==1
            E_BDAS(k,:,nuser) = 10*log10(S_BDAS(1,1,nuser).^2/(NT*sigma2));
            E_BD(k,:,nuser) = 10*log10(S_BD(1,1,nuser).^2/(NT*sigma2));
            E_BMSN_GE_TD(k,:,nuser) = 10*log10((S_BMSN_GE_TD(1,1,nuser).^2)./(RIPGETD(1,nuser)+NT*sigma2));
            E_BMSN_BF(k,:,nuser) = 10*log10((S_BMSN_BF(1,1,nuser).^2)./(RIPBF(1,nuser)+NT*sigma2));
            E_BMSN_GE(k,:,nuser) = 10*log10((S_BMSN_GE(1,1,nuser).^2)./(RIPGE(1,nuser)+NT*sigma2));
        else 
            E_BDAS(k,:,nuser) = 10*log10(diag(S_BDAS(:,:,nuser)).^2/(NT*sigma2));   
            E_BD(k,:,nuser) = 10*log10(diag(S_BD(:,:,nuser)).^2/(NT*sigma2));
            E_BMSN_GE_TD(k,:,nuser) = 10*log10((diag(S_BMSN_GE_TD(:,:,nuser)).^2)./(RIPGETD(:,nuser)+NT*sigma2));
            E_BMSN_BF(k,:,nuser) = 10*log10((diag(S_BMSN_BF(:,:,nuser)).^2)./(RIPBF(:,nuser)+NT*sigma2));
            E_BMSN_GE(k,:,nuser) = 10*log10((diag(S_BMSN_GE(:,:,nuser)).^2)./(RIPGE(:,nuser)+NT*sigma2));
        end
    end

    end

    % ユーザ平均
E_BDASm = mean(E_BDAS,3);
E_BDm = mean(E_BD,3);
E_BMSN_GE_TDm = mean(E_BMSN_GE_TD,3);
E_BMSN_BFm = mean(E_BMSN_BF,3);
E_BMSN_GEm = mean(E_BMSN_GE,3);

E_BDASms(isnr,:) = mean(E_BDASm,1);
E_BDms(isnr,:) = mean(E_BDm,1);
E_BMSN_GE_TDms(isnr,:) = mean(E_BMSN_GE_TDm,1);
E_BMSN_BFms(isnr,:) = mean(E_BMSN_BFm,1);
E_BMSN_GEms(isnr,:) = mean(E_BMSN_GEm,1);

fprintf('SNR = %d dB\n',SNR(isnr));    
end





% CDF of Eigenvalue at Target SNR
rr2(:,1) = (1/SIMU:1/SIMU:1).'*100; 
Y = rr2;
rr3 = zeros(LSNR,NR*5);


for nn=1:NR
    rr3(:,nn) = E_BMSN_BFms(:,nn);
    rr3(:,nn+NR) = E_BMSN_GEms(:,nn);
    rr3(:,nn+2*NR) = E_BMSN_GE_TDms(:,nn);
    rr3(:,nn+3*NR) = E_BDms(:,nn);
end
    rr3(:,13) = E_BDASms(:,1);
    rr3(:,14) = E_BDASms(:,2);
 


%% CDF of EGV at Target SNR 
figure;
mycol = [1 0 1;1 0 1;0 0 0;0 0 0;
      1 0 0;1 0 0;0 0 1;0 0 1;
      1 0 0;1 0 0;
      0 0 0;0 0 0];
set(groot,'defaultAxesColorOrder',mycol)
Holizon_min = round(min(min(rr3(:,11))))-2;
Holizon_max = round(max(max(rr3)))+2;
axis([0 30 Holizon_min Holizon_max]);
grid on;
hold on;
set(gca,'XTick',Holizon_min-5:5:Holizon_max+5,'Fontsize',14,'Fontname','Arial')
xlabel('SINR of eigenvalue [dB]','Fontsize',16,'Fontname','Arial');
ylabel('CDF [%]','Fontsize',16,'Fontname','Arial');
plot(SNR,rr3(:,1),'r-d','Linewidth',2);
plot(SNR,rr3(:,2),'r--d','Linewidth',2);
plot(SNR,rr3(:,3),'r-.d','Linewidth',2);
plot(SNR,rr3(:,4),'b-s','Linewidth',2);
plot(SNR,rr3(:,5),'b--s','Linewidth',2);
plot(SNR,rr3(:,6),'b-.s','Linewidth',2);
plot(SNR,rr3(:,7),'g-o','Linewidth',2);
plot(SNR,rr3(:,8),'g--o','Linewidth',2);
% plot(rr3(:,9),Y,'g-.o',50:100:length(Y),'Linewidth',2);
plot(SNR,rr3(:,10),'c-x','Linewidth',2);
plot(SNR,rr3(:,11),'c--x','Linewidth',2);
plot(SNR,rr3(:,12),'c-.x','Linewidth',2);
plot(SNR,rr3(:,13),'k-<','Linewidth',2);
plot(SNR,rr3(:,14),'k--<','Linewidth',2);
%plot(rr3(:,6),rr3(:,4),'m');
legend('BMSN-BF \lambda_1','BMSN-BF \lambda_2','BMSN-BF \lambda_3',...
    'BMSN-GE \lambda_1','BMSN-GE \lambda_2','BMSN-GE \lambda_3',...
    'BMSN-GE-TD(3) \lambda_1','BMSN-GE-TD(3) \lambda_2',...
    'BD \lambda_1','BD \lambda_2','BD \lambda_3',...
    'BD-AS \lambda_1','BD-AS \lambda_2',...
    'Location','southeast');
%title(target_snr);



