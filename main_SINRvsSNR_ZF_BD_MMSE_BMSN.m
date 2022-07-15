% SINRvsSNR_BD_ZF_GZF_MMSE.m
% Outout INR & Output SNR & Output SINR の計算
% BD, ZF, GZF, and MMSEを比較

clear;%close;

% 出力ファイル名
testfile1 = 'CSV/OINR_16x2x8_BD_2BMSNs_MMSE.csv';
testfile2 = 'CSV/OSNR_16x2x8_BD_2BMSNs_MMSE.csv';
testfile3 = 'CSV/OSINR_16x2x8_BD_2BMSNs_MMSE.csv';
% testfile2 = 'Capacity2x8x4u_BD_CDFSNR20dB.csv';
% testfile3 = 'Eig2x8x4u_BD_CDFSNR20dB.csv';

% パラメータ
% パラメータ条件 NT >= NR*NU
%SN_tar  = 30;        % CDF表示のためのターゲットSNR [dB]
SN_min = 0;
SN_max = 30;         % 最大SNR[dB]
Ntri   = 100;       % 伝搬チャネル行列の発生回数 (通常は1000)
NT     = 24;          % 送信素子数
NR     = 3;          % 受信素子数
NU     = 8;          % ユーザ数 (=2に固定)
I      = eye(NT,NT); % NTxNTの単位行列

% a1 = 1e-2;            % 擬似雑音 for BMSN1
% a2 = 1e-6;            % 擬似雑音 for BMSN2

% 所望のチャネル行列 for BMSN
for nuser = 1:NU
    T(:,:,nuser) = eye(NR,NR);
end

%target_snr=strcat('SNR=',num2str(SN_tar,'%02d'),'dB');

% 伝搬チャネル行列
% Rayleigh
H = (randn(Ntri,NR*NU,NT)+1j*randn(Ntri,NR*NU,NT))/sqrt(2);

SNRs=SN_min:SN_max;
LSNR=length(SNRs);

SP_BD_USER = zeros(1,NU);
SP_ZF_USER = zeros(1,NU);
SP_BMSN_BF_USER = zeros(1,NU);
SP_BMSN_GE_USER = zeros(1,NU);
SP_MMSE_USER = zeros(1,NU);

SP_BDs = zeros(Ntri,LSNR);     
SP_ZFs = zeros(Ntri,LSNR); 
SP_BMSN_BFs = zeros(Ntri,LSNR); 
SP_BMSN_GEs = zeros(Ntri,LSNR);        
SP_MMSEs = zeros(Ntri,LSNR);

RIP_BD_USER = zeros(1,NU);
RIP_ZF_USER = zeros(1,NU);
RIP_BMSN_BF_USER = zeros(1,NU);
RIP_BMSN_GE_USER = zeros(1,NU);
RIP_MMSE_USER = zeros(1,NU);

RIP_BDs = zeros(Ntri,LSNR);     
RIP_ZFs = zeros(Ntri,LSNR); 
RIP_BMSN_BFs = zeros(Ntri,LSNR); 
RIP_BMSN_GEs = zeros(Ntri,LSNR);        
RIP_MMSEs = zeros(Ntri,LSNR);

for isnr=1:LSNR
    SN = SNRs(isnr);
    sigma2 = 1/(10^(SN/10));
    a = sigma2*NT;
    for k = 1:Ntri              % 試行回数のループ

        H0 = squeeze(H(k,:,:)); % k番目の試行回数での伝搬チャネル行列
    
        % BD algorithm
        [W_BD,U_BD,S_BD,RIP_BD,SP_BD] = bd(NT,NR,NU,H0); % function, bd.m を使用
    
        % ZF algorithm
        [W_ZF,U_ZF,S_ZF,RIP_ZF,SP_ZF] = zf(NT,NR,NU,H0); % function, zf.m を使用
    
        % BMSN-BF algorithm
        [W_BMSN_BF,U_BMSN_BF,S_BMSN_BF,RIP_BMSN_BF,SP_BMSN_BF] = bmsn_bf(NT,NR,NU,H0,a,T); % function, bmsn_bf.m を使用
    
        % BMSN_GE algorithm
        [W_BMSN_GE,U_BMSN_GE,S_BMSN_GE,RIP_BMSN_GE,SP_BMSN_GE] = bmsn_gev(NT,NR,NU,H0,a); % function, bmsn_gev.m を使用

        % MMSE algorithm
        [W_MMSE,U_MMSE,S_MMSE,RIP_MMSE,SP_MMSE] = mmse(NT,NR,NU,H0,a); % function, mmse.m を使用
        
        for nuser = 1 : NU
            SP_BD_USER(nuser) = sum(SP_BD(:,nuser))/(NT*sigma2);
            SP_ZF_USER(nuser) = sum(SP_ZF(:,nuser))/(NT*sigma2);
            SP_BMSN_BF_USER(nuser) = sum(SP_BMSN_BF(:,nuser))/(NT*sigma2);
            SP_BMSN_GE_USER(nuser) = sum(SP_BMSN_GE(:,nuser))/(NT*sigma2);
            SP_MMSE_USER(nuser) = sum(SP_MMSE(:,nuser))/(NT*sigma2);
            RIP_BD_USER(nuser) = sum(RIP_BD(:,nuser))/(NT*sigma2);
            RIP_ZF_USER(nuser) = sum(RIP_ZF(:,nuser))/(NT*sigma2);
            RIP_BMSN_BF_USER(nuser) = sum(RIP_BMSN_BF(:,nuser))/(NT*sigma2);
            RIP_BMSN_GE_USER(nuser) = sum(RIP_BMSN_GE(:,nuser))/(NT*sigma2);
            RIP_MMSE_USER(nuser) = sum(RIP_MMSE(:,nuser))/(NT*sigma2);
         end
        SP_BDs(k,isnr) = mean(SP_BD_USER);
        SP_ZFs(k,isnr) = mean(SP_ZF_USER);
        SP_BMSN_BFs(k,isnr) = mean(SP_BMSN_BF_USER);
        SP_BMSN_GEs(k,isnr) = mean(SP_BMSN_GE_USER);
        SP_MMSEs(k,isnr) = mean(SP_MMSE_USER);
        RIP_BDs(k,isnr)   = mean(RIP_BD_USER);      % User average
        RIP_ZFs(k,isnr) = mean(RIP_ZF_USER);  % User average
        RIP_BMSN_BFs(k,isnr) = mean(RIP_BMSN_BF_USER);  % User average
        RIP_BMSN_GEs(k,isnr) = mean(RIP_BMSN_GE_USER);  % User average     
        RIP_MMSEs(k,isnr) = mean(RIP_MMSE_USER);  % User average  
    
    end
    fprintf('SNR = %d dB\n',SN);
end

rr(1,:) = 10*log10(mean(RIP_BMSN_BFs,1).');    % Trial average
rr(2,:) = 10*log10(mean(RIP_BMSN_GEs,1).');    % Trial average
rr(3,:) = 10*log10(mean(RIP_MMSEs,1).');     % Trial average
rr(4,:) = 10*log10(mean(RIP_BDs,1).');       % Trial average
rr(5,:) = 10*log10(mean(RIP_ZFs,1).');    % Trial average
%
rr(6,:) = 10*log10(mean(SP_BMSN_BFs,1).');     % Trial average
rr(7,:) = 10*log10(mean(SP_BMSN_GEs,1).');      % Trial average
rr(8,:) = 10*log10(mean(SP_MMSEs,1).');     % Trial average
rr(9,:) = 10*log10(mean(SP_BDs,1).');        % Trial average
rr(10,:) = 10*log10(mean(SP_ZFs,1).');     % Trial average
%
rr(11,:) = 10*log10(mean(SP_BMSN_BFs,1).'./(mean(RIP_BMSN_BFs,1).'+1));     % Trial average
rr(12,:) = 10*log10(mean(SP_BMSN_GEs,1).'./(mean(RIP_BMSN_GEs,1).'+1));     % Trial average
rr(13,:) = 10*log10(mean(SP_MMSEs,1).'./(mean(RIP_MMSEs,1).'+1));     % Trial average
rr(14,:) = 10*log10(mean(SP_BDs,1).'./(mean(RIP_BDs,1).'+1));     % Trial average
rr(15,:) = 10*log10(mean(SP_ZFs,1).'./(mean(RIP_ZFs,1).'+1));     % Trial average

% csvwrite(testfile1,[SNRs;rr(1:4,:)]);
% csvwrite(testfile2,[SNRs;rr(5:8,:)]);
% csvwrite(testfile3,[SNRs;rr(9:12,:)]);

%% RIP vs SNR 
figure;
% mycol = [0 0 1;
%       0 0 0;
%       0 0.7 0;
%       1 0 1;0.7 0.7 0;
%       1 0 0];
% set(groot,'defaultAxesColorOrder',mycol)
plot(SNRs,rr(1,:),'r-+',SNRs,rr(2,:),'b-o',SNRs,rr(3,:),'g-^',SNRs,rr(4,:),'c-x',...
    SNRs,rr(5,:),'m-v','MarkerIndices',1:5:length(SNRs),'Linewidth',2);
%axis([SN_min SN_max -400 30]);
set(gca,'XTick',SN_min:10:SN_max,'Fontsize',14,'Fontname','Times New Roman')
xlabel('SNR [dB]','Fontsize',16,'Fontname','Times New Roman');
ylabel('Average received INR [dB]','Fontsize',15,'Fontname','Times New Roman');
legend('BMSN-BF','BMSN-GE','MMSE-CI','BD','ZF-CI','Location','Southwest');
grid on;
hold on;

figure;
plot(SNRs,rr(6,:),'r-+',SNRs,rr(7,:),'b-o',SNRs,rr(8,:),'g-^',SNRs,rr(9,:),'c-x',...
    SNRs,rr(10,:),'m-v','MarkerIndices',1:5:length(SNRs),'Linewidth',2);
set(gca,'XTick',SN_min:10:SN_max,'Fontsize',14,'Fontname','Times New Roman')
xlabel('SNR [dB]','Fontsize',16,'Fontname','Times New Roman');
ylabel('Average received SNR [dB]','Fontsize',15,'Fontname','Times New Roman');
legend('BMSN-BF','BMSN-GE','MMSE-CI','BD','ZF-CI','Location','northwest');
grid on;
hold on;

figure;
plot(SNRs,rr(11,:),'r-+',SNRs,rr(12,:),'b-o',SNRs,rr(13,:),'g-^',SNRs,rr(14,:),'c-x',...
    SNRs,rr(15,:),'m-v','MarkerIndices',1:5:length(SNRs),'Linewidth',2);
%axis([SN_min SN_max -400 30]);
set(gca,'XTick',SN_min:10:SN_max,'Fontsize',14,'Fontname','Times New Roman')
xlabel('SNR [dB]','Fontsize',16,'Fontname','Times New Roman');
ylabel('Average received SINR [dB]','Fontsize',15,'Fontname','Times New Roman');
legend('BMSN-BF','BMSN-GE','MMSE-CI','BD','ZF-CI','Location','northwest');
grid on;
hold on;




