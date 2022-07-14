% cap_bd_muser.m
% NU-���[�U�ł�BD�@�̌v�Z
% TDMA, Upper bound�Ɣ�r

clear;%close;
overBD_switch = 0;% 1:on,0:off overBD figure�̗L��

% �p�����[�^���� NT >= NR*NU
K_min = -20;      % �ŏ�K[dB]
K_max = 20;       % �ő�K[dB]

CDF = 10;         % ABR�̃^�[�Q�b�gCDF
SNR_tar = 10;     % �^�[�Q�b�gSNR[dB]
Nt = 16;          % ���M�f�q��
Nr = 2;           % ��M�f�q�� (�A���e�i�I���̏ꍇ1, �����łȂ��ꍇ��2)
Nu = 8;           % ���[�U��
SIMU = 10;      % ���s�񐔁i�ʏ� 1000)

An = 2; % �w�����֐��̌W��(cos�p�^�[���̏ꍇ��2)

I = eye(Nt,Nt);
Nru = Nr*Nu;

d_t = 0.5;   % ���M�A���e�i�Ԋu�iin wavelength)
d_r = 0.5;        % ��M�A���e�i�Ԋu�iin wavelength)
derad = pi/180;   % degree -> rad
K_box=(K_min:5:K_max).'; % figure�̉����̂��߂�K�̔�
LK=length(K_box);        % K�̔��̑傫��

%T ���]�̃`���l���s��
T = zeros(Nr,Nr);
for inu = 1:Nu
    T(:,:,inu) = eye(Nr,Nr);
end
Algorithms = ["BMSN-BF","BMSN-GE","MMSE-CI","GMI1","GMI2","BD","ZF-CI","BMSN-GE3"];
MSt = zeros(SIMU, Nr*Nu, numel(Algorithms)-1);
if Nr > 1
    MSt_stream1 = zeros(SIMU, (Nr-1)*Nu);
else
    MSt_stream1 = zeros(SIMU, Nr*Nu);
end
Q = zeros(SIMU, Nu, numel(Algorithms));
QmC = zeros(LK, numel(Algorithms));
QmCcos = zeros(LK, numel(Algorithms));
% ���̓t�@�C�����
% folder='CSV3/';
% fn1 = 'Eig16x2x8u_BD_SNR';
% fn2 = 'Eig16x2x8u_BMSN1(GEV)_SNR';
% fn3 = 'Eig16x2x8u_BMSN2a_SNR';  % alpha = 1e-2
% fn4 = 'Eig16x2x8u_BMSN2b_SNR';  % alpha = 1e-6
% fn5 = 'Eig16x2x8u_BMSN3_SNR';   % alpha = sigma2Nt equal to MMSE
% fn6 = 'Eig16x2x8u_ZF_SNR';
% fn7 = 'Eig16x2x8u_MMSE_SNR';

% ���̓t�@�C�����
% folder='CSV4/';
% fn1 = 'Eig16x2x8u_ZF_SNR';
% fn2 = 'Eig16x2x8u_BD_SNR';
% %fn3 = 'Eig16x2x8u_BMSN2_SNR'; % 2: alpha = 1e-2, 3: alpha=sigma2xNT
% %fn4 = 'Eig16x2x8u_BDAS_SNR';  % antenna selection
% %fn5 = 'Eig16x2x8u_BMSNAS2_SNR';  % alpha = 1e-2 and antenna selection
% fn3 = 'Eig16x2x8u_BMSN3_SNR';   % alpha = sigma2Nt equal to MMSE
% fn4 = 'Eig16x2x8u_BMSN3(GEV)_SNR';   % alpha = sigma2Nt equal to MMSE
% fn5 = 'Eig16x2x8u_MMSE_SNR';

% �o�̓t�@�C����
%folder= 'CSV/';
%cdfn1 = 'ABRCDFvsSNR_16x2x8u_BD_MMSE_BMSNs_USERave_CDF';
%cdfn2 = 'Eig_CDF_16x2x8u_BMSN_SNR';
%%
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

if d_t < 10
   target_d_t=strcat('d_t=',num2str(d_t,'%.3g'),'\lambda');
else
   target_d_t=strcat('d_t=',num2str(d_t,'%02d'),'\lambda');
end

if d_r < 10
   target_d_r=strcat('d_r=',num2str(d_r,'%.3g'),'\lambda');
else
   target_d_r=strcat('d_r=',num2str(d_r,'%02d'),'\lambda');
end
% �o�̓t�@�C���� with SNR in dB
%cdffile1 = strcat(folder,cdfn1,num2str(CDF,'%02d'),'_1000itr.csv');
%cdffile2 = strcat(folder,cdfn2,num2str(SN_tar,'%02d'),'dB_1000itr.csv');
%%
H_los = zeros(Nu*Nr,Nt); % �`���`���l���s��̒��ڔg����(LOS �`���l��)
sigma2 = 1/(10^(SNR_tar/10)); % noise power
a = sigma2*Nt; % �[���G��

St = zeros(Nr,Nu);
St_GE3 = zeros(Nr-1,Nu);
% MSt_BMSN_GE3 = zeros(SIMU, Nr*Nu);
% MSt_BMSN_BF = zeros(SIMU, Nr*Nu);
% MSt_BMSN_GE = zeros(SIMU, Nr*Nu);
% MSt_MMSE = zeros(SIMU, Nr*Nu);
% MSt_GMI1 = zeros(SIMU, Nr*Nu);
% MSt_GMI2 = zeros(SIMU, Nr*Nu);
% MSt_BD = zeros(SIMU, Nr*Nu);
% MSt_ZF = zeros(SIMU, Nr*Nu);
%%
tic;
for Directivity_switch = 0:1 % 1:on,0:off ����M�f�q�̎w�����l���̗L��
    for ik = 1:LK
    
    K_tar = K_box(ik);

% �o�̓t�@�C���� with alpha and SNR in dB
%evfile1 = strcat(folder,fn1,num2str(snr(isnr),'%02d'),'dB_1000itr.csv');
    
    for isimu = 1:SIMU
        
        % �`���`���l���s��̃}���`�p�X���� (i.i.d.Rayleigh)
        H_iid = (randn(Nr*Nu,Nt)+1j*randn(Nr*Nu,Nt))/sqrt(2);
        % LOS �`���l��
        
        if Directivity_switch == 1
            Theta_t = (rand(1,Nu)-0.5)*180; % ���[�U���̑��M�p �w����:(-90deg - 90deg)
            Theta_r = (rand(1,Nu)-0.5)*180; % ���[�U���̎�M�p �w����:(-90deg - 90deg)
            for n = 1 : Nu
                a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad))*An*cos(Theta_t(1,n)*derad); % ���[�U���̑��M���[�h�x�N�g��
                a_r = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad))*An*cos(Theta_r(1,n)*derad); % ���[�U���̎�M���[�h�x�N�g��
                H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t';                 % ���[�U����LOS�`���l���s��
            end
        else
            Theta_t = (rand(1,Nu)-0.5)*360; % ���[�U���̑��M�p ������:(-180deg - 180deg) 
            Theta_r = (rand(1,Nu)-0.5)*360; % ���[�U���̎�M�p ������:(-180deg - 180deg) 
            for n = 1 : Nu
            a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad)); % ���[�U���̑��M���[�h�x�N�g��
            a_r = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad)); % ���[�U���̎�M���[�h�x�N�g��
            H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t';                 % ���[�U����LOS�`���l���s��
            end
        end
    % K��^�l�ɂ���
    K = 10^(K_tar/10);
      
    % H=[sqrt(K/(K+1))*(LOS �`���l��)]+[sqrt(1/(K+1))*(NLOS �`���l��)]
    H = sqrt(K/(K+1))*H_los + sqrt(1/(K+1))*H_iid;
    
    %% BMSN-BF
    [~,~,STT,RIP,~] = bmsn_bf(Nt,Nr,Nu,H,a,T);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %��M�A���e�i�inr)�̌ŗL�l���z for B-MSN�i���g�������l���j
    % MSt_BMSN_BF(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    MSt(isimu,:,strcmp("BMSN-BF",Algorithms))=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
 
    %% BMSN-GE
    [~,~,STT,RIP,~] = bmsn_gev(Nt,Nr,Nu,H,a);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %��M�A���e�i�inr)�̌ŗL�l���z for B-MSN�i���g�������l���j
    % MSt_BMSN_GE(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    MSt(isimu,:,strcmp("BMSN-GE",Algorithms))=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
% Stge = St
%     RIP
%     reshape(St,[Nru,1])
%     reshape(RIP,[Nru,1])
%     ((reshape(St,[Nru,1]).').^2)
% (reshape(RIP,[Nru,1]).'+Nt*sigma2)
    %% MMSE-CI
    [~,~,STT,RIP,~] = mmse(Nt,Nr,Nu,H,a);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %��M�A���e�i�inr)�̌ŗL�l���z for B-MSN�i���g�������l���j
    % MSt_MMSE(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    MSt(isimu,:,strcmp("MMSE-CI",Algorithms))=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    %% GMI1
    [~,~,STT,RIP,~] = gmmse_m1(Nt,Nr,Nu,H,a);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %��M�A���e�i�inr)�̌ŗL�l���z for B-MSN�i���g�������l���j
    % MSt_GMI1(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    MSt(isimu,:,strcmp("GMI1",Algorithms))=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    %% GMI2
    [~,~,STT,RIP,~] = gmmse_m2(Nt,Nr,Nu,H,a);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %��M�A���e�i�inr)�̌ŗL�l���z for B-MSN�i���g�������l���j
    % MSt_GMI2(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    MSt(isimu,:,strcmp("GMI2",Algorithms))=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    %% BD
    [~,~,STT,RIP,~] = bd(Nt,Nr,Nu,H);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %��M�A���e�i�inr)�̌ŗL�l���z for B-MSN�i���g�������l���j
    % MSt_BD(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    MSt(isimu,:,strcmp("BD",Algorithms))=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);  
    %% ZF-CI
    [~,~,STT,RIP,~] = zf(Nt,Nr,Nu,H);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end    
    %��M�A���e�i�inr)�̌ŗL�l���z for B-MSN�i���g�������l���j
    % MSt_ZF(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    MSt(isimu,:,strcmp("ZF-CI",Algorithms))=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    %% BMSN-GE3
    [~,~,STT,RIP,~] = bmsn_ge3(Nt,Nr,Nu,H,a);
    for inu=1:Nu

% St_3 = St
%     RIP
%     ((St_GE3.').^2)
%     RIP(1,:).'
% RIP(1,:)+Nt*sigma2
    % fprintf('Iteration = %d / %d\n',isimu, SIMU);    
    end % isimu
    
    %% �e�A���S���Y����Channel Capacity
    for nuser=1:Nu
        ns = Nr*(nuser-1)+1:Nr*nuser;
        ns_GE3 = nuser; % BMSN-GE3�p�̃A���e�i��
        
        for n_alg = 1:numel(Algorithms)-1
            Q(:,nuser,n_alg)=sort(sum(log2(1 + MSt(:,ns,n_alg)),2)); % �e�A���S���Y����Channel capacity���v�Z
            % Q(:,nuser,strcmp(" ",Algorithms))=sort(sum(log2(1 + MSt(:,ns,strcmp(" ",Algorithms))),2)); % �e�A���S���Y����Channel capacity���v�Z
        end
        Q(:,nuser,strcmp("BMSN-GE3",Algorithms))=sort(sum(log2(1 + MSt_stream1(:,ns_GE3)),2)); % BMSN-GE3��Channel capacity���v�Z
    end % nuser end

    Qm = mean(Q,2); % Q�̃��[�U�񐔕���
    if Directivity_switch == 0
        QmC(ik,:)=Qm(round(CDF*SIMU/100),:); % Channel capacity
    else
        QmCcos(ik,:)=Qm(round(CDF*SIMU/100),:); % cos Channel capacity
    end
    
    QmC_overBD=QmC(:,:)./QmC(:,strcmp("BD",Algorithms));% C / C_BD

    fprintf('K = %d dB\n',K_box(ik));
    end % ik end
    
    %csvwrite(cdffile1,[K,QmC]);
end % Directivity end
toc;
%% �O���t�\�� figure ABRvsK
figure;
mycol = [1 0 0;0 0 1;0 1 0;0 1 1;1 0 1;0 0 0;1 1 0;0.85 0.325 0.098
         1 0 0;0 0 1;0 1 0;0 1 1;1 0 1;0 0 0;1 1 0;0.85 0.325 0.098]; % �F
set(groot,'defaultAxesColorOrder',mycol)
for n_alg = 1:numel(Algorithms)
    plot(K_box,QmC(:,n_alg),'--','Linewidth',2);
    grid on; hold on;
end
for n_alg = 1:numel(Algorithms)
    plot(K_box,QmCcos(:,n_alg),'-','Linewidth',2);
    grid on; hold on;
end
%plot(K_box,QmC(:,1),'r-',K_box,QmC(:,2),'b-',K_box,QmC(:,3),'g-',K_box,QmC(:,4),'y-',K_box,QmC(:,5),'k-',K_box,QmC(:,6),'c-',K_box,QmC(:,7),'m-','Linewidth',2);

axis([K_min,K_max,0,max(max([QmC QmCcos]))+2])
set(gca,'Fontsize',18,'Fontname','Times New Roman');
lgd = legend;
lgd.NumColumns = numel(Algorithms)/4; % �}��̗񐔂��w��
legend([Algorithms strcat(Algorithms,"cos")],'Location','Northwest');

xlabel('{\it{K}} [dB]','Fontsize',6,'Fontname','Times New Roman');
ylabel('Channel capacity [bits/s/Hz]','Fontsize',6,'Fontname','Times New Roman');
set(gca,'Fontsize',18,'Fontname','Times New Roman');
title(strcat(target_CDF,',',target_SNR,',',target_d_t,',',target_d_r));
grid on; hold on;
%% �O���t�\�� figure Channel Capacity overBD vs K
if overBD_switch == 1
    figure;
    mycol = [1 0 1;0 1 0;1 0 0;0 0 1;0 1 1;0 0 0]; % �F
    set(groot,'defaultAxesColorOrder',mycol)
    %plot(K,QmC(:,1),K,QmC(:,2),K,QmC(:,3),K,QmC(:,4),'Linewidth',2);
    
    plot(K_box,QmC_overBD(:,1),'r-',K_box,QmC_overBD(:,2),'b-',K_box,QmC_overBD(:,3),'g-',K_box,QmC_overBD(:,4),'y-',K_box,QmC_overBD(:,5),'k-',K_box,QmC_overBD(:,7),'m-','Linewidth',2);
    
    axis([K_min,K_max,0,5])
    set(gca,'Fontsize',18,'Fontname','Times New Roman');
    lgd = legend;
    lgd.NumColumns = 3; % �}��̗񐔂��w��
    legend(Algorithms,'Location','Northwest');
    xlabel('{\it{K}} [dB]','Fontsize',6,'Fontname','Times New Roman');
    ylabel('C / C_{BD}','Fontsize',6,'Fontname','Times New Roman');
    set(gca,'Fontsize',18,'Fontname','Times New Roman');
    title(strcat(target_CDF,',',target_SNR,',',target_d_t,',',target_d_r));
    grid on; hold on;
end
% End