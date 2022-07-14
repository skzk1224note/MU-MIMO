% SINR vs. SNR for targer CDF

clear;
Nt = 16;     % ���M�f�q�� (16)
Nr = 2;     % �e���[�U�̎�M�f�q�� (2) (�A���e�i�I���̏ꍇ1, �����łȂ��ꍇ��2)
Nu = 8;     % ���[�U���i8�j

% SNR_tar = 30;  % �^�[�Q�b�gSNR[dB]
SNR_min = 0;        % �ŏ�SNR [dB]
SNR_max = 30;       % �ő�SNR [dB]
SIMU = 1000; %���s�񐔁i�ʏ� 1000�j
CDF = 50;   % SINR��targer CDF

% if SNR_tar < 10
%     tsnr=strcat('SNR= ',num2str(SNR_tar,'%01d'),'dB');
% else
%     tsnr=strcat('SNR=',num2str(SNR_tar,'%02d'),'dB');
% end

I = eye(Nt,Nt);
Nru = Nr*Nu;

%T ���]�̃`���l���s��
for inu = 1:Nu
    T(:,:,inu) = zeros(Nr,Nr);
    T(1,1,inu) = 1;
    T(1,2,inu) = 0;
    T(2,1,inu) = 0;
    T(2,2,inu) = 1;
end

% �o�̓t�@�C����
folder= 'CSV/';
cdfn1 = 'SINRCDFvsSNR_16x2x8u_BD_MMSE_BMSNs_USERave_CDF';

if CDF < 10
    target_CDF=strcat('CDF= ',num2str(CDF,'%01d'),' %');
else
    target_CDF=strcat('CDF= ',num2str(CDF,'%02d'),' %');
end

% �o�̓t�@�C���� with SNR in dB
cdffile1 = strcat(folder,cdfn1,num2str(CDF,'%02d'),'_1000itr.csv');

SNR=(SNR_min:5:SNR_max).';
LSNR=length(SNR);

for isnr = 1:LSNR
    
    SNR_tar = SNR(isnr);

    sigma2 = 1/(10^(SNR_tar/10)); % noise power
    a = sigma2*Nt;
    
for isimu = 1:SIMU

    %H(�`���`���l���s��: iid Rayleigh channel)
    %Hu(���[�Uk�̃`���l���s��Hu(:,:,k))
    H = (randn(Nr*Nu,Nt) + 1j*randn(Nr*Nu,Nt))/sqrt(2);
    Hu = peruser(H,Nu);

    He = zeros((Nu-1)*Nr,Nt,Nu);    % H����1���[�U�̃`���l���s����������s��

    % BMSN_BF
    [~,~,~,RIP,SP] = bmsn_bf(Nt,Nr,Nu,H,a,T);
    MSt_BMSN_BF(isimu,:)=(reshape(SP,[Nru,1]).')./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    
    % BMSN_GE
    [~,~,~,RIP,SP] = bmsn_gev(Nt,Nr,Nu,H,a);
    MSt_BMSN_GE(isimu,:)=(reshape(SP,[Nru,1]).')./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
 
    % MMSE-CI
    [~,~,STT,RIP,SP] = mmse(Nt,Nr,Nu,H,a);
    MSt_MMSE(isimu,:)=(reshape(SP,[Nru,1]).')./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    
    % BD
    [~,~,~,RIP,SP] = bd(Nt,Nr,Nu,H);
    MSt_BD(isimu,:)=(reshape(SP,[Nru,1]).')./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    
    % ZF-CI
    [~,~,~,RIP,SP] = zf(Nt,Nr,Nu,H);
    MSt_ZF(isimu,:)=(reshape(SP,[Nru,1]).')./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
        
    %fprintf('Iteration = %d / %d\n',isimu, SIMU);
    
end % isimu

    % SINR
    for nuser=1:Nu
        Q(:,1,nuser)=sort(mean(MSt_BMSN_BF(:,(nuser-1)*Nr+1:(nuser-1)*Nr+Nr),2));
        Q(:,2,nuser)=sort(mean(MSt_BMSN_GE(:,(nuser-1)*Nr+1:(nuser-1)*Nr+Nr),2));
        Q(:,3,nuser)=sort(mean(MSt_MMSE(:,(nuser-1)*Nr+1:(nuser-1)*Nr+Nr),2));
        Q(:,4,nuser)=sort(mean(MSt_BD(:,(nuser-1)*Nr+1:(nuser-1)*Nr+Nr),2));
        Q(:,5,nuser)=sort(mean(MSt_ZF(:,(nuser-1)*Nr+1:(nuser-1)*Nr+Nr),2));
    end

    Qm = 10*log10(mean(Q,3));

    QmC(isnr,:)=Qm(round(CDF*10),:);

    fprintf('SNR = %d dB\n',SNR(isnr));
    
end % isnr

%csvwrite(cdffile1,[SNR,QmC]);

%Y(:,1) = [0.1:0.1:100].';
%Y(:,1) = [0.01:0.01:100].';

figure;
mycol = [0.7 0.7 0;
      0 1 1;
      0 0.7 0;
      1 0 1;
      0 0 0;
      1 0 0;
      0 0 1];
set(groot,'defaultAxesColorOrder',mycol)
%plot(SNR,QmC(:,7),'y-*',SNR,QmC(:,1),'c-o',SNR,QmC(:,2),'g-s',SNR,QmC(:,3),'m-^',SNR,QmC(:,4),'k-v',SNR,QmC(:,5),'r-x',SNR,QmC(:,6),'b-s','Linewidth',2);
plot(SNR,QmC(:,1),'-d',SNR,QmC(:,2),'-o',SNR,QmC(:,3),'-s',SNR,QmC(:,4),'-^',SNR,QmC(:,5),'-v','Linewidth',2)
%    'MarkerIndices',10:100:length(Y),'Linewidth',2);
%axis([0,15,0,100])
set(gca,'Fontsize',14,'Fontname','Times New Roman');
legend('BMSN-BF','BMSN-GE','MMSE-CI','BD','ZF-CI','Location','Southeast');
xlabel('SNR [dB]','Fontsize',16,'Fontname','Times New Roman');
ylabel('Per-user received SINR [dB]','Fontsize',16,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
title(target_CDF);
grid on;
hold on;

%End


