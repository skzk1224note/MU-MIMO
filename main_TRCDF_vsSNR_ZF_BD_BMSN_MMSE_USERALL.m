clear;

% �p�����[�^���� NT >= NR*NU
SNR_min = 5;        % �ŏ�SNR [dB]
SNR_max = 30;       % �ő�SNR [dB]
%SN_tar  = 25;        % CDF�\���̂��߂̃^�[�Q�b�gSNR [dB]
CDF = 10;   % TR��CDF

NT = 16;       % ���M�f�q��
NR = 2;        % ��M�f�q��/���[�U 
%NA = 1;        % �I�������M�A���e�i��
NU = 8;        % ���[�U��
Ntri = 1000;   % ���s��

Nru = NR*NU;

% ���]�̃`���l���s�� for BMSN
for nuser = 1:NU
    T(:,:,nuser) = eye(NR,NR);
end

SNR=(SNR_min:5:SNR_max).';
LSNR=length(SNR);

if CDF < 10
    target_CDF=strcat('CDF= ',num2str(CDF,'%01d'),' %');
else
    target_CDF=strcat('CDF= ',num2str(CDF,'%02d'),' %');
end

for isnr = 1:LSNR

    SN_tar = SNR(isnr);
    
    % �`���`���l���s��
    % Rayleigh
    H = (randn(Ntri,NR*NU,NT)+1j*randn(Ntri,NR*NU,NT))/sqrt(2);

    E_BD = zeros(Ntri, NR*NU);              % BD�̌ŗL�l
    E_ZF = zeros(Ntri, NR*NU);              % ZF-CI�̌ŗL�l
    E_BMSN_BF = zeros(Ntri, NR*NU);         % BMSN-BF�̌ŗL�l
    E_BMSN_GE = zeros(Ntri, NR*NU);         % BMSN-GE�̌ŗL�l
    E_MMSE  = zeros(Ntri, NR*NU);           % MMSE-CI�̌ŗL�l

    for k = 1:Ntri              % ���s�񐔂̃��[�v

        H0 = squeeze(H(k,:,:)); % k�Ԗڂ̎��s�񐔂ł̓`���`���l���s��
    
        % BD algorithm
        [W_BD,U_BD,S_BD,~,~] = bd(NT,NR,NU,H0); % function, bd.m ���g�p
        
        % ZF-CI algorithm
        [W_ZF,U_ZF,S_ZF,~,~] = zf(NT,NR,NU,H0); % function, zf.m ���g�p
 
        % BMSN-BF algorithm with adapted pseudo-noise power
        a=NT/(10^(SN_tar/10));
        [W_BMSN_BF,U_BMSN_BF,S_BMSN_BF,RIP_BF,~] = bmsn_bf(NT,NR,NU,H0,a,T); % function, bmsn_bf.m ���g�p
    
        % BMSN-GE algorithm with adapted pseudo-noise power
        a=NT/(10^(SN_tar/10));
        [W_BMSN_GE,U_BMSN_GE,S_BMSN_GE,RIP_GE] = bmsn_gev(NT,NR,NU,H0,a); % function, bmsn_gev.m ���g�p
    
        % MMSE-CI algorithm
        a=NT/(10^(SN_tar/10));
        [W_MMSE,U_MMSE,S_MMSE,RIPm] = gmmse_m2(NT,NR,NU,H0,a); % function, mmse.m ���g�p
    
        % ���[�U���̌ŗL�l���z�iRIP���l���j
        snt = 1/(10^(SN_tar/10));
        for nuser=1:NU
            E_BD(k,1+(nuser-1)*NR:NR+(nuser-1)*NR) = 10*log10(diag(S_BD(:,:,nuser)).^2/(NT*snt));  
            E_ZF(k,1+(nuser-1)*NR:NR+(nuser-1)*NR) = 10*log10(diag(S_ZF(:,:,nuser)).^2/(NT*snt));  
            E_BMSN_BF(k,1+(nuser-1)*NR:NR+(nuser-1)*NR) = 10*log10((diag(S_BMSN_BF(:,:,nuser)).^2)./(RIP_BF(:,nuser)+NT*snt));
            E_BMSN_GE(k,1+(nuser-1)*NR:NR+(nuser-1)*NR) = 10*log10((diag(S_BMSN_GE(:,:,nuser)).^2)./(RIP_GE(:,nuser)+NT*snt));
            E_MMSE(k,1+(nuser-1)*NR:NR+(nuser-1)*NR) = 10*log10((diag(S_MMSE(:,:,nuser)).^2)./(RIPm(:,nuser)+NT*snt));
        end
    
    end

    % �ŗL�l�̓���
    Eval_1 = E_BMSN_BF;
    Eval_2 = E_BMSN_GE;
    Eval_3 = E_MMSE;
    Eval_4 = E_BD;
    Eval_5 = E_ZF;

    % �`�����[�g: �ŗL�l����IEEE802.11ac�`�����[�g�ɕϊ�
%     Out_1 = strcat('CSV_TR/TR16x2x8u_BMSNopt_SNR',num2str(SN_tar,'%02d'),'dB.csv');      % �o�̓t�@�C����
%     Out_2 = strcat('CSV_TR/TR16x2x8u_BMSNGopt_SNR',num2str(SN_tar,'%02d'),'dB.csv');      % �o�̓t�@�C����
%     Out_3 = strcat('CSV_TR/TR16x2x8u_MMSE_SNR',num2str(SN_tar,'%02d'),'dB.csv');      % �o�̓t�@�C����
%     Out_4 = strcat('CSV_TR/TR16x2x8u_BD_SNR',num2str(SN_tar,'%02d'),'dB.csv');      % �o�̓t�@�C����

    tr_table; % �`�����[�g��SNR�̊֌W������

    % �`�����[�g�v�Z
    % BMSN-BF
    TR_1 = Eval_to_TR (NR, NU, Ntri,Eval_1,SNRT,TRT);
    %csvwrite(Out_1,TR_1);
    % BMSN-GE
    TR_2 = Eval_to_TR (NR, NU, Ntri,Eval_2,SNRT,TRT);
    %csvwrite(Out_2,TR_2);
    % MMSE-CI
    TR_3 = Eval_to_TR (NR, NU, Ntri,Eval_3,SNRT,TRT);
    %csvwrite(Out_3,TR_3);
    % BD
    TR_4 = Eval_to_TR (NR, NU, Ntri,Eval_4,SNRT,TRT);
    %csvwrite(Out_4,TR_4);
    % ZF-CI
    TR_5 = Eval_to_TR (NR, NU, Ntri,Eval_5,SNRT,TRT);
    % csvwrite(Out_5,TR_5);

    % ��M�A���e�i�̘a
    for nuser=1:NU
        Q(:,1,nuser) = sum(TR_1(:,(nuser-1)*NR+1:(nuser-1)*NR+NR),2);
        Q(:,2,nuser) = sum(TR_2(:,(nuser-1)*NR+1:(nuser-1)*NR+NR),2);
        Q(:,3,nuser) = sum(TR_3(:,(nuser-1)*NR+1:(nuser-1)*NR+NR),2);
        Q(:,4,nuser) = sum(TR_4(:,(nuser-1)*NR+1:(nuser-1)*NR+NR),2);
        Q(:,5,nuser) = sum(TR_5(:,(nuser-1)*NR+1:(nuser-1)*NR+NR),2);
    end

    % ���[�U����
    Qm(:,1) = mean(sort(Q(:,1,nuser),1),3);
    Qm(:,2) = mean(sort(Q(:,2,nuser),1),3);
    Qm(:,3) = mean(sort(Q(:,3,nuser),1),3);
    Qm(:,4) = mean(sort(Q(:,4,nuser),1),3);
    Qm(:,5) = mean(sort(Q(:,5,nuser),1),3);

    QmC(isnr,:)=Qm(round(CDF*10),:);

    fprintf('SNR = %d dB\n',SNR(isnr));
    
end

% �O���t�\��
%% CDF of EGV at Target SNR 
figure;
% mycol = [1 0 0;
%       0 0 1;
%       0 0.7 0;
%       1 0 1;
%       1 1 0;0 0 0];
% set(groot,'defaultAxesColorOrder',mycol)
plot(SNR,QmC(:,1),'r-o',SNR,QmC(:,2),'b-v',SNR,QmC(:,3),'g-^',SNR,QmC(:,4),'c-s',SNR,QmC(:,5),'m-d','Linewidth',2);
%plot(SNR,QmC,'Linewidth',2);
set(gca,'Fontsize',14,'Fontname','Times New Roman');
legend('BMSN-BF','BMSN-GE','MMSE-CI','BD','ZF-CI','Location','Northwest');
xlabel('SNR [dB]','Fontsize',16,'Fontname','Times New Roman');
ylabel('Average transmission rate [Mbps]','Fontsize',16,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
title(target_CDF);
grid on;
hold on;





