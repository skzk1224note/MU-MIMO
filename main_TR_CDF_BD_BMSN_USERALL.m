clear;

% �p�����[�^���� NT >= NR*NU
SN_tar  = 10;        % CDF�\���̂��߂̃^�[�Q�b�gSNR [dB]

NT = 16;       % ���M�f�q��
NR = 2;        % ��M�f�q��/���[�U 
%NA = 1;        % �I�������M�A���e�i��
NU = 8;        % ���[�U��
Ntri = 1000;   % ���s��

Nru = NR*NU;

% ���]�̃`���l���s�� for BMSN
T = zeros(NR,NR,NU);
for nuser = 1:NU
    T(:,:,nuser) = eye(NR,NR);
end

if SN_tar < 10
    target_snr=strcat('SNR= ',num2str(SN_tar,'%01d'),'dB');
else
    target_snr=strcat('SNR=',num2str(SN_tar,'%02d'),'dB');
end

% �`���`���l���s��
% Rayleigh
H = (randn(Ntri,NR*NU,NT)+1j*randn(Ntri,NR*NU,NT))/sqrt(2);

E_BD   = zeros(Ntri, NR*NU);              % BD�̌ŗL�l
% E_BMSN1 = zeros(Ntri, NR*NU);             % BMSN1�̌ŗL�l
% E_BMSN2 = zeros(Ntri, NR*NU);             % BMSN2�̌ŗL�l
E_BMSNBF = zeros(Ntri, NR*NU);             % BMSN3�̌ŗL�l
E_BMSNGE = zeros(Ntri, NR*NU);             % BMSN3�̌ŗL�l
% E_MMSE  = zeros(Ntri, NR*NU);              % MMSE�̌ŗL�l

for k = 1:Ntri              % ���s�񐔂̃��[�v

    H0 = squeeze(H(k,:,:)); % k�Ԗڂ̎��s�񐔂ł̓`���`���l���s��
    
    % BD algorithm
    [W_BD,U_BD,S_BD,RIPBD] = bd(NT,NR,NU,H0); % function, bd.m ���g�p
    
    % BMSN-BF algorithm with adapted pseudo-noise power
    a=NT/(10^(SN_tar/10));
    [W_BMSNBF,U_BMSNBF,S_BMSNBF,RIPBF] = bmsn_bf(NT,NR,NU,H0,a,T); % function, bmsn.m ���g�p
    
    % BMSN-GE algorithm with adapted pseudo-noise power
    a=NT/(10^(SN_tar/10));
    [W_BMSNGE,U_BMSNGE,S_BMSNGE,RIPGE] = bmsn_ge(NT,NR,NU,H0,a); % function, bmsn.m ���g�p
    
%     % MMSE algorithm
%     a=NT/(10^(SN_tar/10));
%     [W_MMSE,U_MMSE,S_MMSE,RIPm] = gmmse_m2(NT,NR,NU,H0,a); % function, bmsn.m ���g�p
    
    % ���[�U���̌ŗL�l���z�iRIP���l���j
    snt = 1/(10^(SN_tar/10));
    for nuser=1:NU
        E_BD(k,1+(nuser-1)*NR:NR+(nuser-1)*NR) = 10*log10((diag(S_BD(:,:,nuser)).^2)./(RIPBD(:,nuser)+NT*snt));   
        E_BMSNBF(k,1+(nuser-1)*NR:NR+(nuser-1)*NR) = 10*log10((diag(S_BMSNBF(:,:,nuser)).^2)./(RIPBF(:,nuser)+NT*snt));
        E_BMSNGE(k,1+(nuser-1)*NR:NR+(nuser-1)*NR) = 10*log10((diag(S_BMSNGE(:,:,nuser)).^2)./(RIPGE(:,nuser)+NT*snt));
%         E_MMSE(k,1+(nuser-1)*NR:NR+(nuser-1)*NR) = 10*log10((diag(S_MMSE(:,:,nuser)).^2)./(RIPm(:,nuser)+NT*snt));
    end
    
end

% �ŗL�l�̓���
Eval_1 = E_BMSNBF;
Eval_2 = E_BMSNGE;
%Eval_3 = E_MMSE;
Eval_3 = E_BD;

% �`�����[�g: �ŗL�l����IEEE802.11ac�`�����[�g�ɕϊ�
% Out_1 = strcat('CSV_TR/TR16x2x8u_BMSNopt_SNR',num2str(SN_tar,'%02d'),'dB.csv');      % �o�̓t�@�C����
% Out_2 = strcat('CSV_TR/TR16x2x8u_BMSNGopt_SNR',num2str(SN_tar,'%02d'),'dB.csv');      % �o�̓t�@�C����
% Out_3 = strcat('CSV_TR/TR16x2x8u_MMSE_SNR',num2str(SN_tar,'%02d'),'dB.csv');      % �o�̓t�@�C����
% Out_4 = strcat('CSV_TR/TR16x2x8u_BD_SNR',num2str(SN_tar,'%02d'),'dB.csv');      % �o�̓t�@�C����

tr_table_40MHz; % �`�����[�g��SNR�̊֌W������

% �`�����[�g�v�Z
% BMSN-BF
TR_1 = Eval_to_TR (NR, NU, Ntri,Eval_1,SNRT,TRT);
%csvwrite(Out_1,TR_1);
% BMSN-GE
TR_2 = Eval_to_TR (NR, NU, Ntri,Eval_2,SNRT,TRT);
%csvwrite(Out_2,TR_2);
% MMSE
%TR_3 = Eval_to_TR (NR, NU, Ntri,Eval_3,SNRT,TRT);
%csvwrite(Out_3,TR_3);
% BD
TR_3 = Eval_to_TR (NR, NU, Ntri,Eval_3,SNRT,TRT);
%csvwrite(Out_4,TR_4);
% BD-AS
% TR_5 = Eval_to_TR (NA, NU, Ntri,Eval_5,SNRT,TRT);
% csvwrite(Out_5,TR_5);
% BMSN-AS
% TR_6 = Eval_to_TR (NA, NU, Ntri,Eval_6,SNRT,TRT);
% csvwrite(Out_6,TR_6);

%nuser=1;
% ��M�A���e�i�̘a
for nuser=1:NU
    Q(:,1,nuser) = sum(TR_1(:,(nuser-1)*NR+1:(nuser-1)*NR+NR),2);
    Q(:,2,nuser) = sum(TR_2(:,(nuser-1)*NR+1:(nuser-1)*NR+NR),2);
    Q(:,3,nuser) = sum(TR_3(:,(nuser-1)*NR+1:(nuser-1)*NR+NR),2);
 %   Q(:,4,nuser) = sum(TR_4(:,(nuser-1)*NR+1:(nuser-1)*NR+NR),2);
end

% ���[�U����
Qm(:,1) = mean(sort(Q(:,1,nuser),1),3);
Qm(:,2) = mean(sort(Q(:,2,nuser),1),3);
Qm(:,3) = mean(sort(Q(:,3,nuser),1),3);
%Qm(:,4) = mean(sort(Q(:,4,nuser),1),3);

% CDF of TR at Target SNR
Y(:,1) = (1/Ntri:1/Ntri:1).'*100; 

% �O���t�\��
%% CDF of EGV at Target SNR 
figure;
mycol = [1 0 0;
      0 0 1;
      0 0.7 0;
      1 0 1;
      0.8 0.6 0;0 0 0];
set(groot,'defaultAxesColorOrder',mycol)
%plot(Q(:,1),Y,'-m',Q(:,2),Y,'-r',Q(:,3),Y,'-b',Q(:,4),Y,'-k',Q(:,5),Y,'-g','Linewidth',2);
plot(Qm(:,1),Y,'-',Qm(:,2),Y,'-',Qm(:,3),Y,'-','Linewidth',2);
set(gca,'Fontsize',14,'Fontname','Times New Roman');
legend('BMSN-BF','BMSN-GE','BD','Location','Southeast');
xlabel('Average transmission rate [Mbps]','Fontsize',16,'Fontname','Times New Roman');
ylabel('CDF [%]','Fontsize',16,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
title(target_snr);
grid on;
hold on;





