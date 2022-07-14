% K_lambda1,lambda2

clear;

% K[dB]�͈̔͐ݒ�
K_min = -20;        % �ŏ�K [dB]
K_max = 20;         % �ő�K [dB]

CDF = 10;   % ����CDF

Nt = 16;    % ���M�f�q��
Nr = 2;     % ��M�f�q�� (�A���e�i�I���̏ꍇ1, �����łȂ��ꍇ��2)
Nu = 8;     % ���[�U��

SIMU = 1000; %���s�񐔁i�ʏ� 1000�j

I = eye(Nt,Nt); % Nt*Nt�̒P�ʍs��
Nru = Nr*Nu;

SNR_tar = 10;   % �^�[�Q�b�gSNR[dB]

% SNR_max = 30; %�ő�SNR[dB]
d_t      = 0.5;      % ���M�A���e�i�Ԋu�iin wavelength)
d_r      = 0.5;      % ��M�A���e�i�Ԋu�iin wavelength)
derad = pi/180;      % degree -> rad

%T ���]�̃`���l���s��
T = zeros(Nr,Nr);
for inu = 1:Nu
    T(:,:,inu) = eye(Nr,Nr);
end

% �G���Ƌ[���G��
snt = 1/(10^(SNR_tar/10));
a = Nt*snt;

if CDF < 10
    target_CDF=strcat('CDF= ',num2str(CDF,'%01d'),'%');
else
    target_CDF=strcat('CDF= ',num2str(CDF,'%02d'),'%');
end

if SNR_tar < 10
    target_SNR=strcat('SNR= ',num2str(SNR_tar,'%01d'),'dB');
else
    target_SNR=strcat('SNR= ',num2str(SNR_tar,'%02d'),'dB');
end

K_box=(K_min:5:K_max).'; % figure�̉����̂��߂�K�̔�
LK=length(K_box);        % K�̔��̑傫��

sigma2 = 1/(10^(SNR_tar/10)); % noise power
a = sigma2*Nt;

% �o�̓t�@�C���� with alpha and SNR in dB
%evfile1 = strcat(folder,fn1,num2str(snr(isnr),'%02d'),'dB_1000itr.csv');

E_ZF = zeros(SIMU, Nr, Nu);             % ZF-CI�̌ŗL�l��SINR
E_BD = zeros(SIMU, Nr, Nu);             % BD�̌ŗL�l��SINR
E_MMSE = zeros(SIMU, Nr, Nu);           % MMSE-CI�̌ŗL�l��SINR
E_GMI1 = zeros(SIMU, Nr, Nu);           % GMI1�̌ŗL�l��SINR
E_GMI2 = zeros(SIMU, Nr, Nu);           % GMI2�̌ŗL�l��SINR
E_BMSN_BF = zeros(SIMU, Nr, Nu);        % BMSN-BF�̌ŗL�l��SINR
E_BMSN_GE = zeros(SIMU, Nr, Nu);        % BMSN-GE�̌ŗL�l��SINR

Eigs_ZF = zeros(SIMU, Nr, Nu);             % ZF-CI�̌ŗL�l
Eigs_BD = zeros(SIMU, Nr, Nu);             % BD�̌ŗL�l
Eigs_MMSE = zeros(SIMU, Nr, Nu);           % MMSE-CI�̌ŗL�l
Eigs_GMI1 = zeros(SIMU, Nr, Nu);           % GMI1�̌ŗL�l��SINR
Eigs_GMI2 = zeros(SIMU, Nr, Nu);           % GMI2�̌ŗL�l��SINR
Eigs_BMSN_BF = zeros(SIMU, Nr, Nu);        % BMSN-BF�̌ŗL�l
Eigs_BMSN_GE = zeros(SIMU, Nr, Nu);        % BMSN-GE�̌ŗL�l

    
for ik = 1:LK
    
    K_tar = K_box(ik);

    for isimu = 1:SIMU
    
    % �`���`���l���s��̃}���`�p�X�g����(NLoS�`���l��)
    H_iid = (randn(Nu*Nr,Nt)+1j*randn(Nu*Nr,Nt))/sqrt(2);
    %�`���`���l���s��̒��ڔg����(LoS�`���l��)
    H_los = zeros(Nu*Nr,Nt);
    
    Theta_t = (rand-0.5)*360;   % ���[�U���̑��M�p (-180deg - 180deg)
    Theta_r = (rand(1,Nu)-0.5)*360; % ���[�U���̎�M�p (-180deg - 180deg)
    a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t*derad));  % ���M���[�h�x�N�g��
    
    for n = 1 : Nu
        a_r = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad)); % ���[�U���̎�M���[�h�x�N�g��
        H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t.';                % ���[�U����LoS�`���l���s��
    end
      
    K = 10^(K_tar/10);
         
    % �`���`���l���s�� H=[sqrt(K/(K+1))*(LOS �`���l��)]...
    %                   .+[sqrt(1/(K+1))*(NLOS �`���l��)]
    H = sqrt(K/(K+1))*H_los + sqrt(1/(K+1))*H_iid;
    
    % ZF-CI algorithm
    [W_ZF,U_ZF,S_ZF,~,~] = zf(Nt,Nr,Nu,H); % function zf.m ���g�p
    
    % BD algorithm
    [W_BD,U_BD,S_BD,~,~] = bd(Nt,Nr,Nu,H); % function bd.m ���g�p
    
    % MMSE-CI algorithm
    [W_MMSE,U_MMSE,S_MMSE,RIPM,~] = mmse(Nt,Nr,Nu,H,a); % function mmse.m ���g�p
    
    % GMI1 algorithm
    [W_GMI1,U_GMI1,S_GMI1,RIPGM1,~] = gmmse_m1(Nt,Nr,Nu,H,a); % function gmmse.m1 ���g�p
    
    % GMI2 algorithm
    [W_GMI2,U_GMI2,S_GMI2,RIPGM2,~] = gmmse_m2(Nt,Nr,Nu,H,a); % function gmmse.m2 ���g�p
    
    % BMSN-BF algorithm
    [W_BMSN_BF,U_BMSN_BF,S_BMSN_BF,RIPBF,~] = bmsn_bf(Nt,Nr,Nu,H,a,T); % function bmsn_bf.m ���g�p
      
    % BMSN-GE algorithm
    [W_BMSN_GE,U_BMSN_GE,S_BMSN_GE,RIPGE,~] = bmsn_gev(Nt,Nr,Nu,H,a); % function bmsn_gev.m ���g�p
    
    
    % ���[�U���̌ŗL�l���z
     for nuser=1:Nu
        if Nr==1
            E_ZF(isimu,:,nuser) = 10*log10(S_ZF(1,1,nuser).^2/(Nt*snt));
            E_BD(isimu,:,nuser) = 10*log10(S_BD(1,1,nuser).^2/(Nt*snt));
            E_MMSE(isimu,:,nuser) = 10*log10((S_MMSE(1,1,nuser).^2)./(RIPM(1,nuser)+Nt*snt));
            E_GMI1(isimu,:,nuser) = 10*log10((S_GMI1(1,1,nuser).^2)./(RIPGM1(1,nuser)+Nt*snt));
            E_GMI2(isimu,:,nuser) = 10*log10((S_GMI2(1,1,nuser).^2)./(RIPGM2(1,nuser)+Nt*snt));
            E_BMSN_BF(isimu,:,nuser) = 10*log10((S_BMSN_BF(1,1,nuser).^2)./(RIPBF(1,nuser)+Nt*snt));
            E_BMSN_GE(isimu,:,nuser) = 10*log10((S_BMSN_GE(1,1,nuser).^2)./(RIPGE(1,nuser)+Nt*snt));
        else 
            E_ZF(isimu,:,nuser) = 10*log10(diag(S_ZF(:,:,nuser)).^2/(Nt*snt));   
            E_BD(isimu,:,nuser) = 10*log10(diag(S_BD(:,:,nuser)).^2/(Nt*snt));
            E_MMSE(isimu,:,nuser) = 10*log10((diag(S_MMSE(:,:,nuser)).^2)./(RIPM(:,nuser)+Nt*snt));
            E_GMI1(isimu,:,nuser) = 10*log10((diag(S_GMI1(:,:,nuser)).^2)./(RIPGM1(:,nuser)+Nt*snt));
            E_GMI2(isimu,:,nuser) = 10*log10((diag(S_GMI2(:,:,nuser)).^2)./(RIPGM2(:,nuser)+Nt*snt));
            E_BMSN_BF(isimu,:,nuser) = 10*log10((diag(S_BMSN_BF(:,:,nuser)).^2)./(RIPBF(:,nuser)+Nt*snt));
            E_BMSN_GE(isimu,:,nuser) = 10*log10((diag(S_BMSN_GE(:,:,nuser)).^2)./(RIPGE(:,nuser)+Nt*snt));
        end % Nr end
    end % nuser end  
    end % isimu end
    
    % �\�[�e�B���O�i�����j
    for nuser=1:Nu
        E_ZF(:,:,nuser) = sort(E_ZF(:,:,nuser),1);
        E_BD(:,:,nuser) = sort(E_BD(:,:,nuser),1);
        E_MMSE(:,:,nuser) = sort(E_MMSE(:,:,nuser),1);
        E_GMI1(:,:,nuser) = sort(E_GMI1(:,:,nuser),1);
        E_GMI2(:,:,nuser) = sort(E_GMI2(:,:,nuser),1);
        E_BMSN_BF(:,:,nuser) = sort(E_BMSN_BF(:,:,nuser),1);
        E_BMSN_GE(:,:,nuser) = sort(E_BMSN_GE(:,:,nuser),1);
    end
    % ���[�U����
    E_ZFm = mean(E_ZF,3);
    E_BDm = mean(E_BD,3);
    E_MMSEm = mean(E_MMSE,3);
    E_GMI1m = mean(E_GMI1,3);
    E_GMI2m = mean(E_GMI2,3);
    E_BMSN_BFm = mean(E_BMSN_BF,3);
    E_BMSN_GEm = mean(E_BMSN_GE,3);
    
    % �^�[�Q�b�gCDF�l�̌ŗL�l�𒊏o
    for nn=1:Nr
        rr3(ik,nn) = E_BMSN_BFm(round(CDF*SIMU/100),nn);     % BMSN_BF��2�̌ŗL�l����������
        rr3(ik,nn+Nr) = E_BMSN_GEm(round(CDF*SIMU/100),nn);  % BMSN_GE��2�̌ŗL�l����������
        rr3(ik,nn+2*Nr) = E_MMSEm(round(CDF*SIMU/100),nn);   % MMSE��2�̌ŗL�l����������
        rr3(ik,nn+3*Nr) = E_GMI1m(round(CDF*SIMU/100),nn);
        rr3(ik,nn+4*Nr) = E_GMI2m(round(CDF*SIMU/100),nn);
        rr3(ik,nn+5*Nr) = E_BDm(round(CDF*SIMU/100),nn);     % BD��2�̌ŗL�l����������
        rr3(ik,nn+6*Nr) = E_ZFm(round(CDF*SIMU/100),nn);     % ZF��2�̌ŗL�l����������
    end
    fprintf('K = %d dB\n',K_box(ik));  % ��͒��̌v�Z�ߒ������邽�߂�K = ??dB�̕\��
end% k_dB end

%% 
figure;
mycol = [0 0 1;
         0 0 1;
         1 0 0;
         1 0 0;
         0 0.7 0;
         0 0.7 0;
         1 0 1;
         1 0 1;
         0 0 0;
         0 0 0]; % �O���t�̐F
set(groot,'defaultAxesColorOrder',mycol) % figure�̏�����ݒ�
axis([K_min K_max SNR_tar-60 SNR_tar]);              % figure�̎��̂Ƃ����ݒ�
grid on;
hold on;
set(gca,'XTick',K_min:5:K_max,'Fontsize',14,'Fontname','Arial') % figure�̏�����ݒ�
plot(K_box,rr3(:,1),'r-d','LineWidth',2);
plot(K_box,rr3(:,2),'r--d','LineWidth',2);
plot(K_box,rr3(:,3),'b-s','LineWidth',2);
plot(K_box,rr3(:,4),'b--s','LineWidth',2);
plot(K_box,rr3(:,5),'g-o','LineWidth',2);
plot(K_box,rr3(:,6),'g--o','LineWidth',2);
plot(K_box,rr3(:,7),'y-+','LineWidth',2);
plot(K_box,rr3(:,8),'y--+','LineWidth',2);
plot(K_box,rr3(:,9),'k-*','LineWidth',2);
plot(K_box,rr3(:,10),'k--*','LineWidth',2);
plot(K_box,rr3(:,11),'c-x','LineWidth',2);
plot(K_box,rr3(:,12),'c--x','LineWidth',2);
plot(K_box,rr3(:,13),'m-^','LineWidth',2);
plot(K_box,rr3(:,14),'m--^','LineWidth',2);

% plot(K_box,rr3(:,1),'r-d','MarkerIndices',5:5:length(rr3(:,1)),'LineWidth',2);
% plot(K_box,rr3(:,2),'r--d','MarkerIndices',10:5:length(rr3(:,2)),'LineWidth',2);
% plot(K_box,rr3(:,3),'b-s','MarkerIndices',5:5:length(rr3(:,3)),'LineWidth',2);
% plot(K_box,rr3(:,4),'b--s','MarkerIndices',10:5:length(rr3(:,4)),'LineWidth',2);
% plot(K_box,rr3(:,5),'g-o','MarkerIndices',5:5:length(rr3(:,5)),'LineWidth',2);
% plot(K_box,rr3(:,6),'g--o','MarkerIndices',10:5:length(rr3(:,6)),'LineWidth',2);
% plot(K_box,rr3(:,7),'c-x','MarkerIndices',5:5:length(rr3(:,7)),'LineWidth',2);
% plot(K_box,rr3(:,8),'c--x','MarkerIndices',10:5:length(rr3(:,8)),'LineWidth',2);
% plot(K_box,rr3(:,9),'m-^','MarkerIndices',5:5:length(rr3(:,9)),'LineWidth',2);
% plot(K_box,rr3(:,10),'m--^','MarkerIndices',10:5:length(rr3(:,10)),'LineWidth',2);

legend('BMSN-BF \lambda_1','BMSN-BF \lambda_2','BMSN-GE \lambda_1','BMSN-GE \lambda_2','MMSE-CI \lambda_1','MMSE-CI \lambda_2',...
    'GMI1 \lambda_1','GMI1 \lambda_2','GMI2 \lambda_1','GMI2 \lambda_2','BD \lambda_1','BD \lambda_2','ZF-CI \lambda_1','ZF-CI \lambda_2','Location','southeast');
xlabel('K [dB]','Fontsize',16,'Fontname','Arial');
ylabel('SINR of eigenvalue [dB]','Fontsize',16,'Fontname','Arial');
set(gca,'Fontsize',16,'Fontname','Arial');
title(strcat(target_CDF,", ",target_SNR));
grid on;
hold on;

%End