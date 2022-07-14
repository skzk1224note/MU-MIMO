% K_lambda1,lambda2
% �ŗL���[�h�`����̓��ْl���ŗL�l�Ƃ��Ă���

clear;

% K[dB]�͈̔͐ݒ�
K_min = -20;        % �ŏ�K [dB]
K_max = 20;         % �ő�K [dB]
weight=0;   % weight figure on=0,off=1

Nt = 16;    % ���M�f�q��
Nr = 1;     % ��M�f�q�� (�A���e�i�I���̏ꍇ1, �����łȂ��ꍇ��2)
Nu = 16;     % ���[�U��

SIMU = 1000; % ���s�񐔁i�ʏ� 1000�j

I = eye(Nt,Nt); % Nt*Nt�̒P�ʍs��

SNR_tar = 10;   % �^�[�Q�b�gSNR[dB]

d_t      = 0.5;      % ���M�A���e�i�Ԋu�iin wavelength)
d_r      = 0.5;      % ��M�A���e�i�Ԋu�iin wavelength)
derad = pi/180;      % degree -> rad

%T ���]�̃`���l���s��
T = zeros(Nr,Nr);
for inu = 1:Nu
    T(:,:,inu) = eye(Nr,Nr);
end
%%
if SNR_tar < 10
    target_SNR=strcat('SNR= ',num2str(SNR_tar,'%01d'),'dB');
else
    target_SNR=strcat('SNR= ',num2str(SNR_tar,'%02d'),'dB');
end
H_isimu=zeros(Nu*Nr,Nt,SIMU);
H_iid_isimu=zeros(Nu*Nr,Nt,SIMU);
H_los_isimu=zeros(Nu*Nr,Nt,SIMU);

Norm_H=zeros(Nu*Nr,1,SIMU);      %���[�U�`���l���̃m����
Norm_H_los=zeros(Nu*Nr,1,SIMU);
Norm_H_iid=zeros(Nu*Nr,1,SIMU);

HW_BMSN_BF=zeros(Nr*Nu,Nr*Nu,SIMU);%HW�̐���
HW_BMSN_GE=zeros(Nr*Nu,Nr*Nu,SIMU);
HW_MMSE=zeros(Nr*Nu,Nr*Nu,SIMU);
HW_BD=zeros(Nr*Nu,Nr*Nu,SIMU);
HW_ZF=zeros(Nr*Nu,Nr*Nu,SIMU);

HW_BMSN_BFn=zeros(1,SIMU);%HW�m�����̐���
HW_BMSN_GEn=zeros(1,SIMU);
HW_MMSEn=zeros(1,SIMU);
HW_BDn=zeros(1,SIMU);
HW_ZFn=zeros(1,SIMU);

Norm_W_ZF_whole=zeros(1,SIMU);%�S��W�m�����̐���
Norm_W_BD_whole=zeros(1,SIMU);
Norm_W_MMSE_whole=zeros(1,SIMU);
Norm_W_BMSN_BF_whole=zeros(1,SIMU);
Norm_W_BMSN_GE_whole=zeros(1,SIMU);

Norm_W_ZF_user=zeros(Nu,SIMU);%���[�UW�m�����̐���
Norm_W_BD_user=zeros(Nu,SIMU);
Norm_W_MMSE_user=zeros(Nu,SIMU);
Norm_W_BMSN_BF_user=zeros(Nu,SIMU);
Norm_W_BMSN_GE_user=zeros(Nu,SIMU);

Norm_W_ZF_user_avsimu=zeros(1,Nu);%���[�UW�m�������s�񐔕��ς̐���
Norm_W_BD_user_avsimu=zeros(1,Nu);
Norm_W_MMSE_user_avsimu=zeros(1,Nu);
Norm_W_BMSN_BF_user_avsimu=zeros(1,Nu);
Norm_W_BMSN_GE_user_avsimu=zeros(1,Nu);
%%
K_box=(K_min:5:K_max).'; % figure�̉����̂��߂�K�̔�
LK=length(K_box);        % K�̔��̑傫��
sigma2 = 1/(10^(SNR_tar/10)); % noise power
a = sigma2*Nt;
 dd3=zeros(LK,5); hh3=zeros(LK,5); ii3=zeros(LK,5); 

% �o�̓t�@�C���� with alpha and SNR in dB
%evfile1 = strcat(folder,fn1,num2str(snr(isnr),'%02d'),'dB_1000itr.csv');

W_ZF_isimu = zeros(Nt,Nr,Nu,SIMU);     %W�̍쐬        
W_BD_isimu = zeros(Nt,Nr,Nu,SIMU);             
W_MMSE_isimu = zeros(Nt,Nr,Nu,SIMU);           
W_BMSN_BF_isimu = zeros(Nt,Nr,Nu,SIMU);       
W_BMSN_GE_isimu = zeros(Nt,Nr,Nu,SIMU);     

W_ZF_whole=zeros(Nt,Nr*Nu,SIMU);            %W�S�̂̍쐬
W_BD_whole=zeros(Nt,Nr*Nu,SIMU);
W_MMSE_whole=zeros(Nt,Nr*Nu,SIMU);
W_BMSN_BF_whole=zeros(Nt,Nr*Nu,SIMU);
W_BMSN_GE_whole=zeros(Nt,Nr*Nu,SIMU);
%%
for ik = 1:LK
    K_tar = K_box(ik);
    
    for isimu = 1:SIMU
    
        H_iid = (randn(Nu*Nr,Nt)+1j*randn(Nu*Nr,Nt))/sqrt(2); % �`���`���l���s��̃}���`�p�X�g����(NLoS�`���l��)
        H_los = zeros(Nu*Nr,Nt);                              % �`���`���l���s��̒��ڔg����(LoS�`���l��)
    
        Theta_t = (rand(1,Nu)-0.5)*360;   % ���[�U���̑��M�p (-180deg - 180deg)
        Theta_r = (rand(1,Nu)-0.5)*360;   % ���[�U���̎�M�p (-180deg - 180deg)  % a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t*derad));  % ���M���[�h�x�N�g��
    
    
        for n = 1 : Nu
            a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad)); % ���[�U���̑��M���[�h�x�N�g��
            a_r = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad)); % ���[�U���̎�M���[�h�x�N�g��
            H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t';                 % ���[�U����LOS�`���l���s��
        end
      
    K = 10^(K_tar/10);
         
    H = sqrt(K/(K+1))*H_los + sqrt(1/(K+1))*H_iid; % �`���`���l���s�� H=[sqrt(K/(K+1))*(LOS �`���l��)]+[sqrt(1/(K+1))*(NLOS �`���l��)]
   
    % ZF-CI algorithm
    [W_ZF,U_ZF,S_ZF,~,~] = zf(Nt,Nr,Nu,H); % function zf.m ���g�p
    
    % BD algorithm
    [W_BD,U_BD,S_BD,~,~] = bd(Nt,Nr,Nu,H); % function bd.m ���g�p
    
    % MMSE-CI algorithm
    [W_MMSE,U_MMSE,S_MMSE,RIPM,~] = mmse(Nt,Nr,Nu,H,a); % function mmse.m ���g�p
    
    % BMSN-BF algorithm
    [W_BMSN_BF,U_BMSN_BF,S_BMSN_BF,RIPBF,~] = bmsn_bf(Nt,Nr,Nu,H,a,T); % function bmsn_bf.m ���g�p
      
    % BMSN-GE algorithm
    [W_BMSN_GE,U_BMSN_GE,S_BMSN_GE,RIPGE,~] = bmsn_gev(Nt,Nr,Nu,H,a); % function bmsn_gev.m ���g�p
    
    W_ZF_isimu(:,:,:,isimu)=W_ZF;          
    W_BD_isimu(:,:,:,isimu) = W_BD;             
    W_MMSE_isimu(:,:,:,isimu) = W_MMSE;           
    W_BMSN_BF_isimu(:,:,:,isimu) =W_BMSN_BF;       
    W_BMSN_GE_isimu(:,:,:,isimu) =W_BMSN_GE;
    
    H_isimu(:,:,isimu)=H;
    H_iid_isimu(:,:,isimu)=H_iid;
    H_los_isimu(:,:,isimu)=H_los;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nuser=1:Nu
    ns = Nr*(nuser-1)+1:Nr*nuser;
    HW_BMSN_BF(ns,ns,isimu)=H_isimu(ns,:,isimu)*W_BMSN_BF_isimu(:,:,nuser,isimu);
    HW_BMSN_GE(ns,ns,isimu)=H_isimu(ns,:,isimu)*W_BMSN_GE_isimu(:,:,nuser,isimu);
    HW_MMSE(ns,ns,isimu)=H_isimu(ns,:,isimu)*W_MMSE_isimu(:,:,nuser,isimu);
    HW_BD(ns,ns,isimu)=H_isimu(ns,:,isimu)*W_BD_isimu(:,:,nuser,isimu);
    HW_ZF(ns,ns,isimu)=H_isimu(ns,:,isimu)*W_ZF_isimu(:,:,nuser,isimu);
end

for n = 1 : Nu*Nr
        Norm_H(n,:,isimu)=norm(H_isimu(n,:,isimu));   % HW�̍쐬
        Norm_H_iid(n,:,isimu)=norm(H_iid_isimu(n,:,isimu));
        Norm_H_los(n,:,isimu)=norm(H_los_isimu(n,:,isimu));
end
    HW_BMSN_BFn(1,isimu)=norm(HW_BMSN_BF(:,:,isimu));  % HW�̃m����
    HW_BMSN_GEn(1,isimu)=norm(HW_BMSN_GE(:,:,isimu));
    HW_MMSEn(1,isimu)=norm(HW_MMSE(:,:,isimu));
    HW_BDn(1,isimu)=norm(HW_BD(:,:,isimu));
    HW_ZFn(1,isimu)=norm(HW_ZF(:,:,isimu));
    
for nuser=1:Nu
    ns = Nr*(nuser-1)+1:Nr*nuser;    
    W_ZF_whole(:,ns,isimu)=W_ZF_isimu(:,:,nuser,isimu);
    W_BD_whole(:,ns,isimu)=W_BD_isimu(:,:,nuser,isimu);
    W_MMSE_whole(:,ns,isimu)=W_MMSE_isimu(:,:,nuser,isimu);
    W_BMSN_BF_whole(:,ns,isimu)=W_BMSN_BF_isimu(:,:,nuser,isimu);
    W_BMSN_GE_whole(:,ns,isimu)=W_BMSN_GE_isimu(:,:,nuser,isimu);
end
    Norm_W_ZF_whole(1,isimu)=norm(W_ZF_whole(:,:,isimu));% W�S�̂̃m����
    Norm_W_BD_whole(1,isimu)=norm(W_BD_whole(:,:,isimu));
    Norm_W_MMSE_whole(1,isimu)=norm(W_MMSE_whole(:,:,isimu));
    Norm_W_BMSN_BF_whole(1,isimu)=norm(W_BMSN_BF_whole(:,:,isimu));
    Norm_W_BMSN_GE_whole(1,isimu)=norm(W_BMSN_GE_whole(:,:,isimu));
    
for nuser=1:Nu    
    Norm_W_ZF_user(nuser,isimu)=norm(W_ZF_isimu(:,:,nuser,isimu));% ���[�U���Ƃ�W�m����
    Norm_W_BD_user(nuser,isimu)=norm(W_BD_isimu(:,:,nuser,isimu));
    Norm_W_MMSE_user(nuser,isimu)=norm(W_MMSE_isimu(:,:,nuser,isimu));
    Norm_W_BMSN_BF_user(nuser,isimu)=norm(W_BMSN_BF_isimu(:,:,nuser,isimu));
    Norm_W_BMSN_GE_user(nuser,isimu)=norm(W_BMSN_GE_isimu(:,:,nuser,isimu));
end
    end   % isimu end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    W_ZF_isimu_av=mean(W_ZF_isimu,4);% ���s�񐔂̕���
    W_BD_isimu_av=mean(W_BD_isimu,4);
    W_MMSE_isimu_av=mean(W_MMSE_isimu,4);
    W_BMSN_BF_isimu_av=mean(W_BMSN_BF_isimu,4);
    W_BMSN_GE_isimu_av=mean(W_BMSN_GE_isimu,4);
    
    W_ZF_user_av=mean(W_ZF_isimu_av,3);% ���[�U�̕���
    W_BD_user_av=mean(W_BD_isimu_av,3);
    W_MMSE_user_av=mean(W_MMSE_isimu_av,3);
    W_BMSN_BF_user_av=mean(W_BMSN_BF_isimu_av,3);
    W_BMSN_GE_user_av=mean(W_BMSN_GE_isimu_av,3);
    
    Norm_W_ZF=mean(Norm_W_ZF_whole,2);% W�S�̂̃m�����̎��s�񐔕���
    Norm_W_BD=mean(Norm_W_BD_whole,2);
    Norm_W_MMSE=mean(Norm_W_MMSE_whole,2);
    Norm_W_BMSN_BF=mean(Norm_W_BMSN_BF_whole,2);
    Norm_W_BMSN_GE=mean(Norm_W_BMSN_GE_whole,2);
      
    Norm_W_ZF_user_avsimu=mean(Norm_W_ZF_user,2);% ���[�U���Ƃ�W�m�����̎��s�񐔕���
    Norm_W_BD_user_avsimu=mean(Norm_W_BD_user,2);
    Norm_W_MMSE_user_avsimu=mean(Norm_W_MMSE_user,2);
    Norm_W_BMSN_BF_user_avsimu=mean(Norm_W_BMSN_BF_user,2);
    Norm_W_BMSN_GE_user_avsimu=mean(Norm_W_BMSN_GE_user,2);
    
    Norm_W_ZF_avsimu_user=mean(Norm_W_ZF_user_avsimu);% ���[�U���Ƃ�W�m�����̎��s�񐔕��ς̂��Ƃ���Ƀ��[�U����
    Norm_W_BD_avsimu_user=mean(Norm_W_BD_user_avsimu);
    Norm_W_MMSE_avsimu_user=mean(Norm_W_MMSE_user_avsimu);
    Norm_W_BMSN_BF_avsimu_user=mean(Norm_W_BMSN_BF_user_avsimu);
    Norm_W_BMSN_GE_avsimu_user=mean(Norm_W_BMSN_GE_user_avsimu);        
 
    Norm_HW_BMSN_BFn=mean(HW_BMSN_BFn,2);
    Norm_HW_BMSN_GEn=mean(HW_BMSN_GEn,2);
    Norm_HW_MMSEn=mean(HW_MMSEn,2);
    Norm_HW_BDn=mean(HW_BDn,2);
    Norm_HW_ZFn=mean(HW_ZFn,2);
    
    Norm_H_av=mean(Norm_H,3);
    Norm_H_iid_av=mean(Norm_H_iid,3);
    Norm_H_los_av=mean(Norm_H_los,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
        dd3(ik,1) = Norm_W_BMSN_GE_avsimu_user;     % ���M�E�G�C�gBMSN_BF�̕���
        dd3(ik,2) = Norm_W_BMSN_BF_avsimu_user;  % ���M�E�G�C�gBMSN_GE�̕���
        dd3(ik,3) = Norm_W_MMSE_avsimu_user;   % ���M�E�G�C�gMMSE�̕���
        dd3(ik,4) = Norm_W_BD_avsimu_user;     % ���M�E�G�C�gBD�̕���
        dd3(ik,5) = Norm_W_ZF_avsimu_user;     % ���M�E�G�C�gZF�̕���            
       
        hh3(ik,1) = Norm_HW_BMSN_BFn;     % BMSN_BF��HW
        hh3(ik,2) = Norm_HW_BMSN_GEn;  % BMSN_GE��HW
        hh3(ik,3) = Norm_HW_MMSEn;   % MMSE��HW
        hh3(ik,4) = Norm_HW_BDn;     % BD��HW
        hh3(ik,5) = Norm_HW_ZFn;     % ZF��HW
        
        ii3(ik,1) = Norm_W_BMSN_BF;     % BMSN_BF��W
        ii3(ik,2) = Norm_W_BMSN_GE;  % BMSN_GE��W
        ii3(ik,3) = Norm_W_MMSE;   % MMSE��W
        ii3(ik,4) = Norm_W_BD;     % BD��W
        ii3(ik,5) = Norm_W_ZF;     % ZF��W
    
    fprintf('K = %d dB\n',K_box(ik));  % ��͒��̌v�Z�ߒ������邽�߂�K = ??dB�̕\��
 end% ik end

if weight==0
%% �O���t�\���@���M�E�G�C�g(���[�U������)
figure;
mycol = [0 0 1;1 0 0;0 0.7 0;1 0 1;0 0 0]; % �O���t�̐F
set(groot,'defaultAxesColorOrder',mycol) % figure�̏�����ݒ�
axis([K_min K_max 0 2]);              % figure�̎��̂Ƃ����ݒ�
grid on;
hold on;
set(gca,'XTick',K_min:5:K_max,'Fontsize',14,'Fontname','Times New Roman') % figure�̏�����ݒ�
plot(K_box,dd3(:,1),'r-d','LineWidth',4);
plot(K_box,dd3(:,2),'b-s','LineWidth',4);
plot(K_box,dd3(:,3),'g-o','LineWidth',4);
plot(K_box,dd3(:,4),'c-x','LineWidth',4);
plot(K_box,dd3(:,5),'m-^','LineWidth',4);

legend('BMSN-BF','BMSN-GE','MMSE-CI','BD','ZF-CI','Location','southeast');
xlabel('{\it{K}} [dB]','Fontsize',36,'Fontname','Times New Roman');
ylabel('norm {\bf{W}}_{user}','Fontsize',36,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
title(target_SNR);
grid on;
hold on;
%% �O���t�\���@���M�E�G�C�g
figure;
mycol = [0 0 1;1 0 0;0 0.7 0;1 0 1;0 0 0]; % �O���t�̐F
set(groot,'defaultAxesColorOrder',mycol) % figure�̏�����ݒ�
axis([K_min K_max 0 3]);              % figure�̎��̂Ƃ����ݒ�
grid on;
hold on;
set(gca,'XTick',K_min:5:K_max,'Fontsize',14,'Fontname','Times New Roman') % figure�̏�����ݒ�
plot(K_box,ii3(:,1),'r-d','LineWidth',4);
plot(K_box,ii3(:,2),'b-s','LineWidth',4);
plot(K_box,ii3(:,3),'g-o','LineWidth',4);
plot(K_box,ii3(:,4),'c-x','LineWidth',4);
plot(K_box,ii3(:,5),'m-^','LineWidth',4);

legend('BMSN-BF','BMSN-GE','MMSE-CI','BD','ZF-CI','Location','southeast');
xlabel('{\it{K}} [dB]','Fontsize',36,'Fontname','Times New Roman');
ylabel('norm \bf{W}','Fontsize',36,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
title(target_SNR);
grid on;
hold on;

end
%% �O���t�\��
figure;
mycol = [0 0 1;1 0 0;0 0.7 0;1 0 1;0 0 0]; % �O���t�̐F
set(groot,'defaultAxesColorOrder',mycol) % figure�̏�����ݒ�
axis([K_min K_max 0 6]);              % figure�̎��̂Ƃ����ݒ�
grid on;
hold on;
set(gca,'XTick',K_min:5:K_max,'Fontsize',14,'Fontname','Times New Roman') % figure�̏�����ݒ�
plot(K_box,hh3(:,1),'r-d','LineWidth',4);
plot(K_box,hh3(:,2),'b-s','LineWidth',4);
plot(K_box,hh3(:,3),'g-o','LineWidth',4);
plot(K_box,hh3(:,4),'c-x','LineWidth',4);
plot(K_box,hh3(:,5),'m-^','LineWidth',4);

legend('BMSN-BF','BMSN-GE','MMSE-CI','BD','ZF-CI','Location','southeast');
xlabel('{\it{K}} [dB]','Fontsize',36,'Fontname','Times New Roman');
ylabel('norm \bf{HW}','Fontsize',36,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
title(target_SNR);
grid on;
hold on;
% End