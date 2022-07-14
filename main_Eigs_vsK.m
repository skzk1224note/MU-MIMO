% K_lambda1,lambda2
% �ŗL���[�h�`����̓��ْl���ŗL�l�Ƃ��Ă���

clear;

% K[dB]�͈̔͐ݒ�
K_min = -20;        % �ŏ�K [dB]
K_max = 20;         % �ő�K [dB]

CDF = 10;   % ����CDF
wantCDF = 1;  % 0�Ȃ�off,1�Ȃ�on(figure�p)

Nt = 16;    % ���M�f�q��
Nr = 2;     % ��M�f�q�� (�A���e�i�I���̏ꍇ1, �����łȂ��ꍇ��2)
Nu = 8;     % ���[�U��

SIMU = 10; % ���s�񐔁i�ʏ� 1000�j
Directivity_switch = 1;% 1:on,0:off ����M�f�q�̎w�����l���̗L��
An = 2; % �w�����֐��̌W��

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
%%
K_box=(K_min:5:K_max).'; % figure�̉����̂��߂�K�̔�
LK=length(K_box);        % K�̔��̑傫��
sigma2 = 1/(10^(SNR_tar/10)); % noise power
a = sigma2*Nt;
aa3=zeros(LK,10); bb3=zeros(LK,10); cc3=zeros(LK,5); dd3=zeros(LK,10); 
ee3=zeros(LK,10); ff3=zeros(LK,5);

% �o�̓t�@�C���� with alpha and SNR in dB
%evfile1 = strcat(folder,fn1,num2str(snr(isnr),'%02d'),'dB_1000itr.csv');

E_ZF = zeros(SIMU, Nr, Nu);             % ZF-CI�̌ŗL�l��SINR
E_BD = zeros(SIMU, Nr, Nu);             % BD�̌ŗL�l��SINR
E_MMSE = zeros(SIMU, Nr, Nu);           % MMSE-CI�̌ŗL�l��SINR
E_BMSN_BF = zeros(SIMU, Nr, Nu);        % BMSN-BF�̌ŗL�l��SINR
E_BMSN_GE = zeros(SIMU, Nr, Nu);        % BMSN-GE�̌ŗL�l��SINR

Eigs_ZF = zeros(SIMU, Nr, Nu);             % ZF-CI�̌ŗL�l
Eigs_BD = zeros(SIMU, Nr, Nu);             % BD�̌ŗL�l
Eigs_MMSE = zeros(SIMU, Nr, Nu);           % MMSE-CI�̌ŗL�l
Eigs_BMSN_BF = zeros(SIMU, Nr, Nu);        % BMSN-BF�̌ŗL�l
Eigs_BMSN_GE = zeros(SIMU, Nr, Nu);        % BMSN-GE�̌ŗL�l

for ik = 1:LK
    K_tar = K_box(ik);
    
    for isimu = 1:SIMU
    
        H_iid = (randn(Nu*Nr,Nt)+1j*randn(Nu*Nr,Nt))/sqrt(2); % �`���`���l���s��̃}���`�p�X�g����(NLoS�`���l��)
        H_los = zeros(Nu*Nr,Nt);                              % �`���`���l���s��̒��ڔg����(LoS�`���l��)
    
        Theta_t = (rand(1,Nu)-0.5)*360;   % ���[�U���̑��M�p (-180deg - 180deg)
        Theta_r = (rand(1,Nu)-0.5)*360;   % ���[�U���̎�M�p (-180deg - 180deg)  % a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t*derad));  % ���M���[�h�x�N�g��
    
    
        if Directivity_switch == 1
            for n = 1 : Nu
                a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad))*An*cos(Theta_t(1,n)*derad); % ���[�U���̑��M���[�h�x�N�g��
                a_r = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad))*An*cos(Theta_r(1,n)*derad); % ���[�U���̎�M���[�h�x�N�g��
                H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t';                 % ���[�U����LOS�`���l���s��
            end
        else
            for n = 1 : Nu
            a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad)); % ���[�U���̑��M���[�h�x�N�g��
            a_r = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad)); % ���[�U���̎�M���[�h�x�N�g��
            H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t';                 % ���[�U����LOS�`���l���s��
            end
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
        
    for nuser=1:Nu         % ���[�U���̌ŗL�l���z
        if Nr==1
            E_ZF(isimu,:,nuser) = 10*log10(S_ZF(1,1,nuser).^2/(Nt*sigma2));
            E_BD(isimu,:,nuser) = 10*log10(S_BD(1,1,nuser).^2/(Nt*sigma2));
            E_MMSE(isimu,:,nuser) = 10*log10((S_MMSE(1,1,nuser).^2)./(Nt*sigma2));
            E_BMSN_BF(isimu,:,nuser) = 10*log10((S_BMSN_BF(1,1,nuser).^2)./(Nt*sigma2));
            E_BMSN_GE(isimu,:,nuser) = 10*log10((S_BMSN_GE(1,1,nuser).^2)./(Nt*sigma2));
        else 
            E_ZF(isimu,:,nuser) = 10*log10(diag(S_ZF(:,:,nuser)).^2/(Nt*sigma2));   
            E_BD(isimu,:,nuser) = 10*log10(diag(S_BD(:,:,nuser)).^2/(Nt*sigma2));
            E_MMSE(isimu,:,nuser) = 10*log10((diag(S_MMSE(:,:,nuser)).^2)./(Nt*sigma2));
            E_BMSN_BF(isimu,:,nuser) = 10*log10((diag(S_BMSN_BF(:,:,nuser)).^2)./(Nt*sigma2));
            E_BMSN_GE(isimu,:,nuser) = 10*log10((diag(S_BMSN_GE(:,:,nuser)).^2)./(Nt*sigma2));
        end % Nr end
    end % nuser end
    
  
    for nuser=1:Nu          %�@��������ŗL�l���̂��̂̒l�����
        if Nr==1
            Eigs_ZF(isimu,:,nuser) = S_ZF(1,1,nuser);
            Eigs_BD(isimu,:,nuser) = S_BD(1,1,nuser);
            Eigs_MMSE(isimu,:,nuser) = S_MMSE(1,1,nuser);
            Eigs_BMSN_BF(isimu,:,nuser) = S_BMSN_BF(1,1,nuser);
            Eigs_BMSN_GE(isimu,:,nuser) = S_BMSN_GE(1,1,nuser);
        else 
            Eigs_ZF(isimu,:,nuser) = diag(S_ZF(:,:,nuser));
            Eigs_BD(isimu,:,nuser) = diag(S_BD(:,:,nuser));
            Eigs_MMSE(isimu,:,nuser) = diag(S_MMSE(:,:,nuser));
            Eigs_BMSN_BF(isimu,:,nuser) = diag(S_BMSN_BF(:,:,nuser));
            Eigs_BMSN_GE(isimu,:,nuser) = diag(S_BMSN_GE(:,:,nuser));
        end % Nr end
    end % nuser end    
    end % isimu end
        
    for nuser=1:Nu    % �\�[�e�B���O�i�����j
        E_ZF(:,:,nuser) = sort(E_ZF(:,:,nuser),1);
        E_BD(:,:,nuser) = sort(E_BD(:,:,nuser),1);
        E_MMSE(:,:,nuser) = sort(E_MMSE(:,:,nuser),1);
        E_BMSN_BF(:,:,nuser) = sort(E_BMSN_BF(:,:,nuser),1);
        E_BMSN_GE(:,:,nuser) = sort(E_BMSN_GE(:,:,nuser),1);% �ŗL�l��SINR
        
        Eigs_ZF(:,:,nuser) = sort( Eigs_ZF(:,:,nuser),1);
        Eigs_BD(:,:,nuser) = sort( Eigs_BD(:,:,nuser),1);
        Eigs_MMSE(:,:,nuser) = sort( Eigs_MMSE(:,:,nuser),1);
        Eigs_BMSN_BF(:,:,nuser) = sort( Eigs_BMSN_BF(:,:,nuser),1);
        Eigs_BMSN_GE(:,:,nuser) = sort( Eigs_BMSN_GE(:,:,nuser),1);% �ŗL�l���̂���
    end
   %%  % ���[�U����
    
    E_ZFm = mean(E_ZF,3);
    E_BDm = mean(E_BD,3);
    E_MMSEm = mean(E_MMSE,3);
    E_BMSN_BFm = mean(E_BMSN_BF,3);
    E_BMSN_GEm = mean(E_BMSN_GE,3);% �ŗL�l��SINR
    
    % �S���s�񐔂̕���
    E_ZF_all = mean(E_ZFm,1);
    E_BD_all = mean(E_BDm,1);
    E_MMSE_all = mean(E_MMSEm,1);
    E_BMSN_BF_all = mean(E_BMSN_BFm,1);
    E_BMSN_GE_all = mean(E_BMSN_GEm,1);
    
    % �ŗL�lSINR��Eigs1��Eigs2�̕���
    E_ZF_whole = mean(E_ZF_all);
    E_BD_whole = mean(E_BD_all);
    E_MMSE_whole = mean(E_MMSE_all);
    E_BMSN_BF_whole = mean(E_BMSN_BF_all);
    E_BMSN_GE_whole = mean(E_BMSN_GE_all);
   %% ���[�U����
    Eigs_ZFm = mean(Eigs_ZF,3);
    Eigs_BDm = mean(Eigs_BD,3);
    Eigs_MMSEm = mean(Eigs_MMSE,3);
    Eigs_BMSN_BFm = mean(Eigs_BMSN_BF,3);
    Eigs_BMSN_GEm = mean(Eigs_BMSN_GE,3);% �ŗL�l���̂���
    
    % �S���s�񐔂̕���
    Eigs_ZF_all = mean(Eigs_ZFm,1);
    Eigs_BD_all = mean(Eigs_BDm,1);
    Eigs_MMSE_all = mean(Eigs_MMSEm,1);
    Eigs_BMSN_BF_all = mean(Eigs_BMSN_BFm,1);
    Eigs_BMSN_GE_all = mean(Eigs_BMSN_GEm,1);
    
    % �ŗL�l���̂��̂�Eigs1��Eigs2�̕���
    Eigs_ZF_whole = mean(Eigs_ZF_all);
    Eigs_BD_whole = mean(Eigs_BD_all);
    Eigs_MMSE_whole = mean(Eigs_MMSE_all);
    Eigs_BMSN_BF_whole = mean(Eigs_BMSN_BF_all);
    Eigs_BMSN_GE_whole = mean(Eigs_BMSN_GE_all);
   %%   �^�[�Q�b�gCDF�l�̌ŗL�l�𒊏o
    
    for nn=1:Nr
        aa3(ik,nn) = Eigs_BMSN_BF_all(:,nn);     % BMSN_BF��2�̌ŗL�l�̎��s�񐔕��ς���������
        aa3(ik,nn+Nr) = Eigs_BMSN_GE_all(:,nn);  % BMSN_GE��2�̌ŗL�l�̎��s�񐔕��ς���������
        aa3(ik,nn+2*Nr) = Eigs_MMSE_all(:,nn);   % MMSE��2�̌ŗL�l�̎��s�񐔕��ς���������
        aa3(ik,nn+3*Nr) = Eigs_BD_all(:,nn);     % BD��2�̌ŗL�l�̎��s�񐔕��ς���������
        aa3(ik,nn+4*Nr) = Eigs_ZF_all(:,nn);     % ZF��2�̌ŗL�l�̎��s�񐔕��ς���������
                
        bb3(ik,nn) = E_BMSN_BF_all(:,nn);     % BMSN_BF��2�̌ŗL�lSINR�̎��s�񐔕��ς���������
        bb3(ik,nn+Nr) = E_BMSN_GE_all(:,nn);  % BMSN_GE��2�̌ŗL�lSINR�̎��s�񐔕��ς���������
        bb3(ik,nn+2*Nr) = E_MMSE_all(:,nn);   % MMSE��2�̌ŗL�lSINR�̎��s�񐔕��ς���������
        bb3(ik,nn+3*Nr) = E_BD_all(:,nn);     % BD��2�̌ŗL�lSINR�̎��s�񐔕��ς���������
        bb3(ik,nn+4*Nr) = E_ZF_all(:,nn);     % ZF��2�̌ŗL�lSINR�̎��s�񐔕��ς���������
        
        cc3(ik,1) = E_BMSN_BF_whole;     % BMSN_BF��2�̌ŗL�lSINR�̕���
        cc3(ik,2) = E_BMSN_GE_whole;  % BMSN_GE��2�̌ŗL�lSINR�̕���
        cc3(ik,3) = E_MMSE_whole;   % MMSE��2�̌ŗL�lSINR�̕���
        cc3(ik,4) = E_BD_whole;     % BD��2�̌ŗL�lSINR�̕���
        cc3(ik,5) = E_ZF_whole;     % ZF��2�̌ŗL�lSINR�̕���
               
        dd3(ik,nn) = Eigs_BMSN_BFm(round(CDF*SIMU/100),nn);     % BMSN_BF��2�̌ŗL�l����������
        dd3(ik,nn+Nr) = Eigs_BMSN_GEm(round(CDF*SIMU/100),nn);  % BMSN_GE��2�̌ŗL�l����������
        dd3(ik,nn+2*Nr) = Eigs_MMSEm(round(CDF*SIMU/100),nn);   % MMSE��2�̌ŗL�l����������
        dd3(ik,nn+3*Nr) = Eigs_BDm(round(CDF*SIMU/100),nn);     % BD��2�̌ŗL�l����������
        dd3(ik,nn+4*Nr) = Eigs_ZFm(round(CDF*SIMU/100),nn);     % ZF��2�̌ŗL�l����������
        
        ee3(ik,nn) = E_BMSN_BFm(round(CDF*SIMU/100),nn);     % BMSN_BF��2�̌ŗL�l��SINR����������
        ee3(ik,nn+Nr) = E_BMSN_GEm(round(CDF*SIMU/100),nn);  % BMSN_GE��2�̌ŗL�l��SINR����������
        ee3(ik,nn+2*Nr) = E_MMSEm(round(CDF*SIMU/100),nn);   % MMSE��2�̌ŗL�l��SINR����������
        ee3(ik,nn+3*Nr) = E_BDm(round(CDF*SIMU/100),nn);     % BD��2�̌ŗL�l��SINR����������
        ee3(ik,nn+4*Nr) = E_ZFm(round(CDF*SIMU/100),nn);     % ZF��2�̌ŗL�l��SINR����������
        
        ff3(ik,1) = Eigs_BMSN_BF_whole;     % BMSN_BF��2�̌ŗL�l�̕���
        ff3(ik,2) = Eigs_BMSN_GE_whole;  % BMSN_GE��2�̌ŗL�l�̕���
        ff3(ik,3) = Eigs_MMSE_whole;   % MMSE��2�̌ŗL�l�̕���
        ff3(ik,4) = Eigs_BD_whole;     % BD��2�̌ŗL�l�̕���
        ff3(ik,5) = Eigs_ZF_whole;     % ZF��2�̌ŗL�l�̕���
    end
    fprintf('K = %d dB\n',K_box(ik));  % ��͒��̌v�Z�ߒ������邽�߂�K = ??dB�̕\��
 end% ik end

%% CDF�ɒ��ڂ����ꍇ�̃O���t
if wantCDF==1
%% K�ɑ΂���ŗL�l[dB]�O���t�iCDF�ɒ��ځj
   figure;
   mycol = [0 0 1;0 0 1;1 0 0;1 0 0;0 0.7 0;0 0.7 0;1 0 1;1 0 1;0 0 0;0 0 0]; % �O���t�̐F
   set(groot,'defaultAxesColorOrder',mycol) % figure�̏�����ݒ�
   axis([K_min K_max SNR_tar-50 SNR_tar+10]);              % figure�̎��̂Ƃ����ݒ�
   grid on;
   hold on;
   set(gca,'XTick',K_min:5:K_max,'Fontsize',14,'Fontname','Times New Roman') % figure�̏�����ݒ�
   if Nr==2
      plot(K_box,ee3(:,1),'r-d','LineWidth',4);
      plot(K_box,ee3(:,2),'r--d','LineWidth',4);
      plot(K_box,ee3(:,3),'b-s','LineWidth',4);
      plot(K_box,ee3(:,4),'b--s','LineWidth',4);
      plot(K_box,ee3(:,5),'g-o','LineWidth',4);
      plot(K_box,ee3(:,6),'g--o','LineWidth',4);
      plot(K_box,ee3(:,7),'c-x','LineWidth',4);
      plot(K_box,ee3(:,8),'c--x','LineWidth',4);
      plot(K_box,ee3(:,9),'m-^','LineWidth',4);
      plot(K_box,ee3(:,10),'m--^','LineWidth',4);
   end
   if Nr==1
      plot(K_box,ee3(:,1),'r-d','LineWidth',4);
      plot(K_box,ee3(:,2),'b-s','LineWidth',4);
      plot(K_box,ee3(:,3),'g-o','LineWidth',4);
      plot(K_box,ee3(:,4),'c-x','LineWidth',4);
      plot(K_box,ee3(:,5),'m-^','LineWidth',4);
   end
   if Nr==1
   legend('BMSN-BF \lambda_1','BMSN-GE \lambda_1','MMSE-CI \lambda_1','BD \lambda_1','ZF-CI \lambda_1','Location','southeast');
   end
   if Nr==2
   legend('BMSN-BF \lambda_1','BMSN-BF \lambda_2','BMSN-GE \lambda_1','BMSN-GE \lambda_2','MMSE-CI \lambda_1','MMSE-CI \lambda_2',...
    'BD \lambda_1','BD \lambda_2','ZF-CI \lambda_1','ZF-CI \lambda_2','Location','southeast');
   end
   xlabel('{\it{K}} [dB]','Fontsize',16,'Fontname','Times New Roman');
   ylabel('Eigenvalues [dB]','Fontsize',16,'Fontname','Times New Roman');
   set(gca,'Fontsize',16,'Fontname','Times New Roman');
   title(strcat(target_CDF,", ",target_SNR));
   grid on;
   hold on;
%% K�ɑ΂���ŗL�l���̂��̂̃O���t�iCDF�ɒ��ځj
   figure;
   mycol = [0 0 1;0 0 1;1 0 0;1 0 0;0 0.7 0;0 0.7 0;1 0 1;1 0 1;0 0 0;0 0 0]; % �O���t�̐F
   set(groot,'defaultAxesColorOrder',mycol) % figure�̏�����ݒ�
   axis([K_min K_max 0 7]);              % figure�̎��̂Ƃ����ݒ�
   grid on;
   hold on;
   set(gca,'XTick',K_min:5:K_max,'Fontsize',14,'Fontname','Times New Roman') % figure�̏�����ݒ�
   if Nr==2
      plot(K_box,dd3(:,1),'r-d','LineWidth',4);
      plot(K_box,dd3(:,2),'r--d','LineWidth',4);
      plot(K_box,dd3(:,3),'b-s','LineWidth',4);
      plot(K_box,dd3(:,4),'b--s','LineWidth',4);
      plot(K_box,dd3(:,5),'g-o','LineWidth',4);
      plot(K_box,dd3(:,6),'g--o','LineWidth',4);
      plot(K_box,dd3(:,7),'c-x','LineWidth',4);
      plot(K_box,dd3(:,8),'c--x','LineWidth',4);
      plot(K_box,dd3(:,9),'m-^','LineWidth',4);
      plot(K_box,dd3(:,10),'m--^','LineWidth',4);
      legend('BMSN-BF \lambda_1','BMSN-BF \lambda_2','BMSN-GE \lambda_1','BMSN-GE \lambda_2','MMSE-CI \lambda_1','MMSE-CI \lambda_2',...
       'BD \lambda_1','BD \lambda_2','ZF-CI \lambda_1','ZF-CI \lambda_2','Location','southeast');
   end
   if Nr==1
      plot(K_box,dd3(:,1),'r-d','LineWidth',4);
      plot(K_box,dd3(:,2),'b-s','LineWidth',4);
      plot(K_box,dd3(:,3),'g-o','LineWidth',4);
      plot(K_box,dd3(:,4),'c-x','LineWidth',4);
      plot(K_box,dd3(:,5),'m-^','LineWidth',4);
   end
   if Nr==1
   legend('BMSN-BF \lambda_1','BMSN-GE \lambda_1','MMSE-CI \lambda_1','BD \lambda_1','ZF-CI \lambda_1','Location','southeast');
   end
   xlabel('{\it{K}} [dB]','Fontsize',16,'Fontname','Times New Roman');
   ylabel('Eigenvalues','Fontsize',16,'Fontname','Times New Roman');
   set(gca,'Fontsize',16,'Fontname','Times New Roman');
   title(strcat(target_CDF,", ",target_SNR));
   grid on;
   hold on;
else
%% K�ɑ΂���Eigs1��Eigs2�̕��ς̃O���t�i�S���s�񐔕��ρj
if Nr==2
   figure;
   mycol = [0 0 1;1 0 0;0 0.7 0;1 0 1;0 0 0]; % �O���t�̐F
   set(groot,'defaultAxesColorOrder',mycol) % figure�̏�����ݒ�
   axis([K_min K_max 0 6]);              % figure�̎��̂Ƃ����ݒ�
   grid on;
   hold on;
   set(gca,'XTick',K_min:5:K_max,'Fontsize',14,'Fontname','Times New Roman') % figure�̏�����ݒ�
   plot(K_box,ff3(:,1),'r-d','LineWidth',4);
   plot(K_box,ff3(:,2),'b-s','LineWidth',4);
   plot(K_box,ff3(:,3),'g-o','LineWidth',4);
   plot(K_box,ff3(:,4),'c-x','LineWidth',4);
   plot(K_box,ff3(:,5),'m-^','LineWidth',4);

   legend('BMSN-BF \lambda_{Av}','BMSN-GE \lambda_{Av}','MMSE-CI \lambda_{Av}','BD \lambda_{Av}','ZF-CI \lambda_{Av}','Location','southeast');
   xlabel('{\it{K}} [dB]','Fontsize',36,'Fontname','Times New Roman');
   ylabel('Eigenvalues','Fontsize',36,'Fontname','Times New Roman');
   set(gca,'Fontsize',16,'Fontname','Times New Roman');
   title(target_SNR);
   grid on;
   hold on;
end
%% K�ɑ΂���Eigs1��Eigs2���ꂼ��̃O���t�i�S���s�񐔕��ρj
figure;
mycol = [0 0 1;1 0 0;0 0.7 0;1 0 1;0 0 0]; % �O���t�̐F
set(groot,'defaultAxesColorOrder',mycol) % figure�̏�����ݒ�
axis([K_min K_max 0 6]);              % figure�̎��̂Ƃ����ݒ�
grid on;
hold on;
set(gca,'XTick',K_min:5:K_max,'Fontsize',14,'Fontname','Times New Roman') % figure�̏�����ݒ�
if Nr==1
   plot(K_box,aa3(:,1),'r-d','LineWidth',4);
   plot(K_box,aa3(:,2),'b-s','LineWidth',4);
   plot(K_box,aa3(:,3),'g-o','LineWidth',4);
   plot(K_box,aa3(:,4),'c-x','LineWidth',4);
   plot(K_box,aa3(:,5),'m-^','LineWidth',4);
end
if Nr==2
   plot(K_box,aa3(:,1),'r-d','LineWidth',4);
   plot(K_box,aa3(:,2),'r--d','LineWidth',4);
   plot(K_box,aa3(:,3),'b-s','LineWidth',4);
   plot(K_box,aa3(:,4),'b--s','LineWidth',4);
   plot(K_box,aa3(:,5),'g-o','LineWidth',4);
   plot(K_box,aa3(:,6),'g--o','LineWidth',4);
   plot(K_box,aa3(:,7),'c-x','LineWidth',4);
   plot(K_box,aa3(:,8),'c--x','LineWidth',4);
   plot(K_box,aa3(:,9),'m-^','LineWidth',4);
   plot(K_box,aa3(:,10),'m--^','LineWidth',4);
end

if Nr==1
   legend('BMSN-BF \lambda_1','BMSN-GE \lambda_1','MMSE-CI \lambda_1','BD \lambda_1','ZF-CI \lambda_1','Location','southeast');
end
if Nr>=2
   legend('BMSN-BF \lambda_1','BMSN-BF \lambda_2','BMSN-GE \lambda_1','BMSN-GE \lambda_2','MMSE-CI \lambda_1','MMSE-CI \lambda_2',...
   'BD \lambda_1','BD \lambda_2','ZF-CI \lambda_1','ZF-CI \lambda_2','Location','southeast');
end
xlabel('{\it{K}} [dB]','Fontsize',36,'Fontname','Times New Roman');
ylabel('Eigenvalues','Fontsize',36,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
title(target_SNR);
grid on;
hold on;
%% K�ɑ΂���Eigs1[dB]��Eigs2[dB]���ꂼ��̃O���t�i�S���s�񐔕��ρj
figure;
mycol = [0 0 1;0 0 1;1 0 0;1 0 0;0 0.7 0;0 0.7 0;1 0 1;1 0 1;0 0 0;0 0 0]; % �O���t�̐F
set(groot,'defaultAxesColorOrder',mycol) % figure�̏�����ݒ�
axis([K_min K_max SNR_tar-50 SNR_tar+10]);              % figure�̎��̂Ƃ����ݒ�
grid on;
hold on;
set(gca,'XTick',K_min:5:K_max,'Fontsize',14,'Fontname','Times New Roman') % figure�̏�����ݒ�
if Nr==1
   plot(K_box,bb3(:,1),'r-d','LineWidth',4);
   plot(K_box,bb3(:,2),'b-s','LineWidth',4);
   plot(K_box,bb3(:,3),'g-o','LineWidth',4);
   plot(K_box,bb3(:,4),'c-x','LineWidth',4);
   plot(K_box,bb3(:,5),'m-^','LineWidth',4);
end
if Nr==2
   plot(K_box,bb3(:,1),'r-d','LineWidth',4);
   plot(K_box,bb3(:,2),'r--d','LineWidth',4);
   plot(K_box,bb3(:,3),'b-s','LineWidth',4);
   plot(K_box,bb3(:,4),'b--s','LineWidth',4);
   plot(K_box,bb3(:,5),'g-o','LineWidth',4);
   plot(K_box,bb3(:,6),'g--o','LineWidth',4);
   plot(K_box,bb3(:,7),'c-x','LineWidth',4);
   plot(K_box,bb3(:,8),'c--x','LineWidth',4);
   plot(K_box,bb3(:,9),'m-^','LineWidth',4);
   plot(K_box,bb3(:,10),'m--^','LineWidth',4);
end

if Nr==1
   legend('BMSN-BF \lambda_1','BMSN-GE \lambda_1','MMSE-CI \lambda_1','BD \lambda_1','ZF-CI \lambda_1','Location','southeast');
end
if Nr>=2
   legend('BMSN-BF \lambda_1','BMSN-BF \lambda_2','BMSN-GE \lambda_1','BMSN-GE \lambda_2','MMSE-CI \lambda_1','MMSE-CI \lambda_2',...
   'BD \lambda_1','BD \lambda_2','ZF-CI \lambda_1','ZF-CI \lambda_2','Location','southeast');
end
xlabel('{\it{K}} [dB]','Fontsize',16,'Fontname','Times New Roman');
ylabel('Eigenvalues [dB]','Fontsize',16,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
title(target_SNR);
grid on;
hold on;
%% K�ɑ΂���Eigs1[dB]��Eigs2[dB]�Ƃ̕��ς̃O���t�i�S���s�񐔕��ρj
figure;
mycol = [0 0 1;0 0 1;1 0 0;1 0 0;0 0.7 0;0 0.7 0;1 0 1;1 0 1;0 0 0;0 0 0]; % �O���t�̐F
set(groot,'defaultAxesColorOrder',mycol) % figure�̏�����ݒ�
axis([K_min K_max SNR_tar-50 SNR_tar+10]);              % figure�̎��̂Ƃ����ݒ�
grid on;
hold on;
set(gca,'XTick',K_min:5:K_max,'Fontsize',14,'Fontname','Times New Roman') % figure�̏�����ݒ�
plot(K_box,cc3(:,1),'r-d','LineWidth',4);
plot(K_box,cc3(:,2),'b-s','LineWidth',4);
plot(K_box,cc3(:,3),'g-o','LineWidth',4);
plot(K_box,cc3(:,4),'c-x','LineWidth',4);
plot(K_box,cc3(:,5),'m-^','LineWidth',4);

legend('BMSN-BF \lambda_{Av}','BMSN-GE \lambda_{Av}','MMSE-CI \lambda_{Av}','BD \lambda_{Av}','ZF-CI \lambda_{Av}','Location','southeast');
xlabel('{\it{K}} [dB]','Fontsize',36,'Fontname','Times New Roman');
ylabel('Eigenvalues[dB]','Fontsize',36,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
title(target_SNR);
grid on;
hold on;
end
%End