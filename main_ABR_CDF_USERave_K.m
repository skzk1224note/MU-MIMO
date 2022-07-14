% Achievable Bit Rate (ABR)

clear;
Nt = 16;    % ���M�f�q��
Nr = 2;     % ��M�f�q�� (�A���e�i�I���̏ꍇ1, �����łȂ��ꍇ��2)
Nu = 8;

SNR_tar = 30;  % �^�[�Q�b�gSNR[dB]
% SNR_max = 30; %�ő�SNR[dB]
SIMU = 1000; %���s�񐔁i�ʏ� 1000�j

% Rice�t�F�[�W���O�̂��߂̕ϐ���`
K_dB   = 0;         % Rician��K�t�@�N�^
K      = 10^(K_dB/10);
% K = 0(K_dB = -inf); % ���C���[�t�F�[�W���O  
d_t      = 0.5;      % ���M�A���e�i�Ԋu�iin wavelength)
d_r      = 0.5;      % ��M�A���e�i�Ԋu�iin wavelength)
derad = pi/180;      % degree -> rad

if SNR_tar < 10
    tsnr=strcat('SNR= ',num2str(SNR_tar,'%01d'),'dB');
else
    tsnr=strcat('SNR=',num2str(SNR_tar,'%02d'),'dB');
end

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

  
sigma2 = 1/(10^(SNR_tar/10)); % noise power
a = sigma2*Nt;

% H (�`���`���l���s��:Rician channel)
% �`���`���l���s��̃}���`�p�X���� (i.i.d. Rayleigh , NLOS �`���l��)
    H_iid = (randn(Nr*Nu,Nt)+1j*randn(Nr*Nu,Nt))/sqrt(2);
% �`���`���l���s��̒��ڔg����(LOS �`���l��)
    H_los = zeros(Nu*Nr,Nt);

% �o�̓t�@�C���� with alpha and SNR in dB
%evfile1 = strcat(folder,fn1,num2str(snr(isnr),'%02d'),'dB_1000itr.csv');
    
for isimu = 1:SIMU

     for dR=1:Nu
        Theta_r(1,dR) = 30*dR-180; % ���[�U���̎�M�p (���x�Ԋu�Ń��[�U��z�u)      
     end
     
     Theta_t = 0;
%    Theta_t = (rand-0.5)*360;   % ���[�U���̑��M�p (-180deg - 180deg)
    %Theta_r = (rand(1,Nu)-0.5)*360; % ���[�U���̎�M�p (-180deg - 180deg)
    a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t*derad));  % ���M���[�h�x�N�g��
    
      for n = 1 : Nu
        a_r = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad)); % ���[�U���̎�M���[�h�x�N�g��
        H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t.';                % ���[�U����LOS�`���l���s��
      end

    % �`���`���l���s��=[sqrt(K/(K+1))*(LOS �`���l��)]...
    %                   .+[sqrt(1/(K+1))*(NLOS �`���l��)]
    H = sqrt(K/(K+1))*H_los + sqrt(1/(K+1))*H_iid;
    Hu = peruser(H,Nu);

    He = zeros((Nu-1)*Nr,Nt,Nu);    % H����1���[�U�̃`���l���s����������s��

    % BMSN3
    [~,~,STT,RIP,~] = bmsn_bf(Nt,Nr,Nu,H,a,T);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %��M�A���e�i�inr)�̌ŗL�l���z for B-MSN�i���g�������l���j
    MSt_BMSN(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
 
    % ZF
    [~,~,STT,RIP,~] = zf(Nt,Nr,Nu,H);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %��M�A���e�i�inr)�̌ŗL�l���z for B-MSN�i���g�������l���j
    MSt_ZF(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    
    % MMSE
    [~,~,STT,RIP,~] = gmmse(Nt,Nr,Nu,H,a);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %��M�A���e�i�inr)�̌ŗL�l���z for B-MSN�i���g�������l���j
    MSt_MMSE(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    
    % BD
    [~,~,STT,RIP,~] = bd(Nt,Nr,Nu,H);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %��M�A���e�i�inr)�̌ŗL�l���z for B-MSN�i���g�������l���j
    MSt_BD(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
        
    %fprintf('Iteration = %d / %d\n',isimu, SIMU);
    
end % isimu

% ABR
for nuser=1:Nu
    Q(:,1,nuser)=sort(sum(log2(1 + MSt_BMSN(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
    Q(:,2,nuser)=sort(sum(log2(1 + MSt_ZF(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
    Q(:,3,nuser)=sort(sum(log2(1 + MSt_MMSE(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
    Q(:,4,nuser)=sort(sum(log2(1 + MSt_BD(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
end

Qm = mean(Q,3);
%ABR(isnr,ia)=Qm(ICDF,ia);

Y(:,1) = [0.1:0.1:100].';

figure;
mycol = [0 0 1;
      1 0 0;
      0 0.7 0;
      1 0 1;
      1 0 0;0 0 0];
set(groot,'defaultAxesColorOrder',mycol)
plot(Qm,Y,'Linewidth',2);
set(gca,'Fontsize',14,'Fontname','Times New Roman');
legend('BMSN','ZF','MMSE','BD','Location','Southeast');
xlabel('Achievable bit rate [bits/s/Hz]','Fontsize',16,'Fontname','Times New Roman');
ylabel('CDF [%]','Fontsize',16,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
title(tsnr);
grid on;
hold on;

%End


