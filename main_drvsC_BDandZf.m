% cap_bd_muser.m
% NU-���[�U�ł�BD�@�̌v�Z
% TDMA, Upper bound�Ɣ�r

clear;%close;


% �p�����[�^���� NT >= NR*NU
K_tar = 20;      % �^�[�Q�b�gK[dB]

CDF = 50;         % ABR�̃^�[�Q�b�gCDF
SNR_tar = 10;     % �^�[�Q�b�gSNR[dB]
Nt = 16;          % ���M�f�q��
Nr = 2;           % ��M�f�q�� (�A���e�i�I���̏ꍇ1, �����łȂ��ꍇ��2)
Nu = 8;           % ���[�U��
SIMU = 1000;      % ���s�񐔁i�ʏ� 1000�j

An = 2; % �w�����֐��̌W��(cos�p�^�[���̏ꍇ��2)

I = eye(Nt,Nt);
Nru = Nr*Nu;


dr_min = 0.5;
dr_max = 1.0;
dr_box = (dr_min:0.05:dr_max); % ���M�A���e�i�Ԋu�iin wavelength)
Ldr=length(dr_box);        % dt_box�̑傫��

d_t = 0.5;        % ��M�A���e�i�Ԋu�iin wavelength)
derad = pi/180;   % degree -> rad



%T ���]�̃`���l���s��
T = zeros(Nr,Nr);
for inu = 1:Nu
    T(:,:,inu) = eye(Nr,Nr);
end


Q = zeros(SIMU,2,Nu);
QmC = zeros(Ldr,2);
QmC_cos = zeros(Ldr,2);

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

if K_tar < 10
   target_K=strcat('{\it{K}}=',num2str(K_tar,'%01d'),'dB');
else
   target_K=strcat('{\it{K}}=',num2str(K_tar,'%02d'),'dB');
end

if d_t < 10
   target_d_t=strcat('d_t=',num2str(d_t,'%.3g'),'\lambda');
else
   target_d_t=strcat('d_t=',num2str(d_t,'%02d'),'\lambda');
end
% �o�̓t�@�C���� with SNR in dB
%cdffile1 = strcat(folder,cdfn1,num2str(CDF,'%02d'),'_1000itr.csv');
%cdffile2 = strcat(folder,cdfn2,num2str(SN_tar,'%02d'),'dB_1000itr.csv');
%%
H_los = zeros(Nu*Nr,Nt); % �`���`���l���s��̒��ڔg����(LOS �`���l��)
sigma2 = 1/(10^(SNR_tar/10)); % noise power
a = sigma2*Nt; % �[���G��

St = zeros(Nr,Nu);
% MSt_BMSN_BF = zeros(SIMU, Nr*Nu);
% MSt_BMSN_GE = zeros(SIMU, Nr*Nu);
% MSt_MMSE = zeros(SIMU, Nr*Nu);
% MSt_GMI1 = zeros(SIMU, Nr*Nu);
% MSt_GMI2 = zeros(SIMU, Nr*Nu);
MSt_BD = zeros(SIMU, Nr*Nu);
MSt_ZF = zeros(SIMU, Nr*Nu);

for Directivity_switch = 0:1 % ����M�f�q�̎w�����l���̗L�� 0:��,1:�L
    for idr = 1:Ldr 
      dr_tar = dr_box(idr);
        for isimu = 1:SIMU      
          % �`���`���l���s��̃}���`�p�X���� (i.i.d.Rayleigh)
          H_iid = (randn(Nr*Nu,Nt)+1j*randn(Nr*Nu,Nt))/sqrt(2);
          % LOS �`���l��      
            if Directivity_switch == 1 % ����M�f�q�̎w�����l���L��
             Theta_t = (rand(1,Nu)-0.5)*180; % ���[�U���̑��M�p �w����:(-90deg - 90deg)
             Theta_r = (rand(1,Nu)-0.5)*180; % ���[�U���̎�M�p �w����:(-90deg - 90deg)
             for n = 1 : Nu
                a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad))*An*cos(Theta_t(1,n)*derad); % ���[�U���̑��M���[�h�x�N�g��
                a_r = exp(-1j*2*pi*dr_tar*(0:Nr-1).'*sin(Theta_r(1,n)*derad))*An*cos(Theta_r(1,n)*derad); % ���[�U���̎�M���[�h�x�N�g��
                H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t';                 % ���[�U����LOS�`���l���s��
             end
            elseif Directivity_switch == 0 % ����M�f�q�̎w�����l������
              Theta_t = (rand(1,Nu)-0.5)*360; % ���[�U���̑��M�p ������:(-180deg - 180deg) 
              Theta_r = (rand(1,Nu)-0.5)*360; % ���[�U���̎�M�p ������:(-180deg - 180deg) 
             for n = 1 : Nu
             a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad)); % ���[�U���̑��M���[�h�x�N�g��
             a_r = exp(-1j*2*pi*dr_tar*(0:Nr-1).'*sin(Theta_r(1,n)*derad)); % ���[�U���̎�M���[�h�x�N�g��
             H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t';                 % ���[�U����LOS�`���l���s��
             end
            end
            % K��^�l�ɂ���
            K = 10^(K_tar/10);
            
            % H=[sqrt(K/(K+1))*(LOS �`���l��)]+[sqrt(1/(K+1))*(NLOS �`���l��)]
            H = sqrt(K/(K+1))*H_los + sqrt(1/(K+1))*H_iid;
            
%             % BMSN-BF
%             [~,~,STT,RIP,~] = bmsn_bf(Nt,Nr,Nu,H,a,T);
%             for inu=1:Nu
%                 St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
%             end
%             %��M�A���e�i�inr)�̌ŗL�l���z for B-MSN�i���g�������l���j
%             MSt_BMSN_BF(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
%             
%             % BMSN-GE
%             [~,~,STT,RIP,~] = bmsn_gev(Nt,Nr,Nu,H,a);
%             for inu=1:Nu
%                 St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
%             end
%             %��M�A���e�i�inr)�̌ŗL�l���z for B-MSN�i���g�������l���j
%             MSt_BMSN_GE(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
%             
%             % MMSE-CI
%             [~,~,STT,RIP,~] = mmse(Nt,Nr,Nu,H,a);
%             for inu=1:Nu
%                 St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
%             end
%             %��M�A���e�i�inr)�̌ŗL�l���z for B-MSN�i���g�������l���j
%             MSt_MMSE(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
%             
%             % GMI1
%             [~,~,STT,RIP,~] = gmmse_m1(Nt,Nr,Nu,H,a);
%             for inu=1:Nu
%                 St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
%             end
%             %��M�A���e�i�inr)�̌ŗL�l���z for B-MSN�i���g�������l���j
%             MSt_GMI1(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
%             
%             % GMI2
%             [~,~,STT,RIP,~] = gmmse_m2(Nt,Nr,Nu,H,a);
%             for inu=1:Nu
%                 St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
%             end
%             %��M�A���e�i�inr)�̌ŗL�l���z for B-MSN�i���g�������l���j
%             MSt_GMI2(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
            
            % BD
                [~,~,STT,RIP,~] = bd(Nt,Nr,Nu,H);
                for inu=1:Nu
                    St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
                end
                %��M�A���e�i�inr)�̌ŗL�l���z for B-MSN�i���g�������l���j
                MSt_BD(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
            
            % ZF-CI
                [~,~,STT,RIP,~] = zf(Nt,Nr,Nu,H);
                for inu=1:Nu
                    St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
                end
                %��M�A���e�i�inr)�̌ŗL�l���z for B-MSN�i���g�������l���j
                MSt_ZF(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
            
            % fprintf('Iteration = %d / %d\n',isimu, SIMU);
            
        end % isimu
        
        % �e�A���S���Y����ABR
        for nuser=1:Nu
            ns = Nr*(nuser-1)+1:Nr*nuser;
%             Q(:,1,nuser)=sort(sum(log2(1 + MSt_BMSN_BF(:,ns)),2));
%             Q(:,2,nuser)=sort(sum(log2(1 + MSt_BMSN_GE(:,ns)),2));
%             Q(:,3,nuser)=sort(sum(log2(1 + MSt_MMSE(:,ns)),2));
%             Q(:,4,nuser)=sort(sum(log2(1 + MSt_GMI1(:,ns)),2));
%             Q(:,5,nuser)=sort(sum(log2(1 + MSt_GMI2(:,ns)),2));
            Q(:,1,nuser)=sort(sum(log2(1 + MSt_BD(:,ns)),2));
            Q(:,2,nuser)=sort(sum(log2(1 + MSt_ZF(:,ns)),2));
        end % nuser end
        
        
        
        Qm = mean(Q,3); % Q�̃��[�U�񐔕���    
        if Directivity_switch == 0 % �w���������̏ꍇ�̍s��
            QmC(idr,:)=Qm(round(CDF*SIMU/100),:);% Channel capacity
        elseif Directivity_switch == 1 % �w�����L��̏ꍇ�̍s��
            QmC_cos(idr,:)=Qm(round(CDF*SIMU/100),:);% Channel capacity
        end
        fprintf('dr = %.3g lambda \n',dr_box(idr));
    end % idr end
    %csvwrite(cdffile1,[dr,QmC]);
end
%% �O���t�\�� figure ABRvsK
figure;
mycol = [1 0 1;0 1 0;1 0 0;0 0 1;0 1 1;0 0 0]; % �F
set(groot,'defaultAxesColorOrder',mycol)
%plot(K,QmC(:,1),K,QmC(:,2),K,QmC(:,3),K,QmC(:,4),'Linewidth',2);

plot(dr_box,QmC(:,1),'c--',dr_box,QmC(:,2),'m--',dr_box,QmC_cos(:,1),'c-',dr_box,QmC_cos(:,2),'m-','Linewidth',2);
axis([dr_min,dr_max,min(min(QmC))-0.05,max(max(QmC_cos))+1])
set(gca,'Fontsize',18,'Fontname','Times New Roman');
lgd = legend;
lgd.NumColumns = 2; % �}��̗񐔂��w��
legend('BD','ZF','BDcos','ZFcos','Location','Northwest');
xlabel('{\it{d_r}} [\lambda]','Fontsize',6,'Fontname','Times New Roman');
ylabel('Channel capacity [bits/s/Hz]','Fontsize',6,'Fontname','Times New Roman');
set(gca,'Fontsize',18,'Fontname','Times New Roman');
title(strcat(target_CDF,',',target_SNR,',',target_K,',',target_d_t));
grid on;
hold on;
% End