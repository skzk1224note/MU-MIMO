% cap_bd_muser.m
% NU-���[�U�ł�BD�@�̌v�Z
% TDMA, Upper bound�Ɣ�r

clear;

% �o�̓t�@�C����
% testfile1 = 'Capacity2x8x4u_BD.csv';
% testfile2 = 'Capacity2x8x4u_BD_CDFSNR20dB.csv';
% testfile3 = 'Eig2x8x4u_BD_CDFSNR20dB.csv';

% �p�����[�^���� NT >= NR*NU

SN_tar = 10;         % CDF�\���̂��߂̃^�[�Q�b�gSNR [dB]
SIMU   = 100;       % ���s��
Nt     = 16;         % ���M�f�q��
Nr     = 2;          % ��M�f�q��(=2�ɌŒ�)
Nu     = 8;          % ���[�U��
NL     = 20;          % �}���`�p�X�g�̑f�g��
I      = eye(Nt,Nt); % NTxNT�̒P�ʍs��

K_dB   = -20;
K = 10^(K_dB/10);    % K�̐^�l
 
d_t = 0.5;      % ���M�A���e�i�Ԋu�iin wavelength)
d_r = 0.5;      % ��M�A���e�i�Ԋu�iin wavelength)
derad = pi/180;      % degree -> rad

% cos�Ƃ�n��in�͐����l�j
n = 0;
if rem(n, 2) == 0 % n�������̂Ƃ�
    a = 2:2:n;
    aa = prod(a);
    b = 1:2:n+1;
    bb = prod(b);
    D = bb/aa;
else              % n����̂Ƃ�
    a = 1:2:n;
    aa = prod(a);
    b = 2:2:n+1;
    bb = prod(b);
    D = (2*bb) / (pi*aa);
end

An = sqrt(D);
%An = 2; % �w�����֐��̌W��
figure_switch = 0;   % figure�\���p�̃X�C�b�` 1:on, 0:off�ion�Ȃ�ŗL�l���̂��̂̃O���t��\���j

a = Nt/(10^(SN_tar/10)); % �[���G��
        
T = zeros(Nr,Nr); % ���]�̃`���l���s�� for BMSN
for nuser = 1:Nu
    T(:,:,nuser) = eye(Nr,Nr);
end
Algorithms = ["BMSN-BF","BMSN-GE","MMSE","GMI1","GMI2","BD","ZF"]; %,"BMSN-GE3"
S = zeros(Nr, Nr, Nu, numel(Algorithms));
E = zeros(SIMU, Nr, Nu, numel(Algorithms));
Eigs = zeros(SIMU, Nr, Nu, numel(Algorithms));
EigsAll = zeros(SIMU, Nr, Nu, numel(Algorithms));
a_t_iid = zeros(Nt, 1, NL);
a_r_iid = zeros(Nr, 1, NL);
X = zeros(Nr, Nt, NL);
if SN_tar < 10
        target_snr=strcat('SNR=',num2str(SN_tar,'%01d'),'dB');
else
        target_snr=strcat('SNR=',num2str(SN_tar,'%02d'),'dB');
end    

if K_dB < 10
   target_K=strcat('K=',num2str(K_dB,'%01d'),'dB');
else
   target_K=strcat('K=',num2str(K_dB,'%02d'),'dB');
end
%�`���`���l���s��̒��ڔg����(LOS �`���l��)
H_los = zeros(Nu*Nr,Nt);

for Directivity_switch = 0:1 % ����M�f�q�̎w�����l���̗L�� 0:��,1:�L
    for isimu = 1:SIMU % ���s�񐔂̃��[�v
        if Directivity_switch == 1 % ����M�f�q�̎w�����l���L��
            Theta_t = (rand(1,Nu)-0.5)*180; % ���[�U���̑��M�p �w����:(-90deg - 90deg)
            Theta_r = (rand(1,Nu)-0.5)*180; % ���[�U���̎�M�p �w����:(-90deg - 90deg)
            % LoS �`���l��
            for n = 1 : Nu
                a_t = exp(1j*Theta_t(1,n)*derad) * exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad)); % ���[�U���̑��M���[�h�x�N�g�� An*cos(Theta_t(1,n)*derad) *
                a_r = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad)); % ���[�U���̎�M���[�h�x�N�g�� An*cos(Theta_r(1,n)*derad) *
                H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t'; % ���[�U����LOS�`���l���s��
                % �􉽊wNLoS�`���l��
                for L = 1:NL
                    Thetat = (rand(1,1)-0.5)*180; % ���[�U���̑��M�p �w����:(-90deg - 90deg)An*cos(Thetat*derad)
                    Thetar = (rand(1,1)-0.5)*180; % ���[�U���̑��M�p �w����:(-90deg - 90deg)An*cos(Thetar*derad)
                    gtheta_t = An * ((cos(Thetat*derad))^(n/2));
                    gtheta_r = An * ((cos(Thetat*derad))^(n/2));
                    a_t_iid(:,:,L) = gtheta_t * sqrt(1/Nt) * exp(1j*Thetat*derad) * exp(-1j*2*pi*d_t*(0:Nt-1).' * sin(Thetat*derad)); % ���[�U����NLoS���M���[�h�x�N�g�� An*cos(Theta_t(1,n)*derad) *
                    a_r_iid(:,:,L) = gtheta_r * sqrt(1/Nr) * exp(-1j*2*pi*d_r*(0:Nr-1).' * sin(Thetar*derad)); % ���[�U����NLoS��M���[�h�x�N�g�� An*cos(Theta_r(1,n)*derad) *
                    X(:,:,L) = (a_r_iid(:,:,L) * a_t_iid(:,:,L)');
                end
                H_iid((n-1)*Nr+1:(n-1)*Nr+Nr,:) = sum(X,3) * sqrt(Nt*Nr/NL);
%                 g_theta_t_iid = 1; % An*cos(Theta_t_iid*derad);
%                 g_theta_r_iid = 1; % An*cos(Theta_r_iid*derad);
%                 a_t_iid = a_t_iid.* g_theta_t_iid';
%                 a_r_iid = a_r_iid.* g_theta_r_iid';
                % H_iid((n-1)*Nr+1:(n-1)*Nr+Nr,:) = (a_r_iid*a_t_iid'); %*sqrt(1/NL); % ���[�U����NLOS�`���l���s��
            end
        elseif Directivity_switch == 0 % ����M�f�q�̎w�����l������
            H_iid = (randn(Nr*Nu,Nt)+1j*randn(Nr*Nu,Nt))/sqrt(2); % �`���`���l���s��̃}���`�p�X���� (i.i.d. Rayleigh , NLOS �`���l��)
            Theta_t = (rand(1,Nu)-0.5)*180; % ���[�U���̑��M�p ������:(-180deg - 180deg) 
            Theta_r = (rand(1,Nu)-0.5)*180; % ���[�U���̎�M�p ������:(-180deg - 180deg) 
            Theta_iid = (rand(Nr*Nu,Nt)-0.5)*180; % ���[�U���̑��M�p �w����:(-90deg - 90deg)
            g_theta_iid = 1; %;An*cos(Theta_iid*derad)
            H_iid = H_iid * g_theta_iid;
            for n = 1 : Nu
                a_t = sqrt(1/Nt) * exp(1j*Theta_t(1,n)*derad) * exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad)); % ���[�U���̑��M���[�h�x�N�g��
                a_r = sqrt(1/Nr) * exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad)); % ���[�U���̎�M���[�h�x�N�g��
                H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t' * sqrt(Nt*Nr);                 % ���[�U����LOS�`���l���s��
            end
            % Normiso = norm(H_iid,'fro');
        end
        
        % �`���`���l���s��̃}���`�p�X���� (i.i.d. Rayleigh, NLOS �`���l��)
%         H_iid = (randn(Nr*Nu,Nt)+1j*randn(Nr*Nu,Nt))/sqrt(2);
        % �`���`���l���s��=[sqrt(K/(K+1))*(LOS �`���l��)]+[sqrt(1/(K+1))*(NLOS �`���l��)]
        H0 = sqrt(K/(K+1))*H_los + sqrt(1/(K+1))*H_iid;
        (norm(H_los,'fro'))^2/(norm(H_iid,'fro'))^2;
        %H0 = squeeze(H(k,:,:)); % k�Ԗڂ̎��s�񐔂ł̓`���`���l���s��
        
        % BMSN-BF algorithm
        [W_BMSN_BF,U_BMSN_BF,S_BMSN_BF,RIPBF,~] = bmsn_bf(Nt,Nr,Nu,H0,a,T); % function bmsn_bf.m ���g�p
        
        % BMSN-GE algorithm
        [W_BMSN_GE,U_BMSN_GE,S_BMSN_GE,RIPGE,~] = bmsn_ge4(Nt,Nr,Nu,H0,a); % function bmsn_ge4.m ���g�p
        
        % MMSE-CI algorithm
        [W_MMSE,U_MMSE,S_MMSE,RIPM,~] = mmse(Nt,Nr,Nu,H0,a); % function mmse.m ���g�p
        
        % GMI1 algorithm
        [W_GMI1,U_GMI1,S_GMI1,RIPGM1,~] = gmmse_m1(Nt,Nr,Nu,H0,a); % function mmse.m ���g�p
        
        % GMI2 algorithm
        [W_GMI2,U_GMI2,S_GMI2,RIPGM2,~] = gmmse_m2(Nt,Nr,Nu,H0,a); % function mmse.m ���g�p
        
        % BD algorithm
        [W_BD,U_BD,S_BD,~,~] = bd(Nt,Nr,Nu,H0); % function bd.m ���g�p
        
        % ZF-CI algorithm
        [W_ZF,U_ZF,S_ZF,~,~] = zf(Nt,Nr,Nu,H0); % function zf.m ���g�p
        
        S(:,:,:,1) = S_BMSN_BF; S(:,:,:,2) = S_BMSN_GE; S(:,:,:,3) = S_MMSE; S(:,:,:,4) = S_GMI1; 
        S(:,:,:,5) = S_GMI2;S(:,:,:,7) = S_ZF;
        if Nr == 1
            S(:,:,:,6) = S_BD(1,1,:);
        else
            S(:,:,:,6) = S_BD;
        end
        % ���[�U���̌ŗL�l���z
        snt = 1/(10^(SN_tar/10));
        for nuser=1:Nu
            if Nr==1
                E(isimu,:,nuser,strcmp("BMSN-BF",Algorithms)) = 10*log10((S_BMSN_BF(1,1,nuser).^2)./(RIPBF(1,nuser)+Nt*snt));
                E(isimu,:,nuser,strcmp("BMSN-GE",Algorithms)) = 10*log10((S_BMSN_GE(1,1,nuser).^2)./(RIPGE(1,nuser)+Nt*snt));
                E(isimu,:,nuser,strcmp("MMSE",Algorithms)) = 10*log10((S_MMSE(1,1,nuser).^2)./(RIPM(1,nuser)+Nt*snt));
                E(isimu,:,nuser,strcmp("GMI1",Algorithms)) = 10*log10((S_GMI1(1,1,nuser).^2)./(RIPGM1(1,nuser)+Nt*snt));
                E(isimu,:,nuser,strcmp("GMI2",Algorithms)) = 10*log10((S_GMI2(1,1,nuser).^2)./(RIPGM2(1,nuser)+Nt*snt));
                E(isimu,:,nuser,strcmp("BD",Algorithms)) = 10*log10(S_BD(1,1,nuser).^2/(Nt*snt));
                E(isimu,:,nuser,strcmp("ZF",Algorithms)) = 10*log10(S_ZF(1,1,nuser).^2/(Nt*snt));    
            else
                E(isimu,:,nuser,strcmp("BMSN-BF",Algorithms)) = 10*log10((diag(S_BMSN_BF(:,:,nuser)).^2)./(RIPBF(:,nuser)+Nt*snt));
                E(isimu,:,nuser,strcmp("BMSN-GE",Algorithms)) = 10*log10((diag(S_BMSN_GE(:,:,nuser)).^2)./(RIPGE(:,nuser)+Nt*snt));
                E(isimu,:,nuser,strcmp("MMSE",Algorithms)) = 10*log10((diag(S_MMSE(:,:,nuser)).^2)./(RIPM(:,nuser)+Nt*snt));
                E(isimu,:,nuser,strcmp("GMI1",Algorithms)) = 10*log10((diag(S_GMI1(:,:,nuser)).^2)./(RIPGM1(:,nuser)+Nt*snt));
                E(isimu,:,nuser,strcmp("GMI2",Algorithms)) = 10*log10((diag(S_GMI2(:,:,nuser)).^2)./(RIPGM2(:,nuser)+Nt*snt));
                E(isimu,:,nuser,strcmp("BD",Algorithms)) = 10*log10(diag(S_BD(:,:,nuser)).^2/(Nt*snt));
                E(isimu,:,nuser,strcmp("ZF",Algorithms)) = 10*log10(diag(S_ZF(:,:,nuser)).^2/(Nt*snt));
            end
        end
        
        % ��������ŗL�l���̂��̂̒l�����
        for nuser=1:Nu
            if Nr==1
                for nalg = 1:numel(Algorithms)
                    Eigs(isimu,:,nuser,nalg) = S(1,1,nuser,nalg);
                end
            else
                for nalg = 1:numel(Algorithms)
                    Eigs(isimu,:,nuser,nalg) = diag(S(:,:,nuser,nalg));
                end
            end % Nr end
        end % nuser end
    end
    for nuser=1:Nu
        E(:,:,nuser,:) = sort(E(:,:,nuser,:),1);
        Eigs(:,:,nuser,:) = sort(Eigs(:,:,nuser,:),1);
    end
    % ���[�U����
    E = mean(E,3);
    Eigs = mean(Eigs,3);
    %�����Eigs1��Eigs2�̕���
    EigsAll = mean(Eigs,2);
    
    % CDF of Eigenvalue at Target SNR
    Y = (1/SIMU:1/SIMU:1).'*100;
    if Directivity_switch == 1 % ����M�f�q�̎w�����l���̗L
        rr3_cos = zeros(SIMU,Nr*numel(Algorithms));
        ee3_cos = zeros(SIMU,Nr*numel(Algorithms));
        for nalg = 1:numel(Algorithms)
            for nn=1:Nr            
                rr3_cos(:,nn+(nalg-1)*Nr) = E(:,nn,nalg);
                ee3_cos(:,nn+(nalg-1)*Nr) = Eigs(:,nn,nalg);
            end
        end
        ss3_cos = EigsAll;
        
    else % ����M�f�q�̎w�����l���̖�
        rr3 = zeros(SIMU,Nr*numel(Algorithms));
        ee3 = zeros(SIMU,Nr*numel(Algorithms)); 
        for nalg = 1:numel(Algorithms)
            for nn=1:Nr            
                rr3(:,nn+(nalg-1)*Nr) = E(:,nn,nalg);
                ee3(:,nn+(nalg-1)*Nr) = Eigs(:,nn,nalg);
            end
        end
        ss3 = EigsAll;
    end
end
%% CDF of EGV at Target SNR
figure;
mycol = [1 0 0;0 0 1;0 1 0;1 1 0;0 0 0;0 1 1;1 0 1
        1 0 0;0 0 1;0 1 0;1 1 0;0 0 0;0 1 1;1 0 1]; % �F
set(groot,'defaultAxesColorOrder',mycol)
Holizon_min = round(min(min(rr3)))-5;
Holizon_max = round(max(max(rr3_cos)))+5;
axis([Holizon_min Holizon_max 0 100]);
grid on;
hold on;
set(gca,'XTick',Holizon_min-5:5:Holizon_max+5,'Fontsize',14,'Fontname','Times New Roman')
xlabel('SINR of eigenvalue [dB]','Fontsize',16,'Fontname','Times New Roman');
ylabel('CDF [%]','Fontsize',16,'Fontname','Times New Roman');
if Nr==1
    plot(rr3(:,1),Y,'r-','Linewidth',2); plot(rr3(:,2),Y,'b-','Linewidth',2);
    plot(rr3(:,3),Y,'g-','Linewidth',2); plot(rr3(:,4),Y,'y-','Linewidth',2);
    plot(rr3(:,5),Y,'k-','Linewidth',2); plot(rr3(:,6),Y,'c-','Linewidth',2);
    plot(rr3(:,7),Y,'m-','Linewidth',2);
    plot(rr3_cos(:,1),Y,'r-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,2),Y,'b-o','MarkerIndices',50:100:length(Y),'Linewidth',2); 
    plot(rr3_cos(:,3),Y,'g-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,4),Y,'y-o','MarkerIndices',50:100:length(Y),'Linewidth',2); 
    plot(rr3_cos(:,5),Y,'k-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,6),Y,'c-o','MarkerIndices',50:100:length(Y),'Linewidth',2);
    plot(rr3_cos(:,7),Y,'m-o','MarkerIndices',10:100:length(Y),'Linewidth',2);
    legend('BMSN-BF \lambda_1','BMSN-GE \lambda_1','MMSE-CI \lambda_1','GMI1 \lambda_1','GMI2 \lambda_1','BD \lambda_1','ZF-CI \lambda_1',...
        'BMSN-BFcos \lambda_1','BMSN-GEcos \lambda_1','MMSE-CIcos \lambda_1','GMI1cos \lambda_1','GMI2cos \lambda_1','BDcos \lambda_1','ZF-CIcos \lambda_1','Location','southeast');
end

if Nr==2
    plot(rr3(:,1),Y,'r-','Linewidth',2); plot(rr3(:,2),Y,'r--','Linewidth',2);
    plot(rr3(:,3),Y,'b-','Linewidth',2); plot(rr3(:,4),Y,'b--','Linewidth',2);
    plot(rr3(:,5),Y,'g-','Linewidth',2); plot(rr3(:,6),Y,'g--','Linewidth',2);
    plot(rr3(:,7),Y,'y-','Linewidth',2); plot(rr3(:,8),Y,'Y--','Linewidth',2);
    plot(rr3(:,9),Y,'k-','Linewidth',2); plot(rr3(:,10),Y,'k--','Linewidth',2);
    plot(rr3(:,11),Y,'c-','Linewidth',2); plot(rr3(:,12),Y,'c--','Linewidth',2);
    plot(rr3(:,13),Y,'m-','Linewidth',2); plot(rr3(:,14),Y,'m--','Linewidth',2);
    
    plot(rr3_cos(:,1),Y,'r-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,2),Y,'r--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
    plot(rr3_cos(:,3),Y,'b-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,4),Y,'b--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
    plot(rr3_cos(:,5),Y,'g-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,6),Y,'g--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
    plot(rr3_cos(:,7),Y,'y-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,8),Y,'Y--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
    plot(rr3_cos(:,9),Y,'k-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,10),Y,'k--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
    plot(rr3_cos(:,11),Y,'c-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,12),Y,'c--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
    plot(rr3_cos(:,13),Y,'m-o','MarkerIndices',10:100:length(Y),'Linewidth',2); plot(rr3_cos(:,14),Y,'m--o','MarkerIndices',50:100:length(Y),'Linewidth',2);
    lgd = legend;
    lgd.NumColumns = 2; % �}��̗񐔂��w��
    legend('BMSN-BF \lambda_1','BMSN-BF \lambda_2','BMSN-GE \lambda_1','BMSN-GE \lambda_2','MMSE-CI \lambda_1','MMSE-CI \lambda_2',...
        'GMI1 \lambda_1','GMI1 \lambda_2','GMI2 \lambda_1','GMI2 \lambda_2','BD \lambda_1','BD \lambda_2','ZF-CI \lambda_1','ZF-CI \lambda_2',...
        'BMSN-BFcos \lambda_1','BMSN-BFcos \lambda_2','BMSN-GEcos \lambda_1','BMSN-GEcos \lambda_2','MMSE-CIcos \lambda_1','MMSE-CIcos \lambda_2',...
        'GMI1cos \lambda_1','GMI1cos \lambda_2','GMI2cos \lambda_1','GMI2cos \lambda_2','BDcos \lambda_1','BDcos \lambda_2',...
        'ZF-CIcos \lambda_1','ZF-CIcos \lambda_2','Location','southeast');
end
title(strcat(target_snr,',',target_K));
%% CDF of EGV_nonSINR at Target SNR 
if figure_switch == 1
    figure;
    if Nr==2
        mycol = [1 0 0;0 0 1;0 1 0;1 1 0;0 0 0;0 1 1;1 0 1
            1 0 0;0 0 1;0 1 0;1 1 0;0 0 0;0 1 1;1 0 1]; % �F
        set(groot,'defaultAxesColorOrder',mycol)        
        for nalg = 1:2:numel(Algorithms)*2-1
            plot(ee3(:,nalg),Y,'-','Linewidth',2);
            grid on; hold on;
        end
        for nalg = 2:2:numel(Algorithms)*2
            plot(ee3(:,nalg),Y,'--','Linewidth',2);
            grid on; hold on;
        end
        lgd = legend;
        lgd.NumColumns = 2; % �}��̗񐔂��w��
        legend([strcat(Algorithms," \lambda_1") strcat(Algorithms," \lambda_2")],'Location','southeast');
    end
    if Nr==1        
        mycol = [1 0 0;0 0 1;0 1 0;1 1 0;0 0 0;0 1 1;1 0 1]; % �F
        set(groot,'defaultAxesColorOrder',mycol)
        for nalg = 1:numel(Algorithms)
            plot(ee3(:,nalg),Y,'-','Linewidth',2);
            grid on; hold on;
        end
        legend(strcat(Algorithms," \lambda_1"),'Location','southeast');
    end
    axis([0 round(max(max(ee3)))+3 0 100]);
    set(gca,'XTick',0:1:round(max(max(ee3)))+3,'Fontsize',14,'Fontname','Times New Roman')
    xlabel('Eigenvalue','Fontsize',16,'Fontname','Times New Roman');
    ylabel('CDF [%]','Fontsize',16,'Fontname','Times New Roman');
    title(strcat(target_snr,',',target_K));
    %% CDF of EGV at Target SNR 
    if Nr==2
        figure;
        mycol = [1 0 0;0 0 1;0 1 0;1 1 0;0 0 0;0 1 1;1 0 1]; % �F
        set(groot,'defaultAxesColorOrder',mycol)        
        axis([0 round(max(max(EigsAll)))+3 0 100]);
        grid on; hold on;
        for nalg = 1:numel(Algorithms)
            plot(EigsAll(:,nalg),Y,'-','Linewidth',2);
        end
        set(gca,'XTick',0:1:round(max(max(EigsAll)))+3,'Fontsize',14,'Fontname','Times New Roman')
        xlabel('Average Eigenvalue','Fontsize',16,'Fontname','Times New Roman');
        ylabel('CDF [%]','Fontsize',16,'Fontname','Times New Roman');     
        legend(strcat(Algorithms,"avr"),'Location','southeast');
        title(strcat(target_snr,',',target_K));
    end
end
% End