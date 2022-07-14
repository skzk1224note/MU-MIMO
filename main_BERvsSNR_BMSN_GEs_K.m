clear;

rng('default');

%% �o�̓t�@�C����
%evfile1 = 'BER_CSV3/BER16x2x8u_BMSN2BF_SNR0-30dB_10000bits_5000itr.csv';

% % evfile1 = 'BER_CSV4/BER16x2x8u_BMSN1GE_SNR0-30dB_10000bits_5000itr.csv';
% % evfile2 = 'BER_CSV4/BER16x2x8u_BMSN2GE_SNR0-30dB_10000bits_5000itr.csv';
% % evfile3 = 'BER_CSV4/BER16x2x8u_BMSN2GEs_SNR0-30dB_10000bits_5000itr.csv';
% % evfile4 = 'BER_CSV4/BER16x2x8u_BMSN2GEls_SNR0-30dB_10000bits_5000itr.csv';
% % evfile5 = 'BER_CSV4/BER16x2x8u_BMSN2GES1_SNR0-30dB_10000bits_5000itr.csv';
% % evfile6 = 'BER_CSV4/BER16x2x8u_BMSN1BF_SNR0-30dB_10000bits_5000itr.csv';
% % evfile7 = 'BER_CSV4/BER16x2x8u_BMSN2BF_SNR0-30dB_10000bits_5000itr.csv';

% evfile8 = 'BER_CSV3/BER16x2x8u_GMMSE2_SNR0-30dB_10000bits_5000itr.csv';
%% �p�����[�^
Nt = 16;          %��n�ǃA���e�i��
Nr = 1;           %���[�U�A���e�i����2(=���M�X�g���[����(not1))
Nu = 8;           %���[�U��
ndata = 10000;    % 1�p�P�b�g������̃V���{�����i�p�P�b�g���j:BER�v�Z�̂��߂�1���s������̃V���{���� 10000
SNR_min = 0;      %�ŏ�SNR[dB]
SNR_max = 30;     %�ő�SNR[dB]
SIMU = 5000;        %���s�񐔁i�ʏ� 5000�j
if Nr == 2
    % pattern = [4,0;2,2;3,1]; %�ϒ��p�^�[��(�擪��[��,0]�������Ă���j
    
    % ���K���ϒ��ł͂Ȃ��C�p�^�[�����Ƃɉ��
    pattern = [4,0];% 2020/4/18 BER�����@��M�A���e�i��1�@4�C0 K=20dB
    % �������܂�
    
    bsr = pattern(1,1);
elseif Nr == 1
    pattern = 4;
    bsr = 4;
end

I = eye(Nt,Nt);

% Rice�t�F�[�W���O�̂��߂̕ϐ���`
K_dB   = 20;         % Rician��K�t�@�N�^
K      = 10^(K_dB/10);
d_t      = 0.5;      % ���M�A���e�i�Ԋu�iin wavelength)
d_r      = 0.5;      % ��M�A���e�i�Ԋu�iin wavelength)
derad = pi/180;      % degree -> rad

An = 2; % �w�����֐��̌W��(cos�p�^�[���̏ꍇ��2)
%%
if K_dB < 10
   target_K=strcat('{\it{K}}=',num2str(K_dB,'%01d'),'dB');
else
   target_K=strcat('{\it{K}}=',num2str(K_dB,'%02d'),'dB');
end

T = zeros(Nr,Nr); % ���]�̃`���l���s��
for inu = 1:Nu
    T(:,:,inu) = eye(Nr,Nr);
end

tic;

snr = SNR_min:SNR_max;  % snr=[SNR_min,....,SNR_max]�̉��x�N�g��

Lsnr = length(snr);     % snr���x�N�g���̑傫��


% �`���`���l���s��̒��ڔg����(LOS �`���l��)
H_los = zeros(Nu*Nr,Nt);
%% simulation
for Directivity_switch = 0:1 % 0,1�ŃI���C�I�t��؂�ւ�
for isnr = 1:Lsnr
    sigma2 = 1/(10^(snr(isnr)/10)); % noise power
    a = sigma2*Nt;
    for isimu = 1:SIMU
        
        % H (�`���`���l���s��:Rician channel)
        % �`���`���l���s��̃}���`�p�X���� (i.i.d. Rayleigh , NLOS �`���l��)
        H_iid = (randn(Nr*Nu,Nt)+1j*randn(Nr*Nu,Nt))/sqrt(2);
        
        % LOS �`���l��      
        if Directivity_switch == 1 % ����M�f�q�̎w�����l���L��
            Theta_t = (rand(1,Nu)-0.5)*180; % ���[�U���̑��M�p �w����:(-90deg - 90deg)
            Theta_r = (rand(1,Nu)-0.5)*180; % ���[�U���̎�M�p �w����:(-90deg - 90deg)
            for n = 1 : Nu
                a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad))*An*cos(Theta_t(1,n)*derad); % ���[�U���̑��M���[�h�x�N�g��
                a_r = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad))*An*cos(Theta_r(1,n)*derad); % ���[�U���̎�M���[�h�x�N�g��
                H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t';                 % ���[�U����LOS�`���l���s��
            end
        elseif Directivity_switch == 0 % ����M�f�q�̎w�����l������
            Theta_t = (rand(1,Nu)-0.5)*360; % ���[�U���̑��M�p ������:(-180deg - 180deg) 
            Theta_r = (rand(1,Nu)-0.5)*360; % ���[�U���̎�M�p ������:(-180deg - 180deg) 
            for n = 1 : Nu
                a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad)); % ���[�U���̑��M���[�h�x�N�g��
                a_r = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad)); % ���[�U���̎�M���[�h�x�N�g��
                H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t';                 % ���[�U����LOS�`���l���s��
            end
        end
        % �`���`���l���s��=[sqrt(K/(K+1))*(LOS �`���l��)]+[sqrt(1/(K+1))*(NLOS �`���l��)]
        H = sqrt(K/(K+1))*H_los + sqrt(1/(K+1))*H_iid;
        %Hu(���[�Uk�̃`���l���s��Hu(:,:,k))
        Hu = peruser(H,Nu);
        HH = H'*H;
        
        %% ZF�E�G�C�g
        %Wzf = H\I;
        %Wzf = H'/(H*H');
        
        %% MMSE�E�G�C�g
        Wmmse = H'/(H*H'+ a*eye(Nr*Nu)); % MMSE
%         
%         for inu = 1:Nu
%             %h = (randn(2,nt) + 1j*randn(2,nt))/sqrt(2);
%             h = H(1+(inu-1)*Nr:2+(inu-1)*Nr,:);
%             S = [svd(h(1,:));svd(h(2,:))];  % nr =2
%             [~,No] = max(S);
%             Hu_as(:,:,inu) = h(No,:);
%         end
        
        %H(�`���`���l���s��)
        %H_as = alluser(Hu_as); %�`���`���l���s��
        He = zeros(Nr*Nu-Nr,Nt,Nu);

        for inu = 1:Nu
            
            %He(H���烆�[�Uk�̃`���l���s����������s��)
            el = 1:Nu;
            el(:,inu) = [];
            He(:,:,inu) = alluser(Hu(:,:,el));
            HeHe = He(:,:,inu)'*He(:,:,inu);
            %He_as(:,:,inu) = alluser(Hu_as(:,:,el));
            
            %% BD
            %Ven(He�̎G��������ԂɑΉ�����ŗL�x�N�g��)
            [~,~,Ve] = svd(He(:,:,inu));
            Ven(:,:,inu) = Ve(:,Nr*(Nu-1)+1:Nt);
            
            %Vts(Hu*Ven�̐M��������ԂɑΉ�����ŗL�x�N�g��)
            %Wtu(���[�Uk�ɑΉ������n�Ǒ��E�G�C�g)
            %St(���ْl���i�[ 2�~nu)
            %Wr(���[�U���E�G�C�g)
            [Ut,Std,Vt] = svd(Hu(:,:,inu)*Ven(:,:,inu));
            St_BD(:,inu) = diag(Std(1:Nr,1:Nr));
            
            Wr_BD(:,:,inu) = Ut';
            Vts_BD(:,:,inu) = Vt(:,1:Nr);
            Wtu_BD(:,:,inu) = Ven(:,:,inu)*Vts_BD(:,:,inu);
            
            %% MMSE
            % inu�̃E�G�C�g�i��x�N�g�����ɐ��K���j
            Wopt(:,:,inu) = Wmmse(:,(inu-1)*Nr+1:(inu-1)*Nr+Nr);
            [Q,~]=qr(Wopt(:,:,inu));
            QMJ(:,:,inu)=Q(:,1:Nr);
            for ij = 1:Nr
                Wopt(:,ij,inu) = Wopt(:,ij,inu)/sqrt(Wopt(:,ij,inu)'*Wopt(:,ij,inu));
            end
            
            % �G��������ԃx�N�g����nuer�̃`���l���s��ɏ�Z
            HTT=Hu(:,:,inu)*Wopt(:,:,inu);      
            % �ϊ��s���SVD
            [Ut,Std,Vt]=svd(HTT);
            St_MMSE(:,inu) = diag(Std(1:Nr,1:Nr));
            Wr_MMSE(:,:,inu) = Ut';
            Vts_MMSE(:,:,inu) = Vt(:,1:Nr);
            Wtu_MMSE(:,:,inu) = Wopt(:,:,inu)*Vts_MMSE(:,:,inu);
            
            %% GMMSE0
            %     % �G��������ԃx�N�g����nuer�̃`���l���s��ɏ�Z                                 
            %     HTT=Hu(:,:,inu)*QMJ(:,:,inu);      
            %     % �ϊ��s���SVD
            %     [Ut,Std,Vt]=svd(HTT);
            %     St_GMMSE0(:,inu) = diag(Std(1:Nr,1:Nr));
            %     Wr_GMMSE0(:,:,inu) = Ut';
            %     Vts_GMMSE0(:,:,inu) = Vt(:,1:Nr);
            %     Wtu_GMMSE0(:,:,inu) = QMJ(:,:,inu)*Vts_GMMSE0(:,:,inu);
                 
            %% GMI1
            TT1(:,:,inu)=(QMJ(:,:,inu)'*HH*QMJ(:,:,inu)+a*eye(Nr))\...
                (QMJ(:,:,inu)'*(Hu(:,:,inu)'*Hu(:,:,inu))*QMJ(:,:,inu));
            bb = 0;
            bb = bb + trace(TT1(:,:,inu)'*TT1(:,:,inu));
            
            Wo=QMJ(:,:,inu)*TT1(:,:,inu);
            HTT=Hu(:,:,inu)*Wo;      
            % �ϊ��s���SVD
            [Ut,Std,Vt]=svd(HTT);
            St_GMI1(:,inu) = diag(Std(1:Nr,1:Nr));
            Wr_GMI1(:,:,inu) = Ut';
            Vts_GMI1(:,:,inu) = Vt(:,1:Nr);
            Wtu_GMI1(:,:,inu) = Wo*Vts_GMI1(:,:,inu);
            
            %% GMI2
            LL=QMJ(:,:,inu)'*HeHe*QMJ(:,:,inu)+a*eye(Nr);
            L=chol((LL+LL')/2);
            TT2=eye(Nr)/L;
            gamma = sqrt(Nr/trace(TT2'*TT2));
            TT2=gamma*TT2;
            
            Wo=QMJ(:,:,inu)*TT2;
            HTT=Hu(:,:,inu)*Wo;      
            % �ϊ��s���SVD
            [Ut,Std,Vt]=svd(HTT);
            St_GMI2(:,inu) = diag(Std(1:Nr,1:Nr));
            Wr_GMI2(:,:,inu) = Ut';
            Vts_GMI2(:,:,inu) = Vt(:,1:Nr);
            Wtu_GMI2(:,:,inu) = Wo*Vts_GMI2(:,:,inu);
            
            
            %% BMSN1-BF
            %MSN�Ɋ�Â����[�U���̃E�G�C�g�̌v�Z
            % MSN�̍œK�E�G�C�gWopt
            Wopt = (HeHe+a*eye(size(He,2)))\Hu(:,:,inu)'*T(:,:,inu);
            %     Wopt(:,:,i) = Wopt(:,:,i)/norm(Wopt(:,:,i),'fro');
            for ij = 1:Nr
                Wopt(:,ij) = Wopt(:,ij)/sqrt(Wopt(:,ij)'*Wopt(:,ij));
            end
            %Vts(Hu*Wopt�̐M��������ԂɑΉ�����ŗL�x�N�g��)
            %Wtu(���[�Uk�ɑΉ������n�Ǒ��E�G�C�g)
            %St(���ْl���i�[ 2�~nu)
            %Wr(���[�U���E�G�C�g)
            [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt);
            St_BMSN_BF(:,inu) = diag(Std(1:Nr,1:Nr));
            
            Wr_BMSN_BF(:,:,inu) = Ut';
            Vts_BMSN_BF(:,:,inu) = Vt(:,1:Nr);
            Wtu_BMSN_BF(:,:,inu) = Wopt*Vts_BMSN_BF(:,:,inu);
            
            %% BMSN-BF2
            %MSN�Ɋ�Â����[�U���̃E�G�C�g�̌v�Z
            % MSN�̍œK�E�G�C�gWopt
            Wo = (HeHe+a*eye(size(He,2)))\Hu(:,:,inu)';
            LL=Wo'*HeHe*Wo+a*(Wo'*Wo);
            L=chol((LL+LL')/2);
            T(:,:,inu)=eye(Nr)/L;
            %T(:,:,nuser)=eye(NR);
            mu=sqrt(Nr/trace(T(:,:,inu)'*T(:,:,inu)));
            T(:,:,inu)=mu*T(:,:,inu);
            gamma = sqrt(Nr/trace(T(:,:,inu)'*(Wo'*Wo)*T(:,:,inu)));
            Wopt = gamma*Wo*T(:,:,inu);

            %     for ij = 1:Nr
            %         Wopt(:,ij,inu) = Wopt(:,ij,inu)/sqrt(Wopt(:,ij,inu)'*Wopt(:,ij,inu));
            %     end
            %Vts(Hu*Wopt�̐M��������ԂɑΉ�����ŗL�x�N�g��)
            %Wtu(���[�Uk�ɑΉ������n�Ǒ��E�G�C�g)
            %St(���ْl���i�[ 2�~nu)
            %Wr(���[�U���E�G�C�g)
            [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt);
            St_BMSN_BF2(:,inu) = diag(Std(1:Nr,1:Nr));
            
            Wr_BMSN_BF2(:,:,inu) = Ut';
            Vts_BMSN_BF2(:,:,inu) = Vt(:,1:Nr);
            Wtu_BMSN_BF2(:,:,inu) = Wopt*Vts_BMSN_BF2(:,:,inu);
            
            %% B-MSN2-M2 ��=sigma2 NT equal to MMSE
            %     Wo =(HeHe+a*eye(size(He,2)))\Hu(:,:,inu)';
            %     
            %     LL=Wo'*HH*Wo+a*(Wo'*Wo)/100;
            %     L=chol((LL+LL')/2);
            %     TT=eye(Nr)/L;
            %     mu=sqrt(Nr/trace(TT'*TT));
            %     TT=mu*TT;
            %     gamma = sqrt(Nr/trace(TT'*(Wo'*Wo)*TT));
            %     Wopt(:,:,inu)=gamma*Wo*TT;
            % 
            %     [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt(:,:,inu));
            %     St_BMSN2M2(:,inu) = diag(Std(1:Nr,1:Nr));
            %     Wr_BMSN2M2(:,:,inu) = Ut';
            %     Vts_BMSN2M2(:,:,inu) = Vt(:,1:Nr);
            %     Wtu_BMSN2M2(:,:,inu) = Wopt(:,:,inu)*Vts_BMSN2M2(:,:,inu);
            
            %% BMSN-GE0 :��ʉ��ŗL�l�����������@(�[���G������j
%             A = HeHe+a*eye(size(He(:,:,inu),2));
%             B = Hu(:,:,inu)'*Hu(:,:,inu);
%             [EW,D] = eig(B,A);
%             [~,IN] = sort(diag(abs(D)).','descend');
%             EW = EW(:,IN);
%             Wopt=EW(:,1:Nr);
%             for ij = 1:Nr
%                 Wopt(:,ij) = Wopt(:,ij)/sqrt(Wopt(:,ij)'*Wopt(:,ij));
%             end
%             %Wopt'*A*Wopt
%             %Vts(HT*W�̐M��������ԂɑΉ�����ŗL�x�N�g��)
%             %Wtu(���[�Uk�ɑΉ������n�Ǒ��E�G�C�g)
%             %St(���ْl���i�[ 2�~nu)
%             %Wr(���[�U���E�G�C�g)
%             [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt);
%             St_BMSN1g(:,inu) = diag(Std(1:Nr,1:Nr));
%             Wr_BMSN1g(:,:,inu) = Ut';
%             Vts_BMSN1g(:,:,inu) = Vt(:,1:Nr);
%             Wtu_BMSN1g(:,:,inu) = Wopt*Vts_BMSN1g(:,:,inu);
            
            %% BMSN-GE1 :��ʉ��ŗL�l�����������@(�[���G������j
            A = HeHe+a*eye(size(He(:,:,inu),2));
            B = Hu(:,:,inu)'*Hu(:,:,inu);
            [EW,D] = eig(B,A);
            [~,IN] = sort(diag(abs(D)).','descend');
            EW = EW(:,IN);
            Wopt=EW(:,1:Nr);
            Wopt = Wopt/norm(Wopt,'fro')*sqrt(Nr);
            %Vts(HT*W�̐M��������ԂɑΉ�����ŗL�x�N�g��)
            %Wtu(���[�Uk�ɑΉ������n�Ǒ��E�G�C�g)
            %St(���ْl���i�[ 2�~nu)
            %Wr(���[�U���E�G�C�g)
            [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt);
            St_BMSN_GE1(:,inu) = diag(Std(1:Nr,1:Nr));
            Wr_BMSN_GE1(:,:,inu) = Ut';
            Vts_BMSN_GE1(:,:,inu) = Vt(:,1:Nr);
            Wtu_BMSN_GE1(:,:,inu) = Wopt*Vts_BMSN_GE1(:,:,inu);
            
            %% BMSN-GE2 ��=sigma2 NT equal to MMSE
            A = HeHe+a*eye(size(He(:,:,inu),2));
            B = Hu(:,:,inu)'*Hu(:,:,inu);
            [EW,D] = eig(B,A);
            [D1,IN] = sort(diag(abs(D)).','descend');
            EW = EW(:,IN);
            D2 = diag(D1(:,1:Nr));
            %Wopt=EW(:,1:Nr)*sqrt(log2(1+D2));
            EWS=EW(:,1:Nr);
            %     for ij = 1:Nr
            %         EWS(:,ij) = EWS(:,ij)/sqrt(EWS(:,ij)'*EWS(:,ij));
            %     end
            Wopt=EWS*sqrt(D2);
            Wopt = Wopt/norm(Wopt,'fro')*sqrt(Nr);
            %     Wopt_k2(:,:,inu) = Wopt_k2(:,:,inu)/norm(Wopt_k2(:,:,inu),'fro')*sqrt(Nr);
            %Vts(HT*W�̐M��������ԂɑΉ�����ŗL�x�N�g��)
            %Wtu(���[�Uk�ɑΉ������n�Ǒ��E�G�C�g)
            %St(���ْl���i�[ 2�~nu)
            %Wr(���[�U���E�G�C�g)  
            [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt);
            St_BMSN_GE2(:,inu) = diag(Std(1:Nr,1:Nr));
            Wr_BMSN_GE2(:,:,inu) = Ut';
            Vts_BMSN_GE2(:,:,inu) = Vt(:,1:Nr);
            Wtu_BMSN_GE2(:,:,inu) = Wopt*Vts_BMSN_GE2(:,:,inu);
            
            %% BMSN-GE3
            A = HeHe+a*eye(size(He(:,:,inu),2));
            B = Hu(:,:,inu)'*Hu(:,:,inu);
            [EW,D] = eig(B,A);
            [D_GE3,IN] = sort(diag(abs(D)).','descend');
            EW = EW(:,IN);
            if D_GE3(:,Nr) < 0.1 && Nr > 1% lambda2��������������lambda1�ɂ���
                EW(:,Nr) = EW(:,Nr-1); % D1(:,Nr) = D1(:,Nr-1)�Ƃ��邽��,�ŗL�x�N�g�����Ή�������
            else
                if Nr ==1
                    EW(:,Nr) = EW(:,Nr);% karini
                %                 fprintf('lambda2 = %.1f \n', D_GE3(:,Nr));
                else
                    EW(:,Nr) = EW(:,Nr-1);% karini
                end
            end
            Wopt = EW(:,1:Nr);
            Wopt = Wopt/norm(Wopt,'fro')*sqrt(Nr);
            
            [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt);            
            St_BMSN_GE3(:,inu) = diag(Std(1,1)); % St(���ْl�͂ЂƂ�)
            Wr_BMSN_GE3(:,:,inu) = Ut(:,1)'; % Wr(���[�U���E�G�C�g)
            Vts_BMSN_GE3(:,:,inu) = Vt(:,1); % Vts(HT*W�̐M��������ԂɑΉ�����ŗL�x�N�g��)
            Wtu_BMSN_GE3(:,:,inu) = Wopt*Vts_BMSN_GE3(:,:,inu); % Wtu(���[�Uk�ɑΉ������n�Ǒ��E�G�C�g)
       
            %% BMSN-GE2-logSINR ��=sigma2 NT equal to MMSE
%             A = HeHe+a*eye(size(He(:,:,inu),2));
%             B = Hu(:,:,inu)'*Hu(:,:,inu);
%             [EW,D] = eig(B,A);
%             [D1,IN] = sort(diag(abs(D)).','descend');
%             EW = EW(:,IN);
%             D2 = diag(D1(:,1:Nr));
%             EWS=EW(:,1:Nr);
%             %     for ij = 1:Nr
%             %         EWS(:,ij) = EWS(:,ij)/sqrt(EWS(:,ij)'*EWS(:,ij));
%             %     end          
%             Wopt=EWS*sqrt(log2(1+D2));
%             %Wopt=EWS(:,1:Nr)*sqrt(D2);
%             Wopt = Wopt/norm(Wopt,'fro')*sqrt(Nr);
%             %     Wopt_k2(:,:,inu) = Wopt_k2(:,:,inu)/norm(Wopt_k2(:,:,inu),'fro')*sqrt(Nr);
%             %Vts(HT*W�̐M��������ԂɑΉ�����ŗL�x�N�g��)
%             %Wtu(���[�Uk�ɑΉ������n�Ǒ��E�G�C�g)
%             %St(���ْl���i�[ 2�~nu)
%             %Wr(���[�U���E�G�C�g)
%             [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt);
%             St_BMSN2gLS(:,inu) = diag(Std(1:Nr,1:Nr));
%             Wr_BMSN2gLS(:,:,inu) = Ut';
%             Vts_BMSN2gLS(:,:,inu) = Vt(:,1:Nr);
%             Wtu_BMSN2gLS(:,:,inu) = Wopt*Vts_BMSN2gLS(:,:,inu);
            
            %% BMSN-GE2-SINR1 ��=sigma2 NT equal to MMSE
            % �����ŗL�x�N�g����p���Ă���_��BMSN-GE3�Ǝ��Ă���i�Ō�̃E�G�C�g�K�i�����Ⴄ�j�����߂Ɏ������ԈႦ����@�Ɠ����H
%             B = HeHe+a*eye(size(He(:,:,inu),2));
%             A = Hu(:,:,inu)'*Hu(:,:,inu);
%             [EW,D] = eig(B,A);
%             [~,IN] = sort(diag(D));
%             EW = EW(:,IN);
%             W_k = EW(:,1)/sqrt(EW(:,1)'*EW(:,1));
%             for ij = 1:Nr
%                 Wopt(:,ij) = W_k;
%             end
%             %Wopt = Wopt/norm(Wopt,'fro')*sqrt(Nr);
%             %Vts(HT*W�̐M��������ԂɑΉ�����ŗL�x�N�g��)
%             %Wtu(���[�Uk�ɑΉ������n�Ǒ��E�G�C�g)
%             %St(���ْl���i�[ 2�~nu)
%             %Wr(���[�U���E�G�C�g)
%             [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt);
%             St_BMSN2gS1(:,inu) = diag(Std(1:Nr,1:Nr));
%             Wr_BMSN2gS1(:,:,inu) = Ut';
%             Vts_BMSN2gS1(:,:,inu) = Vt(:,1:Nr);
%             Wtu_BMSN2gS1(:,:,inu) = Wopt*Vts_BMSN2gS1(:,:,inu);
            
            
            %% ZF(zero-Forcing)
            %     % inu�̃E�G�C�g�i��x�N�g�����ɐ��K���j
            %     Wopt(:,:,inu) = Wzf(:,(inu-1)*Nr+1:(inu-1)*Nr+Nr);
            %     for ij = 1:Nr
            %         Wopt(:,ij,inu) = Wopt(:,ij,inu)/sqrt(Wopt(:,ij,inu)'*Wopt(:,ij,inu));
            %     end
            % 
            %     % �G��������ԃx�N�g����nuer�̃`���l���s��ɏ�Z                                 
            %     HTT=Hu(:,:,inu)*Wopt(:,:,inu);      
            %     % �ϊ��s���SVD
            %     [Ut,Std,Vt]=svd(HTT);
            %     St_ZF(:,inu) = diag(Std(1:Nr,1:Nr));
            %     Wr_ZF(:,:,inu) = Ut';
            %     Vts_ZF(:,:,inu) = Vt(:,1:Nr);
            %     Wtu_ZF(:,:,inu) = Wopt(:,:,inu)*Vts_ZF(:,:,inu);
            %     
            %% GZF(zero-Forcing)
            %     % inu�̃E�G�C�g�i��x�N�g�����ɐ��K���j
            %     Wopt(:,:,inu) = Wzf(:,(inu-1)*Nr+1:(inu-1)*Nr+Nr);
            % %     for ij = 1:Nr
            % %         Wopt(:,ij,inu) = Wopt(:,ij,inu)/sqrt(Wopt(:,ij,inu)'*Wopt(:,ij,inu));
            % %     end
            % 
            %     [Q,~]=qr(Wopt(:,:,inu));
            %     QZJ(:,:,inu)=Q(:,1:Nr);
            % %     for ij = 1:Nr
            % %         Q(:,ij) = Q(:,ij)/sqrt(Q(:,ij)'*Q(:,ij));
            % %     end
            %     
            %     % �G��������ԃx�N�g����nuer�̃`���l���s��ɏ�Z                                 
            %     HTT=Hu(:,:,inu)*QZJ(:,:,inu);      
            %     % �ϊ��s���SVD
            %     [Ut,Std,Vt]=svd(HTT);
            %     St_GZF(:,inu) = diag(Std(1:Nr,1:Nr));
            %     Wr_GZF(:,:,inu) = Ut';
            %     Vts_GZF(:,:,inu) = Vt(:,1:Nr);
            %     Wtu_GZF(:,:,inu) = QZJ(:,:,inu)*Vts_GZF(:,:,inu);
            
        end
        
        %Wt(��n�Ǒ��E�G�C�g)
        %Wt = reshape(Wtu,[nt,nt]);
        
        amp2 = 10^(snr(isnr)/10)/Nt; % amp(���M�d��)
        
        % �ȍ~ nr=2 �̏ꍇ�̂�
        if Nr == 2
            %data(���[�Uk�ւ̑��M�f�[�^)
            for inu = 1:Nu
                for ii = 1:size(pattern,1)
                    data1 = randi([0 2^pattern(ii,1)-1],1,ndata);
                    data2 = randi([0 2^pattern(ii,2)-1],1,ndata);
                    data(:,:,ii,inu) = [data1;data2];
                end
            end
            %s(���[�Uk�ւ̑��M�X�g���[��)
            for inu = 1:Nu
                for ii = 1:size(pattern,1)
                    s1 = Mapping(data(1,:,ii,inu),pattern(ii,1));
                    s2 = Mapping(data(2,:,ii,inu),pattern(ii,2));
                    s(:,:,ii,inu) = [s1;s2];                    
                    s_GE3(:,:,ii,inu) = s1; % BMSN-GE3��p�̑��M�X�g���[��
                end
            end
        elseif Nr==1
            %data(���[�Uk�ւ̑��M�f�[�^)
            for inu = 1:Nu
                for ii = 1:size(pattern,1)
                    data1 = randi([0 2^pattern(ii,1)-1],1,ndata);
                    data(:,:,ii,inu) = data1;
                end
            end
            %s(���[�Uk�ւ̑��M�X�g���[��)
            for inu = 1:Nu
                for ii = 1:size(pattern,1)
                    s1 = Mapping(data(1,:,ii,inu),pattern(ii,1));
                    s(:,:,ii,inu) = s1;
                    s_GE3(:,:,ii,inu) = s1; % BMSN-GE3��p�̑��M�X�g���[��
                end
            end
        end
        % %data(���M�f�[�^) for antenna selection (AS)
        % data_as = randi([0 2^bsr-1],Nu,ndata);
        % 
        % %s(���M�X�g���[��) for antenna selection (AS)
        % s_as = Mapping(data_as,bsr);
        
        
        %y(���[�Uk�̓��͐M��)
        %Y(���[�Uk�̕����M��)
        for inu = 1:Nu
            n = (randn(Nr,ndata) + 1j*randn(Nr,ndata))/sqrt(2); % n(�M�G��)
            
            for ii = 1:size(pattern,1)
                Hs_BD = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
                % Hs_ZF = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
                % Hs_GZF = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
                Hs_GMI1 = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
                Hs_GMI2 = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
                Hs_BMSN_BF = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
                Hs_BMSN_BF2 = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
                %Hs_BMSN1g = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
                Hs_BMSN_GE1 = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
                Hs_BMSN_GE2 = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
                Hs_BMSN_GE3 = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
%                 Hs_BMSN2gLS = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
%                 Hs_BMSN2gS1 = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
                for ij = 1:Nu
                    Hs_BD = Hs_BD + Hu(:,:,inu)*Wtu_BD(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2); % amp2 = 10^(snr(isnr)/10)/Nt (���M�d��)
                    %             Hs_ZF = Hs_ZF + Hu(:,:,inu)*Wtu_ZF(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
                    %             Hs_GZF = Hs_GZF + Hu(:,:,inu)*Wtu_GZF(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
                    Hs_GMI1 = Hs_GMI1 + Hu(:,:,inu)*Wtu_GMI1(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
                    Hs_GMI2 = Hs_GMI2 + Hu(:,:,inu)*Wtu_GMI2(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
                    Hs_BMSN_BF = Hs_BMSN_BF + Hu(:,:,inu)*Wtu_BMSN_BF(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
                    Hs_BMSN_BF2 = Hs_BMSN_BF2 + Hu(:,:,inu)*Wtu_BMSN_BF2(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
                    % Hs_BMSN1g = Hs_BMSN1g + Hu(:,:,inu)*Wtu_BMSN1g(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
                    Hs_BMSN_GE1 = Hs_BMSN_GE1 + Hu(:,:,inu)*Wtu_BMSN_GE1(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
                    Hs_BMSN_GE2 = Hs_BMSN_GE2 + Hu(:,:,inu)*Wtu_BMSN_GE2(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
                    Hs_BMSN_GE3 = Hs_BMSN_GE3 + Hu(:,:,inu)*Wtu_BMSN_GE3(:,:,ij)*s_GE3(:,:,ii,ij)*sqrt(amp2);
%                     Hs_BMSN2gLS = Hs_BMSN2gLS + Hu(:,:,inu)*Wtu_BMSN2gLS(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
%                     Hs_BMSN2gS1 = Hs_BMSN2gS1 + Hu(:,:,inu)*Wtu_BMSN2gS1(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
                end
                yu_BD(:,:,ii,inu) = Wr_BD(:,:,inu)*(Hs_BD+n);
                %         yu_ZF(:,:,ii,inu) = Wr_ZF(:,:,inu)*(Hs_ZF+n);
                %         yu_GZF(:,:,ii,inu) = Wr_GZF(:,:,inu)*(Hs_GZF+n);
                yu_GMI1(:,:,ii,inu) = Wr_GMI1(:,:,inu)*(Hs_GMI1+n);
                yu_GMI2(:,:,ii,inu) = Wr_GMI2(:,:,inu)*(Hs_GMI2+n);
                yu_BMSN_BF(:,:,ii,inu) = Wr_BMSN_BF(:,:,inu)*(Hs_BMSN_BF+n);
                yu_BMSN_BF2(:,:,ii,inu) = Wr_BMSN_BF2(:,:,inu)*(Hs_BMSN_BF2+n);
                % yu_BMSN1g(:,:,ii,inu) = Wr_BMSN1g(:,:,inu)*(Hs_BMSN1g+n);
                yu_BMSN_GE1(:,:,ii,inu) = Wr_BMSN_GE1(:,:,inu)*(Hs_BMSN_GE1+n);
                yu_BMSN_GE2(:,:,ii,inu) = Wr_BMSN_GE2(:,:,inu)*(Hs_BMSN_GE2+n);
                yu_BMSN_GE3(:,:,ii,inu) = Wr_BMSN_GE3(:,:,inu)*(Hs_BMSN_GE3+n);
                % yu_BMSN2gLS(:,:,ii,inu) = Wr_BMSN2gLS(:,:,inu)*(Hs_BMSN2gLS+n);
                % yu_BMSN2gS1(:,:,ii,inu) = Wr_BMSN2gS1(:,:,inu)*(Hs_BMSN2gS1+n);
                
                
                
                %     
                %         yu_ZF(1,:,ii,inu) = yu_ZF(1,:,ii,inu)/St_ZF(1,inu)/sqrt(amp2);
                %         yu_ZF(2,:,ii,inu) = yu_ZF(2,:,ii,inu)/St_ZF(2,inu)/sqrt(amp2);
                %         Yu_ZF(1,:,ii,inu) = Decode(yu_ZF(1,:,ii,inu),pattern(ii,1));
                %         Yu_ZF(2,:,ii,inu) = Decode(yu_ZF(2,:,ii,inu),pattern(ii,2));
                % 
                %         yu_GZF(1,:,ii,inu) = yu_GZF(1,:,ii,inu)/St_GZF(1,inu)/sqrt(amp2);
                %         yu_GZF(2,:,ii,inu) = yu_GZF(2,:,ii,inu)/St_GZF(2,inu)/sqrt(amp2);
                %         Yu_GZF(1,:,ii,inu) = Decode(yu_GZF(1,:,ii,inu),pattern(ii,1));
                %         Yu_GZF(2,:,ii,inu) = Decode(yu_GZF(2,:,ii,inu),pattern(ii,2));
                
                yu_BD(1,:,ii,inu) = yu_BD(1,:,ii,inu)/St_BD(1,inu)/sqrt(amp2);
                Yu_BD(1,:,ii,inu) = Decode(yu_BD(1,:,ii,inu),pattern(ii,1));
                yu_GMI1(1,:,ii,inu) = yu_GMI1(1,:,ii,inu)/St_GMI1(1,inu)/sqrt(amp2);
                Yu_GMI1(1,:,ii,inu) = Decode(yu_GMI1(1,:,ii,inu),pattern(ii,1));
                yu_GMI2(1,:,ii,inu) = yu_GMI2(1,:,ii,inu)/St_GMI2(1,inu)/sqrt(amp2);
                Yu_GMI2(1,:,ii,inu) = Decode(yu_GMI2(1,:,ii,inu),pattern(ii,1));
                yu_BMSN_BF(1,:,ii,inu) = yu_BMSN_BF(1,:,ii,inu)/St_BMSN_BF(1,inu)/sqrt(amp2);
                Yu_BMSN_BF(1,:,ii,inu) = Decode(yu_BMSN_BF(1,:,ii,inu),pattern(ii,1));
                yu_BMSN_BF2(1,:,ii,inu) = yu_BMSN_BF2(1,:,ii,inu)/St_BMSN_BF2(1,inu)/sqrt(amp2);
                Yu_BMSN_BF2(1,:,ii,inu) = Decode(yu_BMSN_BF2(1,:,ii,inu),pattern(ii,1));
%                 yu_BMSN1g(1,:,ii,inu) = yu_BMSN1g(1,:,ii,inu)/St_BMSN1g(1,inu)/sqrt(amp2);
%                 Yu_BMSN1g(1,:,ii,inu) = Decode(yu_BMSN1g(1,:,ii,inu),pattern(ii,1));
                yu_BMSN_GE1(1,:,ii,inu) = yu_BMSN_GE1(1,:,ii,inu)/St_BMSN_GE1(1,inu)/sqrt(amp2);
                Yu_BMSN_GE1(1,:,ii,inu) = Decode(yu_BMSN_GE1(1,:,ii,inu),pattern(ii,1));
                yu_BMSN_GE2(1,:,ii,inu) = yu_BMSN_GE2(1,:,ii,inu)/St_BMSN_GE2(1,inu)/sqrt(amp2);
                Yu_BMSN_GE2(1,:,ii,inu) = Decode(yu_BMSN_GE2(1,:,ii,inu),pattern(ii,1));
                yu_BMSN_GE3(1,:,ii,inu) = yu_BMSN_GE3(1,:,ii,inu)/St_BMSN_GE3(1,inu)/sqrt(amp2);
                Yu_BMSN_GE3(1,:,ii,inu) = Decode(yu_BMSN_GE3(1,:,ii,inu),pattern(ii,1));                
%                 yu_BMSN2gLS(1,:,ii,inu) = yu_BMSN2gLS(1,:,ii,inu)/St_BMSN2gLS(1,inu)/sqrt(amp2);
%                 Yu_BMSN2gLS(1,:,ii,inu) = Decode(yu_BMSN2gLS(1,:,ii,inu),pattern(ii,1));
%                 yu_BMSN2gS1(1,:,ii,inu) = yu_BMSN2gS1(1,:,ii,inu)/St_BMSN2gS1(1,inu)/sqrt(amp2);
%                 Yu_BMSN2gS1(1,:,ii,inu) = Decode(yu_BMSN2gS1(1,:,ii,inu),pattern(ii,1));
                
                if Nr == 2
                    yu_BD(2,:,ii,inu) = yu_BD(2,:,ii,inu)/St_BD(2,inu)/sqrt(amp2);
                    Yu_BD(2,:,ii,inu) = Decode(yu_BD(2,:,ii,inu),pattern(ii,2));
                    yu_GMI1(2,:,ii,inu) = yu_GMI1(2,:,ii,inu)/St_GMI1(2,inu)/sqrt(amp2);
                    Yu_GMI1(2,:,ii,inu) = Decode(yu_GMI1(2,:,ii,inu),pattern(ii,2));
                    yu_GMI2(2,:,ii,inu) = yu_GMI2(2,:,ii,inu)/St_GMI2(2,inu)/sqrt(amp2);
                    Yu_GMI2(2,:,ii,inu) = Decode(yu_GMI2(2,:,ii,inu),pattern(ii,2));
                    yu_BMSN_BF(2,:,ii,inu) = yu_BMSN_BF(2,:,ii,inu)/St_BMSN_BF(2,inu)/sqrt(amp2);%                
                    Yu_BMSN_BF(2,:,ii,inu) = Decode(yu_BMSN_BF(2,:,ii,inu),pattern(ii,2));%                
                    yu_BMSN_BF2(2,:,ii,inu) = yu_BMSN_BF2(2,:,ii,inu)/St_BMSN_BF2(2,inu)/sqrt(amp2);%                
                    Yu_BMSN_BF2(2,:,ii,inu) = Decode(yu_BMSN_BF2(2,:,ii,inu),pattern(ii,2));%                                
%                     yu_BMSN1g(2,:,ii,inu) = yu_BMSN1g(2,:,ii,inu)/St_BMSN1g(2,inu)/sqrt(amp2);%                
%                     Yu_BMSN1g(2,:,ii,inu) = Decode(yu_BMSN1g(2,:,ii,inu),pattern(ii,2));%                               
                    yu_BMSN_GE1(2,:,ii,inu) = yu_BMSN_GE1(2,:,ii,inu)/St_BMSN_GE1(2,inu)/sqrt(amp2);%                
                    Yu_BMSN_GE1(2,:,ii,inu) = Decode(yu_BMSN_GE1(2,:,ii,inu),pattern(ii,2));%                                
                    yu_BMSN_GE2(2,:,ii,inu) = yu_BMSN_GE2(2,:,ii,inu)/St_BMSN_GE2(2,inu)/sqrt(amp2);%            
                    Yu_BMSN_GE2(2,:,ii,inu) = Decode(yu_BMSN_GE2(2,:,ii,inu),pattern(ii,2));%
                    
%                     yu_BMSN_GE3(2,:,ii,inu) = yu_BMSN_GE3(2,:,ii,inu)/St_BMSN_GE3(2,inu)/sqrt(amp2); % BMSN-GE3�͑��M�X�g���[���ЂƂ�
%                     Yu_BMSN_GE3(2,:,ii,inu) = Decode(yu_BMSN_GE3(2,:,ii,inu),pattern(ii,2));         % BMSN-GE3�͑��M�X�g���[���ЂƂ�
                    
%                     yu_BMSN2gLS(2,:,ii,inu) = yu_BMSN2gLS(2,:,ii,inu)/St_BMSN2gLS(2,inu)/sqrt(amp2);%                
%                     Yu_BMSN2gLS(2,:,ii,inu) = Decode(yu_BMSN2gLS(2,:,ii,inu),pattern(ii,2));%
%                     yu_BMSN2gS1(2,:,ii,inu) = yu_BMSN2gS1(2,:,ii,inu)/St_BMSN2gS1(2,inu)/sqrt(amp2);%                
%                     Yu_BMSN2gS1(2,:,ii,inu) = Decode(yu_BMSN2gS1(2,:,ii,inu),pattern(ii,2));%
                end
                
            end
        end
        
        %ber(�ϒ��p�^�[������BER)  for BD and B-MSN
        for inu = 1:Nu
            for ii = 1:size(pattern,1)
                if ii == 1
                    ber_BD(ii,inu) = BER2(BER1(Yu_BD(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
                    ber_GMI1(ii,inu) = BER2(BER1(Yu_GMI1(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
                    ber_GMI2(ii,inu) = BER2(BER1(Yu_GMI2(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
                    %     ber_ZF(ii,inu) = BER2(BER1(Yu_ZF(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
                    %     ber_GZF(ii,inu) = BER2(BER1(Yu_GZF(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
                    ber_BMSN_BF(ii,inu) = BER2(BER1(Yu_BMSN_BF(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
                    ber_BMSN_BF2(ii,inu) = BER2(BER1(Yu_BMSN_BF2(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
                    % ber_BMSN1g(ii,inu) = BER2(BER1(Yu_BMSN1g(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
                    ber_BMSN_GE1(ii,inu) = BER2(BER1(Yu_BMSN_GE1(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
                    ber_BMSN_GE2(ii,inu) = BER2(BER1(Yu_BMSN_GE2(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
                    
                    if size(Wtu_BMSN_GE3,2) == 1
                        ber_BMSN_GE3(ii,inu) = BER2(BER1(Yu_BMSN_GE3(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
                    end                     
                    
%                     ber_BMSN2gLS(ii,inu) = BER2(BER1(Yu_BMSN2gLS(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
%                     ber_BMSN2gS1(ii,inu) = BER2(BER1(Yu_BMSN2gS1(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
                else
                    ber_BD(ii,inu) = BER2(BER1(Yu_BD(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
                        + BER2(BER1(Yu_BD(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
                    ber_GMI1(ii,inu) = BER2(BER1(Yu_GMI1(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
                        + BER2(BER1(Yu_GMI1(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
                    ber_GMI2(ii,inu) = BER2(BER1(Yu_GMI2(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
                        + BER2(BER1(Yu_GMI2(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
                    %     ber_ZF(ii,inu) = BER2(BER1(Yu_ZF(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
                    %         + BER2(BER1(Yu_ZF(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
                    %     ber_GZF(ii,inu) = BER2(BER1(Yu_GZF(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
                    %         + BER2(BER1(Yu_GZF(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
                    ber_BMSN_BF(ii,inu) = BER2(BER1(Yu_BMSN_BF(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
                        + BER2(BER1(Yu_BMSN_BF(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
                    ber_BMSN_BF2(ii,inu) = BER2(BER1(Yu_BMSN_BF2(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
                        + BER2(BER1(Yu_BMSN_BF2(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
%                     ber_BMSN1g(ii,inu) = BER2(BER1(Yu_BMSN1g(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
%                         + BER2(BER1(Yu_BMSN1g(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
                    ber_BMSN_GE1(ii,inu) = BER2(BER1(Yu_BMSN_GE1(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
                        + BER2(BER1(Yu_BMSN_GE1(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
                    ber_BMSN_GE2(ii,inu) = BER2(BER1(Yu_BMSN_GE2(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
                        + BER2(BER1(Yu_BMSN_GE2(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
                    
                    ber_BMSN_GE3(ii,inu) = BER2(BER1(Yu_BMSN_GE3(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1));
                    
%                     ber_BMSN2gLS(ii,inu) = BER2(BER1(Yu_BMSN2gLS(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
%                         + BER2(BER1(Yu_BMSN2gLS(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
%                     ber_BMSN2gS1(ii,inu) = BER2(BER1(Yu_BMSN2gS1(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
%                         + BER2(BER1(Yu_BMSN2gS1(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
                end
            end
        end
        
        %ber_min(���ŏ��̕ϒ��p�^�[����BER) for BD and B-MSN
        ber_min_BD(isimu,isnr) = mean(min(ber_BD));
        ber_min_GMI1(isimu,isnr) = mean(min(ber_GMI1));
        ber_min_GMI2(isimu,isnr) = mean(min(ber_GMI2));
        % ber_min_ZF(isimu,isnr) = mean(min(ber_ZF));
        % ber_min_GZF(isimu,isnr) = mean(min(ber_GZF));
        ber_min_BMSN_BF(isimu,isnr) = mean(min(ber_BMSN_BF));
        ber_min_BMSN_BF2(isimu,isnr) = mean(min(ber_BMSN_BF2));
        % ber_min_BMSN1g(isimu,isnr) = mean(min(ber_BMSN1g));
        ber_min_BMSN_GE1(isimu,isnr) = mean(min(ber_BMSN_GE1));
        ber_min_BMSN_GE2(isimu,isnr) = mean(min(ber_BMSN_GE2));        
        ber_min_BMSN_GE3(isimu,isnr) = mean(min(ber_BMSN_GE3));        
%         ber_min_BMSN2gLS(isimu,isnr) = mean(min(ber_BMSN2gLS));
%         ber_min_BMSN2gS1(isimu,isnr) = mean(min(ber_BMSN2gS1));
        %fprintf('Iteration = %d / %d\n',isimu, SIMU);
        
    end % isimu
    
    fprintf('SNR = %d dB\n',snr(isnr));
end % isnr

%% �f�[�^
%BER(���s�񐔕���)
if Directivity_switch == 0
    BER_BD = mean(ber_min_BD);
    BER_GMI1 = mean(ber_min_GMI1);
    BER_GMI2 = mean(ber_min_GMI2);
    % BER_ZF = mean(ber_min_ZF);
    % BER_GZF = mean(ber_min_GZF);
    BER_BMSN_BF = mean(ber_min_BMSN_BF);
    BER_BMSN_BF2 = mean(ber_min_BMSN_BF2);
    % BER_BMSN1g = mean(ber_min_BMSN1g);
    BER_BMSN_GE1 = mean(ber_min_BMSN_GE1);
    BER_BMSN_GE2 = mean(ber_min_BMSN_GE2);
    BER_BMSN_GE3 = mean(ber_min_BMSN_GE3);
%     BER_BMSN2gLS = mean(ber_min_BMSN2gLS);
%     BER_BMSN2gS1 = mean(ber_min_BMSN2gS1);
elseif Directivity_switch == 1
    BER_BDcos = mean(ber_min_BD);
    BER_GMI1cos = mean(ber_min_GMI1);
    BER_GMI2cos = mean(ber_min_GMI2);
    % BER_ZFcos = mean(ber_min_ZF);
    % BER_GZFcos = mean(ber_min_GZF);
    BER_BMSN_BFcos = mean(ber_min_BMSN_BF);
    BER_BMSN_BF2cos = mean(ber_min_BMSN_BF2);
    % BER_BMSN1gcos = mean(ber_min_BMSN1g);
    BER_BMSN_GE1cos = mean(ber_min_BMSN_GE1);
    BER_BMSN_GE2cos = mean(ber_min_BMSN_GE2);
    BER_BMSN_GE3cos = mean(ber_min_BMSN_GE3);
%     BER_BMSN2gLScos = mean(ber_min_BMSN2gLS);
%     BER_BMSN2gS1cos = mean(ber_min_BMSN2gS1);
end

end

% % csvwrite(evfile1,[snr;BER_BD]);
% % csvwrite(evfile2,[snr;BER_ZF]);
% % csvwrite(evfile3,[snr;BER_GZF]);
% % csvwrite(evfile6,[snr;BER_BMSN2g2]);
% 
% csvwrite(evfile1,[snr;BER_BMSN1g]);
% csvwrite(evfile2,[snr;BER_BMSN2g]);
% csvwrite(evfile3,[snr;BER_BMSN2gS]);
% csvwrite(evfile4,[snr;BER_BMSN2gLS]);
% csvwrite(evfile5,[snr;BER_BMSN2gS1]);
% csvwrite(evfile6,[snr;BER_BMSN1]);
% csvwrite(evfile7,[snr;BER_BMSN2]);


toc;
%%
%�`��
figure;
semilogy(snr,BER_BMSN_GE1,'--b',snr,BER_BMSN_GE2,'--k',snr,BER_BMSN_GE3,'--r',...
    snr,BER_BMSN_GE1cos,'-b',snr,BER_BMSN_GE2cos,'-k',snr,BER_BMSN_GE3cos,'-ro','Linewidth',2);
axis([SNR_min SNR_max 1e-4 1e0]);
set(gca,'XTick',SNR_min:5:SNR_max,'YTick',[1e-4, 1e-3, 1e-2, 1e-1 1e0],'Fontsize',14,'Fontname','Times New Roman')
lgd = legend;
lgd.NumColumns = 2; % �}��̗񐔂��w��
legend('BMSN-GE1','BMSN-GE2','BMSN-GE3',...
    'BMSN-GE1cos','BMSN-GE2cos','BMSN-GE3cos','Location','southwest')
xlabel('SNR [dB]','Fontsize',16,'Fontname','Times New Roman');
ylabel('Average BER','Fontsize',16,'Fontname','Times New Roman');
title(target_K);
grid on;
hold on;
%% ���ׂĕ`��
figure;
semilogy(snr,BER_BD,'--b',snr,BER_GMI1,'--g',snr,BER_GMI2,'--k',snr,BER_BMSN_BF,'--c',snr,BER_BMSN_BF2,'--m',snr,BER_BMSN_GE3,'--r',...
    snr,BER_BDcos,'-b',snr,BER_GMI1cos,'-g',snr,BER_GMI2cos,'-k',snr,BER_BMSN_BFcos,'-c',snr,BER_BMSN_BF2cos,'-m',snr,BER_BMSN_GE3cos,'-ro','Linewidth',2);
axis([SNR_min SNR_max 1e-4 1e0]);
set(gca,'XTick',SNR_min:5:SNR_max,'YTick',[1e-4, 1e-3, 1e-2, 1e-1 1e0],'Fontsize',14,'Fontname','Times New Roman')
lgd = legend;
lgd.NumColumns = 2; % �}��̗񐔂��w��
legend('BD','GMI1','GMI2','BMSN-BF','BMSN-BF2','BMSN-GE3','BDcos','GMI1cos','GMI2cos','BMSN-BFcos','BMSN-BF2cos','BMSN-GE3cos','Location','southwest')
xlabel('SNR [dB]','Fontsize',16,'Fontname','Times New Roman');
ylabel('Average BER','Fontsize',16,'Fontname','Times New Roman');
title(target_K);
grid on;
hold on;

