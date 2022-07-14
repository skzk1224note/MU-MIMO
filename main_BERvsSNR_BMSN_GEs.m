clear;

rng('default');

% �o�̓t�@�C����
%evfile1 = 'BER_CSV3/BER16x2x8u_BMSN2BF_SNR0-30dB_10000bits_5000itr.csv';
evfile1 = 'BER_CSV4/BER16x2x8u_BMSN1GE_SNR0-30dB_10000bits_5000itr.csv';
evfile2 = 'BER_CSV4/BER16x2x8u_BMSN2GE_SNR0-30dB_10000bits_5000itr.csv';
evfile3 = 'BER_CSV4/BER16x2x8u_BMSN2GEs_SNR0-30dB_10000bits_5000itr.csv';
evfile4 = 'BER_CSV4/BER16x2x8u_BMSN2GEls_SNR0-30dB_10000bits_5000itr.csv';
evfile5 = 'BER_CSV4/BER16x2x8u_BMSN2GES1_SNR0-30dB_10000bits_5000itr.csv';
evfile6 = 'BER_CSV4/BER16x2x8u_BMSN1BF_SNR0-30dB_10000bits_5000itr.csv';
evfile7 = 'BER_CSV4/BER16x2x8u_BMSN2BF_SNR0-30dB_10000bits_5000itr.csv';
% evfile8 = 'BER_CSV3/BER16x2x8u_GMMSE2_SNR0-30dB_10000bits_5000itr.csv';

Nt = 16; %��n�ǃA���e�i��
Nr = 2; %���[�U�A���e�i����2(=���M�X�g���[����(not1))
Nu = 8; %���[�U��
ndata = 100; % 1�p�P�b�g������̃V���{�����i�p�P�b�g���j:BER�v�Z�̂��߂�1���s������̃V���{���� 10000
SNR_min = 0;  %�ŏ�SNR[dB]
SNR_max = 30; %�ő�SNR[dB]
SIMU = 50; %���s�񐔁i�ʏ� 1000,5000�j
pattern = [4,0;2,2;3,1]; %�ϒ��p�^�[��(�擪��[��,0]�������Ă���j
bsr = pattern(1,1);
% a = 1e-2; % BMSN�̋[���G��
I = eye(Nt,Nt);

%T ���]�̃`���l���s��
T = zeros(Nr,Nr,Nu);

for inu = 1:Nu
%     T(:,:,inu) = zeros(Nr,Nr);
    T(1,1,inu) = 1;
    T(1,2,inu) = 0;
    T(2,1,inu) = 0;
    T(2,2,inu) = 1;
%    T(:,:,inu) = eye(Nr,Nr);
end

tic;

snr = SNR_min:SNR_max;  % snr=[SNR_min,....,SNR_max]�̉��x�N�g��

Lsnr = length(snr);     % snr���x�N�g���̑傫��
    
for isnr = 1:Lsnr
    sigma2 = 1/(10^(snr(isnr)/10)); % noise power
    a = sigma2*Nt;
    %a2 = sigma2*Nt;
    
for isimu = 1:SIMU

%H(�`���`���l���s��: iid Rayleigh channel)

%Hu(���[�Uk�̃`���l���s��Hu(:,:,k))

H = (randn(Nr*Nu,Nt) + 1j*randn(Nr*Nu,Nt))/sqrt(2);
Hu = peruser(H,Nu);
HH = H'*H;

% ZF�E�G�C�g
%Wzf = H\I;
%Wzf = H'/(H*H');

% MMSE�E�G�C�g
%Wmmse = H'/(H*H'+ a*eye(Nr*Nu)); % MMSE

% for inu = 1:Nu
%     %h = (randn(2,nt) + 1j*randn(2,nt))/sqrt(2);
%     h = H(1+(inu-1)*Nr:2+(inu-1)*Nr,:);
%     S = [svd(h(1,:));svd(h(2,:))];  % nr =2
%     [~,No] = max(S);
%     Hu_as(:,:,inu) = h(No,:);
% end

%H(�`���`���l���s��)
%H_as = alluser(Hu_as); %�`���`���l���s��

He = zeros(Nt-Nr,Nt,Nu);

for inu = 1:Nu
    
    %He(H���烆�[�Uk�̃`���l���s����������s��)
    el = 1:Nu;
    el(:,inu) = [];
    He(:,:,inu) = alluser(Hu(:,:,el));
    HeHe = He(:,:,inu)'*He(:,:,inu);
    %He_as(:,:,inu) = alluser(Hu_as(:,:,el));

%     % BD
%     %Ven(He�̎G��������ԂɑΉ�����ŗL�x�N�g��)
%     [~,~,Ve] = svd(He(:,:,inu));
%     Ven(:,:,inu) = Ve(:,Nr*(Nu-1)+1:Nt);
% 
%     %Vts(Hu*Ven�̐M��������ԂɑΉ�����ŗL�x�N�g��)
%     %Wtu(���[�Uk�ɑΉ������n�Ǒ��E�G�C�g)
%     %St(���ْl���i�[ 2�~nu)
%     %Wr(���[�U���E�G�C�g)
%     [Ut,Std,Vt] = svd(Hu(:,:,inu)*Ven(:,:,inu));
%     St_BD(:,inu) = diag(Std(1:Nr,1:Nr));
%     
%     Wr_BD(:,:,inu) = Ut';
%     Vts_BD(:,:,inu) = Vt(:,1:Nr);
%     Wtu_BD(:,:,inu) = Ven(:,:,inu)*Vts_BD(:,:,inu);

%     % MMSE
%     % inu�̃E�G�C�g�i��x�N�g�����ɐ��K���j
%     Wopt(:,:,inu) = Wmmse(:,(inu-1)*Nr+1:(inu-1)*Nr+Nr);
%     [Q,~]=qr(Wopt(:,:,inu));
%     QMJ(:,:,inu)=Q(:,1:Nr);
%     for ij = 1:Nr
%         Wopt(:,ij,inu) = Wopt(:,ij,inu)/sqrt(Wopt(:,ij,inu)'*Wopt(:,ij,inu));
%     end
% 
%     % �G��������ԃx�N�g����nuer�̃`���l���s��ɏ�Z                                 
%     HTT=Hu(:,:,inu)*Wopt(:,:,inu);      
%     % �ϊ��s���SVD
%     [Ut,Std,Vt]=svd(HTT);
%     St_MMSE(:,inu) = diag(Std(1:Nr,1:Nr));
%     Wr_MMSE(:,:,inu) = Ut';
%     Vts_MMSE(:,:,inu) = Vt(:,1:Nr);
%     Wtu_MMSE(:,:,inu) = Wopt(:,:,inu)*Vts_MMSE(:,:,inu);

%     % GMMSE0
%     % �G��������ԃx�N�g����nuer�̃`���l���s��ɏ�Z                                 
%     HTT=Hu(:,:,inu)*QMJ(:,:,inu);      
%     % �ϊ��s���SVD
%     [Ut,Std,Vt]=svd(HTT);
%     St_GMMSE0(:,inu) = diag(Std(1:Nr,1:Nr));
%     Wr_GMMSE0(:,:,inu) = Ut';
%     Vts_GMMSE0(:,:,inu) = Vt(:,1:Nr);
%     Wtu_GMMSE0(:,:,inu) = QMJ(:,:,inu)*Vts_GMMSE0(:,:,inu);
%     
%     % GMMSE1
%     TT1(:,:,inu)=(QMJ(:,:,inu)'*HH*QMJ(:,:,inu)+a*eye(Nr))\...
%         (QMJ(:,:,inu)'*(Hu(:,:,inu)'*Hu(:,:,inu))*QMJ(:,:,inu));
%     bb = bb + trace(TT1(:,:,inu)'*TT1(:,:,inu));
%     
%     % GMMSE2
%     LL=QMJ(:,:,inu)'*HeHe*QMJ(:,:,inu)+a*eye(Nr);
%     L=chol((LL+LL')/2);
%     TT2=eye(Nr)/L;
%     gamma = sqrt(Nr/trace(TT2'*TT2));
%     TT2=gamma*TT2;
%     
%     Wo=QMJ(:,:,inu)*TT2;
%     HTT=Hu(:,:,inu)*Wo;      
%     % �ϊ��s���SVD
%     [Ut,Std,Vt]=svd(HTT);
%     St_GMMSE2(:,inu) = diag(Std(1:Nr,1:Nr));
%     Wr_GMMSE2(:,:,inu) = Ut';
%     Vts_GMMSE2(:,:,inu) = Vt(:,1:Nr);
%     Wtu_GMMSE2(:,:,inu) = Wo*Vts_GMMSE2(:,:,inu);
    
   
    % BMSN1-BF
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
    St_BMSN1(:,inu) = diag(Std(1:Nr,1:Nr));
    
    Wr_BMSN1(:,:,inu) = Ut';
    Vts_BMSN1(:,:,inu) = Vt(:,1:Nr);
    Wtu_BMSN1(:,:,inu) = Wopt*Vts_BMSN1(:,:,inu);

    % BMSN2-BF
    %MSN�Ɋ�Â����[�U���̃E�G�C�g�̌v�Z
    % MSN�̍œK�E�G�C�gWopt
    Wopt = (HeHe+a*eye(size(He,2)))\Hu(:,:,inu)'*T(:,:,inu);
    Wopt = Wopt/norm(Wopt,'fro')*sqrt(Nr);
%     for ij = 1:Nr
%         Wopt(:,ij,inu) = Wopt(:,ij,inu)/sqrt(Wopt(:,ij,inu)'*Wopt(:,ij,inu));
%     end
    %Vts(Hu*Wopt�̐M��������ԂɑΉ�����ŗL�x�N�g��)
    %Wtu(���[�Uk�ɑΉ������n�Ǒ��E�G�C�g)
    %St(���ْl���i�[ 2�~nu)
    %Wr(���[�U���E�G�C�g)
    [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt);
    St_BMSN2(:,inu) = diag(Std(1:Nr,1:Nr));
    
    Wr_BMSN2(:,:,inu) = Ut';
    Vts_BMSN2(:,:,inu) = Vt(:,1:Nr);
    Wtu_BMSN2(:,:,inu) = Wopt*Vts_BMSN2(:,:,inu);
   
    % B-MSN2-M2 ��=sigma2 NT equal to MMSE
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

    % BMSN1-GE :��ʉ��ŗL�l�����������@(�[���G������j
    A = HeHe+a*eye(size(He(:,:,inu),2));
    B = Hu(:,:,inu)'*Hu(:,:,inu);
    [EW,D] = eig(B,A);
    [~,IN] = sort(diag(abs(D)).','descend');
    EW = EW(:,IN);
    Wopt=EW(:,1:Nr);
    for ij = 1:Nr
        Wopt(:,ij) = Wopt(:,ij)/sqrt(Wopt(:,ij)'*Wopt(:,ij));
    end
    %Wopt'*A*Wopt
    %Vts(HT*W�̐M��������ԂɑΉ�����ŗL�x�N�g��)
    %Wtu(���[�Uk�ɑΉ������n�Ǒ��E�G�C�g)
    %St(���ْl���i�[ 2�~nu)
    %Wr(���[�U���E�G�C�g)
    [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt);
    St_BMSN1g(:,inu) = diag(Std(1:Nr,1:Nr));
    Wr_BMSN1g(:,:,inu) = Ut';
    Vts_BMSN1g(:,:,inu) = Vt(:,1:Nr);
    Wtu_BMSN1g(:,:,inu) = Wopt*Vts_BMSN1g(:,:,inu);
    
    % BMSN2-GE :��ʉ��ŗL�l�����������@(�[���G������j
    A = HeHe+a*eye(size(He(:,:,inu),2));
    B = Hu(:,:,inu)'*Hu(:,:,inu);
    [EW,D] = eig(B,A);
    %EW'*A*EW
    [~,IN] = sort(diag(abs(D)).','descend');
    EW = EW(:,IN);
    Wopt=EW(:,1:Nr);
%     for ij = 1:Nr
%         Wopt(:,ij) = Wopt(:,ij)/sqrt(Wopt(:,ij)'*Wopt(:,ij));
%     end             
    Wopt = Wopt/norm(Wopt,'fro')*sqrt(Nr);
    %Vts(HT*W�̐M��������ԂɑΉ�����ŗL�x�N�g��)
    %Wtu(���[�Uk�ɑΉ������n�Ǒ��E�G�C�g)
    %St(���ْl���i�[ 2�~nu)
    %Wr(���[�U���E�G�C�g)
    [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt);
    St_BMSN2g(:,inu) = diag(Std(1:Nr,1:Nr));
    Wr_BMSN2g(:,:,inu) = Ut';
    Vts_BMSN2g(:,:,inu) = Vt(:,1:Nr);
    Wtu_BMSN2g(:,:,inu) = Wopt*Vts_BMSN2g(:,:,inu);

    % BMSN2-GE-SINR ��=sigma2 NT equal to MMSE
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
    St_BMSN2gS(:,inu) = diag(Std(1:Nr,1:Nr));
    Wr_BMSN2gS(:,:,inu) = Ut';
    Vts_BMSN2gS(:,:,inu) = Vt(:,1:Nr);
    Wtu_BMSN2gS(:,:,inu) = Wopt*Vts_BMSN2gS(:,:,inu);

    % BMSN2-GE-logSINR ��=sigma2 NT equal to MMSE
    A = HeHe+a*eye(size(He(:,:,inu),2));
    B = Hu(:,:,inu)'*Hu(:,:,inu);
    [EW,D] = eig(B,A);
    [D1,IN] = sort(diag(abs(D)).','descend');
    EW = EW(:,IN);
    D2 = diag(D1(:,1:Nr));
    EWS=EW(:,1:Nr);
%     for ij = 1:Nr
%         EWS(:,ij) = EWS(:,ij)/sqrt(EWS(:,ij)'*EWS(:,ij));
%     end          
    Wopt=EWS*sqrt(log2(1+D2));
    %Wopt=EWS(:,1:Nr)*sqrt(D2);
    Wopt = Wopt/norm(Wopt,'fro')*sqrt(Nr);
%     Wopt_k2(:,:,inu) = Wopt_k2(:,:,inu)/norm(Wopt_k2(:,:,inu),'fro')*sqrt(Nr);
    %Vts(HT*W�̐M��������ԂɑΉ�����ŗL�x�N�g��)
    %Wtu(���[�Uk�ɑΉ������n�Ǒ��E�G�C�g)
    %St(���ْl���i�[ 2�~nu)
    %Wr(���[�U���E�G�C�g)
    [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt);
    St_BMSN2gLS(:,inu) = diag(Std(1:Nr,1:Nr));
    Wr_BMSN2gLS(:,:,inu) = Ut';
    Vts_BMSN2gLS(:,:,inu) = Vt(:,1:Nr);
    Wtu_BMSN2gLS(:,:,inu) = Wopt*Vts_BMSN2gLS(:,:,inu);
    
    % BMSN2-GE-SINR1 ��=sigma2 NT equal to MMSE
    B = HeHe+a*eye(size(He(:,:,inu),2));
    A = Hu(:,:,inu)'*Hu(:,:,inu);
    [EW,D] = eig(B,A);
    [~,IN] = sort(diag(D));
    EW = EW(:,IN);
    W_k = EW(:,1)/sqrt(EW(:,1)'*EW(:,1));
    for ij = 1:Nr
         Wopt(:,ij) = W_k;
    end             
    %Wopt = Wopt/norm(Wopt,'fro')*sqrt(Nr);
    %Vts(HT*W�̐M��������ԂɑΉ�����ŗL�x�N�g��)
    %Wtu(���[�Uk�ɑΉ������n�Ǒ��E�G�C�g)
    %St(���ْl���i�[ 2�~nu)
    %Wr(���[�U���E�G�C�g)
    [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt);
    St_BMSN2gS1(:,inu) = diag(Std(1:Nr,1:Nr));
    Wr_BMSN2gS1(:,:,inu) = Ut';
    Vts_BMSN2gS1(:,:,inu) = Vt(:,1:Nr);
    Wtu_BMSN2gS1(:,:,inu) = Wopt*Vts_BMSN2gS1(:,:,inu);
    
    
%     % ZF(zero-Forcing)
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
%     % GZF(zero-Forcing)
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

%amp(���M�d��)
amp2 = 10^(snr(isnr)/10)/Nt;

% �ȍ~ nr=2 �̏ꍇ�̂�
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
end
end

% %data(���M�f�[�^) for antenna selection (AS)
% data_as = randi([0 2^bsr-1],Nu,ndata);
% 
% %s(���M�X�g���[��) for antenna selection (AS)
% s_as = Mapping(data_as,bsr);


%n(�M�G��)
%y(���[�Uk�̓��͐M��)
%Y(���[�Uk�̕����M��)
for inu = 1:Nu
    n = (randn(Nr,ndata) + 1j*randn(Nr,ndata))/sqrt(2);
    
    for ii = 1:size(pattern,1)
%         Hs_BD = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
%         Hs_ZF = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
%         Hs_GZF = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
        Hs_BMSN1 = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
        Hs_BMSN2 = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
        Hs_BMSN1g = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
        Hs_BMSN2g = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
        Hs_BMSN2gS = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
        Hs_BMSN2gLS = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
        Hs_BMSN2gS1 = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
        for ij = 1:Nu
%             Hs_BD = Hs_BD + Hu(:,:,inu)*Wtu_BD(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
%             Hs_ZF = Hs_ZF + Hu(:,:,inu)*Wtu_ZF(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
%             Hs_GZF = Hs_GZF + Hu(:,:,inu)*Wtu_GZF(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
            Hs_BMSN1 = Hs_BMSN1 + Hu(:,:,inu)*Wtu_BMSN1(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
            Hs_BMSN2 = Hs_BMSN2 + Hu(:,:,inu)*Wtu_BMSN2(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
            Hs_BMSN1g = Hs_BMSN1g + Hu(:,:,inu)*Wtu_BMSN1g(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
            Hs_BMSN2g = Hs_BMSN2g + Hu(:,:,inu)*Wtu_BMSN2g(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
            Hs_BMSN2gS = Hs_BMSN2gS + Hu(:,:,inu)*Wtu_BMSN2gS(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
            Hs_BMSN2gLS = Hs_BMSN2gLS + Hu(:,:,inu)*Wtu_BMSN2gLS(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
            Hs_BMSN2gS1 = Hs_BMSN2gS1 + Hu(:,:,inu)*Wtu_BMSN2gS1(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
        end
%         yu_BD(:,:,ii,inu) = Wr_BD(:,:,inu)*(Hs_BD+n);
%         yu_ZF(:,:,ii,inu) = Wr_ZF(:,:,inu)*(Hs_ZF+n);
%         yu_GZF(:,:,ii,inu) = Wr_GZF(:,:,inu)*(Hs_GZF+n);
        yu_BMSN1(:,:,ii,inu) = Wr_BMSN1(:,:,inu)*(Hs_BMSN1+n);
        yu_BMSN2(:,:,ii,inu) = Wr_BMSN2(:,:,inu)*(Hs_BMSN2+n);
        yu_BMSN1g(:,:,ii,inu) = Wr_BMSN1g(:,:,inu)*(Hs_BMSN1g+n);
        yu_BMSN2g(:,:,ii,inu) = Wr_BMSN2g(:,:,inu)*(Hs_BMSN2g+n);
        yu_BMSN2gS(:,:,ii,inu) = Wr_BMSN2gS(:,:,inu)*(Hs_BMSN2gS+n);
        yu_BMSN2gLS(:,:,ii,inu) = Wr_BMSN2gLS(:,:,inu)*(Hs_BMSN2gLS+n);
        yu_BMSN2gS1(:,:,ii,inu) = Wr_BMSN2gS1(:,:,inu)*(Hs_BMSN2gS1+n);
        
%         yu_BD(1,:,ii,inu) = yu_BD(1,:,ii,inu)/St_BD(1,inu)/sqrt(amp2);
%         yu_BD(2,:,ii,inu) = yu_BD(2,:,ii,inu)/St_BD(2,inu)/sqrt(amp2);
%         Yu_BD(1,:,ii,inu) = Decode(yu_BD(1,:,ii,inu),pattern(ii,1));
%         Yu_BD(2,:,ii,inu) = Decode(yu_BD(2,:,ii,inu),pattern(ii,2));
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
 
        yu_BMSN1(1,:,ii,inu) = yu_BMSN1(1,:,ii,inu)/St_BMSN1(1,inu)/sqrt(amp2);
        yu_BMSN1(2,:,ii,inu) = yu_BMSN1(2,:,ii,inu)/St_BMSN1(2,inu)/sqrt(amp2);
        Yu_BMSN1(1,:,ii,inu) = Decode(yu_BMSN1(1,:,ii,inu),pattern(ii,1));
        Yu_BMSN1(2,:,ii,inu) = Decode(yu_BMSN1(2,:,ii,inu),pattern(ii,2));
        
        yu_BMSN2(1,:,ii,inu) = yu_BMSN2(1,:,ii,inu)/St_BMSN2(1,inu)/sqrt(amp2);
        yu_BMSN2(2,:,ii,inu) = yu_BMSN2(2,:,ii,inu)/St_BMSN2(2,inu)/sqrt(amp2);
        Yu_BMSN2(1,:,ii,inu) = Decode(yu_BMSN2(1,:,ii,inu),pattern(ii,1));
        Yu_BMSN2(2,:,ii,inu) = Decode(yu_BMSN2(2,:,ii,inu),pattern(ii,2));
        
        yu_BMSN1g(1,:,ii,inu) = yu_BMSN1g(1,:,ii,inu)/St_BMSN1g(1,inu)/sqrt(amp2);
        yu_BMSN1g(2,:,ii,inu) = yu_BMSN1g(2,:,ii,inu)/St_BMSN1g(2,inu)/sqrt(amp2);
        Yu_BMSN1g(1,:,ii,inu) = Decode(yu_BMSN1g(1,:,ii,inu),pattern(ii,1));
        Yu_BMSN1g(2,:,ii,inu) = Decode(yu_BMSN1g(2,:,ii,inu),pattern(ii,2));

        yu_BMSN2g(1,:,ii,inu) = yu_BMSN2g(1,:,ii,inu)/St_BMSN2g(1,inu)/sqrt(amp2);
        yu_BMSN2g(2,:,ii,inu) = yu_BMSN2g(2,:,ii,inu)/St_BMSN2g(2,inu)/sqrt(amp2);
        Yu_BMSN2g(1,:,ii,inu) = Decode(yu_BMSN2g(1,:,ii,inu),pattern(ii,1));
        Yu_BMSN2g(2,:,ii,inu) = Decode(yu_BMSN2g(2,:,ii,inu),pattern(ii,2));
        
        yu_BMSN2gS(1,:,ii,inu) = yu_BMSN2gS(1,:,ii,inu)/St_BMSN2gS(1,inu)/sqrt(amp2);
        yu_BMSN2gS(2,:,ii,inu) = yu_BMSN2gS(2,:,ii,inu)/St_BMSN2gS(2,inu)/sqrt(amp2);
        Yu_BMSN2gS(1,:,ii,inu) = Decode(yu_BMSN2gS(1,:,ii,inu),pattern(ii,1));
        Yu_BMSN2gS(2,:,ii,inu) = Decode(yu_BMSN2gS(2,:,ii,inu),pattern(ii,2));
        
        yu_BMSN2gLS(1,:,ii,inu) = yu_BMSN2gLS(1,:,ii,inu)/St_BMSN2gLS(1,inu)/sqrt(amp2);
        yu_BMSN2gLS(2,:,ii,inu) = yu_BMSN2gLS(2,:,ii,inu)/St_BMSN2gLS(2,inu)/sqrt(amp2);
        Yu_BMSN2gLS(1,:,ii,inu) = Decode(yu_BMSN2gLS(1,:,ii,inu),pattern(ii,1));
        Yu_BMSN2gLS(2,:,ii,inu) = Decode(yu_BMSN2gLS(2,:,ii,inu),pattern(ii,2));
        
        yu_BMSN2gS1(1,:,ii,inu) = yu_BMSN2gS1(1,:,ii,inu)/St_BMSN2gS1(1,inu)/sqrt(amp2);
        yu_BMSN2gS1(2,:,ii,inu) = yu_BMSN2gS1(2,:,ii,inu)/St_BMSN2gS1(2,inu)/sqrt(amp2);
        Yu_BMSN2gS1(1,:,ii,inu) = Decode(yu_BMSN2gS1(1,:,ii,inu),pattern(ii,1));
        Yu_BMSN2gS1(2,:,ii,inu) = Decode(yu_BMSN2gS1(2,:,ii,inu),pattern(ii,2));
        
    end
end

%ber(�ϒ��p�^�[������BER)  for BD and B-MSN
for inu = 1:Nu
for ii = 1:size(pattern,1)
if ii == 1
%     ber_BD(ii,inu) = BER2(BER1(Yu_BD(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
%     ber_ZF(ii,inu) = BER2(BER1(Yu_ZF(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
%     ber_GZF(ii,inu) = BER2(BER1(Yu_GZF(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
    ber_BMSN1(ii,inu) = BER2(BER1(Yu_BMSN1(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
    ber_BMSN2(ii,inu) = BER2(BER1(Yu_BMSN2(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
    ber_BMSN1g(ii,inu) = BER2(BER1(Yu_BMSN1g(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
    ber_BMSN2g(ii,inu) = BER2(BER1(Yu_BMSN2g(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
    ber_BMSN2gS(ii,inu) = BER2(BER1(Yu_BMSN2gS(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
    ber_BMSN2gLS(ii,inu) = BER2(BER1(Yu_BMSN2gLS(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
    ber_BMSN2gS1(ii,inu) = BER2(BER1(Yu_BMSN2gS1(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
else
%     ber_BD(ii,inu) = BER2(BER1(Yu_BD(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
%         + BER2(BER1(Yu_BD(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
%     ber_ZF(ii,inu) = BER2(BER1(Yu_ZF(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
%         + BER2(BER1(Yu_ZF(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
%     ber_GZF(ii,inu) = BER2(BER1(Yu_GZF(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
%         + BER2(BER1(Yu_GZF(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
    ber_BMSN1(ii,inu) = BER2(BER1(Yu_BMSN1(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
        + BER2(BER1(Yu_BMSN1(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
    ber_BMSN2(ii,inu) = BER2(BER1(Yu_BMSN2(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
        + BER2(BER1(Yu_BMSN2(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
    ber_BMSN1g(ii,inu) = BER2(BER1(Yu_BMSN1g(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
        + BER2(BER1(Yu_BMSN1g(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
    ber_BMSN2g(ii,inu) = BER2(BER1(Yu_BMSN2g(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
        + BER2(BER1(Yu_BMSN2g(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
    ber_BMSN2gS(ii,inu) = BER2(BER1(Yu_BMSN2gS(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
        + BER2(BER1(Yu_BMSN2gS(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
    ber_BMSN2gLS(ii,inu) = BER2(BER1(Yu_BMSN2gLS(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
        + BER2(BER1(Yu_BMSN2gLS(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
    ber_BMSN2gS1(ii,inu) = BER2(BER1(Yu_BMSN2gS1(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
        + BER2(BER1(Yu_BMSN2gS1(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
end
end
end

%ber_min(���ŏ��̕ϒ��p�^�[����BER) for BD and B-MSN
% ber_min_BD(isimu,isnr) = mean(min(ber_BD));
% ber_min_ZF(isimu,isnr) = mean(min(ber_ZF));
% ber_min_GZF(isimu,isnr) = mean(min(ber_GZF));
ber_min_BMSN1(isimu,isnr) = mean(min(ber_BMSN1));
ber_min_BMSN2(isimu,isnr) = mean(min(ber_BMSN2));
ber_min_BMSN1g(isimu,isnr) = mean(min(ber_BMSN1g));
ber_min_BMSN2g(isimu,isnr) = mean(min(ber_BMSN2g));
ber_min_BMSN2gS(isimu,isnr) = mean(min(ber_BMSN2gS));
ber_min_BMSN2gLS(isimu,isnr) = mean(min(ber_BMSN2gLS));
ber_min_BMSN2gS1(isimu,isnr) = mean(min(ber_BMSN2gS1));
%fprintf('Iteration = %d / %d\n',isimu, SIMU);

end % isimu

fprintf('SNR = %d dB\n',snr(isnr));
end % isnr


%BER(���s�񐔕���)
% BER_BD = mean(ber_min_BD);
% BER_ZF = mean(ber_min_ZF);
% BER_GZF = mean(ber_min_GZF);
BER_BMSN1 = mean(ber_min_BMSN1);
BER_BMSN2 = mean(ber_min_BMSN2);
BER_BMSN1g = mean(ber_min_BMSN1g);
BER_BMSN2g = mean(ber_min_BMSN2g);
BER_BMSN2gS = mean(ber_min_BMSN2gS);
BER_BMSN2gLS = mean(ber_min_BMSN2gLS);
BER_BMSN2gS1 = mean(ber_min_BMSN2gS1);

% % csvwrite(evfile1,[snr;BER_BD]);
% % csvwrite(evfile2,[snr;BER_ZF]);
% % csvwrite(evfile3,[snr;BER_GZF]);
% % csvwrite(evfile6,[snr;BER_BMSN2g2]);
csvwrite(evfile1,[snr;BER_BMSN1g]);
csvwrite(evfile2,[snr;BER_BMSN2g]);
csvwrite(evfile3,[snr;BER_BMSN2gS]);
csvwrite(evfile4,[snr;BER_BMSN2gLS]);
csvwrite(evfile5,[snr;BER_BMSN2gS1]);
csvwrite(evfile6,[snr;BER_BMSN1]);
csvwrite(evfile7,[snr;BER_BMSN2]);

toc;

%�`��
figure;
semilogy(snr,BER_BMSN1g,'-m',snr,BER_BMSN2g,'-b',snr,BER_BMSN2gS,'-k',snr,BER_BMSN2gLS,'-g',snr,BER_BMSN2gS1,'-r','Linewidth',2);
axis([SNR_min SNR_max 1e-4 1e0]);
set(gca,'XTick',SNR_min:5:SNR_max,'YTick',[1e-4, 1e-3, 1e-2, 1e-1 1e0],'Fontsize',14,'Fontname','Times New Roman')
legend('BMSN1-GE','BMSN2-GE','BMSN2-GES','BMSN2-GELS','BMSN2-GES1','Location','southwest')
xlabel('SNR [dB]','Fontsize',16,'Fontname','Times New Roman');
ylabel('Average BER','Fontsize',16,'Fontname','Times New Roman');
%title('iid Rayleigh channel')
grid on;
hold on;


