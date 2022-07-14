clear;

rng('default');

% �o�̓t�@�C����
folder = 'CSV_BER/';
evfile1 = 'BER16x2x8u_BD_SNR0-30dB_10000bits_5000itr.csv';
evfile2 = 'BER16x2x8u_BDAS_SNR0-30dB_10000bits_5000itr.csv';
evfile3 = 'BER16x2x8u_BMSN2_SNR0-30dB_10000bits_5000itr.csv';
evfile4 = 'BER16x2x8u_BMSNAS2_SNR0-30dB_10000bits_5000itr.csv';
evfile5 = 'BER16x2x8u_ZF_SNR0-30dB_10000bits_5000itr.csv';
evfile6 = 'BER16x2x8u_BMSN1b_SNR0-30dB_10000bits_5000itr.csv';

Nt = 16; %��n�ǃA���e�i��
Nr = 2; %���[�U�A���e�i����2(=���M�X�g���[����(not1))
Nu = 8; %���[�U��
ndata = 10; % 1�p�P�b�g������̃V���{�����i�p�P�b�g���j:BER�v�Z�̂��߂�1���s������̃V���{���� 10000
SNR_min = 0;  %�ŏ�SNR[dB]
SNR_max = 30; %�ő�SNR[dB]
SIMU = 10; %���s�񐔁i�ʏ� 1000,5000�j
pattern = [4,0;2,2;3,1]; %�ϒ��p�^�[��(�擪��[��,0]�������Ă���j
bsr = pattern(1,1);
a = 1e-2; % BMSN�̋[���G��
I = eye(Nt,Nt);

%T ���]�̃`���l���s��
for inu = 1:Nu
    T(:,:,inu) = zeros(Nr,Nr);
    T(1,1,inu) = 1;
    T(1,2,inu) = 0;
    T(2,1,inu) = 0;
    T(2,2,inu) = 1;
end

tic;

snr = SNR_min:SNR_max;
Lsnr = length(snr);
% MSt_BD1=zeros(Nr,Lsnr);
% MSt_BD2=zeros(Nr-1,Lsnr);
% MSt_MSN1=zeros(Nr,Lsnr);
% MSt_MSN2=zeros(Nr-1,Lsnr);
% MSt_ZF=zeros(Nr,Lsnr);
%MSt_MMSE=zeros(Nr,Lsnr);
    
for isnr = 1:Lsnr
    sigma2 = 1/(10^(snr(isnr)/10)); % noise power
    ao = sigma2*Nt; % adaptive pseudo-noise power for BMSN

for isimu = 1:SIMU

%H(�`���`���l���s��: iid Rayleigh channel)
%Hu(���[�Uk�̃`���l���s��Hu(:,:,k))
H = (randn(Nr*Nu,Nt) + 1j*randn(Nr*Nu,Nt))/sqrt(2);
Hu = peruser(H,Nu);

% ZF�E�G�C�g
Wzf = H\I;
%Wzf = H'/(H*H');

% MMSE�E�G�C�g
%Wmmse = H'/(H*H'+ sigma2*Nu*eye(Nr*Nu)); % MMSE

for inu = 1:Nu
    %h = (randn(2,nt) + 1j*randn(2,nt))/sqrt(2);
    h = H(1+(inu-1)*Nr:2+(inu-1)*Nr,:);
    S = [svd(h(1,:));svd(h(2,:))];  % nr =2
    [~,No] = max(S);
    Hu_as(:,:,inu) = h(No,:);
end

%H(�`���`���l���s��)
H_as = alluser(Hu_as); %�`���`���l���s��

He = zeros(14,16,8);

for inu = 1:Nu
    
    %He(H���烆�[�Uk�̃`���l���s����������s��)
    el = 1:Nu;
    el(:,inu) = [];
    He(:,:,inu) = alluser(Hu(:,:,el));
    He_as(:,:,inu) = alluser(Hu_as(:,:,el));

    % BD
    %Ven(He�̎G��������ԂɑΉ�����ŗL�x�N�g��)
    [~,~,Ve] = svd(He(:,:,inu));
    Ven1(:,:,inu) = Ve(:,Nr*(Nu-1)+1:Nt);

    %Vts(Hu*Ven�̐M��������ԂɑΉ�����ŗL�x�N�g��)
    %Wtu(���[�Uk�ɑΉ������n�Ǒ��E�G�C�g)
    %St(���ْl���i�[ 2�~nu)
    %Wr(���[�U���E�G�C�g)
    [Ut,Std,Vt] = svd(Hu(:,:,inu)*Ven1(:,:,inu));
    St_BD1(:,inu) = diag(Std(1:Nr,1:Nr));
    
    Wr_BD1(:,:,inu) = Ut';
    Vts_BD1(:,:,inu) = Vt(:,1:Nr);
    Wtu_BD1(:,:,inu) = Ven1(:,:,inu)*Vts_BD1(:,:,inu);
    
    % B-MSN
    %MSN�Ɋ�Â����[�U���̃E�G�C�g�̌v�Z
    % MSN�̍œK�E�G�C�gWopt
    Wopt(:,:,inu) = (He(:,:,inu)'*He(:,:,inu)+a*eye(size(He,2)))\Hu(:,:,inu)'*T(:,:,inu);
%     Wopt(:,:,i) = Wopt(:,:,i)/norm(Wopt(:,:,i),'fro');
    for ij = 1:Nr
        Wopt(:,ij,inu) = Wopt(:,ij,inu)/sqrt(Wopt(:,ij,inu)'*Wopt(:,ij,inu));
    end
    
    %Vts(Hu*Wopt�̐M��������ԂɑΉ�����ŗL�x�N�g��)
    %Wtu(���[�Uk�ɑΉ������n�Ǒ��E�G�C�g)
    %St(���ْl���i�[ 2�~nu)
    %Wr(���[�U���E�G�C�g)
    [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt(:,:,inu));
    St_BMSN1(:,inu) = diag(Std(1:Nr,1:Nr));
    
    Wr_BMSN1(:,:,inu) = Ut';
    Vts_BMSN1(:,:,inu) = Vt(:,1:Nr);
    Wtu_BMSN1(:,:,inu) = Wopt(:,:,inu)*Vts_BMSN1(:,:,inu);

    % BD with antenna selection
    %Ven(He�̎G��������ԂɑΉ�����ŗL�x�N�g��)
    [~,~,Ve] = svd(He_as(:,:,inu));
    Ven2(:,:,inu) = Ve(:,Nu:Nt);

    %Vts(Hu*Ven�̐M��������ԂɑΉ�����ŗL�x�N�g��)
    %Wtu(���[�Uk�ɑΉ������n�Ǒ��E�G�C�g)
    %St(���ْl���i�[ 1�~nu)
    %Wr(���[�U���E�G�C�g)
    [Ut,Std,Vt] = svd(Hu_as(:,:,inu)*Ven2(:,:,inu));
    St_BD2(:,inu) = Std(1,1);
    
    Wr_BD2(:,:,inu) = Ut';
    Vts_BD2(:,:,inu) = Vt(:,1);
    Wtu_BD2(:,:,inu) = Ven2(:,:,inu)*Vts_BD2(:,:,inu);
    
    % B-MSN with antenna selection
    %MSN�Ɋ�Â����[�U���̃E�G�C�g�̌v�Z
    % MSN�̍œK�E�G�C�gWopt
    Wopt2(:,:,inu) = (He_as(:,:,inu)'*He_as(:,:,inu)+ao*eye(size(He_as,2)))\Hu_as(:,:,inu)'*T(1,1,inu);
%     Wopt(:,:,i) = Wopt(:,:,i)/norm(Wopt(:,:,i),'fro');
    Wopt2(:,1,inu) = Wopt2(:,1,inu)/sqrt(Wopt2(:,1,inu)'*Wopt2(:,1,inu));

    %Vts(Hu*Ven�̐M��������ԂɑΉ�����ŗL�x�N�g��)
    %Wtu(���[�Uk�ɑΉ������n�Ǒ��E�G�C�g)
    %St(���ْl���i�[ 1�~nu)
    %Wr(���[�U���E�G�C�g)
    [Ut,Std,Vt] = svd(Hu_as(:,:,inu)*Wopt2(:,:,inu));
    St_BMSN2(:,inu) = Std(1,1);
    
    Wr_BMSN2(:,:,inu) = Ut';
    Vts_BMSN2(:,:,inu) = Vt(:,1);
    Wtu_BMSN2(:,:,inu) = Wopt2(:,:,inu)*Vts_BMSN2(:,:,inu);

    % BMSN :��ʉ��ŗL�l�����������@
    B(:,:,inu) = He(:,:,inu)'*He(:,:,inu);
    A(:,:,inu) = Hu(:,:,inu)'*Hu(:,:,inu);
    [Wopt_k(:,:,inu),D] = eig(B(:,:,inu),A(:,:,inu));
    [D1,IN] = sort(diag(D));
    D2(:,:,inu) = D1.';
    Wopt_k(:,:,inu) = Wopt_k(:,IN,inu);
    Wopt_k2(:,:,inu)=Wopt_k(:,1:Nr,inu);
    for ij = 1:Nr
        Wopt_k2(:,ij,inu) = Wopt_k2(:,ij,inu)/sqrt(Wopt_k2(:,ij,inu)'*Wopt_k2(:,ij,inu));
    end             
    % nu�̃E�G�C�g(�M��������Ԃ𗘗p)                                 
    %Vts(HT*W�̐M��������ԂɑΉ�����ŗL�x�N�g��)
    %Wtu(���[�Uk�ɑΉ������n�Ǒ��E�G�C�g)
    %St(���ْl���i�[ 2�~nu)
    %Wr(���[�U���E�G�C�g)
    [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt_k2(:,:,inu));
    St_BMSN1b(:,inu) = diag(Std(1:Nr,1:Nr));
    Wr_BMSN1b(:,:,inu) = Ut';
    Vts_BMSN1b(:,:,inu) = Vt(:,1:Nr);
    Wtu_BMSN1b(:,:,inu) = Wopt_k2(:,:,inu)*Vts_BMSN1b(:,:,inu);
    
    % ZF(zero-Forcing)
    % inu�̃E�G�C�g�i��x�N�g�����ɐ��K���j
    Wopt(:,:,inu) = Wzf(:,(inu-1)*Nr+1:(inu-1)*Nr+Nr);
    for ij = 1:Nr
        Wopt(:,ij,inu) = Wopt(:,ij,inu)/sqrt(Wopt(:,ij,inu)'*Wopt(:,ij,inu));
    end

    % �G��������ԃx�N�g����nuer�̃`���l���s��ɏ�Z                                 
    HTT=Hu(:,:,inu)*Wopt(:,:,inu);      
    % �ϊ��s���SVD
    [Ut,Std,Vt]=svd(HTT);
    St_ZF(:,inu) = diag(Std(1:Nr,1:Nr));
    Wr_ZF(:,:,inu) = Ut';
    Vts_ZF(:,:,inu) = Vt(:,1:Nr);
    Wtu_ZF(:,:,inu) = Wopt(:,:,inu)*Vts_ZF(:,:,inu);
               
    % MMSE
    % inu�̃E�G�C�g�i��x�N�g�����ɐ��K���j
%     Wopt(:,:,inu) = Wmmse(:,(inu-1)*Nr+1:(inu-1)*Nr+Nr);
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
    
end

%Wt(��n�Ǒ��E�G�C�g)
%Wt = reshape(Wtu,[nt,nt]);

%amp(���M�d��)
amp = 10^(snr(isnr)/10)/Nt;

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

%data(���M�f�[�^) for antenna selection (AS)
data_as = randi([0 2^bsr-1],Nu,ndata);

%s(���M�X�g���[��) for antenna selection (AS)
s_as = Mapping(data_as,bsr);


%n(�M�G��)
%y(���[�Uk�̓��͐M��)
%Y(���[�Uk�̕����M��)
for inu = 1:Nu
    n = (randn(Nr,ndata) + 1j*randn(Nr,ndata))/sqrt(2);
    
    %  for BD and BMSN with antenna selection (BD-AS, BMSN-AS)
    Hs_BD2 = zeros(1,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
    Hs_BMSN2 = zeros(1,ndata);
    for ij = 1:Nu
        Hs_BD2 = Hs_BD2 + Hu_as(:,:,inu)*Wtu_BD2(:,:,ij)*s_as(ij,:)*sqrt(amp);
        Hs_BMSN2 = Hs_BMSN2 + Hu_as(:,:,inu)*Wtu_BMSN2(:,:,ij)*s_as(ij,:)*sqrt(amp);
    end
    yu_BD2(:,:,inu) = Wr_BD2(:,:,inu)*(Hs_BD2+n(1,:))/sqrt(amp)/St_BD2(1,inu);
    yu_BMSN2(:,:,inu) = Wr_BMSN2(:,:,inu)*(Hs_BMSN2+n(1,:))/sqrt(amp)/St_BMSN2(1,inu);

    for ii = 1:size(pattern,1)
        Hs_BD1 = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
        Hs_BMSN1 = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
        Hs_ZF = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
        %Hs_MMSE = zeros(Nr,ndata);    % ���[�Ui�̎�M�M���i��M�d�ݕt���O�j
        for ij = 1:Nu
            Hs_BD1 = Hs_BD1 + Hu(:,:,inu)*Wtu_BD1(:,:,ij)*s(:,:,ii,ij)*sqrt(amp);
            Hs_BMSN1 = Hs_BMSN1 + Hu(:,:,inu)*Wtu_BMSN1(:,:,ij)*s(:,:,ii,ij)*sqrt(amp);
            Hs_ZF = Hs_ZF + Hu(:,:,inu)*Wtu_ZF(:,:,ij)*s(:,:,ii,ij)*sqrt(amp);
            %Hs_MMSE = Hs_MMSE + Hu(:,:,inu)*Wtu_MMSE(:,:,ij)*s(:,:,ii,ij)*sqrt(amp);
        end
        yu_BD1(:,:,ii,inu) = Wr_BD1(:,:,inu)*(Hs_BD1+n);
        yu_BMSN1(:,:,ii,inu) = Wr_BMSN1(:,:,inu)*(Hs_BMSN1+n);
        yu_ZF(:,:,ii,inu) = Wr_ZF(:,:,inu)*(Hs_ZF+n);
        %yu_MMSE(:,:,ii,inu) = Wr_MMSE(:,:,inu)*(Hs_MMSE+n);

        yu_BD1(1,:,ii,inu) = yu_BD1(1,:,ii,inu)/St_BD1(1,inu)/sqrt(amp);
        yu_BD1(2,:,ii,inu) = yu_BD1(2,:,ii,inu)/St_BD1(2,inu)/sqrt(amp);
        Yu_BD1(1,:,ii,inu) = Decode(yu_BD1(1,:,ii,inu),pattern(ii,1));
        Yu_BD1(2,:,ii,inu) = Decode(yu_BD1(2,:,ii,inu),pattern(ii,2));
    
        yu_BMSN1(1,:,ii,inu) = yu_BMSN1(1,:,ii,inu)/St_BMSN1(1,inu)/sqrt(amp);
        yu_BMSN1(2,:,ii,inu) = yu_BMSN1(2,:,ii,inu)/St_BMSN1(2,inu)/sqrt(amp);
        Yu_BMSN1(1,:,ii,inu) = Decode(yu_BMSN1(1,:,ii,inu),pattern(ii,1));
        Yu_BMSN1(2,:,ii,inu) = Decode(yu_BMSN1(2,:,ii,inu),pattern(ii,2));
        
        yu_ZF(1,:,ii,inu) = yu_ZF(1,:,ii,inu)/St_ZF(1,inu)/sqrt(amp);
        yu_ZF(2,:,ii,inu) = yu_ZF(2,:,ii,inu)/St_ZF(2,inu)/sqrt(amp);
        Yu_ZF(1,:,ii,inu) = Decode(yu_ZF(1,:,ii,inu),pattern(ii,1));
        Yu_ZF(2,:,ii,inu) = Decode(yu_ZF(2,:,ii,inu),pattern(ii,2));
        
%         yu_MMSE(1,:,ii,inu) = yu_MMSE(1,:,ii,inu)/St_MMSE(1,inu)/sqrt(amp);
%         yu_MMSE(2,:,ii,inu) = yu_MMSE(2,:,ii,inu)/St_MMSE(2,inu)/sqrt(amp);
%         Yu_MMSE(1,:,ii,inu) = Decode(yu_MMSE(1,:,ii,inu),pattern(ii,1));
%         Yu_MMSE(2,:,ii,inu) = Decode(yu_MMSE(2,:,ii,inu),pattern(ii,2));
        
    end
end

% for BD and BMSN with antenna selection (BD-AS, BMSN-AS)
y_bdas = alluser(yu_BD2);
Y_bdas = Decode(y_bdas,bsr);
y_msnas = alluser(yu_BMSN2);
Y_msnas = Decode(y_msnas,bsr);
%ber
ber_BD2(isimu,isnr) = BER2(BER1(Y_bdas,data_as,bsr),bsr);
ber_MSN2(isimu,isnr) = BER2(BER1(Y_msnas,data_as,bsr),bsr);

%ber(�ϒ��p�^�[������BER)  for BD and B-MSN
for inu = 1:Nu
for ii = 1:size(pattern,1)
if ii == 1
    ber_BD1(ii,inu) = BER2(BER1(Yu_BD1(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
    ber_BMSN1(ii,inu) = BER2(BER1(Yu_BMSN1(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
    ber_ZF(ii,inu) = BER2(BER1(Yu_ZF(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
    %ber_MMSE(ii,inu) = BER2(BER1(Yu_MMSE(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
else
    ber_BD1(ii,inu) = BER2(BER1(Yu_BD1(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
        + BER2(BER1(Yu_BD1(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
    ber_BMSN1(ii,inu) = BER2(BER1(Yu_BMSN1(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
        + BER2(BER1(Yu_BMSN1(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
    ber_ZF(ii,inu) = BER2(BER1(Yu_ZF(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
        + BER2(BER1(Yu_ZF(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
%     ber_MMSE(ii,inu) = BER2(BER1(Yu_MMSE(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
%         + BER2(BER1(Yu_MMSE(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
end
end
end

%ber_min(���ŏ��̕ϒ��p�^�[����BER) for BD and B-MSN
ber_min_BD1(isimu,isnr) = mean(min(ber_BD1));
ber_min_BMSN1(isimu,isnr) = mean(min(ber_BMSN1));
ber_min_ZF(isimu,isnr) = mean(min(ber_ZF));
%ber_min_MMSE(isimu,isnr) = mean(min(ber_MMSE));

%��M�A���e�i�inr)�̕��ϓ��ْl for BD and B-MSN
%10*log10(MSt_BD1(:,inudiag(S_BD(:,:,nuser)).^2/(NT*snt));   
% MSt_BD1(:,isnr)=MSt_BD1(:,isnr)+mean(St_BD1(:,inu),2);
% MSt_BD2(:,isnr)=MSt_BD2(:,isnr)+mean(St_BD2(:,inu),2);
% MSt_MSN1(:,isnr)=MSt_MSN1(:,isnr)+mean(St_BMSN1(:,inu),2);
% MSt_MSN2(:,isnr)=MSt_MSN2(:,isnr)+mean(St_BMSN2(:,inu),2);
% MSt_ZF(:,isnr)=MSt_ZF(:,isnr)+mean(St_ZF(:,inu),2);
%MSt_MMSE(:,isnr)=MSt_MMSE(:,isnr)+mean(St_MMSE(:,inu),2);

%fprintf('Iteration = %d / %d\n',isimu, SIMU);

end % isimu
% MSt_BD1(:,isnr)=10*log10((MSt_BD1(:,isnr)/SIMU).^2/(Nt*sigma2));
% MSt_BD2(:,isnr)=10*log10((MSt_BD2(:,isnr)/SIMU).^2/(Nt*sigma2));
% MSt_MSN1(:,isnr)=10*log10((MSt_MSN1(:,isnr)/SIMU).^2/(Nt*sigma2));
% MSt_MSN2(:,isnr)=10*log10((MSt_MSN2(:,isnr)/SIMU).^2/(Nt*sigma2));
% MSt_ZF(:,isnr)=10*log10((MSt_ZF(:,isnr)/SIMU).^2/(Nt*sigma2));
%MSt_MMSE(:,isnr)=10*log10((MSt_MMSE(:,isnr)/SIMU).^2/(Nt*sigma2));

fprintf('SNR = %d dB\n',snr(isnr));
end % isnr


%BER(���s�񐔕���)
BER_BD1 = mean(ber_min_BD1);
BER_BMSN1 = mean(ber_min_BMSN1);
BER_BD2 = mean(ber_BD2);
BER_BMSN2 = mean(ber_MSN2);
BER_ZF = mean(ber_min_ZF);
%BER_MMSE = mean(ber_min_MMSE);

% csvwrite(strcat(folder,evfile1),[snr;BER_BD1]);
% csvwrite(strcat(folder,evfile2),[snr;BER_BD2]);
% csvwrite(strcat(folder,evfile3),[snr;BER_BMSN1]);
% csvwrite(strcat(folder,evfile4),[snr;BER_BMSN2]);
% csvwrite(strcat(folder,evfile5),[snr;BER_ZF]);
% csvwrite(strcat(folder,evfile6),[snr;BER_MMSE]);

toc;

%�`��
figure;
semilogy(snr,BER_BD1,'-b',snr,BER_BMSN1,'-r',snr,BER_BD2,'-g',snr,...
    BER_BMSN2,'-k',snr,BER_ZF,'-m','Linewidth',2);
axis([SNR_min SNR_max 1e-4 1e0]);
set(gca,'XTick',SNR_min:5:SNR_max,'YTick',[1e-4, 1e-3, 1e-2, 1e-1 1e0],'Fontsize',14,'Fontname','Times New Roman')
legend('BD','BMSN','BD-AS','BMSN-AS','ZF','Location','southwest')
xlabel('SNR [dB]','Fontsize',16,'Fontname','Times New Roman');
ylabel('Average BER','Fontsize',16,'Fontname','Times New Roman');
%title('iid Rayleigh channel')
grid on;
hold on;


