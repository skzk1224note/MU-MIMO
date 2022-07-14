
% Achievable Bit Rate (ABR)

clear;
Nt = 16;    % ���M�f�q��(16)
Nr = 2;     % ��M�f�q�� (�A���e�i�I���̏ꍇ1, �����łȂ��ꍇ��2)
Nu = 8;     % ���[�U��(8)

SNR_min = 5;  %�ŏ�SNR[dB]
SNR_max = 30; %�ő�SNR[dB]
snr = SNR_min:5:SNR_max;
Lsnr = length(snr);

SIMU = 1000; %���s�񐔁i�ʏ� 1000�j

% �[���G��
a_p=[-6:1:-3 -2:0.1:2];
%a_p=-6:1:2;
alpha = 10.^a_p;
Lal=length(alpha);

% CDF
CDF=50; % CDF (%)
ICDF=round(CDF*10);

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

ABR=zeros(Lsnr,Lal);
a_opt=zeros(1,Lsnr);
St_BMSN=zeros(Nr,Nu);
MSt_BMSN=zeros(SIMU,Nru);

for ia = 1:Lal
    a=alpha(ia); % noise power
    
    for isnr = 1:Lsnr
        sigma2 = 1/(10^(snr(isnr)/10)); % noise power
        a_opt(isnr)=Nt*sigma2; % pseudo-noise
        %a_opt(isnr)=Nr*Nt*sigma2; % pseudo-noise for BMSN2
        
        for isimu = 1:SIMU

            %H(�`���`���l���s��: iid Rayleigh channel)
            %Hu(���[�Uk�̃`���l���s��Hu(:,:,k))
            H = (randn(Nr*Nu,Nt) + 1j*randn(Nr*Nu,Nt))/sqrt(2);
            Hu = peruser(H,Nu);

            He = zeros((Nu-1)*Nr,Nt,Nu);    % H����1���[�U�̃`���l���s����������s��

            [W,UTT,STT,RIP,~] = bmsn_bf(Nt,Nr,Nu,H,a,T);
            for inu=1:Nu
                St_BMSN(:,inu) = diag(STT(1:Nr,1:Nr,inu));
            end
            %��M�A���e�i�inr)�̌ŗL�l���z for B-MSN�i���g�������l���j
            MSt_BMSN(isimu,:)=((reshape(St_BMSN,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
 
            %fprintf('Iteration = %d / %d\n',isimu, SIMU);
            
        end %isimu
        
        % ABR
        for nuser=1:Nu
            Q(:,nuser) = sort(sum(log2(1 + MSt_BMSN(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
        end
        Qm = mean(Q(:,:),2);
        ABR(isnr,ia)=mean(Qm);
        %fprintf('SNR = %d dB\n',snr(isnr));
    
    end % isnr
    fprintf('alpha = %e\n',alpha(ia));
end % ia

figure;
mycol = [1 0 0;
      1 0 1;
      0 1 1;
      0 0.8 0;
      0 0 1;
      0 0 0];
set(groot,'defaultAxesColorOrder',mycol)

semilogx(alpha,ABR,'-','Linewidth',2);
axis([1e-6 1e2 0 14]);hold on
for isnr=1:Lsnr
    semilogx([a_opt(isnr),a_opt(isnr)],ylim,'--','Color',mycol(isnr,:),'Linewidth',1.2);hold on
end
set(gca,'Fontsize',14,'Fontname','Times New Roman');
xlabel('Pseudo-noise in BMSN','Fontsize',16,'Fontname','Times New Roman');
ylabel('Average achievable bit rate [bits/s/Hz]','Fontsize',16,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
legend('SNR=5dB','SNR=10dB','SNR=15dB','SNR=20dB','SNR=25dB','SNR=30dB','Location','Southwest');
title('BMSN_BF');
grid on;
hold on;

%End


