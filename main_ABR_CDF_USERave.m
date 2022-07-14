% Achievable Bit Rate (ABR)

clear;
Nt = 16;    % 送信素子数
Nr = 2;     % 受信素子数 (アンテナ選択の場合1, そうでない場合は2)
Nu = 8;

SNR_tar = 30;  % ターゲットSNR[dB]
% SNR_max = 30; %最大SNR[dB]
SIMU = 1000; %試行回数（通常 1000）

if SNR_tar < 10
    tsnr=strcat('SNR= ',num2str(SNR_tar,'%01d'),'dB');
else
    tsnr=strcat('SNR=',num2str(SNR_tar,'%02d'),'dB');
end

I = eye(Nt,Nt);
Nru = Nr*Nu;

%T 所望のチャネル行列
for inu = 1:Nu
    T(:,:,inu) = zeros(Nr,Nr);
    T(1,1,inu) = 1;
    T(1,2,inu) = 0;
    T(2,1,inu) = 0;
    T(2,2,inu) = 1;
end

  
sigma2 = 1/(10^(SNR_tar/10)); % noise power
a = sigma2*Nt;

% 出力ファイル名 with alpha and SNR in dB
%evfile1 = strcat(folder,fn1,num2str(snr(isnr),'%02d'),'dB_1000itr.csv');
    
for isimu = 1:SIMU

    %H(伝搬チャネル行列: iid Rayleigh channel)
    %Hu(ユーザkのチャネル行列Hu(:,:,k))
    H = (randn(Nr*Nu,Nt) + 1j*randn(Nr*Nu,Nt))/sqrt(2);
    Hu = peruser(H,Nu);

    He = zeros((Nu-1)*Nr,Nt,Nu);    % Hから1ユーザのチャネル行列を除いた行列

    % BMSN3
    [~,~,STT,RIP,~] = bmsn_bf(Nt,Nr,Nu,H,a,T);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    MSt_BMSN(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
 
    % ZF
    [~,~,STT,RIP,~] = zf(Nt,Nr,Nu,H);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    MSt_ZF(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    
    % MMSE
    [~,~,STT,RIP,~] = gmmse(Nt,Nr,Nu,H,a);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    MSt_MMSE(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    
    % BD
    [~,~,STT,RIP,~] = bd(Nt,Nr,Nu,H);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
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


