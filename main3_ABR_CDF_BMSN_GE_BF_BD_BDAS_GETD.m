% Achievable Bit Rate (ABR)

clear;
Nt = 24;    % 送信素子数
Nr = 3;     % 受信素子数
Nu = 8;

Nrs = Nr-1; % Antenna Selection

SNR_tar = 40;  % ターゲットSNR[dB]
% SNR_max = 30; %最大SNR[dB]
SIMU = 1000; %試行回数（通常 1000）

if SNR_tar < 10
    tsnr=strcat('SNR= ',num2str(SNR_tar,'%01d'),'dB');
else
    tsnr=strcat('SNR=',num2str(SNR_tar,'%02d'),'dB');
end

I = eye(Nt,Nt);
Nru = Nr*Nu;
Nru1 = Nrs*Nu;

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

for isimu = 1:SIMU

    %H(伝搬チャネル行列: iid Rayleigh channel)
    %Hu(ユーザkのチャネル行列Hu(:,:,k))
    H = (randn(Nr*Nu,Nt) + 1j*randn(Nr*Nu,Nt))/sqrt(2);
    Hu = peruser(H,Nu);

    %He = zeros((Nu-1)*Nr,Nt,Nu);    % Hから1ユーザのチャネル行列を除いた行列

%     for inu = 1:Nu
%         h = H(1+(inu-1)*Nr:2+(inu-1)*Nr,:);
%         S = [svd(h(1,:));svd(h(2,:))];  % nr =2
%         [~,No] = max(S);
%         Hu_as(:,:,inu) = h(No,:);
%     end
% 
%     %H_as(伝搬チャネル行列)
%     H_as = alluser(Hu_as); %伝搬チャネル行列

    %% BD
    [~,~,STT,RIP,~] = bd(Nt,Nr,Nu,H);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    MSt_BD(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    clear St;
    
    %% BMSN-GE-TD
    [~,~,STT,RIP,~] = bmsn_getd_3(Nt,Nr,Nu,H,a);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    MSt_BMSN_GETD_3(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
     clear St;
     
    %% BMSN-GE
    [~,~,STT,RIP,~] = bmsn_gev(Nt,Nr,Nu,H,a);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    MSt_BMSN_GE(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
     clear St;  
     
    %% BMSN-BF
    [~,~,STT,RIP,~] = bmsn_bf(Nt,Nr,Nu,H,a,T);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    MSt_BMSN_BF(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
     clear St;
    
    %% BD_AS
    [~,~,STT,RIP,~] = bd_as(Nt,Nr,Nrs,Nu,H);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nrs,1:Nrs,inu)); % STT(1,1,inu)
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    MSt_BDAS(isimu,:)=((reshape(St,[Nru1,1]).').^2)./(reshape(RIP,[Nru1,1]).'+Nt*sigma2); % ((St.').^2)./(RIP.'+Nt*sigma2);
     clear St;
    
     %% BMSN-GETD-CW
%     [~,~,STT,RIP,~] = bmsn_getd_3_CW(Nt,Nr,Nu,H,a);
%     for inu=1:Nu
%         St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
%     end
%     %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
%     MSt_BMSN_GETD_3_CW(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
%      clear St;  
end % isimu

% ABR
for nuser=1:Nu
    Q(:,1,nuser)=sort(sum(log2(1 + MSt_BMSN_BF(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
    Q(:,2,nuser)=sort(sum(log2(1 + MSt_BMSN_GE(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
    Q(:,3,nuser)=sort(sum(log2(1 + MSt_BMSN_GETD_3(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
    Q(:,4,nuser)=sort(sum(log2(1 + MSt_BD(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
    Q(:,5,nuser)=sort(sum(log2(1 + MSt_BDAS(:,(nuser-1)*Nrs+1:(nuser-1)*Nrs+Nrs)),2));
%     Q(:,6,nuser)=sort(sum(log2(1 + MSt_BMSN_GETD_3_CW(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
end
Qm = mean(Q,3);
%ABR(isnr,ia)=Qm(ICDF,ia);

Y(:,1) = (1/SIMU:1/SIMU:1).'*100;

figure;
plot(Qm,Y,'Linewidth',2);
set(gca,'Fontsize',14,'Fontname','Times New Roman');
legend('BMSN-BF',...
    'BMSN-GE',...
    'BMSN-GE-TD(3)',...
    'BD',...
    'BD-AS',...
    'Location','Southeast');
xlabel('Average Capacity [bits/s/Hz]','Fontsize',16,'Fontname','Times New Roman');
ylabel('CDF [%]','Fontsize',16,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
title(tsnr);
grid on;
hold on;