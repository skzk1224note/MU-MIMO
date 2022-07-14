% Achievable Bit Rate (ABR)

clear;
Nt = 16;    % 送信素子数
Nr = 2;     % 受信素子数 (アンテナ選択の場合1, そうでない場合は2)
Nu = 8;

Nrs = Nr-1; % Antenna Selection

SNR_tar = 20;  % ターゲットSNR[dB]

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

% 出力ファイル名 with alpha and SNR in dB
%evfile1 = strcat(folder,fn1,num2str(snr(isnr),'%02d'),'dB_1000itr.csv');
    
for isimu = 1:SIMU

    %H(伝搬チャネル行列: iid Rayleigh channel)
    %Hu(ユーザkのチャネル行列Hu(:,:,k))
    H = (randn(Nr*Nu,Nt) + 1j*randn(Nr*Nu,Nt))/sqrt(2);
    Hu = peruser(H,Nu);

%     for inu = 1:Nu
%         h = H(1+(inu-1)*Nr:2+(inu-1)*Nr,:);
%         S = [svd(h(1,:));svd(h(2,:))];  % nr =2
%         [~,No] = max(S);
%         Hu_as(:,:,inu) = h(No,:);
%     end
% 
%     %H_as(伝搬チャネル行列)
%     H_as = alluser(Hu_as); %伝搬チャネル行列

    %He = zeros((Nu-1)*Nr,Nt,Nu);    % Hから1ユーザのチャネル行列を除いた行列

    %% BMSN-BF1
    [~,~,STT,RIP,~] = bmsn_bf1(Nt,Nr,Nu,H,a,T);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布(SINR) for B-MSN（干渉波成分を考慮）
    MSt_BMSN1(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    clear St;
    
      %% BMSN-BF2
%     [~,~,STT,RIP,~] = bmsn_bf2(Nt,Nr,Nu,H,a);
%     for inu=1:Nu
%         St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
%     end
%     %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
%     MSt_BMSN2(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    
    %% BMSN-GE1
    [~,~,STT,RIP,~] = bmsn_ge1(Nt,Nr,Nu,H,a);
    %Nr1 = Nr-1;
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    MSt_BMSN2(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    clear St;
    
      %% BMSN-GE2
%     [~,~,STT,RIP,~] = bmsn_ge2(Nt,Nr,Nu,H,a);
%     %Nr1 = Nr-1;
%     for inu=1:Nu
%         St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
%     end
%     %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
%     MSt_BMSN4(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    
      %% BMSN-GE3
    [~,~,STT,RIP,~] = bmsn_getd_3(Nt,Nr,Nu,H,a);
%     %Nr1 = Nr-1;
     %i=SNR_tar;
    for inu=1:Nu
       St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
     
    end
   
%     %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）(SINR)
 MSt_BMSN5(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    
      %% MMSE-CI
%     [~,~,STT,RIP,~] = mmse(Nt,Nr,Nu,H,a);
%     %Nr1 = Nr-1;
%     for inu=1:Nu
%         St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
%     end
%     %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
%     MSt_MMSE(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    
      %% GMMSE-0
%     [~,~,STT,RIP,~] = gmmse(Nt,Nr,Nu,H,a);
%     for inu=1:Nu
%         St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
%     end
%     %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
%     MSt_GMMSE0(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    
      %% GMMSE-1
%     [~,~,STT,RIP,~] = gmmse_m1(Nt,Nr,Nu,H,a);
%     for inu=1:Nu
%         St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
%     end
%     %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
%     MSt_GMMSE1(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
  
      %% GMMSE-2
%     [~,~,STT,RIP,~] = gmmse_m2(Nt,Nr,Nu,H,a);
%     for inu=1:Nu
%         St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
%     end
%     %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
%     MSt_GMMSE2(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    
    %% BD
    [~,~,STT,RIP,~] = bd(Nt,Nr,Nu,H);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    MSt_BD(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    clear St;
    
    %% BD_AS
    [~,~,STT,RIP,~] = bd_as(Nt,Nr,Nrs,Nu,H);
    for inu=1:Nu
        St(:,inu) = STT(1,1,inu);
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    MSt_BDAS(isimu,:)=((St.').^2)./(RIP.'+Nt*sigma2);
    clear St;
    
    %% ZF
    [~,~,STT,RIP,~] = zf(Nt,Nr,Nu,H);
    for inu=1:Nu
        St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
    end
    %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
    MSt_ZF(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
    clear St;
        
    %fprintf('Iteration = %d / %d\n',isimu, SIMU);
    
end % isimu

% ABR
for nuser=1:Nu
%     Q(:,1,nuser)=sort(sum(log2(1 + MSt_MMSE(:,(nuser-1)*Nr+1:(nuser-1)*Nr+Nr)),2));
%     Q(:,2,nuser)=sort(sum(log2(1 + MSt_GMMSE0(:,(nuser-1)*Nr+1:(nuser-1)*Nr+Nr)),2));
%     Q(:,3,nuser)=sort(sum(log2(1 + MSt_GMMSE1(:,(nuser-1)*Nr+1:(nuser-1)*Nr+Nr)),2));
%     Q(:,4,nuser)=sort(sum(log2(1 + MSt_GMMSE2(:,(nuser-1)*Nr+1:(nuser-1)*Nr+Nr)),2));
    Q(:,1,nuser)=sort(sum(log2(1 + MSt_BMSN1(:,(nuser-1)*Nr+1:(nuser-1)*Nr+Nr)),2));
    Q(:,2,nuser)=sort(sum(log2(1 + MSt_BMSN2(:,(nuser-1)*Nr+1:(nuser-1)*Nr+Nr)),2));
%     Q(:,7,nuser)=sort(sum(log2(1 + MSt_BMSN3(:,(nuser-1)*Nr+1:(nuser-1)*Nr+Nr)),2));
%     Q(:,8,nuser)=sort(sum(log2(1 + MSt_BMSN4(:,(nuser-1)*Nr+1:(nuser-1)*Nr+Nr)),2));
     Q(:,9,nuser)=sort(sum(log2(1 + MSt_BMSN5(:,(nuser-1)*Nr+1:(nuser-1)*Nr+Nr)),2));
    Q(:,3,nuser)=sort(sum(log2(1 + MSt_BD(:,(nuser-1)*Nr+1:(nuser-1)*Nr+Nr)),2));
    Q(:,4,nuser)=sort(sum(log2(1 + MSt_BDAS(:,(nuser-1)*Nrs+1:(nuser-1)*Nrs+Nrs)),2));
    Q(:,5,nuser)=sort(sum(log2(1 + MSt_ZF(:,(nuser-1)*Nr+1:(nuser-1)*Nr+Nr)),2));
end
%三番目の次元に沿って、平均値を求める(1000*9*8のうち、1000*9のそれぞれのページの平均値)
%今はユーザー数が8なので8ページ分の平均値
Qm = mean(Q,3);
%ABR(isnr,ia)=Qm(ICDF,ia);

Y(:,1) = [0.1:0.1:100].';

figure;
mycol = [1 0 1;
      0 1 0;
      1 0 0;
      0 0 1;
      0 1 1;
      0 0 0]; %0.7 0.7 0;
set(groot,'defaultAxesColorOrder',mycol)
plot(Qm(:,1),Y,'r-o',Qm(:,2),Y,'b-s',Qm(:,3),Y,'m-^',Qm(:,4),Y,'c-v',Qm(:,5),Y,'g-v',Qm(:,9),Y,'y-v',...
    'MarkerIndices',10:100:length(Y),'Linewidth',2);
axis([0,15,0,100]);
%axis([0,20,0,100]);
% plot(Qm(:,5),Y,'b-+',...
%      Qm(:,6),Y,'g-o',...
%      Qm(:,7),Y,'r-^',...
%      Qm(:,8),Y,'m-v',...
%      Qm(:,9),Y,'c-s',...
%      Qm(:,10),Y,'k-x',...
%     'MarkerIndices',1:100:length(Y),'Linewidth',2);
set(gca,'Fontsize',14,'Fontname','Times New Roman');
legend('BMSN-BF','BMSN-GE','BD','BD-AS','ZF','GE-TD','Location','Southeast');
xlabel('Channel capacity [bits/s/Hz]','Fontsize',16,'Fontname','Times New Roman');
ylabel('CDF [%]','Fontsize',16,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
title(tsnr);
grid on;
hold on;

% figure;
% mycol = [0 0 0;
%       0 1 0;
%       0 0 1;
%       1 0 1;
%       1 0 0;
%       0.7 0.7 0];
% set(groot,'defaultAxesColorOrder',mycol)
% plot(Qm(:,1),Y,'k-s',Qm(:,2),Y,'g--^',Qm(:,3),Y,'b-v',Qm(:,4),Y,'m-x',Qm(:,7),Y,'r--^',...
%     'MarkerIndices',10:100:length(Y),'Linewidth',2);
% %axis([2,7,0,100]);
% axis([4,12,0,100]);
% %axis([6,18,0,100]);
% % plot(Qm(:,5),Y,'b-+',...
% %      Qm(:,6),Y,'g-o',...
% %      Qm(:,7),Y,'r-^',...
% %      Qm(:,8),Y,'m-v',...
% %      Qm(:,9),Y,'c-s',...
% %      Qm(:,10),Y,'k-x',...
% %     'MarkerIndices',1:100:length(Y),'Linewidth',2);
% set(gca,'Fontsize',14,'Fontname','Times New Roman');
% legend('MMSE-CI','GMI-0','GMI-1','GMI-2','BMSN-GE1','Location','Northwest');
% xlabel('Channel capacity [bits/s/Hz]','Fontsize',16,'Fontname','Times New Roman');
% ylabel('CDF [%]','Fontsize',16,'Fontname','Times New Roman');
% set(gca,'Fontsize',16,'Fontname','Times New Roman');
% title(tsnr);
% grid on;
% hold on;

%End


