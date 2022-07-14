clear;%close;

% パラメータ
% パラメータ条件 NT >= NR*NU
SNR_min = 0;        % 最小SNR [dB]
SNR_max = 30;       % 最大SNR [dB]

CDF = 50;   % ABRのCDF
Nt = 16;    % 送信素子数
Nr = 2;     % 受信素子数 (アンテナ選択の場合1, そうでない場合は2)
Nu = 8;
SIMU = 1000; %試行回数（通常 1000）

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

% % 入力ファイル基幹名
% folder='CSV3/';
% fn1 = 'Eig16x2x8u_BD_SNR';
% fn2 = 'Eig16x2x8u_BMSN1(GEV)_SNR';
% fn3 = 'Eig16x2x8u_BMSN2a_SNR';  % alpha = 1e-2
% fn4 = 'Eig16x2x8u_BMSN2b_SNR';  % alpha = 1e-6
% fn5 = 'Eig16x2x8u_BMSN3_SNR';   % alpha = sigma2Nt equal to MMSE
% fn6 = 'Eig16x2x8u_ZF_SNR';
% fn7 = 'Eig16x2x8u_MMSE_SNR';

% 入力ファイル基幹名
% folder='CSV4/';
% fn1 = 'Eig16x2x8u_ZF_SNR';
% fn2 = 'Eig16x2x8u_BD_SNR';
% %fn3 = 'Eig16x2x8u_BMSN2_SNR'; % 2: alpha = 1e-2, 3: alpha=sigma2xNT
% %fn4 = 'Eig16x2x8u_BDAS_SNR';  % antenna selection
% %fn5 = 'Eig16x2x8u_BMSNAS2_SNR';  % alpha = 1e-2 and antenna selection
% fn3 = 'Eig16x2x8u_BMSN3_SNR';   % alpha = sigma2Nt equal to MMSE
% fn4 = 'Eig16x2x8u_BMSN3(GEV)_SNR';   % alpha = sigma2Nt equal to MMSE
% fn5 = 'Eig16x2x8u_MMSE_SNR';

% 出力ファイル名
folder= 'CSV/';
cdfn1 = 'ABRCDFvsSNR_16x2x8u_BD_MMSE_BMSNs_USERave_CDF';
%cdfn2 = 'Eig_CDF_16x2x8u_BMSN_SNR';

if CDF < 10
    target_CDF=strcat('CDF= ',num2str(CDF,'%01d'),' %');
else
    target_CDF=strcat('CDF= ',num2str(CDF,'%02d'),' %');
end

% 出力ファイル名 with SNR in dB
cdffile1 = strcat(folder,cdfn1,num2str(CDF,'%02d'),'_1000itr.csv');
%cdffile2 = strcat(folder,cdfn2,num2str(SN_tar,'%02d'),'dB_1000itr.csv');

SNR=(SNR_min:5:SNR_max).';
LSNR=length(SNR);

for isnr = 1:LSNR

    SNR_tar = SNR(isnr);

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

    % BMSN-GE
        [~,~,STT,RIP,~] = bmsn_gev(Nt,Nr,Nu,H,a);
        for inu=1:Nu
            St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
        end
        %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
        MSt_BMSN_GE(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
 
    % BMSN-GE4
        [~,~,STT,RIP,~] = bmsn_ge4(Nt,Nr,Nu,H,a);
        for inu=1:Nu
            St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
        end
        %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
        MSt_BMSN_GE4(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
 
    % BMSN-GE5
        [~,~,STT,RIP,~] = bmsn_ge5(Nt,Nr,Nu,H,a);
        for inu=1:Nu
            St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
        end
        %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
        MSt_BMSN_GE5(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
        
    % BMSN-GE6
        [~,~,STT,RIP,~] = bmsn_ge6(Nt,Nr,Nu,H,a);
        for inu=1:Nu
            St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
        end
        %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
        MSt_BMSN_GE6(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
        
    % BMSN-GE3
        [~,~,STT,RIP,~] = bmsn_ge3(Nt,Nr,Nu,H,a);
        for inu=1:Nu
            St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
        end
        %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
        MSt_BMSN_GE3(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
        
    % BMSN-GE3EX
        [~,~,STT,RIP,~] = bmsn_ge3EX(Nt,Nr,Nu,H,a);
        for inu=1:Nu
            St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
        end
        %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
        MSt_BMSN_GE3EX(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);    
    
    % BD
        [~,~,STT,RIP,~] = bd(Nt,Nr,Nu,H);
        for inu=1:Nu
            St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
        end
        %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
        MSt_BD(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
        
    % ZF-CI
%         [~,~,STT,RIP,~] = zf(Nt,Nr,Nu,H);
%         for inu=1:Nu
%             St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
%         end
        %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
%         MSt_ZF(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
        
        
    %fprintf('Iteration = %d / %d\n',isimu, SIMU);
    
    end % isimu

    % ABR
    for nuser=1:Nu
        Q(:,1,nuser)=sort(sum(log2(1 + MSt_BMSN_GE(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
        Q(:,2,nuser)=sort(sum(log2(1 + MSt_BMSN_GE4(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
        Q(:,3,nuser)=sort(sum(log2(1 + MSt_BMSN_GE5(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
        Q(:,4,nuser)=sort(sum(log2(1 + MSt_BMSN_GE6(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
        Q(:,5,nuser)=sort(sum(log2(1 + MSt_BMSN_GE3(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
        Q(:,6,nuser)=sort(sum(log2(1 + MSt_BMSN_GE3EX(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
        Q(:,7,nuser)=sort(sum(log2(1 + MSt_BD(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
%         Q(:,4,nuser)=sort(sum(log2(1 + MSt_BMSN1(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
%         Q(:,5,nuser)=sort(sum(log2(1 + MSt_BMSN2(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));

    end

    Qm = mean(Q,3);

    QmC(isnr,:)=Qm(round(CDF*10),:);

    fprintf('SNR = %d dB\n',SNR(isnr));
end

%csvwrite(cdffile1,[SNR,QmC]);

% グラフ表示
%% CDF of EGV at Target SNR 
figure;
mycol = [1 0 0;
      0 0 1;
      0 1 0;
      1 0 1;
      1 0.65 0;
      0 0.5 0.5;
      0.5 0 0.5]; %0.7 0.7 0;
set(groot,'defaultAxesColorOrder',mycol)
plot(SNR,QmC,'Linewidth',1.5);
axis([5,30,0,12])
set(gca,'Fontsize',14,'Fontname','Arial');
legend('BMSN-GE','BMSN-GE (log2(1+λ))','BMSN-GE (√λ)','BMSN-GE (log2(1+λ))^2','BMSN-GEcopy','BMSN-GEcopy \lambda_2<0.1','BD','Location','Northwest');
xlabel('SNR [dB]','Fontsize',16,'Fontname','Arial');
ylabel('Average Capacity [bits/s/Hz]','Fontsize',16,'Fontname','Arial');
set(gca,'Fontsize',16,'Fontname','Arial');
title(target_CDF);
grid on;
hold on;

% figure;
% mycol = [0 1 1;
%       0 1 0;
%       1 0 1;
%       0 0 1;
%       1 0 0;
%       0 0 0]; %0.7 0.7 0;
% set(groot,'defaultAxesColorOrder',mycol)
% %plot(SNR,QmC(:,1),SNR,QmC(:,2),SNR,QmC(:,3),SNR,QmC(:,4),'Linewidth',2);
% plot(SNR,QmC(:,10),'c-d',SNR,QmC(:,1),'k-s',SNR,QmC(:,2),'g-o',SNR,QmC(:,3),'b--x',SNR,QmC(:,4),'m-v',SNR,QmC(:,7),'r--^',...
%     'Linewidth',2);
% axis([5,30,0,12])
% set(gca,'Fontsize',14,'Fontname','Times New Roman');
% legend('BD','MMSE-CI','GMI-0','GMI-1','GMI-2','BMSN-GE1','Location','Northwest');
% xlabel('SNR [dB]','Fontsize',16,'Fontname','Times New Roman');
% ylabel('Achievable bit rate [bits/s/Hz]','Fontsize',16,'Fontname','Times New Roman');
% set(gca,'Fontsize',16,'Fontname','Times New Roman');
% title(target_CDF);
% grid on;
% hold on;