clear;%close;

% パラメータ
% パラメータ条件 NT >= NR*NU
SNR_min = 0;        % 最小SNR [dB]
SNR_max = 30;       % 最大SNR [dB]

CDF = 10;   % ABRのCDF
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

if CDF < 10
    target_CDF=strcat('CDF= ',num2str(CDF,'%01d'),' %');
else
    target_CDF=strcat('CDF= ',num2str(CDF,'%02d'),' %');
end

SNR=(SNR_min:1:SNR_max).';
LSNR=length(SNR);

for isnr = 1:LSNR

    SNR_tar = SNR(isnr);

    sigma2 = 1/(10^(SNR_tar/10)); % noise power
    a = sigma2*Nt;
    
    for isimu = 1:SIMU

    %H(伝搬チャネル行列: iid Rayleigh channel)
    %Hu(ユーザkのチャネル行列Hu(:,:,k))
        H = (randn(Nr*Nu,Nt) + 1j*randn(Nr*Nu,Nt))/sqrt(2);
        Hu = peruser(H,Nu);

        He = zeros((Nu-1)*Nr,Nt,Nu);    % Hから1ユーザのチャネル行列を除いた行列

    %% BMSN-GE3EX_0.01
        [~,~,STT,RIP,~] = bmsn_ge3EX_0_01(Nt,Nr,Nu,H,a);
        for inu=1:Nu
            St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
        end
        %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
        MSt_BMSN_GE3EX_0_01(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
 
    %% BMSN-GE3EX_0.02
        [~,~,STT,RIP,~] = bmsn_ge3EX_0_02(Nt,Nr,Nu,H,a);
        for inu=1:Nu
            St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
        end
        %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
        MSt_BMSN_GE3EX_0_02(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
 
    %% BMSN-GE3EX_0.03
        [~,~,STT,RIP,~] = bmsn_ge3EX_0_03(Nt,Nr,Nu,H,a);
        for inu=1:Nu
            St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
        end
        %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
        MSt_BMSN_GE3EX_0_03(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
        
    %% BMSN-GE3EX_0.04
        [~,~,STT,RIP,~] = bmsn_ge3EX_0_04(Nt,Nr,Nu,H,a);
        for inu=1:Nu
            St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
        end
        %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
        MSt_BMSN_GE3EX_0_04(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
        
    %% BMSN-GE3EX_0.05
        [~,~,STT,RIP,~] = bmsn_ge3EX_0_05(Nt,Nr,Nu,H,a);
        for inu=1:Nu
            St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
        end
        %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
        MSt_BMSN_GE3EX_0_05(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
        
    %% BMSN-GE3EX_0.1
        [~,~,STT,RIP,~] = bmsn_ge3EX_0_1(Nt,Nr,Nu,H,a);
        for inu=1:Nu
            St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
        end
        %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
        MSt_BMSN_GE3EX_0_1(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);    
    
    %% BMSN-GE3EX_0.15
        [~,~,STT,RIP,~] = bmsn_ge3EX_0_15(Nt,Nr,Nu,H,a);
        for inu=1:Nu
            St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
        end
        %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
        MSt_BMSN_GE3EX_0_15(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
        
    %% BMSN-GE3EX_0.2
        [~,~,STT,RIP,~] = bmsn_ge3EX_0_2(Nt,Nr,Nu,H,a);
        for inu=1:Nu
            St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
        end
        %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
        MSt_BMSN_GE3EX_0_2(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
        
    %% BMSN-GE3EX_0.5
        [~,~,STT,RIP,~] = bmsn_ge3EX_0_5(Nt,Nr,Nu,H,a);
        for inu=1:Nu
            St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
        end
        %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
        MSt_BMSN_GE3EX_0_5(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
        
    %% BMSN-GE3EX_1
        [~,~,STT,RIP,~] = bmsn_ge3EX_1(Nt,Nr,Nu,H,a);
        for inu=1:Nu
            St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
        end
        %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
        MSt_BMSN_GE3EX_1(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
        
        
    %fprintf('Iteration = %d / %d\n',isimu, SIMU);
    
    end % isimu

    % ABR
    for nuser=1:Nu
        Q(:,1,nuser)=sort(sum(log2(1 + MSt_BMSN_GE3EX_0_01(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
        Q(:,2,nuser)=sort(sum(log2(1 + MSt_BMSN_GE3EX_0_02(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
        Q(:,3,nuser)=sort(sum(log2(1 + MSt_BMSN_GE3EX_0_03(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
        Q(:,4,nuser)=sort(sum(log2(1 + MSt_BMSN_GE3EX_0_04(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
        Q(:,5,nuser)=sort(sum(log2(1 + MSt_BMSN_GE3EX_0_05(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
        Q(:,6,nuser)=sort(sum(log2(1 + MSt_BMSN_GE3EX_0_1(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
        Q(:,7,nuser)=sort(sum(log2(1 + MSt_BMSN_GE3EX_0_15(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
        Q(:,8,nuser)=sort(sum(log2(1 + MSt_BMSN_GE3EX_0_2(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
        Q(:,9,nuser)=sort(sum(log2(1 + MSt_BMSN_GE3EX_0_5(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
        Q(:,10,nuser)=sort(sum(log2(1 + MSt_BMSN_GE3EX_1(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));

    end

    Qm = mean(Q,3);

    QmC(isnr,:)=Qm(round(CDF*10),:);

    fprintf('SNR = %d dB\n',SNR(isnr));
end

%csvwrite(cdffile1,[SNR,QmC]);

% グラフ表示
%% CDF of EGV at Target SNR 
figure;
plot(SNR,QmC,'Linewidth',1.5);
% axis([0,30,0,12])
set(gca,'Fontsize',14,'Fontname','Arial');
legend('\lambda_2/\lambda_1 < 0.01','\lambda_2/\lambda_1 < 0.02','\lambda_2/\lambda_1 < 0.03','\lambda_2/\lambda_1 < 0.04',...
    '\lambda_2/\lambda_1 < 0.05','\lambda_2/\lambda_1 < 0.1','\lambda_2/\lambda_1 < 0.15','\lambda_2/\lambda_1 < 0.2',...
    '\lambda_2/\lambda_1 < 0.5','\lambda_2/\lambda_1 < 1','Location','Northwest');
xlabel('SNR [dB]','Fontsize',16,'Fontname','Arial');
ylabel('Average Capacity [bits/s/Hz]','Fontsize',16,'Fontname','Arial');
set(gca,'Fontsize',16,'Fontname','Arial');
title(target_CDF);
grid on;
hold on;


% %plot(SNR,QmC(:,1),SNR,QmC(:,2),SNR,QmC(:,3),SNR,QmC(:,4),'Linewidth',2);
% plot(SNR,QmC(:,10),'c-d',SNR,QmC(:,1),'k-s',SNR,QmC(:,2),'g-o',SNR,QmC(:,3),'b--x',SNR,QmC(:,4),'m-v',SNR,QmC(:,7),'r--^',...
%     'Linewidth',2);
