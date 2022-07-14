clear;%close;

% パラメータ
% パラメータ条件 NT >= NR*NU
lamr_min = 0;        % 第一固有値と第二固有値の比率の最小 [％]
lamr_max = 1;       % 第一固有値と第二固有値の比率の最大 [％]

SNR = 30;    % 入力SNR [dB]
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

sigma2 = 1/(10^(SNR/10)); % noise power
a = sigma2*Nt;

lamr=(lamr_min:0.01:lamr_max).';
Llamr=length(lamr);

for ilamr = 1:Llamr
    
    lamr_tar = lamr(ilamr);
    
    for isimu = 1:SIMU

    %H(伝搬チャネル行列: iid Rayleigh channel)5
    %Hu(ユーザkのチャネル行列Hu(:,:,k))
        H = (randn(Nr*Nu,Nt) + 1j*randn(Nr*Nu,Nt))/sqrt(2);
        Hu = peruser(H,Nu);

        He = zeros((Nu-1)*Nr,Nt,Nu);    % Hから1ユーザのチャネル行列を除いた行列
        
    % BMSN-GE7
        [~,~,STT,RIP,~] = bmsn_ge7(Nt,Nr,Nu,H,a,lamr_tar);
        for inu=1:Nu
            St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
        end
        %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
        MSt_BMSN_GE7(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);    
    
    end % isimu
    
    % ABR
    for nuser=1:Nu
        Q(:,1,nuser)=sort(sum(log2(1 + MSt_BMSN_GE7(:,(nuser-1)*Nr+1:(nuser-1)*Nr+2)),2));
    end
    
    Qm = mean(Q,3);

    QmC(ilamr,:)=Qm(round(CDF*SIMU/100),:);

    fprintf('lamr = %f \n',lamr(ilamr));
end

% グラフ表示
%% CDF of EGV at Target SNR 
figure;
plot(lamr,QmC,'Linewidth',1.5);
set(gca,'Fontsize',14,'Fontname','Arial');
xlabel('\lambda_2/\lambda_1','Fontsize',16,'Fontname','Arial');
ylabel('Average Capacity [bits/s/Hz]','Fontsize',16,'Fontname','Arial');
set(gca,'Fontsize',16,'Fontname','Arial');
title('CDF = 50%, SNR = 30dB');

grid on;
hold on;
