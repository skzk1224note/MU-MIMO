clear;

% パラメータ条件 NT >= NR*NU
SN_tar  = 10;        % CDF表示のためのターゲットSNR [dB]

NT = 16;       % 送信素子数
NR = 2;        % 受信素子数/ユーザ 
%NA = 1;        % 選択する受信アンテナ数
NU = 8;        % ユーザ数
Ntri = 1000;   % 試行回数

Nru = NR*NU;

% 所望のチャネル行列 for BMSN
T = zeros(NR,NR,NU);
for nuser = 1:NU
    T(:,:,nuser) = eye(NR,NR);
end

if SN_tar < 10
    target_snr=strcat('SNR= ',num2str(SN_tar,'%01d'),'dB');
else
    target_snr=strcat('SNR=',num2str(SN_tar,'%02d'),'dB');
end

% 伝搬チャネル行列
% Rayleigh
H = (randn(Ntri,NR*NU,NT)+1j*randn(Ntri,NR*NU,NT))/sqrt(2);

E_BD   = zeros(Ntri, NR*NU);              % BDの固有値
% E_BMSN1 = zeros(Ntri, NR*NU);             % BMSN1の固有値
% E_BMSN2 = zeros(Ntri, NR*NU);             % BMSN2の固有値
E_BMSNBF = zeros(Ntri, NR*NU);             % BMSN3の固有値
E_BMSNGE = zeros(Ntri, NR*NU);             % BMSN3の固有値
% E_MMSE  = zeros(Ntri, NR*NU);              % MMSEの固有値

for k = 1:Ntri              % 試行回数のループ

    H0 = squeeze(H(k,:,:)); % k番目の試行回数での伝搬チャネル行列
    
    % BD algorithm
    [W_BD,U_BD,S_BD,RIPBD] = bd(NT,NR,NU,H0); % function, bd.m を使用
    
    % BMSN-BF algorithm with adapted pseudo-noise power
    a=NT/(10^(SN_tar/10));
    [W_BMSNBF,U_BMSNBF,S_BMSNBF,RIPBF] = bmsn_bf(NT,NR,NU,H0,a,T); % function, bmsn.m を使用
    
    % BMSN-GE algorithm with adapted pseudo-noise power
    a=NT/(10^(SN_tar/10));
    [W_BMSNGE,U_BMSNGE,S_BMSNGE,RIPGE] = bmsn_ge(NT,NR,NU,H0,a); % function, bmsn.m を使用
    
%     % MMSE algorithm
%     a=NT/(10^(SN_tar/10));
%     [W_MMSE,U_MMSE,S_MMSE,RIPm] = gmmse_m2(NT,NR,NU,H0,a); % function, bmsn.m を使用
    
    % ユーザ毎の固有値分布（RIPを考慮）
    snt = 1/(10^(SN_tar/10));
    for nuser=1:NU
        E_BD(k,1+(nuser-1)*NR:NR+(nuser-1)*NR) = 10*log10((diag(S_BD(:,:,nuser)).^2)./(RIPBD(:,nuser)+NT*snt));   
        E_BMSNBF(k,1+(nuser-1)*NR:NR+(nuser-1)*NR) = 10*log10((diag(S_BMSNBF(:,:,nuser)).^2)./(RIPBF(:,nuser)+NT*snt));
        E_BMSNGE(k,1+(nuser-1)*NR:NR+(nuser-1)*NR) = 10*log10((diag(S_BMSNGE(:,:,nuser)).^2)./(RIPGE(:,nuser)+NT*snt));
%         E_MMSE(k,1+(nuser-1)*NR:NR+(nuser-1)*NR) = 10*log10((diag(S_MMSE(:,:,nuser)).^2)./(RIPm(:,nuser)+NT*snt));
    end
    
end

% 固有値の入力
Eval_1 = E_BMSNBF;
Eval_2 = E_BMSNGE;
%Eval_3 = E_MMSE;
Eval_3 = E_BD;

% 伝送レート: 固有値からIEEE802.11ac伝送レートに変換
% Out_1 = strcat('CSV_TR/TR16x2x8u_BMSNopt_SNR',num2str(SN_tar,'%02d'),'dB.csv');      % 出力ファイル名
% Out_2 = strcat('CSV_TR/TR16x2x8u_BMSNGopt_SNR',num2str(SN_tar,'%02d'),'dB.csv');      % 出力ファイル名
% Out_3 = strcat('CSV_TR/TR16x2x8u_MMSE_SNR',num2str(SN_tar,'%02d'),'dB.csv');      % 出力ファイル名
% Out_4 = strcat('CSV_TR/TR16x2x8u_BD_SNR',num2str(SN_tar,'%02d'),'dB.csv');      % 出力ファイル名

tr_table_40MHz; % 伝送レートとSNRの関係を入れる

% 伝送レート計算
% BMSN-BF
TR_1 = Eval_to_TR (NR, NU, Ntri,Eval_1,SNRT,TRT);
%csvwrite(Out_1,TR_1);
% BMSN-GE
TR_2 = Eval_to_TR (NR, NU, Ntri,Eval_2,SNRT,TRT);
%csvwrite(Out_2,TR_2);
% MMSE
%TR_3 = Eval_to_TR (NR, NU, Ntri,Eval_3,SNRT,TRT);
%csvwrite(Out_3,TR_3);
% BD
TR_3 = Eval_to_TR (NR, NU, Ntri,Eval_3,SNRT,TRT);
%csvwrite(Out_4,TR_4);
% BD-AS
% TR_5 = Eval_to_TR (NA, NU, Ntri,Eval_5,SNRT,TRT);
% csvwrite(Out_5,TR_5);
% BMSN-AS
% TR_6 = Eval_to_TR (NA, NU, Ntri,Eval_6,SNRT,TRT);
% csvwrite(Out_6,TR_6);

%nuser=1;
% 受信アンテナの和
for nuser=1:NU
    Q(:,1,nuser) = sum(TR_1(:,(nuser-1)*NR+1:(nuser-1)*NR+NR),2);
    Q(:,2,nuser) = sum(TR_2(:,(nuser-1)*NR+1:(nuser-1)*NR+NR),2);
    Q(:,3,nuser) = sum(TR_3(:,(nuser-1)*NR+1:(nuser-1)*NR+NR),2);
 %   Q(:,4,nuser) = sum(TR_4(:,(nuser-1)*NR+1:(nuser-1)*NR+NR),2);
end

% ユーザ平均
Qm(:,1) = mean(sort(Q(:,1,nuser),1),3);
Qm(:,2) = mean(sort(Q(:,2,nuser),1),3);
Qm(:,3) = mean(sort(Q(:,3,nuser),1),3);
%Qm(:,4) = mean(sort(Q(:,4,nuser),1),3);

% CDF of TR at Target SNR
Y(:,1) = (1/Ntri:1/Ntri:1).'*100; 

% グラフ表示
%% CDF of EGV at Target SNR 
figure;
mycol = [1 0 0;
      0 0 1;
      0 0.7 0;
      1 0 1;
      0.8 0.6 0;0 0 0];
set(groot,'defaultAxesColorOrder',mycol)
%plot(Q(:,1),Y,'-m',Q(:,2),Y,'-r',Q(:,3),Y,'-b',Q(:,4),Y,'-k',Q(:,5),Y,'-g','Linewidth',2);
plot(Qm(:,1),Y,'-',Qm(:,2),Y,'-',Qm(:,3),Y,'-','Linewidth',2);
set(gca,'Fontsize',14,'Fontname','Times New Roman');
legend('BMSN-BF','BMSN-GE','BD','Location','Southeast');
xlabel('Average transmission rate [Mbps]','Fontsize',16,'Fontname','Times New Roman');
ylabel('CDF [%]','Fontsize',16,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
title(target_snr);
grid on;
hold on;





