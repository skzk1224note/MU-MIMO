clear;

% パラメータ条件 NT >= NR*NU
SNR_min = 0;        % 最小SNR [dB]
SNR_max = 30;       % 最大SNR [dB]
%SN_tar  = 25;        % CDF表示のためのターゲットSNR [dB]
%CDF = 50;   % TRのCDF

NT = 24;       % 送信素子数
NR =3;        % 受信素子数/ユーザ 
NRs = NR-1;        % 選択する受信アンテナ数
NU = 8;        % ユーザ数
SINU = 1000;   % 試行回数

Nru = NR*NU;    % 受信アンテナ総数
Nrus = NRs*NU;  % 選択する受信アンテナの総数

% 所望のチャネル行列 for BMSN-BF
T = zeros(NR,NR,NU);
for nuser = 1:NU
    T(:,:,nuser) = eye(NR,NR);
end

SNR=(SNR_min:5:SNR_max).';
LSNR=length(SNR);

% if CDF < 10
%     target_CDF=strcat('CDF= ',num2str(CDF,'%01d'),' %');
% else
%     target_CDF=strcat('CDF= ',num2str(CDF,'%02d'),' %');
% end

for isnr = 1:LSNR

    SN_tar = SNR(isnr);
    
    % 伝搬チャネル行列
    % Rayleigh
    H = (randn(SINU,NR*NU,NT)+1j*randn(SINU,NR*NU,NT))/sqrt(2);

    E_BD   = zeros(SINU, NR*NU);              % BDの固有値
    % E_BMSN1 = zeros(Ntri, NR*NU);             % BMSN1の固有値
    % E_BMSN2 = zeros(Ntri, NR*NU);             % BMSN2の固有値
    E_BMSNBF = zeros(SINU, NR*NU);             % BMSN-BFの固有値
    E_BMSNGE = zeros(SINU, NR*NU);              % BMSN-GEの固有値
    E_BMSNGEc = zeros(SINU, NR*NU);              % BMSN-GEcopyの固有値
%     E_ZF  = zeros(SINU, NR*NU);                % ZFの固有値
    E_BDAS  = zeros(SINU, NRs*NU);                % ZFの固有値

    for k = 1:SINU              % 試行回数のループ

        H0 = squeeze(H(k,:,:)); % k番目の試行回数での伝搬チャネル行列
  
%         % ZF algorithm
%         [W_ZF,U_ZF,S_ZF,RIPZF] = zf(NT,NR,NU,H0); % function, zf.m を使用
        a=NT/(10^(SN_tar/10));
        % BD algorithm
        [W_BD,U_BD,S_BD,RIPBD] = bmsn_getd_1_2_1(NT,NR,NU,H0,a); % function, bd.m を使用
    
        % BMSN-BF algorithm with adapted pseudo-noise power
        a=NT/(10^(SN_tar/10));
        [W_BMSNBF,U_BMSNBF,S_BMSNBF,RIPBF] = bmsn_getd_1_1_1(NT,NR,NU,H0,a); % function, bmsn_bf.m を使用
    
        % BMSN-GE algorithm with adapted pseudo-noise power
        a=NT/(10^(SN_tar/10));
        [W_BMSNGE,U_BMSNGE,S_BMSNGE,RIPGE] = bmsn_ge(NT,NR,NU,H0,a); % function, bmsn_gev.m を使用
        
        % BMSN-GEcopy algorithm with adapted pseudo-noise power
        a=NT/(10^(SN_tar/10));
        [W_BMSNGEc,U_BMSNGEc,S_BMSNGEc,RIPGEc] = bmsn_geatd_1(NT,NR,NU,H0,a); % function, bmsn_getd_3.m を使用
        a=NT/(10^(SN_tar/10));
        % BD-AS algorithm
        [W_BDAS,U_BDAS,S_BDAS,RIPBDAS] = bmsn_getd_1_2_2(NT,NR,NU,H0,a); % function, bd_as.m を使用
        
        % BMSN-GE-AS algorithm
%         [W_BMSNAS,U_BMSNAS,S_BMSNAS,RIPBMSNAS] = bmsn_ge_as(NT,NR,NRs,NU,H0,a); % function, bmsn_ge_as.m を使用
    
%         % MMSE algorithm
%         a=NT/(10^(SN_tar/10));
%         [W_MMSE,U_MMSE,S_MMSE,RIPm] = gmmse_m2(NT,NR,NU,H0,a); % function, bmsn.m を使用
    
        % ユーザ毎の固有値分布（RIPを考慮）
        snt = 1/(10^(SN_tar/10));
        for nuser=1:NU
%             E_ZF(k,1+(nuser-1)*NR:NR+(nuser-1)*NR) = 10*log10((diag(S_ZF(1:NR,1:NR,nuser)).^2)./(RIPZF(:,nuser)+NT*snt));    
            E_BD(k,1+(nuser-1)*NR:NR+(nuser-1)*NR) = 10*log10((diag(S_BD(1:NR,1:NR,nuser)).^2)./(RIPBD(:,nuser)+NT*snt));    
            E_BMSNBF(k,1+(nuser-1)*NR:NR+(nuser-1)*NR) = 10*log10((diag(S_BMSNBF(1:NR,1:NR,nuser)).^2)./(RIPBF(:,nuser)+NT*snt));
            E_BMSNGE(k,1+(nuser-1)*NR:NR+(nuser-1)*NR) = 10*log10((diag(S_BMSNGE(1:NR,1:NR,nuser)).^2)./(RIPGE(:,nuser)+NT*snt));
            E_BMSNGEc(k,1+(nuser-1)*NR:NR+(nuser-1)*NR) = 10*log10((diag(S_BMSNGEc(1:NR,1:NR,nuser)).^2)./(RIPGEc(:,nuser)+NT*snt));
            E_BDAS(k,1+(nuser-1)*NR:NR+(nuser-1)*NR) = 10*log10((diag(S_BDAS(1:NR,1:NR,nuser)).^2)./(RIPBDAS(:,nuser)+NT*snt));
%             E_BMSNAS(k,1+(nuser-1)*NRs:NRs+(nuser-1)*NRs) = 10*log10((diag(S_BMSNAS(1:NRs,1:NRs,nuser)).^2)./(RIPBMSNAS(:,nuser)+NT*snt));
        end
    
    end

    % 固有値の入力
    Eval_1 = E_BMSNBF;
    Eval_2 = E_BMSNGE;
%     Eval_3 = E_MMSE;
    Eval_3 = E_BMSNGEc;
    Eval_4 = E_BD;
%     Eval_5 = E_BMSNAS;
    Eval_5 = E_BDAS;

    % 伝送レート: 固有値からIEEE802.11ac伝送レートに変換
%     Out_1 = strcat('CSV_TR/TR16x2x8u_BMSNopt_SNR',num2str(SN_tar,'%02d'),'dB.csv');      % 出力ファイル名
%     Out_2 = strcat('CSV_TR/TR16x2x8u_BMSNGopt_SNR',num2str(SN_tar,'%02d'),'dB.csv');      % 出力ファイル名
%     Out_3 = strcat('CSV_TR/TR16x2x8u_MMSE_SNR',num2str(SN_tar,'%02d'),'dB.csv');      % 出力ファイル名
%     Out_4 = strcat('CSV_TR/TR16x2x8u_BD_SNR',num2str(SN_tar,'%02d'),'dB.csv');      % 出力ファイル名

    tr_table_40MHz; % 伝送レートとSNRの関係を入れる（20MHz or 40MHz）

    % 伝送レート計算
    % BMSN-BF
    TR_1 = Eval_to_TR (NR, NU, SINU,Eval_1,SNRT,TRT);
    %csvwrite(Out_1,TR_1);
    % BMSN-GE
    TR_2 = Eval_to_TR (NR, NU, SINU,Eval_2,SNRT,TRT);
    %csvwrite(Out_2,TR_2);
%     % MMSE
%     TR_3 = Eval_to_TR (NR, NU, Ntri,Eval_3,SNRT,TRT);
%     %csvwrite(Out_3,TR_3);
    % BD
    TR_3 = Eval_to_TR (NR, NU, SINU,Eval_3,SNRT,TRT);
    % BD-AS
    TR_4 = Eval_to_TR (NR, NU, SINU,Eval_4,SNRT,TRT);
    %csvwrite(Out_4,TR_4);
    % BMSN-AS
%     TR_5 = Eval_to_TR (NRs, NU, SINU,Eval_5,SNRT,TRT);
    % csvwrite(Out_5,TR_5);
    % BMSN-GEcopy
    TR_5 = Eval_to_TR (NR, NU, SINU,Eval_5,SNRT,TRT);
    % csvwrite(Out_5,TR_5);
    

    % 受信アンテナの和
    for nuser=1:NU
        Q(:,1,nuser) = sum(TR_1(:,(nuser-1)*NR+1:(nuser-1)*NR+NR),2);
        Q(:,2,nuser) = sum(TR_2(:,(nuser-1)*NR+1:(nuser-1)*NR+NR),2);
        Q(:,3,nuser) = sum(TR_3(:,(nuser-1)*NR+1:(nuser-1)*NR+NR),2);
        Q(:,4,nuser) = sum(TR_4(:,(nuser-1)*NR+1:(nuser-1)*NR+NR),2);
%         Q(:,5,nuser) = sum(TR_5(:,(nuser-1)*NRs+1:(nuser-1)*NRs+NRs),2);
        Q(:,5,nuser) = sum(TR_5(:,(nuser-1)*NR+1:(nuser-1)*NR+NR),2);
    end

    % ユーザ平均
    Qm(:,1) = mean(sort(Q(:,1,nuser),1),3);
    Qm(:,2) = mean(sort(Q(:,2,nuser),1),3);
    Qm(:,3) = mean(sort(Q(:,3,nuser),1),3);
    Qm(:,4) = mean(sort(Q(:,4,nuser),1),3);
%     Qm(:,5) = mean(sort(Q(:,5,nuser),1),3);
    Qm(:,5) = mean(sort(Q(:,5,nuser),1),3);

%    QmC(isnr,:)=Qm(round(CDF*10),:);
    % 試行平均
    QmC(isnr,:)=mean(Qm,1);

    fprintf('SNR = %d dB\n',SNR(isnr));
    
end

% グラフ表示
%% CDF of EGV at Target SNR 
figure;
mycol = [1 0 0;
      0 0 1;
      0 0.7 0;
      1 0 1;
      0.8 0.6 0;0 0 0];
set(groot,'defaultAxesColorOrder',mycol)
plot(SNR,QmC(:,1),'-o',SNR,QmC(:,2),'-v',SNR,QmC(:,3),'-^',SNR,QmC(:,4),'-s',SNR,QmC(:,5),'-x','Linewidth',2);
%plot(SNR,QmC,'Linewidth',2);
set(gca,'Fontsize',14,'Fontname','Times New Roman');
legend('GETD-1-1-1','BMSN-GE(1-2-3)','BMSN-GE-TD(1-1-3)','GETD-1-2-1','GETD-1-2-2','Location','Northwest');
xlabel('SNR [dB]','Fontsize',16,'Fontname','Times New Roman');
ylabel('Average transmission rate [Mbps]','Fontsize',16,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
%title(target_CDF);
grid on;
hold on;