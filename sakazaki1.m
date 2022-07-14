clear;

rng('default');

% 出力ファイル名
evfile1 = 'BER_CSV5/BER16x2x8u_BD_SNR0-30dB_10000bits_5000itr.csv';
evfile2 = 'BER_CSV5/BER16x2x8u_BMSN2_SNR0-30dB_10000bits_5000itr.csv';
evfile3 = 'BER_CSV5/BER16x2x8u_BMSN2GE_SNR0-30dB_10000bits_5000itr.csv';
evfile4 = 'BER_CSV5/BER16x2x8u_BMSN2GEc_SNR0-30dB_10000bits_5000itr.csv';
evfile5 = 'BER_CSV5/BER16x2x8u_BDAS_SNR0-30dB_10000bits_5000itr.csv';

Nt = 16; %基地局アンテナ数
Nr = 2; %ユーザアンテナ数は2(=送信ストリーム数(not1))
Nu = 8; %ユーザ数
ndata = 10000; % 1パケット当たりのシンボル数（パケット長）:BER計算のための1試行当たりのシンボル数 10000
SNR_min = 0;  %最小SNR[dB]
SNR_max = 30; %最大SNR[dB]
SIMU = 1000; %試行回数（通常 1000,5000）
pattern = [4,0;2,2;3,1]; %変調パターン(先頭に[●,0]をもってくる）
bsr = pattern(1,1);
% a = 1e-2; % BMSNの擬似雑音
I = eye(Nt,Nt);

%T 所望のチャネル行列
for inu = 1:Nu
    T(:,:,inu) = zeros(Nr,Nr);
    T(1,1,inu) = 1;
    T(1,2,inu) = 0;
    T(2,1,inu) = 0;
    T(2,2,inu) = 1;
end

tic;

snr = SNR_min:SNR_max;  % snr=[SNR_min,....,SNR_max]の横ベクトル
Lsnr = length(snr);     % snr横ベクトルの大きさ
    
for isnr = 1:Lsnr
    sigma2 = 1/(10^(snr(isnr)/10)); % noise power
    a = sigma2*Nt;
    %a2 = sigma2*Nt;
    
for isimu = 1:SIMU

%H(伝搬チャネル行列: iid Rayleigh channel)
%Hu(ユーザkのチャネル行列Hu(:,:,k))

H = (randn(Nr*Nu,Nt) + 1j*randn(Nr*Nu,Nt))/sqrt(2);
Hu = peruser(H,Nu);
HH = H'*H;

% ZFウエイト
%Wzf = H\I;
%Wzf = H'/(H*H');

% MMSEウエイト
%Wmmse = H'/(H*H'+ a*eye(Nr*Nu)); % MMSE

for inu = 1:Nu
    %h = (randn(2,nt) + 1j*randn(2,nt))/sqrt(2);
    h = H(1+(inu-1)*Nr:2+(inu-1)*Nr,:);
    S = [svd(h(1,:));svd(h(2,:))];  % nr =2
    [~,No] = max(S);
    Hu_as(:,:,inu) = h(No,:);
end

%H(伝搬チャネル行列)
H_as = alluser(Hu_as); %伝搬チャネル行列

He = zeros(Nt-Nr,Nt,Nu);

for inu = 1:Nu
    
    %He(Hからユーザkのチャネル行列を除いた行列)
    el = 1:Nu;
    el(:,inu) = [];
    He(:,:,inu) = alluser(Hu(:,:,el));
    HeHe = He(:,:,inu)'*He(:,:,inu);
    He_as(:,:,inu) = alluser(Hu_as(:,:,el));

    %% BD
    %Ven(Heの雑音部分空間に対応する固有ベクトル)
    [~,~,Ve] = svd(He(:,:,inu));
    Ven(:,:,inu) = Ve(:,Nr*(Nu-1)+1:Nt);

    %Vts(Hu*Venの信号部分空間に対応する固有ベクトル)
    %Wtu(ユーザkに対応する基地局側ウエイト)
    %St(特異値を格納 2×nu)
    %Wr(ユーザ側ウエイト)
    [Ut,Std,Vt] = svd(Hu(:,:,inu)*Ven(:,:,inu));
    St_BD(:,inu) = diag(Std(1:Nr,1:Nr));
  
    Wr_BD(:,:,inu) = Ut';
    Vts_BD(:,:,inu) = Vt(:,1:Nr);
    Wtu_BD(:,:,inu) = Ven(:,:,inu)*Vts_BD(:,:,inu); 

    %% BMSN2-BF (bmsn_bf.m)
    %MSNに基づくユーザ毎のウエイトの計算
    % MSNの最適ウエイトWopt
    Wopt = (HeHe+a*eye(size(He,2)))\Hu(:,:,inu)'*T(:,:,inu);
    Wopt = Wopt/norm(Wopt,'fro')*sqrt(Nr);
%     for ij = 1:Nr
%         Wopt(:,ij,inu) = Wopt(:,ij,inu)/sqrt(Wopt(:,ij,inu)'*Wopt(:,ij,inu));
%     end
    %Vts(Hu*Woptの信号部分空間に対応する固有ベクトル)
    %Wtu(ユーザkに対応する基地局側ウエイト)
    %St(特異値を格納 2×nu)
    %Wr(ユーザ側ウエイト)
    [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt);
    St_BMSN2(:,inu) = diag(Std(1:Nr,1:Nr));
    
    Wr_BMSN2(:,:,inu) = Ut';
    Vts_BMSN2(:,:,inu) = Vt(:,1:Nr);
    Wtu_BMSN2(:,:,inu) = Wopt*Vts_BMSN2(:,:,inu);
    
     %% BMSN2-GE α=sigma2 NT equal to MMSE
    A = HeHe+a*eye(size(He(:,:,inu),2));
    B = Hu(:,:,inu)'*Hu(:,:,inu);
    [EW,D] = eig(B,A);
    [~,IN] = sort(diag(abs(D)).','descend');
    EW = EW(:,IN);
    EWS=EW(:,1:Nr);
%     for ij = 1:Nr
%         EWS(:,ij) = EWS(:,ij)/sqrt(EWS(:,ij)'*EWS(:,ij));
%     end          
    Wopt=EWS;
%     Wopt=EWS(:,1:Nr)*sqrt(D2);
    Wopt = Wopt/norm(Wopt,'fro')*sqrt(Nr);
%     Wopt_k2(:,:,inu) = Wopt_k2(:,:,inu)/norm(Wopt_k2(:,:,inu),'fro')*sqrt(Nr);
    %Vts(HT*Wの信号部分空間に対応する固有ベクトル)
    %Wtu(ユーザkに対応する基地局側ウエイト)
    %St(特異値を格納 2×nu)
    %Wr(ユーザ側ウエイト)
    [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt);
    St_BMSN2g(:,inu) = diag(Std(1:Nr,1:Nr));
    Wr_BMSN2g(:,:,inu) = Ut';
    Vts_BMSN2g(:,:,inu) = Vt(:,1:Nr);
    Wtu_BMSN2g(:,:,inu) = Wopt*Vts_BMSN2g(:,:,inu);
    
    %% BMSN2-GE-TD α=sigma2 NT equal to MMSE
    A = HeHe+a*eye(size(He(:,:,inu),2));
    B = Hu(:,:,inu)'*Hu(:,:,inu);
    [EW,D] = eig(B,A);
    [D1,IN] = sort(diag(abs(D)).','descend');
    EW = EW(:,IN);
    if Nr > 1
            EW(:,Nr) = EW(:,Nr-1);
    end
    EWS=EW(:,1:Nr);
%     for ij = 1:Nr
%         EWS(:,ij) = EWS(:,ij)/sqrt(EWS(:,ij)'*EWS(:,ij));
%     end          
    Wopt=EWS;
%     Wopt=EWS(:,1:Nr)*sqrt(D2);
    Wopt = Wopt/norm(Wopt,'fro')*sqrt(Nr);
%     Wopt_k2(:,:,inu) = Wopt_k2(:,:,inu)/norm(Wopt_k2(:,:,inu),'fro')*sqrt(Nr);
    %Vts(HT*Wの信号部分空間に対応する固有ベクトル)
    %Wtu(ユーザkに対応する基地局側ウエイト)
    %St(特異値を格納 2×nu)
    %Wr(ユーザ側ウエイト)
    [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt);
    St_BMSN2gc(:,inu) = diag(Std(1:Nr,1:Nr));
    Wr_BMSN2gc(:,:,inu) = Ut';
    Vts_BMSN2gc(:,:,inu) = Vt(:,1:Nr);
    Wtu_BMSN2gc(:,:,inu) = Wopt*Vts_BMSN2gc(:,:,inu);
    
    %% BD with antenna selection
    %Ven(Heの雑音部分空間に対応する固有ベクトル)
    [~,~,Ve] = svd(He_as(:,:,inu));
    Ven2(:,:,inu) = Ve(:,Nu:Nt);

    %Vts(Hu*Venの信号部分空間に対応する固有ベクトル)
    %Wtu(ユーザkに対応する基地局側ウエイト)
    %St(特異値を格納 1×nu)
    %Wr(ユーザ側ウエイト)
    [Ut,Std,Vt] = svd(Hu_as(:,:,inu)*Ven2(:,:,inu));
    St_BD2(:,inu) = Std(1,1);
    
    Wr_BD2(:,:,inu) = Ut';
    Vts_BD2(:,:,inu) = Vt(:,1);
    Wtu_BD2(:,:,inu) = Ven2(:,:,inu)*Vts_BD2(:,:,inu);

end

%Wt(基地局側ウエイト)
%Wt = reshape(Wtu,[nt,nt]);

%amp(送信電力)
amp2 = 10^(snr(isnr)/10)/Nt;

% 以降 nr=2 の場合のみ
%data(ユーザkへの送信データ)
for inu = 1:Nu
for ii = 1:size(pattern,1)
    data1 = randi([0 2^pattern(ii,1)-1],1,ndata);
    data2 = randi([0 2^pattern(ii,2)-1],1,ndata);
    data(:,:,ii,inu) = [data1;data2];
end
end

%s(ユーザkへの送信ストリーム)
for inu = 1:Nu
for ii = 1:size(pattern,1)
    s1 = Mapping(data(1,:,ii,inu),pattern(ii,1));
    s2 = Mapping(data(2,:,ii,inu),pattern(ii,2));
    s(:,:,ii,inu) = [s1;s2];
end
end

% %data(送信データ) for antenna selection (AS)
data_as = randi([0 2^bsr-1],Nu,ndata);

% %s(送信ストリーム) for antenna selection (AS)
s_as = Mapping(data_as,bsr);


%n(熱雑音)
%y(ユーザkの入力信号)
%Y(ユーザkの復調信号)
for inu = 1:Nu
    n = (randn(Nr,ndata) + 1j*randn(Nr,ndata))/sqrt(2);
    
    %  for BD with antenna selection (BD-AS, BMSN-AS)
    Hs_BD2 = zeros(1,ndata);    % ユーザiの受信信号（受信重み付け前）
    for ij = 1:Nu
        Hs_BD2 = Hs_BD2 + Hu_as(:,:,inu)*Wtu_BD2(:,:,ij)*s_as(ij,:)*sqrt(amp2);
    end
    yu_BD2(:,:,inu) = Wr_BD2(:,:,inu)*(Hs_BD2+n(1,:))/sqrt(amp2)/St_BD2(1,inu);
    
    for ii = 1:size(pattern,1)
        Hs_BD = zeros(Nr,ndata);    % ユーザiの受信信号（受信重み付け前）
        Hs_BMSN2 = zeros(Nr,ndata);    % ユーザiの受信信号（受信重み付け前）
        Hs_BMSN2g = zeros(Nr,ndata);    % ユーザiの受信信号（受信重み付け前）
        Hs_BMSN2gc = zeros(Nr,ndata);    % ユーザiの受信信号（受信重み付け前）
        for ij = 1:Nu
            Hs_BD = Hs_BD + Hu(:,:,inu)*Wtu_BD(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
            Hs_BMSN2 = Hs_BMSN2 + Hu(:,:,inu)*Wtu_BMSN2(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
            Hs_BMSN2g = Hs_BMSN2g + Hu(:,:,inu)*Wtu_BMSN2g(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
            Hs_BMSN2gc = Hs_BMSN2gc + Hu(:,:,inu)*Wtu_BMSN2gc(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
        end
        yu_BD(:,:,ii,inu) = Wr_BD(:,:,inu)*(Hs_BD+n);
        yu_BMSN2(:,:,ii,inu) = Wr_BMSN2(:,:,inu)*(Hs_BMSN2+n);
        yu_BMSN2g(:,:,ii,inu) = Wr_BMSN2g(:,:,inu)*(Hs_BMSN2g+n);
        yu_BMSN2gc(:,:,ii,inu) = Wr_BMSN2gc(:,:,inu)*(Hs_BMSN2gc+n);
        
        yu_BD(1,:,ii,inu) = yu_BD(1,:,ii,inu)/St_BD(1,inu)/sqrt(amp2);
        yu_BD(2,:,ii,inu) = yu_BD(2,:,ii,inu)/St_BD(2,inu)/sqrt(amp2);
        Yu_BD(1,:,ii,inu) = Decode(yu_BD(1,:,ii,inu),pattern(ii,1));
        Yu_BD(2,:,ii,inu) = Decode(yu_BD(2,:,ii,inu),pattern(ii,2));
        
        yu_BMSN2(1,:,ii,inu) = yu_BMSN2(1,:,ii,inu)/St_BMSN2(1,inu)/sqrt(amp2);
        yu_BMSN2(2,:,ii,inu) = yu_BMSN2(2,:,ii,inu)/St_BMSN2(2,inu)/sqrt(amp2);
        Yu_BMSN2(1,:,ii,inu) = Decode(yu_BMSN2(1,:,ii,inu),pattern(ii,1));
        Yu_BMSN2(2,:,ii,inu) = Decode(yu_BMSN2(2,:,ii,inu),pattern(ii,2));
        
        yu_BMSN2g(1,:,ii,inu) = yu_BMSN2g(1,:,ii,inu)/St_BMSN2g(1,inu)/sqrt(amp2);
        yu_BMSN2g(2,:,ii,inu) = yu_BMSN2g(2,:,ii,inu)/St_BMSN2g(2,inu)/sqrt(amp2);
        Yu_BMSN2g(1,:,ii,inu) = Decode(yu_BMSN2g(1,:,ii,inu),pattern(ii,1));
        Yu_BMSN2g(2,:,ii,inu) = Decode(yu_BMSN2g(2,:,ii,inu),pattern(ii,2));
        
        yu_BMSN2gc(1,:,ii,inu) = yu_BMSN2gc(1,:,ii,inu)/St_BMSN2gc(1,inu)/sqrt(amp2);
        yu_BMSN2gc(2,:,ii,inu) = yu_BMSN2gc(2,:,ii,inu)/St_BMSN2gc(2,inu)/sqrt(amp2);
        Yu_BMSN2gc(1,:,ii,inu) = Decode(yu_BMSN2gc(1,:,ii,inu),pattern(ii,1));
        Yu_BMSN2gc(2,:,ii,inu) = Decode(yu_BMSN2gc(2,:,ii,inu),pattern(ii,2));
        
    end
end
end
end

% for BD and BMSN with antenna selection (BD-AS, BMSN-AS)
y_bdas = alluser(yu_BD2);
Y_bdas = Decode(y_bdas,bsr);


TR_1 = Eval_to_TR (NR, NU, SINU,Eval_1,SNRT,TRT);
    %csvwrite(Out_1,TR_1);
    % BMSN-GE
    TR_2 = Eval_to_TR (NR, NU, SINU,Eval_2,SNRT,TRT);
    %csvwrite(Out_2,TR_2);
%     % MMSE
%     TR_3 = Eval_to_TR (NR, NU, Ntri,Eval_3,SNRT,TRT);
%     %csvwrite(Out_3,TR_3);
    % BMSN-GE-TD
    TR_3 = Eval_to_TR (NR, NU, SINU,Eval_3,SNRT,TRT);
    % BD
    TR_4 = Eval_to_TR (NR, NU, SINU,Eval_4,SNRT,TRT);
    %csvwrite(Out_4,TR_4);
    % BMSN-AS
%     TR_5 = Eval_to_TR (NRs, NU, SINU,Eval_5,SNRT,TRT);
    % csvwrite(Out_5,TR_5);
    % BD-AS
    TR_5 = Eval_to_TR (NR, NU, SINU,Eval_5,SNRT,TRT);
    % csvwrite(Out_5,TR_5);
    

    % 受信アンテナの和
    for nuser=1:NU
        Q(:,1,nuser) = sum(TR_1(:,(nuser-1)*NR+1:(nuser-1)*NR+NR),2);
        Q(:,2,nuser) = sum(TR_2(:,(nuser-1)*NR+1:(nuser-1)*NR+NR),2);
        Q(:,3,nuser) = sum(TR_3(:,(nuser-1)*NR+1:(nuser-1)*NR+NR),2);
        Q(:,4,nuser) = sum(TR_4(:,(nuser-1)*NR+1:(nuser-1)*NR+NR),2);
%         Q(:,5,nuser) = sum(TR_5(:,(nuser-1)*NRs+1:(nuser-1)*NRs+NRs),2);
        Q(:,5,nuser) = sum(TR_5(:,(nuser-1)*NRs+1:(nuser-1)*NRs+NRs),2);
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
plot(SNR,QmC(:,1),'-o',SNR,QmC(:,2),'-v',SNR,QmC(:,3),'-^',SNR,QmC(:,4),'-s',SNR,QmC(:,5),'k--x','Linewidth',2);
%plot(SNR,QmC,'Linewidth',2);
set(gca,'Fontsize',14,'Fontname','Times New Roman');
legend('BMSN-BF','BMSN-GE','BMSN-GE-TD','BD','BD-AS','Location','Northwest');
xlabel('SNR [dB]','Fontsize',16,'Fontname','Times New Roman');
ylabel('Average transmission rate [Mbps]','Fontsize',16,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
%title(target_CDF);
grid on;
hold on;