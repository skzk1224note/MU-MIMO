clear;

rng('default');

% 出力ファイル名
% evfile1 = 'Eig16x2x8u_BD_SNR0-30dB_100itr.csv';
% evfile2 = 'Eig16x2x8u_BDAS_SNR0-30dB_100itr.csv';
% evfile3 = 'Eig16x2x8u_BMSN2_SNR0-30dB_100itr.csv';
% evfile4 = 'Eig16x2x8u_BMSNAS_SNR0-30dB_100itr.csv';
% evfile5 = 'Eig16x2x8u_ZF_SNR0-30dB_100itr.csv';
% evfile6 = 'Eig16x2x8u_BMSN1b_SNR0-30dB_100itr.csv';

Nt = 16; %基地局アンテナ数
Nr = 2; %ユーザアンテナ数は2(=送信ストリーム数(not1))
Nu = 8; %ユーザ数

ndata = 10000; % 1パケット当たりのシンボル数（パケット長）:BER計算のための1試行当たりのシンボル数 10000
SIMU = 5000; %試行回数（通常 5000）

SNR_min = 0;  %最小SNR[dB]
SNR_max = 30; %最大SNR[dB]

Angle_interval = 0.1;   % 受信角の間隔[deg]

pattern = [4,0;2,2;3,1]; %変調パターン(先頭に[●,0]をもってくる）
bsr = pattern(1,1);
a = 1e-2; % BMSNの擬似雑音
I = eye(Nt,Nt);

% Riceフェージングのための変数定義
K_dB   = 20;         % ライスファクタK
K      = 10^(K_dB/10);
% K = 0(K_dB = -20dB); % レイリーフェージング 

d_t      = 0.5;      % 送信アンテナ間隔（in wavelength)
d_r      = 0.5;      % 受信アンテナ間隔（in wavelength)
derad = pi/180;      % degree -> rad

if K_dB < 10
   target_K=strcat('{\it{K}}=',num2str(K_dB,'%01d'),'dB');
else
   target_K=strcat('{\it{K}}=',num2str(K_dB,'%02d'),'dB');
end


% T 所望のチャネル行列
for inu = 1:Nu
    T(:,:,inu) = zeros(Nr,Nr);
    T(1,1,inu) = 1;
    T(1,2,inu) = 0;
    T(2,1,inu) = 0;
    T(2,2,inu) = 1;
end

tic;

snr = SNR_min:SNR_max;
Lsnr = length(snr);
MSt_BD=zeros(Nr,Lsnr);
%MSt_ZF=zeros(Nr-1,Lsnr);
MSt_BMSN1=zeros(Nr,Lsnr);
MSt_BMSN2=zeros(Nr,Lsnr);
MSt_BMSN3=zeros(Nr,Lsnr);
MSt_MMSE=zeros(Nr,Lsnr);
MSt_BMSN1b=zeros(Nr,Lsnr);

% H (伝搬チャネル行列)
% 伝搬チャネル行列のマルチパス成分 (i.i.d. Rayleigh , NLOS チャネル)
    H_iid = (randn(Nr*Nu,Nt)+1j*randn(Nr*Nu,Nt))/sqrt(2);
% 伝搬チャネル行列の直接波成分(LOS チャネル)
    H_los = zeros(Nu*Nr,Nt);

    
for isnr = 1:Lsnr
    sigma2 = 1/(10^(snr(isnr)/10)); % noise power

for isimu = 1:SIMU
    
    for RR=1:Nu
        Theta_r(1,RR) = -180+Angle_interval*RR; % ユーザ毎の受信角 (5度間隔でユーザを配置)
    end
    
    Theta_t = (rand-0.5)*360;
   % Theta_t = 0;   % ユーザ毎の送信角固定0
   % Theta_r = (rand(1,Nu)-0.5)*360; % ユーザ毎の受信角 (-180deg - 180deg)
    a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t*derad));  % 送信モードベクトル
    
      for n = 1 : Nu
        a_r = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad)); % ユーザ毎の受信モードベクトル
        H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t.';                % ユーザ毎のLOSチャネル行列
      end

    % 伝搬チャネル行列=[sqrt(K/(K+1))*(LOS チャネル)]...
    %                   .+[sqrt(1/(K+1))*(NLOS チャネル)]
    H = sqrt(K/(K+1))*H_los + sqrt(1/(K+1))*H_iid;
    Hu = peruser(H,Nu);

% ZFウエイト
%Wzf = H\I;
%Wzf = H'/(H*H');

% MMSEウエイト
Wmmse = H'/(H*H'+ sigma2*Nt*eye(Nt)); % MMSE

% for inu = 1:Nu
%     %h = (randn(2,nt) +
%     1j*randn(2,nt))/sqrt(2);rrrrrrrrqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq
%     h = H(1+(inu-1)*Nr:2+(inu-1)*Nr,:);
%     S = [svd(h(1,:));svd(h(2,:))];  % nr =2
%     [~,No] = max(S);
%     Hu_as(:,:,inu) = h(No,:);
% end

%H(伝搬チャネル行列)
%H_as = alluser(Hu_as); %伝搬チャネル行列

He = zeros(Nt-Nr,Nt,Nu);

for inu = 1:Nu
    
    %He(Hからユーザkのチャネル行列を除いた行列)
    el = 1:Nu;
    el(:,inu) = [];
    He(:,:,inu) = alluser(Hu(:,:,el));
    %He_as(:,:,inu) = alluser(Hu_as(:,:,el));

    % BD
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
    
    % B-MSN α=1e-2
    % MSNに基づくユーザ毎のウエイトの計算
    % MSNの最適ウエイトWopt
    a1=1e-2;
    Wopt(:,:,inu) = (He(:,:,inu)'*He(:,:,inu)+a1*eye(size(He,2)))\Hu(:,:,inu)'*T(:,:,inu);
%     Wopt(:,:,i) = Wopt(:,:,i)/norm(Wopt(:,:,i),'fro');
    for ij = 1:Nr
        Wopt(:,ij,inu) = Wopt(:,ij,inu)/sqrt(Wopt(:,ij,inu)'*Wopt(:,ij,inu));
    end
    
    %Vts(Hu*Woptの信号部分空間に対応する固有ベクトル)
    %Wtu(ユーザkに対応する基地局側ウエイト)
    %St(特異値を格納 2×nu)
    %Wr(ユーザ側ウエイト)
    [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt(:,:,inu));
    St_BMSN1(:,inu) = diag(Std(1:Nr,1:Nr));
    
    Wr_BMSN1(:,:,inu) = Ut';
    Vts_BMSN1(:,:,inu) = Vt(:,1:Nr);
    Wtu_BMSN1(:,:,inu) = Wopt(:,:,inu)*Vts_BMSN1(:,:,inu);

    % B-MSN2 α=1e-6
    %MSNに基づくユーザ毎のウエイトの計算
    % MSNの最適ウエイトWopt
    a2=1e-6;
    Wopt(:,:,inu) = (He(:,:,inu)'*He(:,:,inu)+a2*eye(size(He,2)))\Hu(:,:,inu)'*T(:,:,inu);
%     Wopt(:,:,i) = Wopt(:,:,i)/norm(Wopt(:,:,i),'fro');
    for ij = 1:Nr
        Wopt(:,ij,inu) = Wopt(:,ij,inu)/sqrt(Wopt(:,ij,inu)'*Wopt(:,ij,inu));
    end
    
    %Vts(Hu*Woptの信号部分空間に対応する固有ベクトル)
    %Wtu(ユーザkに対応する基地局側ウエイト)
    %St(特異値を格納 2×nu)
    %Wr(ユーザ側ウエイト)
    [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt(:,:,inu));
    St_BMSN2(:,inu) = diag(Std(1:Nr,1:Nr));
    
    Wr_BMSN2(:,:,inu) = Ut';
    Vts_BMSN2(:,:,inu) = Vt(:,1:Nr);
    Wtu_BMSN2(:,:,inu) = Wopt(:,:,inu)*Vts_BMSN2(:,:,inu);
   
    % BMSN-BF法 α=sigma2 NT equal to MMSE
    %MSNに基づくユーザ毎のウエイトの計算
    % MSNの最適ウエイトWopt
    a3=sigma2*Nt;
    Wopt(:,:,inu) = (He(:,:,inu)'*He(:,:,inu)+a3*eye(size(He,2)))\Hu(:,:,inu)'*T(:,:,inu);
%     Wopt(:,:,i) = Wopt(:,:,i)/norm(Wopt(:,:,i),'fro');
    for ij = 1:Nr
        Wopt(:,ij,inu) = Wopt(:,ij,inu)/sqrt(Wopt(:,ij,inu)'*Wopt(:,ij,inu));
    end
    
    %Vts(Hu*Woptの信号部分空間に対応する固有ベクトル)
    %Wtu(ユーザkに対応する基地局側ウエイト)
    %St(特異値を格納 2×nu)
    %Wr(ユーザ側ウエイト)
    [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt(:,:,inu));
    St_BMSN3(:,inu) = diag(Std(1:Nr,1:Nr));
    
    Wr_BMSN3(:,:,inu) = Ut';
    Vts_BMSN3(:,:,inu) = Vt(:,1:Nr);
    Wtu_BMSN3(:,:,inu) = Wopt(:,:,inu)*Vts_BMSN3(:,:,inu);
    
    % BMSN-GE :一般化固有値問題を解く方法
    B(:,:,inu) = He(:,:,inu)'*He(:,:,inu)+a3*eye(size(He(:,:,inu),2));
    A(:,:,inu) = Hu(:,:,inu)'*Hu(:,:,inu);
    [Wopt_k(:,:,inu),D] = eig(B(:,:,inu),A(:,:,inu));
    [D1,IN] = sort(diag(D));
    D2(:,:,inu) = D1.';
    Wopt_k(:,:,inu) = Wopt_k(:,IN,inu);
    Wopt_k2(:,:,inu)=Wopt_k(:,1:Nr,inu);
    for ij = 1:Nr
        Wopt_k2(:,ij,inu) = Wopt_k2(:,ij,inu)/sqrt(Wopt_k2(:,ij,inu)'*Wopt_k2(:,ij,inu));
    end             
    % nuのウエイト(信号部分空間を利用)                                 
    %Vts(HT*Wの信号部分空間に対応する固有ベクトル)
    %Wtu(ユーザkに対応する基地局側ウエイト)
    %St(特異値を格納 2×nu)
    %Wr(ユーザ側ウエイト)
    [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt_k2(:,:,inu));
    St_BMSN1b(:,inu) = diag(Std(1:Nr,1:Nr));
    Wr_BMSN1b(:,:,inu) = Ut';
    Vts_BMSN1b(:,:,inu) = Vt(:,1:Nr);
    Wtu_BMSN1b(:,:,inu) = Wopt_k2(:,:,inu)*Vts_BMSN1b(:,:,inu);
    
    % ZF(zero-Forcing)
    % inuのウエイト（列ベクトル毎に正規化）
%     Wopt(:,:,inu) = Wmmse(:,(inu-1)*Nr+1:(inu-1)*Nr+Nr);
%     for ij = 1:Nr
%         Wopt(:,ij,inu) = Wopt(:,ij,inu)/sqrt(Wopt(:,ij,inu)'*Wopt(:,ij,inu));
%     end
% 
%     % 雑音部分空間ベクトルをnuerのチャネル行列に乗算                                 
%     HTT=Hu(:,:,inu)*Wopt(:,:,inu);      
%     % 変換行列をSVD
%     [Ut,Std,Vt]=svd(HTT);
%     St_ZF(:,inu) = diag(Std(1:Nr,1:Nr));
%     Wr_ZF(:,:,inu) = Ut';
%     Vts_ZF(:,:,inu) = Vt(:,1:Nr);
%     Wtu_ZF(:,:,inu) = Wopt(:,:,inu)*Vts_ZF(:,:,inu);
               
    % MMSE
    % inuのウエイト（列ベクトル毎に正規化）
    Wopt(:,:,inu) = Wmmse(:,(inu-1)*Nr+1:(inu-1)*Nr+Nr);
    for ij = 1:Nr
        Wopt(:,ij,inu) = Wopt(:,ij,inu)/sqrt(Wopt(:,ij,inu)'*Wopt(:,ij,inu));
    end

    % 雑音部分空間ベクトルをnuerのチャネル行列に乗算                                 
    HTT=Hu(:,:,inu)*Wopt(:,:,inu);      
    % 変換行列をSVD
    [Ut,Std,Vt]=svd(HTT);
    St_MMSE(:,inu) = diag(Std(1:Nr,1:Nr));
    Wr_MMSE(:,:,inu) = Ut';
    Vts_MMSE(:,:,inu) = Vt(:,1:Nr);
    Wtu_MMSE(:,:,inu) = Wopt(:,:,inu)*Vts_MMSE(:,:,inu);
    
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
% data_as = randi([0 2^bsr-1],Nu,ndata);
% 
% %s(送信ストリーム) for antenna selection (AS)
% s_as = Mapping(data_as,bsr);


%n(熱雑音)
%y(ユーザkの入力信号)
%Y(ユーザkの復調信号)
for inu = 1:Nu
    n = (randn(Nr,ndata) + 1j*randn(Nr,ndata))/sqrt(2);
    
    for ii = 1:size(pattern,1)
        Hs_BD = zeros(Nr,ndata);    % ユーザiの受信信号（受信重み付け前）
        Hs_BMSN1 = zeros(Nr,ndata);    % ユーザiの受信信号（受信重み付け前）
        Hs_BMSN2 = zeros(Nr,ndata);    % ユーザiの受信信号（受信重み付け前）
        Hs_BMSN3 = zeros(Nr,ndata);    % ユーザiの受信信号（受信重み付け前）
        Hs_MMSE = zeros(Nr,ndata);    % ユーザiの受信信号（受信重み付け前）
        Hs_BMSN1b = zeros(Nr,ndata);    % ユーザiの受信信号（受信重み付け前）
        for ij = 1:Nu
            Hs_BD = Hs_BD + Hu(:,:,inu)*Wtu_BD(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
            Hs_BMSN1 = Hs_BMSN1 + Hu(:,:,inu)*Wtu_BMSN1(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
            Hs_BMSN2 = Hs_BMSN2 + Hu(:,:,inu)*Wtu_BMSN2(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
            Hs_BMSN3 = Hs_BMSN3 + Hu(:,:,inu)*Wtu_BMSN3(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
            Hs_MMSE = Hs_MMSE + Hu(:,:,inu)*Wtu_MMSE(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
            Hs_BMSN1b = Hs_BMSN1b + Hu(:,:,inu)*Wtu_BMSN1b(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
        end
        yu_BD(:,:,ii,inu) = Wr_BD(:,:,inu)*(Hs_BD+n);
        yu_BMSN1(:,:,ii,inu) = Wr_BMSN1(:,:,inu)*(Hs_BMSN1+n);
        yu_BMSN2(:,:,ii,inu) = Wr_BMSN2(:,:,inu)*(Hs_BMSN2+n);
        yu_BMSN3(:,:,ii,inu) = Wr_BMSN3(:,:,inu)*(Hs_BMSN3+n);
        yu_MMSE(:,:,ii,inu) = Wr_MMSE(:,:,inu)*(Hs_MMSE+n);
        yu_BMSN1b(:,:,ii,inu) = Wr_BMSN1b(:,:,inu)*(Hs_BMSN1b+n);

        yu_BD(1,:,ii,inu) = yu_BD(1,:,ii,inu)/St_BD(1,inu)/sqrt(amp2);
        yu_BD(2,:,ii,inu) = yu_BD(2,:,ii,inu)/St_BD(2,inu)/sqrt(amp2);
        Yu_BD(1,:,ii,inu) = Decode(yu_BD(1,:,ii,inu),pattern(ii,1));
        Yu_BD(2,:,ii,inu) = Decode(yu_BD(2,:,ii,inu),pattern(ii,2));
    
        yu_BMSN1(1,:,ii,inu) = yu_BMSN1(1,:,ii,inu)/St_BMSN1(1,inu)/sqrt(amp2);
        yu_BMSN1(2,:,ii,inu) = yu_BMSN1(2,:,ii,inu)/St_BMSN1(2,inu)/sqrt(amp2);
        Yu_BMSN1(1,:,ii,inu) = Decode(yu_BMSN1(1,:,ii,inu),pattern(ii,1));
        Yu_BMSN1(2,:,ii,inu) = Decode(yu_BMSN1(2,:,ii,inu),pattern(ii,2));

        yu_BMSN2(1,:,ii,inu) = yu_BMSN2(1,:,ii,inu)/St_BMSN2(1,inu)/sqrt(amp2);
        yu_BMSN2(2,:,ii,inu) = yu_BMSN2(2,:,ii,inu)/St_BMSN2(2,inu)/sqrt(amp2);
        Yu_BMSN2(1,:,ii,inu) = Decode(yu_BMSN2(1,:,ii,inu),pattern(ii,1));
        Yu_BMSN2(2,:,ii,inu) = Decode(yu_BMSN2(2,:,ii,inu),pattern(ii,2));
 
        yu_BMSN3(1,:,ii,inu) = yu_BMSN3(1,:,ii,inu)/St_BMSN3(1,inu)/sqrt(amp2);
        yu_BMSN3(2,:,ii,inu) = yu_BMSN3(2,:,ii,inu)/St_BMSN3(2,inu)/sqrt(amp2);
        Yu_BMSN3(1,:,ii,inu) = Decode(yu_BMSN3(1,:,ii,inu),pattern(ii,1));
        Yu_BMSN3(2,:,ii,inu) = Decode(yu_BMSN3(2,:,ii,inu),pattern(ii,2));
        
        yu_MMSE(1,:,ii,inu) = yu_MMSE(1,:,ii,inu)/St_MMSE(1,inu)/sqrt(amp2);
        yu_MMSE(2,:,ii,inu) = yu_MMSE(2,:,ii,inu)/St_MMSE(2,inu)/sqrt(amp2);
        Yu_MMSE(1,:,ii,inu) = Decode(yu_MMSE(1,:,ii,inu),pattern(ii,1));
        Yu_MMSE(2,:,ii,inu) = Decode(yu_MMSE(2,:,ii,inu),pattern(ii,2));
        
        yu_BMSN1b(1,:,ii,inu) = yu_BMSN1b(1,:,ii,inu)/St_BMSN1b(1,inu)/sqrt(amp2);
        yu_BMSN1b(2,:,ii,inu) = yu_BMSN1b(2,:,ii,inu)/St_BMSN1b(2,inu)/sqrt(amp2);
        Yu_BMSN1b(1,:,ii,inu) = Decode(yu_BMSN1b(1,:,ii,inu),pattern(ii,1));
        Yu_BMSN1b(2,:,ii,inu) = Decode(yu_BMSN1b(2,:,ii,inu),pattern(ii,2));
        
    end
end

%ber(変調パターン毎のBER)  for BD and B-MSN
for inu = 1:Nu
for ii = 1:size(pattern,1)
if ii == 1
    ber_BD(ii,inu) = BER2(BER1(Yu_BD(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
    ber_BMSN1(ii,inu) = BER2(BER1(Yu_BMSN1(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
    ber_BMSN2(ii,inu) = BER2(BER1(Yu_BMSN2(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
    ber_BMSN3(ii,inu) = BER2(BER1(Yu_BMSN3(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
    ber_MMSE(ii,inu) = BER2(BER1(Yu_MMSE(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
    ber_BMSN1b(ii,inu) = BER2(BER1(Yu_BMSN1b(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
else
    ber_BD(ii,inu) = BER2(BER1(Yu_BD(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
        + BER2(BER1(Yu_BD(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
    ber_BMSN1(ii,inu) = BER2(BER1(Yu_BMSN1(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
        + BER2(BER1(Yu_BMSN1(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
    ber_BMSN2(ii,inu) = BER2(BER1(Yu_BMSN2(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
        + BER2(BER1(Yu_BMSN2(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
    ber_BMSN3(ii,inu) = BER2(BER1(Yu_BMSN3(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
        + BER2(BER1(Yu_BMSN3(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
    ber_MMSE(ii,inu) = BER2(BER1(Yu_MMSE(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
        + BER2(BER1(Yu_MMSE(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
    ber_BMSN1b(ii,inu) = BER2(BER1(Yu_BMSN1b(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
        + BER2(BER1(Yu_BMSN1b(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
end
end
end

%ber_min(誤り最小の変調パターンのBER) for BD and B-MSN
ber_min_BD(isimu,isnr) = mean(min(ber_BD));
ber_min_BMSN1(isimu,isnr) = mean(min(ber_BMSN1));
ber_min_BMSN2(isimu,isnr) = mean(min(ber_BMSN2));
ber_min_BMSN3(isimu,isnr) = mean(min(ber_BMSN3));
ber_min_MMSE(isimu,isnr) = mean(min(ber_MMSE));
ber_min_BMSN1b(isimu,isnr) = mean(min(ber_BMSN1b));

%受信アンテナ（nr)の平均特異値 for BD and B-MSN
%10*log10(MSt_BD1(:,inudiag(S_BD(:,:,nuser)).^2/(NT*snt));   
MSt_BD(:,isnr)=MSt_BD(:,isnr)+mean(St_BD(:,inu),2);
MSt_BMSN1(:,isnr)=MSt_BMSN1(:,isnr)+mean(St_BMSN1(:,inu),2);
MSt_BMSN2(:,isnr)=MSt_BMSN2(:,isnr)+mean(St_BMSN2(:,inu),2);
MSt_BMSN3(:,isnr)=MSt_BMSN3(:,isnr)+mean(St_BMSN3(:,inu),2);
MSt_MMSE(:,isnr)=MSt_MMSE(:,isnr)+mean(St_MMSE(:,inu),2);
MSt_BMSN1b(:,isnr)=MSt_BMSN1b(:,isnr)+mean(St_BMSN1b(:,inu),2);

%fprintf('Iteration = %d / %d\n',isimu, SIMU);

end % isimu
MSt_BD(:,isnr)=10*log10((MSt_BD(:,isnr)/SIMU).^2/(Nt*sigma2));
MSt_BMSN1(:,isnr)=10*log10((MSt_BMSN1(:,isnr)/SIMU).^2/(Nt*sigma2));
MSt_BMSN2(:,isnr)=10*log10((MSt_BMSN2(:,isnr)/SIMU).^2/(Nt*sigma2));
MSt_BMSN3(:,isnr)=10*log10((MSt_BMSN3(:,isnr)/SIMU).^2/(Nt*sigma2));
MSt_MMSE(:,isnr)=10*log10((MSt_MMSE(:,isnr)/SIMU).^2/(Nt*sigma2));
MSt_BMSN1b(:,isnr)=10*log10((MSt_BMSN1b(:,isnr)/SIMU).^2/(Nt*sigma2));
fprintf('SNR = %d dB\n',snr(isnr));
end % isnr

% csvwrite(evfile1,[snr;MSt_BD1]);
% csvwrite(evfile2,[snr;MSt_BD2]);
% csvwrite(evfile3,[snr;MSt_MSN1]);
% csvwrite(evfile4,[snr;MSt_MSN2]);
% csvwrite(evfile5,[snr;MSt_ZF]);
% csvwrite(evfile6,[snr;MSt_MMSE]);

%BER(試行回数平均)
BER_BD = mean(ber_min_BD);
BER_BMSN1 = mean(ber_min_BMSN1);
BER_BMSN2 = mean(ber_min_BMSN2);
BER_BMSN3 = mean(ber_min_BMSN3);
BER_MMSE = mean(ber_min_MMSE);
BER_BMSN1b = mean(ber_min_BMSN1b);

toc;

%描写
figure;
mycol = [1 0 1;
         0 1 0;
         1 0 0;
         0 0 1;
         0 1 1;
         0 0 0]; % 色
set(groot,'defaultAxesColorOrder',mycol)
semilogy(snr,BER_BD,'-',snr,BER_MMSE,'-',snr,BER_BMSN3,snr,BER_BMSN1b,'-','Linewidth',2);
axis([SNR_min SNR_max 1e-1 1e0]);
set(gca,'XTick',SNR_min:5:SNR_max,'YTick',[1e-1 1e0],'Fontsize',14,'Fontname','Times New Roman')
legend('BD','MMSE-CI','BMSN-BF','BMSN-GE','Location','southwest')
xlabel('SNR [dB]','Fontsize',16,'Fontname','Times New Roman');
ylabel('Average BER','Fontsize',16,'Fontname','Times New Roman');
title(target_K);
grid on;
hold on;

% End