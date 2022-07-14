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
Nr = 1; %ユーザアンテナ数は2(=送信ストリーム数(not1))
Nu = 8; %ユーザ数

K_min = -20;      % 最小K[dB]
K_max = 20;       % 最大K[dB]
SNR_tar = 15;     % ターゲットSNR[dB]
ndata = 10000; % 1パケット当たりのシンボル数（パケット長）:BER計算のための1試行当たりのシンボル数 10000
SIMU = 10 ;%試行回数（通常 5000）
if Nr==2
    pattern = [4,0;2,2;3,1]; %変調パターン(先頭に[●,0]をもってくる）
end
if Nr==1
    pattern = 4;
end
bsr = pattern(1,1);
sigma2 = 1/(10^(SNR_tar/10)); % noise power

I = eye(Nt,Nt);

d_t      = 0.5;      % 送信アンテナ間隔（in wavelength)
d_r      = 0.5;      % 受信アンテナ間隔（in wavelength)
derad = pi/180;      % degree -> rad

if SNR_tar < 10
   target_SNR=strcat('SNR=',num2str(SNR_tar,'%01d'),'dB');
else
   target_SNR=strcat('SNR=',num2str(SNR_tar,'%02d'),'dB');
end

%T 所望のチャネル行列
T = zeros(Nr,Nr);
for inu = 1:Nu
    T(:,:,inu) = eye(Nr,Nr);
end


tic;
K_box=(K_min:5:K_max).'; % figureの横軸のためのKの箱
LK=length(K_box);        % Kの箱の大きさ

MSt_BD=zeros(Nr,LK);
% MSt_ZF=zeros(Nr-1,LK);
% MSt_BMSN1=zeros(Nr,LK);
% MSt_BMSN2=zeros(Nr,LK);
MSt_BMSN_BF=zeros(Nr,LK);
MSt_MMSE=zeros(Nr,LK);
MSt_BMSN_GE=zeros(Nr,LK);

% H (伝搬チャネル行列)
% 伝搬チャネル行列のマルチパス成分 (i.i.d. Rayleigh , NLOS チャネル)
    H_iid = (randn(Nr*Nu,Nt)+1j*randn(Nr*Nu,Nt))/sqrt(2);
% 伝搬チャネル行列の直接波成分(LOS チャネル)
    H_los = zeros(Nu*Nr,Nt);

    K_box=(K_min:5:K_max).'; % figureの横軸のためのKの箱
    LK=length(K_box);        % Kの箱の大きさ

for ik = 1:LK
    
    K_tar = K_box(ik);

for isimu = 1:SIMU

    Theta_t = (rand(1,Nu)-0.5)*360; % ユーザ毎の送信角 (-180deg - 180deg)
    Theta_r = (rand(1,Nu)-0.5)*360; % ユーザ毎の受信角 (-180deg - 180deg)
    a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t*derad));  % 送信モードベクトル
    
      for n = 1 : Nu
          a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad));  % 送信モードベクトル
          a_r = exp(-1j*2*pi*d_r*(0:Nr-1).'*sin(Theta_r(1,n)*derad)); % ユーザ毎の受信モードベクトル
          H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t.';                % ユーザ毎のLOSチャネル行列
      end

      
      % K真値
      K = 10^(K_tar/10);
    
    % 伝搬チャネル行列=[sqrt(K/(K+1))*(LOS チャネル)]...
    %                   .+[sqrt(1/(K+1))*(NLOS チャネル)]
    H = sqrt(K/(K+1))*H_los + sqrt(1/(K+1))*H_iid;
    Hu = peruser(H,Nu);

% ZFウエイト
%Wzf = H\I;
%Wzf = H'/(H*H');

% MMSEウエイト
Wmmse = H'/(H*H'+ sigma2*Nt*eye(Nr*Nu,Nr*Nu)); % MMSE

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

He = zeros(Nr*Nu-Nr,Nt,Nu);

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
%     a1=1e-2;
%     Wopt(:,:,inu) = (He(:,:,inu)'*He(:,:,inu)+a1*eye(size(He,2)))\Hu(:,:,inu)'*T(:,:,inu);
% %     Wopt(:,:,i) = Wopt(:,:,i)/norm(Wopt(:,:,i),'fro');
%     for ij = 1:Nr
%         Wopt(:,ij,inu) = Wopt(:,ij,inu)/sqrt(Wopt(:,ij,inu)'*Wopt(:,ij,inu));
%     end
%     
%     %Vts(Hu*Woptの信号部分空間に対応する固有ベクトル)
%     %Wtu(ユーザkに対応する基地局側ウエイト)
%     %St(特異値を格納 2×nu)
%     %Wr(ユーザ側ウエイト)
%     [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt(:,:,inu));
%     St_BMSN1(:,inu) = diag(Std(1:Nr,1:Nr));
%     
%     Wr_BMSN1(:,:,inu) = Ut';
%     Vts_BMSN1(:,:,inu) = Vt(:,1:Nr);
%     Wtu_BMSN1(:,:,inu) = Wopt(:,:,inu)*Vts_BMSN1(:,:,inu);

    % B-MSN2 α=1e-6
    %MSNに基づくユーザ毎のウエイトの計算
    % MSNの最適ウエイトWopt
%     a2=1e-6;
%     Wopt(:,:,inu) = (He(:,:,inu)'*He(:,:,inu)+a2*eye(size(He,2)))\Hu(:,:,inu)'*T(:,:,inu);
% %     Wopt(:,:,i) = Wopt(:,:,i)/norm(Wopt(:,:,i),'fro');
%     for ij = 1:Nr
%         Wopt(:,ij,inu) = Wopt(:,ij,inu)/sqrt(Wopt(:,ij,inu)'*Wopt(:,ij,inu));
%     end
%     
%     %Vts(Hu*Woptの信号部分空間に対応する固有ベクトル)
%     %Wtu(ユーザkに対応する基地局側ウエイト)
%     %St(特異値を格納 2×nu)
%     %Wr(ユーザ側ウエイト)
%     [Ut,Std,Vt] = svd(Hu(:,:,inu)*Wopt(:,:,inu));
%     St_BMSN2(:,inu) = diag(Std(1:Nr,1:Nr));
%     
%     Wr_BMSN2(:,:,inu) = Ut';
%     Vts_BMSN2(:,:,inu) = Vt(:,1:Nr);
%     Wtu_BMSN2(:,:,inu) = Wopt(:,:,inu)*Vts_BMSN2(:,:,inu);
   
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
    St_BMSN_BF(:,inu) = diag(Std(1:Nr,1:Nr));    
    Wr_BMSN_BF(:,:,inu) = Ut';
    Vts_BMSN_BF(:,:,inu) = Vt(:,1:Nr);
    Wtu_BMSN_BF(:,:,inu) = Wopt(:,:,inu)*Vts_BMSN_BF(:,:,inu);
    
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
    St_BMSN_GE(:,inu) = diag(Std(1:Nr,1:Nr));
    Wr_BMSN_GE(:,:,inu) = Ut';
    Vts_BMSN_GE(:,:,inu) = Vt(:,1:Nr);
    Wtu_BMSN_GE(:,:,inu) = Wopt_k2(:,:,inu)*Vts_BMSN_GE(:,:,inu);
    
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
amp2 = 10^(SNR_tar/10)/Nt;

% 以降 nr=2 の場合のみ
if Nr==2
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
end

if Nr==1
    %data(ユーザkへの送信データ)
    for inu = 1:Nu
        for ii = 1:size(pattern,1)
            data1 = randi([0 2^pattern(ii,1)-1],1,ndata);
            data(:,:,ii,inu) = [data1];
        end
    end
    %s(ユーザkへの送信ストリーム)
    for inu = 1:Nu
        for ii = 1:size(pattern,1)
            s1 = Mapping(data(1,:,ii,inu),pattern(ii,1));
            s(:,:,ii,inu) = [s1];
        end
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
        % Hs_BMSN1 = zeros(Nr,ndata);    % ユーザiの受信信号（受信重み付け前）
        % Hs_BMSN2 = zeros(Nr,ndata);    % ユーザiの受信信号（受信重み付け前）
        Hs_BMSN_BF = zeros(Nr,ndata);    % ユーザiの受信信号（受信重み付け前）
        Hs_MMSE = zeros(Nr,ndata);    % ユーザiの受信信号（受信重み付け前）
        Hs_BMSN_GE = zeros(Nr,ndata);    % ユーザiの受信信号（受信重み付け前）
        for ij = 1:Nu
            Hs_BD = Hs_BD + Hu(:,:,inu)*Wtu_BD(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
            % Hs_BMSN1 = Hs_BMSN1 + Hu(:,:,inu)*Wtu_BMSN1(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
            % Hs_BMSN2 = Hs_BMSN2 + Hu(:,:,inu)*Wtu_BMSN2(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
            Hs_BMSN_BF = Hs_BMSN_BF + Hu(:,:,inu)*Wtu_BMSN_BF(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
            Hs_MMSE = Hs_MMSE + Hu(:,:,inu)*Wtu_MMSE(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
            Hs_BMSN_GE = Hs_BMSN_GE + Hu(:,:,inu)*Wtu_BMSN_GE(:,:,ij)*s(:,:,ii,ij)*sqrt(amp2);
        end
        yu_BD(:,:,ii,inu) = Wr_BD(:,:,inu)*(Hs_BD+n);
        % yu_BMSN1(:,:,ii,inu) = Wr_BMSN1(:,:,inu)*(Hs_BMSN1+n);
        % yu_BMSN2(:,:,ii,inu) = Wr_BMSN2(:,:,inu)*(Hs_BMSN2+n);
        yu_BMSN_BF(:,:,ii,inu) = Wr_BMSN_BF(:,:,inu)*(Hs_BMSN_BF+n);
        yu_MMSE(:,:,ii,inu) = Wr_MMSE(:,:,inu)*(Hs_MMSE+n);
        yu_BMSN_GE(:,:,ii,inu) = Wr_BMSN_GE(:,:,inu)*(Hs_BMSN_GE+n);

        yu_BD(1,:,ii,inu) = yu_BD(1,:,ii,inu)/St_BD(1,inu)/sqrt(amp2);
        Yu_BD(1,:,ii,inu) = Decode(yu_BD(1,:,ii,inu),pattern(ii,1));            
        % yu_BMSN1(1,:,ii,inu) = yu_BMSN1(1,:,ii,inu)/St_BMSN1(1,inu)/sqrt(amp2);
        % Yu_BMSN1(1,:,ii,inu) = Decode(yu_BMSN1(1,:,ii,inu),pattern(ii,1));        
        % yu_BMSN2(1,:,ii,inu) = yu_BMSN2(1,:,ii,inu)/St_BMSN2(1,inu)/sqrt(amp2);        
        % Yu_BMSN2(1,:,ii,inu) = Decode(yu_BMSN2(1,:,ii,inu),pattern(ii,1));
        yu_BMSN_BF(1,:,ii,inu) = yu_BMSN_BF(1,:,ii,inu)/St_BMSN_BF(1,inu)/sqrt(amp2);        
        Yu_BMSN_BF(1,:,ii,inu) = Decode(yu_BMSN_BF(1,:,ii,inu),pattern(ii,1));       
        yu_MMSE(1,:,ii,inu) = yu_MMSE(1,:,ii,inu)/St_MMSE(1,inu)/sqrt(amp2);        
        Yu_MMSE(1,:,ii,inu) = Decode(yu_MMSE(1,:,ii,inu),pattern(ii,1));       
        yu_BMSN_GE(1,:,ii,inu) = yu_BMSN_GE(1,:,ii,inu)/St_BMSN_GE(1,inu)/sqrt(amp2);        
        Yu_BMSN_GE(1,:,ii,inu) = Decode(yu_BMSN_GE(1,:,ii,inu),pattern(ii,1));
        
        if Nr==2
            yu_BD(2,:,ii,inu) = yu_BD(2,:,ii,inu)/St_BD(2,inu)/sqrt(amp2);
            Yu_BD(2,:,ii,inu) = Decode(yu_BD(2,:,ii,inu),pattern(ii,2));
            % yu_BMSN1(2,:,ii,inu) = yu_BMSN1(2,:,ii,inu)/St_BMSN1(2,inu)/sqrt(amp2);
            % Yu_BMSN1(2,:,ii,inu) = Decode(yu_BMSN1(2,:,ii,inu),pattern(ii,2));
            % yu_BMSN2(2,:,ii,inu) = yu_BMSN2(2,:,ii,inu)/St_BMSN2(2,inu)/sqrt(amp2);
            % Yu_BMSN2(2,:,ii,inu) = Decode(yu_BMSN2(2,:,ii,inu),pattern(ii,2));
            yu_BMSN_BF(2,:,ii,inu) = yu_BMSN_BF(2,:,ii,inu)/St_BMSN_BF(2,inu)/sqrt(amp2);
            Yu_BMSN_BF(2,:,ii,inu) = Decode(yu_BMSN_BF(2,:,ii,inu),pattern(ii,2));
            yu_MMSE(2,:,ii,inu) = yu_MMSE(2,:,ii,inu)/St_MMSE(2,inu)/sqrt(amp2);
            Yu_MMSE(2,:,ii,inu) = Decode(yu_MMSE(2,:,ii,inu),pattern(ii,2));
            yu_BMSN_GE(2,:,ii,inu) = yu_BMSN_GE(2,:,ii,inu)/St_BMSN_GE(2,inu)/sqrt(amp2);
            Yu_BMSN_GE(2,:,ii,inu) = Decode(yu_BMSN_GE(2,:,ii,inu),pattern(ii,2));
        end
   end
end


%ber(変調パターン毎のBER)  for BD and B-MSN
for inu = 1:Nu
for ii = 1:size(pattern,1)
if ii == 1
    ber_BD(ii,inu) = BER2(BER1(Yu_BD(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
    % ber_BMSN1(ii,inu) = BER2(BER1(Yu_BMSN1(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
    % ber_BMSN2(ii,inu) = BER2(BER1(Yu_BMSN2(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
    ber_BMSN_BF(ii,inu) = BER2(BER1(Yu_BMSN_BF(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
    ber_MMSE(ii,inu) = BER2(BER1(Yu_MMSE(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
    ber_BMSN_GE(ii,inu) = BER2(BER1(Yu_BMSN_GE(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(1,1));
else
    ber_BD(ii,inu) = BER2(BER1(Yu_BD(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
        + BER2(BER1(Yu_BD(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
    % ber_BMSN1(ii,inu) = BER2(BER1(Yu_BMSN1(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
    %    + BER2(BER1(Yu_BMSN1(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
    % ber_BMSN2(ii,inu) = BER2(BER1(Yu_BMSN2(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
    %    + BER2(BER1(Yu_BMSN2(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
    ber_BMSN_BF(ii,inu) = BER2(BER1(Yu_BMSN_BF(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
        + BER2(BER1(Yu_BMSN_BF(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
    ber_MMSE(ii,inu) = BER2(BER1(Yu_MMSE(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
        + BER2(BER1(Yu_MMSE(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
    ber_BMSN_GE(ii,inu) = BER2(BER1(Yu_BMSN_GE(1,:,ii,inu),data(1,:,ii,inu),pattern(ii,1)),pattern(ii,1)) ...
        + BER2(BER1(Yu_BMSN_GE(2,:,ii,inu),data(2,:,ii,inu),pattern(ii,2)),pattern(ii,2));
end
end
end

%ber_min(誤り最小の変調パターンのBER) for BD and B-MSN
ber_min_BD(isimu,ik) = mean(min(ber_BD));
% ber_min_BMSN1(isimu,ik) = mean(min(ber_BMSN1));
% ber_min_BMSN2(isimu,ik) = mean(min(ber_BMSN2));
ber_min_BMSN_BF(isimu,ik) = mean(min(ber_BMSN_BF));
ber_min_MMSE(isimu,ik) = mean(min(ber_MMSE));
ber_min_BMSN_GE(isimu,ik) = mean(min(ber_BMSN_GE));

%受信アンテナ（nr)の平均特異値 for BD and B-MSN
%10*log10(MSt_BD1(:,inudiag(S_BD(:,:,nuser)).^2/(NT*snt));   
MSt_BD(:,ik)=MSt_BD(:,ik)+mean(St_BD(:,inu),2);
% MSt_BMSN1(:,ik)=MSt_BMSN1(:,ik)+mean(St_BMSN1(:,inu),2);
% MSt_BMSN2(:,ik)=MSt_BMSN2(:,ik)+mean(St_BMSN2(:,inu),2);
MSt_BMSN_BF(:,ik)=MSt_BMSN_BF(:,ik)+mean(St_BMSN_BF(:,inu),2);
MSt_MMSE(:,ik)=MSt_MMSE(:,ik)+mean(St_MMSE(:,inu),2);
MSt_BMSN_GE(:,ik)=MSt_BMSN_GE(:,ik)+mean(St_BMSN_GE(:,inu),2);

%fprintf('Iteration = %d / %d\n',isimu, SIMU);

end % isimu
MSt_BD(:,ik)=10*log10((MSt_BD(:,ik)/SIMU).^2/(Nt*sigma2));
% MSt_BMSN1(:,ik)=10*log10((MSt_BMSN1(:,ik)/SIMU).^2/(Nt*sigma2));
% MSt_BMSN2(:,ik)=10*log10((MSt_BMSN2(:,ik)/SIMU).^2/(Nt*sigma2));
MSt_BMSN_BF(:,ik)=10*log10((MSt_BMSN_BF(:,ik)/SIMU).^2/(Nt*sigma2));
MSt_MMSE(:,ik)=10*log10((MSt_MMSE(:,ik)/SIMU).^2/(Nt*sigma2));
MSt_BMSN_GE(:,ik)=10*log10((MSt_BMSN_GE(:,ik)/SIMU).^2/(Nt*sigma2));
fprintf('K = %d dB\n',K_box(ik));
end % ik

% csvwrite(evfile1,[snr;MSt_BD1]);
% csvwrite(evfile2,[snr;MSt_BD2]);
% csvwrite(evfile3,[snr;MSt_MSN1]);
% csvwrite(evfile4,[snr;MSt_MSN2]);
% csvwrite(evfile5,[snr;MSt_ZF]);
% csvwrite(evfile6,[snr;MSt_MMSE]);

%BER(試行回数平均)
BER_BD = mean(ber_min_BD);
% BER_BMSN1 = mean(ber_min_BMSN1);
% BER_BMSN2 = mean(ber_min_BMSN2);
BER_BMSN_BF = mean(ber_min_BMSN_BF);
BER_MMSE = mean(ber_min_MMSE);
BER_BMSN_GE = mean(ber_min_BMSN_GE);

toc;

%描写
figure;
mycol = [0 0 1;1 0 0;0 0.7 0;1 0 1;0 0 0]; % グラフの色
set(groot,'defaultAxesColorOrder',mycol)
semilogy(K_box,BER_BMSN_BF,'r-',K_box,BER_BMSN_GE,'b-',K_box,BER_MMSE,'g-',K_box,BER_BD,'c-','Linewidth',2);
axis([K_min K_max 1e-4 1e0]);
set(gca,'XTick',K_min:5:K_max,'YTick',[1e-4, 1e-3, 1e-2, 1e-1 1e0],'Fontsize',14,'Fontname','Times New Roman')
legend('BMSN-BF','BMSN-GE','MMSE-CI','BD','Location','southwest')
xlabel('{\it{K}} [dB]','Fontsize',16,'Fontname','Times New Roman');
ylabel('Average BER','Fontsize',16,'Fontname','Times New Roman');
title(target_SNR);
grid on;
hold on;

% End