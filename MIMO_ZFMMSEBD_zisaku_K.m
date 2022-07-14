clear;

nt = 4; %基地局アンテナ数
%ユーザアンテナ数は2
nu = 2; %ユーザ数
ndata = 250; %パケット数
SNR = 30; %試行SNR[dB]
SIMU = 100; %試行回数
bsr = 4; %bits/symbol/user


% Riceフェージングのための変数定義
K_dB   = -10000000000000;         % RicianのKファクタ
K      = 10^(K_dB/10);
% K = 0(K_dB = -inf); % レイリーフェージング  
d_t      = 0.5;      % 送信アンテナ間隔（in wavelength)
d_r      = 0.5;      % 受信アンテナ間隔（in wavelength)
derad = pi/180;      % degree -> rad

ber_zf = zeros(1,1);
ber_mmse = zeros(1,1);
ber_bd = zeros(1,1);

% H (伝搬チャネル行列:Rician channel)
% 伝搬チャネル行列のマルチパス成分 (i.i.d. Rayleigh , NLOS チャネル)
    H_iid = (randn(2*nu,nt)+1j*randn(2*nu,nt))/sqrt(2);
% 伝搬チャネル行列の直接波成分(LOS チャネル)
    H_los = zeros(nu*2,nt);

for snr = 0:SNR
for simu = 1:SIMU

 Theta_t = (rand-0.5)*360;   % ユーザ毎の送信角 (-180deg - 180deg)
    Theta_r = (rand(1,nu)-0.5)*360; % ユーザ毎の受信角 (-180deg - 180deg)
    a_t = exp(-1j*2*pi*d_t*(0:nt-1).'*sin(Theta_t*derad));  % 送信モードベクトル
    
      for n = 1 : nu
        a_r = exp(-1j*2*pi*d_r*(0:2-1).'*sin(Theta_r(1,n)*derad)); % ユーザ毎の受信モードベクトル
        H_los((n-1)*2+1:(n-1)*2+2,:) = a_r*a_t.';                % ユーザ毎のLOSチャネル行列
      end

    % 伝搬チャネル行列=[sqrt(K/(K+1))*(LOS チャネル)]...
    %                   .+[sqrt(1/(K+1))*(NLOS チャネル)]
    H = sqrt(K/(K+1))*H_los + sqrt(1/(K+1))*H_iid;
Hu = peruser(H,nu);
%Hu(アンテナ選択後のユーザkのチャネル行列)
S = zeros(1,1);
No = zeros(1,1);

for i = 1:nu
    h = Hu(:,:,i)'*Hu(:,:,i);
    for ii = 1:nt
        S(ii,1) = svd(h(ii,:));
    end
   
    [~,No(1,i)] = max(S);
    beam(i,:) = Hu(:,No(1,i),i)';
end

%Hu_beam
for i = 1:nu
    Hu_beam(:,:,i) = beam(i,:)*Hu(:,:,i);
end

%He(Hからユーザkのチャネル行列を除いた行列)
for i = 1:nu
    e = 1:nu;    e(:,i) = [];
    He(:,:,i) = alluser(Hu_beam(:,:,e));
end

%Ven(Heの雑音部分空間に対応する固有ベクトル)
for i = 1:nu
    [~,~,Ve] = svd(He(:,:,i));
    Ven(:,:,i) = Ve(:,nu:nt);
end

%Vts(Hu*Venの信号部分空間に対応する固有ベクトル)
%Wtu(ユーザkに対応する基地局側ウエイト)
%St(特異値を格納 1×nu)
%Wr(ユーザ側ウエイト)
for i = 1:nu
    [Ut,~,Vt] = svd(Hu_beam(:,:,i)*Ven(:,:,i));
    St(:,i) = svd(Hu_beam(:,:,i)*Ven(:,:,i));
    Wr(:,:,i) = Ut';
    Vts(:,:,i) = Vt(:,1);
    Wtu(:,:,i) = Ven(:,:,i)*Vts(:,:,i);
end

%amp(送信電力)
amp = 10^(snr/10)/nt;

%data(送信データ)
%data_bd(BDの送信データ)
data_bd = randi([0 2^bsr-1],nu,ndata);
data = randi([0 2^(bsr/2)-1],2*nu,ndata);
datau = peruser(data,nu);

%s(送信ストリーム)
%s_bd(BD法の送信ストリーム)
s_bd = Mapping(data_bd,bsr);
s = Mapping(data,bsr/2);

% Wt = pinv(H)/norm(pinv(H),'fro');
W_bd = reshape(Wtu,[nt,nu]);
W_mmse = H'*pinv(H*H' + nt*4*eye(nt)/10^(snr/10))/norm(H'*pinv(H*H' + nt*4*eye(nt)/10^(snr/10)),'fro');
W_zf = pinv(H)/norm(pinv(H),'fro') ;  

%n_bd(BDの熱雑音)
%x(受信信号)
%yu(ユーザkの入力信号)
n = (randn(2*nu,ndata) + 1i*randn(2*nu,ndata))/sqrt(2);
x_bd = H*W_bd*s_bd*sqrt(amp) + n;
x_zf = H*W_zf*s*sqrt(amp) + n;
x_mmse = H*W_mmse*s*sqrt(amp) + n;
for i = 1:nu
    yu(:,:,i) = Wr(:,:,i)*beam(i,:)*x_bd(i*2-1:i*2,:)/sqrt(amp)/St(1,i);
end

%y(入力信号)
%Y(復調信号)
y_bd = alluser(yu);
Y_bd = Decode(y_bd,bsr);
y_zf = alluser(x_zf);
Y_zf = Decode(y_zf,bsr/2);
Yu_zf = peruser(Y_zf,nu);
y_mmse = alluser(x_mmse);
Y_mmse = Decode(y_mmse,bsr/2);
Yu_mmse = peruser(Y_mmse,nu);

%ber
Ber_zf = zeros(1,1);
Ber_mmse = zeros(1,1);
for i = 1:nu
    Ber_zf(1,i) = BER2(BER1(Yu_zf(1,:,i),datau(1,:,i),bsr/2),bsr) + BER2(BER1(Yu_zf(2,:,i),datau(2,:,i),bsr/2),bsr);
    Ber_mmse(1,i) = BER2(BER1(Yu_mmse(1,:,i),datau(1,:,i),bsr/2),bsr) + BER2(BER1(Yu_mmse(2,:,i),datau(2,:,i),bsr/2),bsr);
end

ber_zf(simu,snr+1) = mean(Ber_zf);
ber_mmse(simu,snr+1) = mean(Ber_mmse);
ber_bd(simu,snr+1) = BER2(BER1(Y_bd,data_bd,bsr),bsr);
end
snr


            
%BER(試行回数平均)
BER_zf = mean(ber_zf);
BER_mmse = mean(ber_mmse);
BER_bd = mean(ber_bd);
end
%描写
figure
semilogy(0:SNR,BER_zf(1:SNR+1),'b'); %青→BER(ZF)
hold on
set(gca,'FontName','Times New Roman','FontSize',15);
xlabel('SNR [dB]','fontsize',20);
ylabel('BER [dB]','fontsize',20);
semilogy(0:SNR,BER_mmse(1:SNR+1),'r'); %赤→BER(MMSE)
hold on
semilogy(0:SNR,BER_bd(1:SNR+1),'k'); %黒→BER(BD)
legend('ZF','MMSE','BD');
ylim([10^(-5) 1]);
grid on;
hold off;