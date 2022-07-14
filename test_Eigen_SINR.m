Nu=8;
%NR=2;
NR=3;
NT=24;
SINU=1000;%試行回数1000回
SNR_min=0;
SNR_max=10;
SNR=(SNR_min:1:SNR_max).';
LSNR=length(SNR);
SNR1=zeros(Nu,LSNR);
SNR2=zeros(Nu,LSNR);
SNR3=zeros(Nu,LSNR);
SNR1ave=zeros(SINU,LSNR);
SNR2ave=zeros(SINU,LSNR);
SNR3ave=zeros(SINU,LSNR);
%SNR_tar=
    
   
for isinu=1:SINU
H = (randn(NR*Nu,NT) + 1j*randn(NR*Nu,NT))/sqrt(2);
%H = (randn(SINU,NR*Nu,NT) + 1j*randn(SINU,NR*Nu,NT))/sqrt(2);
    Hu = peruser(H,Nu);

for isnr = 1:LSNR
  
    SNR_tar = SNR(isnr);
    sigma2 = 1/(10^(SNR_tar/10)); % noise power
    a = sigma2*NT;
    %D12=zeros(1,isnr);

for nuser=1:Nu
    ns = NR*(nuser-1)+1:NR*nuser; 
    % 全ユーザを統合したチャネル行列 (NU*NR) * NT   
    HT=H;
    
    % HT行列の内のns行を全て抽出
    Hu=HT(ns,:);
    % nuserのチャネル行列を抜き取り
    HT(ns,:)=[];
    A = HT'*HT+a*eye(size(HT,2));
    B = Hu'*Hu;
    %固有ベクトルと固有値の行列を返す
    [EW,D] = eig(B,A);
    % sort(diag(abs(D)).','descend')=Db
    [D1,IN] = sort(diag(abs(D)).','descend');    %~ = D1（下のセクション）
    %D12(nuser,isnr)=D1(:,1)/D1(:,2);%固有値比を行列の要素に入れ込む
    %D13(nuser,isnr)=D1(:,1)/D1(:,3);
    %D23(nuser,isnr)=D1(:,2)/D1(:,3);
    %固有値のSNRを計算する

    SNR1(nuser,isnr)=D1(:,1)./a;
    SNR2(nuser,isnr)=D1(:,2)./a;
    SNR3(nuser,isnr)=D1(:,3)./a;
     
end
end
%D12ave=mean(D12);
%D12ave(isinu,:)=mean(D12);%それぞれのSNRに対しての、固有値比の平均をとる
%D13ave(isinu,:)=mean(D13);
%D23ave(isinu,:)=mean(D23);
SNR1ave(isinu,:)=mean(SNR1);
SNR2ave(isinu,:)=mean(SNR2);
SNR3ave(isinu,:)=mean(SNR3);
end
%D12con=mean(D12ave);
%D13con=mean(D13ave);
%D23con=mean(D23ave);%試行回数1000回の平均を取る
SNR1con=mean(SNR1ave);
SNR2con=mean(SNR2ave);
SNR3con=mean(SNR3ave);
%plot(SNR,D12con,'-o');
%plot(SNR,D12con,'-o',SNR,D13con,'-v',SNR,D23con,'-^');
%legend('lamda1to2','lamda1to3','lamda2to3','Location','Northwest');
%xlabel('SNR [dB]','Fontsize',16,'Fontname','Times New Roman');
%ylabel('lamda','Fontsize',16,'Fontname','Times New Roman');
plot(SNR,SNR1con,'-o',SNR,SNR2con,'-v',SNR,SNR3con,'-x');
legend('Eigen1 SNR','Eigen2 SNR','Eigen3 SNR','Location','Northwest');
xlabel('SNR [dB]','Fontsize',16,'Fontname','Times New Roman');
ylabel('Eigen SNR[dB]','Fontsize',16,'Fontname','Times New Roman');

