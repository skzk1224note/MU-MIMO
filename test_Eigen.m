
Nu=8;
NR=2;
NT=16;
SNR_tar=20;
D2=zeros(NR,NR,Nu);
H = (randn(NR*Nu,NT) + 1j*randn(NR*Nu,NT))/sqrt(2);
    Hu = peruser(H,Nu);
    sigma2 = 1/(10^(SNR_tar/10)); % noise power
a = sigma2*NT;
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
    D3=sqrt(D1) ;
    %D12=D1(:,1)/D1(:,2)
    SNR1=D1(:,1)./a;
    SNR2=D1(:,2)./a;
    EWs = EW(:,IN);
   D2(:,:,nuser) = diag(D3(:,1:NR))
    %% Eigen SNR2が使えない状態であると判断されるdB以下で第一固有値の道へ変更
    if SNR2<20
       % D1(:,NR) < 0.1
        EWs(:,NR) = EWs(:,NR-1);
    end
    %BFへの切り替え
    if SNR2<50
    Wopt(:,:,nuser) = EWs(:,1:NR);
    else
    Wopt(:,:,nuser) = EWs(:,1:NR)*(D2(:,:,nuser));
    end
    Wopt(:,:,nuser) = Wopt(:,:,nuser)/norm(Wopt(:,:,nuser),'fro')*sqrt(NR);
end
