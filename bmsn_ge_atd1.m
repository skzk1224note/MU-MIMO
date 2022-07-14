% bmsn3_gev_SINR.m
% BMSN(GEV) algorithm for NU-user
% Condition : NT >= NR*NU

function [W,UTT,STT,RIP,SP] = bmsn_ge_atd1(NT,NR,NU,H,a)
Wopt=zeros(NT,NR,NU);W=zeros(NT,NR,NU);
UTT=zeros(NR,NR,NU);STT=zeros(NR,NR,NU);VTT=zeros(NR,NR,NU);
% ne = NR*(NU-1)+1:NT;

for nuser=1:NU
    % nuserにおける受信アンテナ番号
    ns = NR*(nuser-1)+1:NR*nuser; 
    % 全ユーザを統合したチャネル行列 (NU*NR) * NT   
    HT=H;
    
    % HT行列の内のns行を全て抽出
    Hu=HT(ns,:);
    % nuserのチャネル行列を抜き取り
    HT(ns,:)=[];
    
    % BMSN with GEV
    %一般固有値展開
    A = HT'*HT+a*eye(size(HT,2));
    B = Hu'*Hu;
    %固有ベクトルと固有値の行列を返す
    [EW,D] = eig(B,A);
    % sort(diag(abs(D)).','descend')=Db
    [D1,IN] = sort(diag(abs(D)).','descend');    %~ = D1（下のセクション）
    %largeD=diag(abs(D));
    %largeD()
    D_bf=squatm(D1);
    %固有値1と2の比率
    D12=D1(:,1)./D1(:,2);
    %sigma2=a/NT;
    %列のIN個をそれぞれ縦に取り出す
    %固有値の大きい順に対応する固有ベクトルを並べる
    EW = EW(:,IN);
    %% lambda2が小さかったらlambda1とし，lambda1に対応する固有ベクトルのみをふたつ使用
    if  NR>1
    %if D12 > 3
            EW(:,NR) = EW(:,NR-1);
    %end
%NRとNR-1の通路での固有ベクトルを同じにする
        %EW(:,NR) = EW(:,NR-1);
    end
    %%
    %一般固有値問題の解を送信ウエイトにする
    Wopt(:,:,nuser) = EW(:,1:NR);
    %GEからBFへの切り替え条件重み付け
    
    Wopt(:,:,nuser) = Wopt(:,:,nuser)/norm(Wopt(:,:,nuser),'fro')*sqrt(NR); % Normalization
    
    %固有モード伝送
    % Woptをnuerのチャネル行列に乗算 -> Block channel matrix: HTT                                 
    HTT=Hu*Wopt(:,:,nuser);
    
    % HTT3 = HTT
    % Block channel matrixをSVD
    [UTT(:,:,nuser),STT(:,:,nuser),VTT(:,:,nuser)]=svd(HTT);
    % nuserのウエイト(信号部分空間を利用)
    if NR > 1
        W(:,:,nuser)=Wopt(:,:,nuser)*VTT(:,1:NR,nuser); %lambda1に対応する固有ベクトルのみをふたつ使用
    end
    % 確認用：現在はコメント
    %UTT(:,:,nuser)'*H(ns,:)*W(:,:,nuser)
    %STT(:,:,nuser)
    %W2=[W2,W(:,:,nuser)];
%     NormU3 = norm(UTT(:,:,nuser),'fro')
%     NormV3 = norm(VTT(:,:,nuser),'fro')
%     NormW3 = norm(Wopt(:,:,nuser),'fro')
end
%% 固有モード伝送後のウエイト確認
%D2_bmsnge3 = D2
% HTT3 = HTT
% Wopt3 = Wopt
%  UTT_bmsnge3 = UTT
% STT_bmsnge3 = STT
%VTT_bmsnge3 = VTT

%%
% 所望波＆干渉波電力の計算
SP = zeros(NR,NU);
RIP = zeros(NR,NU);
for nuser=1:NU
    % nuserにおける受信アンテナ番号
    ns = NR*(nuser-1)+1:NR*nuser; 
    nuser2=1:NU;
    nuser2(nuser)=[];
    YI = zeros(NR,NR);
    %1からNUまで繰り返し
    for  nn=nuser2
        YI=YI+UTT(:,1:NR,nuser)'*H(ns,:)*W(:,:,nn);
   
    end
    
    
   
    %sumでYIのそれぞれの要素の絶対値の各行の要素を足し合わせる
    RIP(:,nuser) = sum(abs(YI).^2,2); % 干渉波電力
   
    YS=UTT(:,1:NR,nuser)'*H(ns,:)*W(:,:,nuser);
    
    SP(:,nuser) = sum(abs(YS).^2,2); % 所望波電力
    
end
% YI3 = YI
% YS3 = YS
% RIP3 = RIP
% 確認用：現在はコメント
% abs(H*W2)