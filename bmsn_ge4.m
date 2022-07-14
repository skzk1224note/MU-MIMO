% bmsn3_gev_SINR.m
% BMSN(GEV) algorithm for NU-user
% Condition : NT >= NR*NU

function [W,UTT,STT,RIP,SP] = bmsn_ge4(NT,NR,NU,H,a)
Wopt=zeros(NT,NR,NU);W=zeros(NT,NR,NU);
UTT=zeros(NR,NR,NU);STT=zeros(NR,NR,NU);VTT=zeros(NR,NR,NU);

D2=zeros(NR,NR,NU);
% ne = NR*(NU-1)+1:NT;

for nuser=1:NU
    % nuserにおける受信アンテナ番号
    ns = NR*(nuser-1)+1:NR*nuser; 
    % 全ユーザを統合したチャネル行列 (NU*NR) * NT   
    HT=H;
    % nuserのチャネル行列を抜き取り
    Hu=HT(ns,:);
    HT(ns,:)=[];
    % BMSN with GEV
    A = HT'*HT+a*eye(size(HT,2));
    B = Hu'*Hu;
    [EW,D] = eig(B,A);
    [D1,IN] = sort(diag(abs(D)).','descend');
    EW = EW(:,IN);
%     EW(:,NR) = EW(:,1);
%     if D1(:,NR) < 0.1 % lambda2が小さかったら0にする
%         D1(:,NR) = 0;
%         EW(:,NR) = EW(:,1);
%     end
    D2(:,:,nuser) = diag(D1(:,1:NR)); 
    D3(:,:,nuser)=sqrt(log2(D2(:,:,nuser)+1));%　重み付け3　
    Wopt(:,:,nuser)=EW(:,1:NR)*D3(:,:,nuser); % Weighting by sqrt of SINR　　（※ここを変える）
    Wopt(:,:,nuser) = Wopt(:,:,nuser)/norm(Wopt(:,:,nuser),'fro')*sqrt(NR); % Normalization

    % Woptをnuerのチャネル行列に乗算 -> Block channel matrix: HTT                                 
    HTT=Hu*Wopt(:,:,nuser);      
    % Block channel matrixをSVD
    [UTT(:,:,nuser),STT(:,:,nuser),VTT(:,:,nuser)]=svd(HTT);
    % nuserのウエイト(信号部分空間を利用)                                 
    W(:,:,nuser)=Wopt(:,:,nuser)*VTT(:,1:NR,nuser);
        
    % 確認用：現在はコメント
    %UTT(:,:,nuser)'*H(ns,:)*W(:,:,nuser)
    %STT(:,:,nuser)
    %W2=[W2,W(:,:,nuser)];
%     NormU4 = norm(UTT(:,:,nuser),'fro')
%     NormV4 = norm(VTT(:,:,nuser),'fro')
end
%% 固有モード伝送後のウエイト確認
% D2_bmsnge4 = D2
% 
% UTT_bmsnge4 = UTT
% STT_bmsnge4 = STT
% VTT_bmsnge4 = VTT

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
    for nn=nuser2
        YI=YI+UTT(:,1:NR,nuser)'*H(ns,:)*W(:,:,nn);
    end
    RIP(:,nuser) = sum(abs(YI).^2,2); % 干渉波電力
    
    YS=UTT(:,1:NR,nuser)'*H(ns,:)*W(:,:,nuser);
    SP(:,nuser) = sum(abs(YS).^2,2); % 所望波電力
    
end

% 確認用：現在はコメント
% abs(H*W2)