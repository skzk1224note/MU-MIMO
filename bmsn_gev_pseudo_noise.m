% bmsnb.m
% BMSN(GEV) algorithm for NU-user
% Condition : NT >= NR*NU

function [W,UTT,STT,RIP,SP] = bmsn_gev_pseudo_noise(NT,NR,NU,H)
W=zeros(NT,NR,NU);
% W2=[];
% ne = NR*(NU-1)+1:NT;

for nuser=1:NU
    % nuserにおける受信アンテナ番号
    ns = NR*(nuser-1)+1:NR*nuser; 
    % 全ユーザを統合したチャネル行列 (NU*NR) * NT   
    HT=H;
    % nuserのチャネル行列を抜き取り
    Hu=HT(ns,:);
    HT(ns,:)=[];  
    % 擬似雑音をチャネルから演算
    a = sqrt(norm(HT,'fro'))/(NU-1);
    
    % BMSN with GEV
    A = HT'*HT+a*eye(size(HT,2)); % H(k)の項に当たる
    B = Hu'*Hu;                   % H(k)からユーザkを抜き取った項に当たる
    [EW,D] = eig(B,A);
    [~,IN] = sort(diag(abs(D)).','descend'); % 対角行列において固有値を昇順に並べる
    %D2(:,:,nuser) = D1.';
    EW = EW(:,IN);
    Wopt(:,:,nuser)=EW(:,1:NR);
    Wopt(:,:,nuser) = Wopt(:,:,nuser)/norm(Wopt(:,:,nuser),'fro')*sqrt(NR);
%     for ij = 1:NR
%         Wopt(:,ij,nuser) = Wopt(:,ij,nuser)/sqrt(Wopt(:,ij,nuser)'*Wopt(:,ij,nuser));
%     end
    %[~,~,VT(:,:,nuser)]=svd(HT); 
    % 雑音部分空間ベクトルをnuerのチャネル行列に乗算                                 
    HTT=Hu*Wopt(:,:,nuser);      
    % 変換行列をSVD
    [UTT(:,:,nuser),STT(:,:,nuser),VTT(:,:,nuser)]=svd(HTT); 
    % nuserのウエイト(信号部分空間を利用)                                 
    W(:,:,nuser)=Wopt(:,:,nuser)*VTT(:,1:NR,nuser);
        
    % 確認用：現在はコメント
    %UTT(:,:,nuser)'*H(ns,:)*W(:,:,nuser)
    %STT(:,:,nuser)
    %W2=[W2,W(:,:,nuser)];
    %VTTge = VTT
end

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
% YIge = YI
% RIPge = RIP
% 確認用：現在はコメント
% abs(H*W2)