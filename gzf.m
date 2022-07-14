% zf.m
% ZF+(EM-BF) algorithm for NU-user
% Condition : NT >= NR*NU

function [W,UTT,STT,RIP,SP] = gzf(NT,NR,NU,H)
W=zeros(NT,NR,NU);
% W2=[];
% ne = NR*(NU-1)+1:NT;

Wzf = H'/(H*H');

for nuser=1:NU
    % nuserにおける受信アンテナ番号
    ns = NR*(nuser-1)+1:NR*nuser; 

    % nuserのチャネル行列を抜き取り
    Hu=H(ns,:);

    % nuserのウエイト（列ベクトル毎に正規化）
    Wopt = Wzf(:,ns);
%     for ij = 1:NR
%         Wopt(:,ij,nuser) = Wopt(:,ij,nuser)/sqrt(Wopt(:,ij,nuser)'*Wopt(:,ij,nuser));
%     end

    [Q,~]=qr(Wopt);
    QJ(:,:,nuser)=Q(:,1:NR);
%     for ins=1:NU
%         ins
%         inss = [NR*(ins-1)+1:NR*ins]; 
%         H(inss,:)*Q
%     end
%     for ij = 1:NR
%         Q(:,ij) = Q(:,ij)/sqrt(Q(:,ij)'*Q(:,ij));
%     end
    % 雑音部分空間ベクトルをnuerのチャネル行列に乗算                                 
    HTT=Hu*QJ(:,:,nuser);
    % 変換行列をSVD
    [UTT(:,:,nuser),STT(:,:,nuser),VTT(:,:,nuser)]=svd(HTT);
    % nuserのウエイト(信号部分空間を利用)                                 
    W(:,:,nuser)=QJ(:,:,nuser)*VTT(:,1:NR,nuser); 
               
    % 確認用：現在はコメント
    %UTT(:,:,nuser)'*H(ns,:)*W(:,:,nuser)
    %STT(:,:,nuser)
    %W2=[W2,W(:,:,nuser)];
    
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
% 確認用：現在はコメント
% abs(H*W2)