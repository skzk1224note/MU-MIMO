% bmsn.m
% B-MSN algorithm for NU-user
% Condition : NT >= NR*NU

function [W,UTT,STT,RIP,SP] = bmsn_bf_pseudo_noise(NT,NR,NU,H,T)
W=zeros(NT,NR,NU);
% W2=[];
% ne = NR*(NU-1)+1:NT;
%a = 1e-2; % 擬似雑音

%T 所望のチャネル行列

if nargin < 6
    for nuser = 1:NU
        T(:,:,nuser) = eye(NR);
    end
end

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
    
    % nuser以外のチャネル行列をSVD
    Wopt(:,:,nuser) = (HT'*HT+a*eye(size(HT,2)))\Hu'*T(:,:,nuser);
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