% bmsn_as.m 
% B-MSN algorithm for NU-user
% Condition : NT >= NR*NU
% Antenna Selection : NRs

function [W,UTT,STT,RIP,SP] = bmsn_bf_as(NT,NR,NRs,NU,H,a,T)
% Example
% NR = 2;
% NRs = 1;
W=zeros(NT,NRs,NU);
% W2=[];
ne = NRs*(NU-1)+1:NT;
%a: 擬似雑音

Hu_as = zeros(NRs,NT,NU);
Tas = zeros(NRs,NRs,NU);
for inu = 1:NU
    h = H(1+(inu-1)*NR:NR+(inu-1)*NR,:);
    h_norm=sqrt(sum(abs(h).^2,2));
    [~,Ind]=sort(h_norm,'descend');
%     [~,No] = max(h_norm);
    No = Ind(1:NRs);
    Hu_as(:,:,inu) = h(No,:);
    Tas(:,:,inu) = T(No,No,inu);
end
%H_as(伝搬チャネル行列)
H_as = alluser(Hu_as); %伝搬チャネル行列

for nuser=1:NU
        
    % nuserにおける受信アンテナ番号
    ns = NRs*(nuser-1)+1:NRs*nuser; 
    % 全ユーザを統合したチャネル行列 (NU*NR) * NT   
    HT=H_as;
    % nuserのチャネル行列を抜き取り
    Hu=HT(ns,:);
    HT(ns,:)=[];  
    % nuser以外のチャネル行列をSVD
    Wopt(:,:,nuser) = (HT'*HT+a*eye(size(HT,2)))\Hu'*Tas(:,:,nuser);
    Wopt(:,:,nuser) = Wopt(:,:,nuser)/norm(Wopt(:,:,nuser),'fro')*sqrt(NRs);
    % 雑音部分空間ベクトルをnuerのチャネル行列に乗算                                 
    HTT=Hu*Wopt(:,:,nuser);      
    % 変換行列をSVD
    [UTT(:,:,nuser),STT(:,:,nuser),VTT(:,:,nuser)]=svd(HTT,'econ'); 
    % nuserのウエイト(信号部分空間を利用)                                 
    W(:,:,nuser)=Wopt(:,:,nuser)*VTT(:,1:NRs,nuser);
        
    % 確認用：現在はコメント
    %UTT(:,:,nuser)'*H(ns,:)*W(:,:,nuser)
    %STT(:,:,nuser)
    %W2=[W2,W(:,:,nuser)];
    
end

% 所望波＆干渉波電力の計算
SP = zeros(NRs,NU);
RIP = zeros(NRs,NU);
for nuser=1:NU
    % nuserにおける受信アンテナ番号
    ns = NRs*(nuser-1)+1:NRs*nuser; 
        
    nuser2=1:NU;
    nuser2(nuser)=[];
    YI = zeros(NRs,NRs);
    for nn=nuser2
        YI=YI+UTT(:,1:NRs,nuser)'*H_as(ns,:)*W(:,:,nn);
    end
    RIP(:,nuser) = sum(abs(YI).^2,2); % 干渉波電力
    
    YS=UTT(:,1:NRs,nuser)'*H_as(ns,:)*W(:,:,nuser);
    SP(:,nuser) = sum(abs(YS).^2,2); % 所望波電力
    
end
    
    % 確認用：現在はコメント
    %UTT(:,:,nuser)'*H(ns,:)*W(:,:,nuser)
    %STT(:,:,nuser)
    %W2=[W2,W(:,:,nuser)];
    
end
% 確認用：現在はコメント
% abs(H*W2)