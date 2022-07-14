% regularized bd.m
% BD algorithm for NU-user
% Condition : NT >= NR*NU

function [W,UTT,STT,RIP,SP] = rbd(NT,NR,NU,H,a)
W=zeros(NT,NR,NU);
% W2=[];
%ne = NR*(NU-1)+1:NT;
for nuser=1:NU
    % nuserにおける受信アンテナ番号
    ns = NR*(nuser-1)+1:NR*nuser; 
    % 全ユーザを統合したチャネル行列 (NU*NR) * NT   
    HT=H;
    % nuserのチャネル行列を抜き取り                        
    HT(ns,:)=[];  
    % nuser以外のチャネル行列をSVD                  
    %[UT(:,:,nuser),ST(:,:,nuser),VT(:,:,nuser)]=svd(HT); 
    [~,ST(:,:,nuser),VT(:,:,nuser)]=svd(HT); 
    % 雑音部分空間ベクトルをnuerのチャネル行列に乗算
    D(:,:,nuser)=sqrt(eye(NT)/(ST(:,:,nuser).'*ST(:,:,nuser)+a*eye(NT)));
    F(:,:,nuser)=VT(:,:,nuser)*D(:,:,nuser);
    F(:,:,nuser) = F(:,:,nuser)/norm(F(:,:,nuser),'fro')*sqrt(NR);
    HTT=H(ns,:)*F(:,:,nuser);      
    % 変換行列をSVD
    [UTT(:,:,nuser),STT(:,:,nuser),VTT(:,:,nuser)]=svd(HTT); 
    % nuserのウエイト(信号部分空間を利用)                                 
    W(:,:,nuser)=F(:,:,nuser)*VTT(:,1:NR,nuser); 
                                     
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