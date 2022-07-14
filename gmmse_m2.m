% GMMSE-CI Method2
% 
% Condition : NT >= NR*NU

function [W,UTT,STT,RIP,SP] = gmmse_m2(NT,NR,NU,H,a)
W=zeros(NT,NR,NU);
% W2=[];
% ne = NR*(NU-1)+1:NT;

%Wzf = H\eye(NT,NT);
%Wmmse = H'/(H*H'+ a*eye(NT,NT)); % MMSE
Wmmse = H'/(H*H'+ a*eye(NR*NU,NR*NU)); % MMSE

for nuser=1:NU
    % nuserにおける受信アンテナ番号
    ns = NR*(nuser-1)+1:NR*nuser; 

    % nuserのチャネル行列を抜き取り 
    HT=H;
    HT(ns,:)=[];  
    HH = HT'*HT;
    
    % nuserのウエイトをQR分解し，Qを抽出
    Wopt = Wmmse(:,ns);
    [Q,~] = qr(Wopt);
    QJ=Q(:,1:NR);
    
    % 各ユーザのT行列の計算
    LL=QJ'*HH*QJ+a*eye(NR);
    L=chol((LL+LL')/2);
    T=eye(NR)/L;
    gamma = sqrt(NR/trace(T'*T));
%    T=gamma*T/sqrt(NR*NU);
    T=gamma*T;
    
    Wo=QJ*T;
    Hu=H(ns,:);
    HTT=Hu*Wo;
    % 変換行列をSVD
    [UTT(:,:,nuser),STT(:,:,nuser),VTT(:,:,nuser)]=svd(HTT);
    % nuserのウエイト(信号部分空間を利用)                                 
    W(:,:,nuser)=Wo*VTT(:,1:NR,nuser);     
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
