% GMMSE-CI Method1
% 
% Condition : NT >= NR*NU

function [W,UTT,STT,RIP,SP] = gmmse_m1(NT,NR,NU,H,a)
W=zeros(NT,NR,NU);
% W2=[];
% ne = NR*(NU-1)+1:NT;

%Wzf = H\eye(NT,NT);
%Wmmse = H'/(H*H'+ a*eye(NT,NT)); % MMSE
Wmmse = H'/(H*H'+ a*eye(NR*NU,NR*NU)); % MMSE
%
HH = H'*H;

bb = 0;
for nuser=1:NU
    % nuserにおける受信アンテナ番号
    ns = NR*(nuser-1)+1:NR*nuser; 

    % nuserのウエイトをQR分解し，Qを抽出
    Wopt = Wmmse(:,ns);
    [Q,~] = qr(Wopt);
    QJ(:,:,nuser)=Q(:,1:NR);
    
    % 各ユーザのT行列の計算
    Hu=H(ns,:);     % nuserのチャネル行列
    T(:,:,nuser)=(QJ(:,:,nuser)'*HH*QJ(:,:,nuser)+a*eye(NR))\...
        (QJ(:,:,nuser)'*(Hu'*Hu)*QJ(:,:,nuser));
    bb = bb + trace(T(:,:,nuser)'*T(:,:,nuser));
end
%
beta = sqrt(NT)/sqrt(bb);
T = beta*T;
%
for nuser=1:NU
    ns = NR*(nuser-1)+1:NR*nuser; 
    Hu=H(ns,:);
    Wo=QJ(:,:,nuser)*T(:,:,nuser);
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
