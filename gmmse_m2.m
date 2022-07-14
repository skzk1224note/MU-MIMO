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
    % nuser�ɂ������M�A���e�i�ԍ�
    ns = NR*(nuser-1)+1:NR*nuser; 

    % nuser�̃`���l���s��𔲂���� 
    HT=H;
    HT(ns,:)=[];  
    HH = HT'*HT;
    
    % nuser�̃E�G�C�g��QR�������CQ�𒊏o
    Wopt = Wmmse(:,ns);
    [Q,~] = qr(Wopt);
    QJ=Q(:,1:NR);
    
    % �e���[�U��T�s��̌v�Z
    LL=QJ'*HH*QJ+a*eye(NR);
    L=chol((LL+LL')/2);
    T=eye(NR)/L;
    gamma = sqrt(NR/trace(T'*T));
%    T=gamma*T/sqrt(NR*NU);
    T=gamma*T;
    
    Wo=QJ*T;
    Hu=H(ns,:);
    HTT=Hu*Wo;
    % �ϊ��s���SVD
    [UTT(:,:,nuser),STT(:,:,nuser),VTT(:,:,nuser)]=svd(HTT);
    % nuser�̃E�G�C�g(�M��������Ԃ𗘗p)                                 
    W(:,:,nuser)=Wo*VTT(:,1:NR,nuser);     
end

% ���]�g�����g�d�͂̌v�Z
SP = zeros(NR,NU);
RIP = zeros(NR,NU);
for nuser=1:NU
    % nuser�ɂ������M�A���e�i�ԍ�
    ns = NR*(nuser-1)+1:NR*nuser; 
        
    nuser2=1:NU;
    nuser2(nuser)=[];
    YI = zeros(NR,NR);
    for nn=nuser2
        YI=YI+UTT(:,1:NR,nuser)'*H(ns,:)*W(:,:,nn);
    end
    RIP(:,nuser) = sum(abs(YI).^2,2); % ���g�d��
    
    YS=UTT(:,1:NR,nuser)'*H(ns,:)*W(:,:,nuser);
    SP(:,nuser) = sum(abs(YS).^2,2); % ���]�g�d��
    
end
