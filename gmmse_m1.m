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
    % nuser�ɂ������M�A���e�i�ԍ�
    ns = NR*(nuser-1)+1:NR*nuser; 

    % nuser�̃E�G�C�g��QR�������CQ�𒊏o
    Wopt = Wmmse(:,ns);
    [Q,~] = qr(Wopt);
    QJ(:,:,nuser)=Q(:,1:NR);
    
    % �e���[�U��T�s��̌v�Z
    Hu=H(ns,:);     % nuser�̃`���l���s��
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
