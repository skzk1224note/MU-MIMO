% zf.m
% ZF+(EM-BF) algorithm for NU-user
% Condition : NT >= NR*NU

function [W,UTT,STT,RIP,SP] = gzf(NT,NR,NU,H)
W=zeros(NT,NR,NU);
% W2=[];
% ne = NR*(NU-1)+1:NT;

Wzf = H'/(H*H');

for nuser=1:NU
    % nuser�ɂ������M�A���e�i�ԍ�
    ns = NR*(nuser-1)+1:NR*nuser; 

    % nuser�̃`���l���s��𔲂����
    Hu=H(ns,:);

    % nuser�̃E�G�C�g�i��x�N�g�����ɐ��K���j
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
    % �G��������ԃx�N�g����nuer�̃`���l���s��ɏ�Z                                 
    HTT=Hu*QJ(:,:,nuser);
    % �ϊ��s���SVD
    [UTT(:,:,nuser),STT(:,:,nuser),VTT(:,:,nuser)]=svd(HTT);
    % nuser�̃E�G�C�g(�M��������Ԃ𗘗p)                                 
    W(:,:,nuser)=QJ(:,:,nuser)*VTT(:,1:NR,nuser); 
               
    % �m�F�p�F���݂̓R�����g
    %UTT(:,:,nuser)'*H(ns,:)*W(:,:,nuser)
    %STT(:,:,nuser)
    %W2=[W2,W(:,:,nuser)];
    
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
% �m�F�p�F���݂̓R�����g
% abs(H*W2)