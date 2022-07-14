% bmsn with beamforming condition (T=I)
% B-MSN algorithm for NU-user
% Condition : NT >= NR*NU

function [W,UTT,STT,RIP,SP] = bmsn_bf1(NT,NR,NU,H,a,T)
% W2=[];
% ne = NR*(NU-1)+1:NT;
%a = 1e-2; % �[���G��

%T ���]�̃`���l���s��

if nargin < 6
    for nuser = 1:NU
        T(:,:,nuser) = eye(NR);
    end
end

W=zeros(NT,NR,NU);

for nuser=1:NU
    % nuser�ɂ������M�A���e�i�ԍ�
    ns = NR*(nuser-1)+1:NR*nuser; 
    % �S���[�U�𓝍������`���l���s�� (NU*NR) * NT   
    HT=H;
    % nuser�̃`���l���s��𔲂����
    Hu=HT(ns,:);
    HT(ns,:)=[];  
    % nuser�ȊO�̃`���l���s���SVD
    Wopt(:,:,nuser) = (HT'*HT+a*eye(size(HT,2)))\Hu'*T(:,:,nuser);
    Wopt(:,:,nuser) = Wopt(:,:,nuser)/norm(Wopt(:,:,nuser),'fro')*sqrt(NR);
%     for ij = 1:NR
%         Wopt(:,ij,nuser) = Wopt(:,ij,nuser)/sqrt(Wopt(:,ij,nuser)'*Wopt(:,ij,nuser));
%     end
    %[~,~,VT(:,:,nuser)]=svd(HT); 
    % �G��������ԃx�N�g����nuer�̃`���l���s��ɏ�Z                                 
    HTT=Hu*Wopt(:,:,nuser);      
    % �ϊ��s���SVD
    [UTT(:,:,nuser),STT(:,:,nuser),VTT(:,:,nuser)]=svd(HTT); 
    % nuser�̃E�G�C�g(�M��������Ԃ𗘗p)                                 
    W(:,:,nuser)=Wopt(:,:,nuser)*VTT(:,1:NR,nuser);
        
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