% regularized bd.m
% BD algorithm for NU-user
% Condition : NT >= NR*NU

function [W,UTT,STT,RIP,SP] = rbd(NT,NR,NU,H,a)
W=zeros(NT,NR,NU);
% W2=[];
%ne = NR*(NU-1)+1:NT;
for nuser=1:NU
    % nuser�ɂ������M�A���e�i�ԍ�
    ns = NR*(nuser-1)+1:NR*nuser; 
    % �S���[�U�𓝍������`���l���s�� (NU*NR) * NT   
    HT=H;
    % nuser�̃`���l���s��𔲂����                        
    HT(ns,:)=[];  
    % nuser�ȊO�̃`���l���s���SVD                  
    %[UT(:,:,nuser),ST(:,:,nuser),VT(:,:,nuser)]=svd(HT); 
    [~,ST(:,:,nuser),VT(:,:,nuser)]=svd(HT); 
    % �G��������ԃx�N�g����nuer�̃`���l���s��ɏ�Z
    D(:,:,nuser)=sqrt(eye(NT)/(ST(:,:,nuser).'*ST(:,:,nuser)+a*eye(NT)));
    F(:,:,nuser)=VT(:,:,nuser)*D(:,:,nuser);
    F(:,:,nuser) = F(:,:,nuser)/norm(F(:,:,nuser),'fro')*sqrt(NR);
    HTT=H(ns,:)*F(:,:,nuser);      
    % �ϊ��s���SVD
    [UTT(:,:,nuser),STT(:,:,nuser),VTT(:,:,nuser)]=svd(HTT); 
    % nuser�̃E�G�C�g(�M��������Ԃ𗘗p)                                 
    W(:,:,nuser)=F(:,:,nuser)*VTT(:,1:NR,nuser); 
                                     
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