% bd.m
% BD algorithm for NU-user�iBD�@�j
% Condition : NT >= NR*NU

function [W,UTT,STT,RIP,SP] = bd(NT,NR,NU,H)
W=zeros(NT,NR,NU);
% W2=[];
ne = NR*(NU-1)+1:NT;
for nuser=1:NU
    % nuser�ɂ������M�A���e�i�ԍ�
    ns = NR*(nuser-1)+1:NR*nuser; 
    % �S���[�U�𓝍������`���l���s�� (NU*NR) * NT   
    HT=H;
    % nuser�̃`���l���s��𔲂����i��s���}������nuser�̃`���l���s����폜�j                        
    HT(ns,:)=[];   
    % nuser�ȊO�̃`���l���s���SVD�i���ْl�����j                  
    %[UT(:,:,nuser),ST(:,:,nuser),VT(:,:,nuser)]=svd(HT); 
    [~,~,VT(:,:,nuser)]=svd(HT); % ~ = �j�������o��
    % �G��������ԃx�N�g����nuer�̃`���l���s��ɏ�Z                                 
    HTT=H(ns,:)*VT(:,ne,nuser);      
    % �ϊ��s���SVD
    [UTT(:,:,nuser),STT(:,:,nuser),VTT(:,:,nuser)]=svd(HTT); 
    % nuser�̃E�G�C�g(�M��������Ԃ𗘗p)                                 
    W(:,:,nuser)=VT(:,ne,nuser)*VTT(:,1:NR,nuser); 
                                     
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
        YI=YI+UTT(:,1:NR,nuser)'*H(ns,:)*W(:,:,nn); % ' = ���f����]�u�@.' = �]�u
    end
    RIP(:,nuser) = sum(abs(YI).^2,2); % ���g�d��
    % abs(a)�͔z��a�̊e�v�f�̐�Βl��Ԃ�
    YS=UTT(:,1:NR,nuser)'*H(ns,:)*W(:,:,nuser);
    SP(:,nuser) = sum(abs(YS).^2,2); % ���]�g�d��
    
end

% �m�F�p�F���݂̓R�����g
% abs(H*W2)