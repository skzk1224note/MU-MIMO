% bmsnb.m
% BMSN(GEV) algorithm for NU-user
% Condition : NT >= NR*NU

function [W,UTT,STT,RIP,SP] = bmsn_ge1(NT,NR,NU,H,a)
W=zeros(NT,NR,NU);
% W2=[];
% ne = NR*(NU-1)+1:NT;

for nuser=1:NU
    % nuser�ɂ������M�A���e�i�ԍ�
    ns = NR*(nuser-1)+1:NR*nuser; 
    % �S���[�U�𓝍������`���l���s�� (NU*NR) * NT   
    HT=H;
    % nuser�̃`���l���s��𔲂����
    Hu=HT(ns,:);
    HT(ns,:)=[];  
    % BMSN with GEV
    A = HT'*HT+a*eye(size(HT,2));
    B = Hu'*Hu;
    [EW,D] = eig(B,A); % [EW,D] = eig(A,B)�́A��ʉ��ŗL�l����Ȃ�Ίp�s��D�ƁA�Ή�����E�ŗL�x�N�g�����ɂ���X�p�[�X�s�� EW(�܂� A*EW = B*EW*D) ��Ԃ��܂��B
    [~,IN] = sort(diag(abs(D)).','descend'); % B = sort(A,'descend') ��A�̗���~���ɕ��בւ���B�@IN = �C���f�b�N�X�s��
    %D2(:,:,nuser) = D1.';
    EW = EW(:,IN); % �ŗL�l�ɑΉ������ŗL�x�N�g���ɕ��ёւ�
    Wopt(:,:,nuser)=EW(:,1:NR); % ��(3.10) �i�Q�l�F�ĒÏC�_�j
    Wopt(:,:,nuser) = Wopt(:,:,nuser)/norm(Wopt(:,:,nuser),'fro')*sqrt(NR); % �K�i��
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