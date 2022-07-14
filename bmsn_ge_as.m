% bmsn_as.m 
% B-MSN algorithm for NU-user
% Condition : NT >= NR*NU
% Antenna Selection : NRs

function [W,UTT,STT,RIP,SP] = bmsn_ge_as(NT,NR,NRs,NU,H,a)
% Example
% NR = 2;
% NRs = 1;
W=zeros(NT,NRs,NU);
% W2=[];
ne = NRs*(NU-1)+1:NT;
%a: �[���G��

Hu_as = zeros(NRs,NT,NU);
for inu = 1:NU
    h = H(1+(inu-1)*NR:NR+(inu-1)*NR,:);
    h_norm=sqrt(sum(abs(h).^2,2));
    [~,Ind]=sort(h_norm,'descend');
%     [~,No] = max(h_norm);
    No = Ind(1:NRs);
    Hu_as(:,:,inu) = h(No,:);
end
%H_as(�`���`���l���s��)
H_as = alluser(Hu_as); %�`���`���l���s��

for nuser=1:NU
        
    % nuser�ɂ������M�A���e�i�ԍ�
    ns = NRs*(nuser-1)+1:NRs*nuser; 
    % �S���[�U�𓝍������`���l���s�� (NU*NR) * NT   
    HT=H_as;
    % nuser�̃`���l���s��𔲂����
    Hu=HT(ns,:);
    HT(ns,:)=[];  
    % BMSN with GEV
    A = HT'*HT+a*eye(size(HT,2));
    B = Hu'*Hu;
    [EW,D] = eig(B,A);
    [~,IN] = sort(diag(abs(D)).','descend');
    %D2(:,:,nuser) = D1.';
    EW = EW(:,IN);
    Wopt(:,:,nuser)=EW(:,1:NRs);
    Wopt(:,:,nuser) = Wopt(:,:,nuser)/norm(Wopt(:,:,nuser),'fro')*sqrt(NRs);
    %[~,~,VT(:,:,nuser)]=svd(HT); 
    % �G��������ԃx�N�g����nuer�̃`���l���s��ɏ�Z                                 
    HTT=Hu*Wopt(:,:,nuser);      
    % �ϊ��s���SVD
    [UTT(:,:,nuser),STT(:,:,nuser),VTT(:,:,nuser)]=svd(HTT,'econ'); 
    % nuser�̃E�G�C�g(�M��������Ԃ𗘗p)                                 
    W(:,:,nuser)=Wopt(:,:,nuser)*VTT(:,1:NRs,nuser);
        
    % �m�F�p�F���݂̓R�����g
    %UTT(:,:,nuser)'*H(ns,:)*W(:,:,nuser)
    %STT(:,:,nuser)
    %W2=[W2,W(:,:,nuser)];
    
end

% ���]�g�����g�d�͂̌v�Z
SP = zeros(NRs,NU);
RIP = zeros(NRs,NU);
for nuser=1:NU
    % nuser�ɂ������M�A���e�i�ԍ�
    ns = NRs*(nuser-1)+1:NRs*nuser; 
        
    nuser2=1:NU;
    nuser2(nuser)=[];
    YI = zeros(NRs,NRs);
    for nn=nuser2
        YI=YI+UTT(:,1:NRs,nuser)'*H_as(ns,:)*W(:,:,nn);
    end
    RIP(:,nuser) = sum(abs(YI).^2,2); % ���g�d��
    
    YS=UTT(:,1:NRs,nuser)'*H_as(ns,:)*W(:,:,nuser);
    SP(:,nuser) = sum(abs(YS).^2,2); % ���]�g�d��
    
end
    
    % �m�F�p�F���݂̓R�����g
    %UTT(:,:,nuser)'*H(ns,:)*W(:,:,nuser)
    %STT(:,:,nuser)
    %W2=[W2,W(:,:,nuser)];
    
end
% �m�F�p�F���݂̓R�����g
% abs(H*W2)