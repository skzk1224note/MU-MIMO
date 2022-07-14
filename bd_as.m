% bd_as.m
% BD algorithm for NU-user
% Condition : NT >= NR*NU
% Antenna Selection : NRs

function [W,UTT,STT,RIP,SP] = bd_as(NT,NR,NRs,NU,H)
%
% NR = 2;
% NRs = 1;
W=zeros(NT,NRs,NU);
% W2=[];
ne = NRs*(NU-1)+1:NT;

Hu_as = zeros(NRs,NT,NU);
for inu = 1:NU
    h = H(1+(inu-1)*NR:NR+(inu-1)*NR,:);
    h_norm = sqrt(sum(abs(h).^2,2));
    [~,No] = sort(h_norm,'descend');    % [~,No] = max(h_norm)
    h = h(No,:);
    Hu_as(:,:,inu) = h(1:NRs,:);   % Hu_as(:,:,inu) = h(No,:);
end

%H_as(�`���`���l���s��)
H_as = alluser(Hu_as); %�`���`���l���s��
    
for nuser=1:NU
    % nuser�ɂ������M�A���e�i�ԍ�
    ns = [NRs*(nuser-1)+1:NRs*nuser]; 
    % �S���[�U�𓝍������`���l���s�� (NU*NR) * NT   
    HT=H_as;
    % nuser�̃`���l���s��𔲂����                        
    HT(ns,:)=[];  
    % nuser�ȊO�̃`���l���s���SVD                  
    %[UT(:,:,nuser),ST(:,:,nuser),VT(:,:,nuser)]=svd(HT); 
    [~,~,VT(:,:,nuser)]=svd(HT); 
    % �G��������ԃx�N�g����nuer�̃`���l���s��ɏ�Z                                 
    HTT=H_as(ns,:)*VT(:,ne,nuser);      
    % �ϊ��s���SVD
    [UTT(:,:,nuser),STT(:,:,nuser),VTT(:,:,nuser)]=svd(HTT); 
    % nuser�̃E�G�C�g(�M��������Ԃ𗘗p)                                 
    W(:,:,nuser)=VT(:,ne,nuser)*VTT(:,1:NRs,nuser); 
                                     
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