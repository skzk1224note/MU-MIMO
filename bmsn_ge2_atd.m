% bmsn2_gev_SINR.m
% BMSN(GEV) algorithm for NU-user
% Condition : NT >= NR*NU

function [W,UTT,STT,RIP,SP] = bmsn_ge2_atd(NT,NR,NU,H,a)
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
    [EW,D] = eig(B,A);
    [D1,IN] = sort(diag(abs(D)).','descend');
    D12=D1(:,1)/D1(:,2);
    D2(:,:,nuser) = diag(D1(:,1:NR));
    EW = EW(:,IN);
    if D12 < 2.5
    Wopt(:,:,nuser)=EW(:,1:NR)*sqrt(D2(:,:,nuser)); % Weighting by sqrt of SINR 
    else
    Wopt(:,:,nuser)=EW(:,1:NR);
    end
    Wopt(:,:,nuser) = Wopt(:,:,nuser)/norm(Wopt(:,:,nuser),'fro')*sqrt(NR); % Normalization
    

    % Wopt��nuer�̃`���l���s��ɏ�Z -> Block channel matrix: HTT                                 
    HTT=Hu*Wopt(:,:,nuser);
%     Hu2 = Hu
%     HTT2 = HTT
    % Block channel matrix��SVD
    [UTT(:,:,nuser),STT(:,:,nuser),VTT(:,:,nuser)]=svd(HTT); 
    % nuser�̃E�G�C�g(�M��������Ԃ𗘗p)     
    W(:,:,nuser)=Wopt(:,:,nuser)*VTT(:,1:NR,nuser);
    % �m�F�p�F���݂̓R�����g
    %UTT(:,:,nuser)'*H(ns,:)*W(:,:,nuser)
    %STT(:,:,nuser)
    %W2=[W2,W(:,:,nuser)];
    %NormU2 = norm(UTT(:,:,nuser),'fro')
    %NormV2 = norm(VTT(:,:,nuser),'fro')
    %NormW2 = norm(Wopt(:,:,nuser),'fro')
end
% %% �ŗL���[�h�`����̃E�G�C�g�m�F
% D2_bmsnge2 = D2
% 
% Wopt2 = Wopt
%  UTT_bmsnge2 = UTT
% % STT_bmsnge2 = STT
%  VTT_bmsnge2 = VTT

%%
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