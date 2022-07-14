% bmsn3_gev_SINR.m
% BMSN(GEV) algorithm for NU-user
% Condition : NT >= NR*NU

function [W,UTT,STT,RIP,SP] = bmsn_ge_atd1(NT,NR,NU,H,a)
Wopt=zeros(NT,NR,NU);W=zeros(NT,NR,NU);
UTT=zeros(NR,NR,NU);STT=zeros(NR,NR,NU);VTT=zeros(NR,NR,NU);
% ne = NR*(NU-1)+1:NT;

for nuser=1:NU
    % nuser�ɂ������M�A���e�i�ԍ�
    ns = NR*(nuser-1)+1:NR*nuser; 
    % �S���[�U�𓝍������`���l���s�� (NU*NR) * NT   
    HT=H;
    
    % HT�s��̓���ns�s��S�Ē��o
    Hu=HT(ns,:);
    % nuser�̃`���l���s��𔲂����
    HT(ns,:)=[];
    
    % BMSN with GEV
    %��ʌŗL�l�W�J
    A = HT'*HT+a*eye(size(HT,2));
    B = Hu'*Hu;
    %�ŗL�x�N�g���ƌŗL�l�̍s���Ԃ�
    [EW,D] = eig(B,A);
    % sort(diag(abs(D)).','descend')=Db
    [D1,IN] = sort(diag(abs(D)).','descend');    %~ = D1�i���̃Z�N�V�����j
    %largeD=diag(abs(D));
    %largeD()
    D_bf=squatm(D1);
    %�ŗL�l1��2�̔䗦
    D12=D1(:,1)./D1(:,2);
    %sigma2=a/NT;
    %���IN�����ꂼ��c�Ɏ��o��
    %�ŗL�l�̑傫�����ɑΉ�����ŗL�x�N�g������ׂ�
    EW = EW(:,IN);
    %% lambda2��������������lambda1�Ƃ��Clambda1�ɑΉ�����ŗL�x�N�g���݂̂��ӂ��g�p
    if  NR>1
    %if D12 > 3
            EW(:,NR) = EW(:,NR-1);
    %end
%NR��NR-1�̒ʘH�ł̌ŗL�x�N�g���𓯂��ɂ���
        %EW(:,NR) = EW(:,NR-1);
    end
    %%
    %��ʌŗL�l���̉��𑗐M�E�G�C�g�ɂ���
    Wopt(:,:,nuser) = EW(:,1:NR);
    %GE����BF�ւ̐؂�ւ������d�ݕt��
    
    Wopt(:,:,nuser) = Wopt(:,:,nuser)/norm(Wopt(:,:,nuser),'fro')*sqrt(NR); % Normalization
    
    %�ŗL���[�h�`��
    % Wopt��nuer�̃`���l���s��ɏ�Z -> Block channel matrix: HTT                                 
    HTT=Hu*Wopt(:,:,nuser);
    
    % HTT3 = HTT
    % Block channel matrix��SVD
    [UTT(:,:,nuser),STT(:,:,nuser),VTT(:,:,nuser)]=svd(HTT);
    % nuser�̃E�G�C�g(�M��������Ԃ𗘗p)
    if NR > 1
        W(:,:,nuser)=Wopt(:,:,nuser)*VTT(:,1:NR,nuser); %lambda1�ɑΉ�����ŗL�x�N�g���݂̂��ӂ��g�p
    end
    % �m�F�p�F���݂̓R�����g
    %UTT(:,:,nuser)'*H(ns,:)*W(:,:,nuser)
    %STT(:,:,nuser)
    %W2=[W2,W(:,:,nuser)];
%     NormU3 = norm(UTT(:,:,nuser),'fro')
%     NormV3 = norm(VTT(:,:,nuser),'fro')
%     NormW3 = norm(Wopt(:,:,nuser),'fro')
end
%% �ŗL���[�h�`����̃E�G�C�g�m�F
%D2_bmsnge3 = D2
% HTT3 = HTT
% Wopt3 = Wopt
%  UTT_bmsnge3 = UTT
% STT_bmsnge3 = STT
%VTT_bmsnge3 = VTT

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
    %1����NU�܂ŌJ��Ԃ�
    for  nn=nuser2
        YI=YI+UTT(:,1:NR,nuser)'*H(ns,:)*W(:,:,nn);
   
    end
    
    
   
    %sum��YI�̂��ꂼ��̗v�f�̐�Βl�̊e�s�̗v�f�𑫂����킹��
    RIP(:,nuser) = sum(abs(YI).^2,2); % ���g�d��
   
    YS=UTT(:,1:NR,nuser)'*H(ns,:)*W(:,:,nuser);
    
    SP(:,nuser) = sum(abs(YS).^2,2); % ���]�g�d��
    
end
% YI3 = YI
% YS3 = YS
% RIP3 = RIP
% �m�F�p�F���݂̓R�����g
% abs(H*W2)