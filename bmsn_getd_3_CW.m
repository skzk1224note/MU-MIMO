% bmsn3_gev_SINR.m
% BMSN(GEV) algorithm for NU-user
% Condition : NT >= NR*NU

function [W,UTT,STT,RIP,SP] = bmsn_getd_3_CW(NT,NR,NU,H,a)
% Wopt=zeros(NT,NR,NU);
% W=zeros(NT,NR-1,NU);
% UTT=zeros(NR,NR,NU);STT=zeros(NR,NR,NU);VTT=zeros(NR,NR,NU);
% ne = NR*(NU-1)+1:NT;

for nuser=1:NU
    % nuserにおける受信アンテナ番号
    ns = NR*(nuser-1)+1:NR*nuser; 
    % 全ユーザを統合したチャネル行列 (NU*NR) * NT   
    HT=H;
    % nuserのチャネル行列を抜き取り
    Hu=HT(ns,:);
    HT(ns,:)=[];
    % BMSN with GEV
    A = HT'*HT+a*eye(size(HT,2));
    B = Hu'*Hu;
    [EW,D] = eig(B,A);
    [D1,IN] = sort(diag(abs(D)).','descend');
    EW = EW(:,IN);
    %EW(:,NR) = EW(:,NR-1);
    %EW(:,NR-1) = EW(:,NR-2);
    
    De(:,:,nuser) = diag(D1(:,1:NR))
    %%
    Wopt(:,:,nuser) = EW(:,1:NR)*(log2(De(:,:,nuser)+1)); % Weight
    Wopt(:,:,nuser) = Wopt(:,:,nuser)/norm(Wopt(:,:,nuser),'fro')*sqrt(NR); % Normalization
    
    
    % Woptをnuerのチャネル行列に乗算 -> Block channel matrix: HTT                                 
    HTT=Hu*Wopt(:,:,nuser);
    
    % HTT3 = HTT
    % Block channel matrixをSVD
    [UTT(:,:,nuser),STT(:,:,nuser),VTT(:,:,nuser)]=svd(HTT,'econ');
    % nuserのウエイト(信号部分空間を利用)
    if NR > 1
        W(:,:,nuser)=Wopt(:,:,nuser)*VTT(:,1:NR,nuser);
    else
        W(:,:,nuser)=Wopt(:,:,nuser)*VTT(:,1:NR,nuser);
    end
    % 確認用：現在はコメント
    %UTT(:,:,nuser)'*H(ns,:)*W(:,:,nuser)
    %STT(:,:,nuser)
    %W2=[W2,W(:,:,nuser)];
%     NormU3 = norm(UTT(:,:,nuser),'fro')
%     NormV3 = norm(VTT(:,:,nuser),'fro')
%     NormW3 = norm(Wopt(:,:,nuser),'fro')
end
%% 固有モード伝送後のウエイト確認
%D2_bmsnge3 = D2
% HTT3 = HTT
% Wopt3 = Wopt
%  UTT_bmsnge3 = UTT
% STT_bmsnge3 = STT
%VTT_bmsnge3 = VTT

%%
% 所望波＆干渉波電力の計算
% SP = zeros(NR,NU);
% RIP = zeros(NR,NU);
for nuser=1:NU
    % nuserにおける受信アンテナ番号
    ns = NR*(nuser-1)+1:NR*nuser; 
        
    nuser2=1:NU;
    nuser2(nuser)=[];
    
   if NR > 1
        YI = zeros(NR-1,NR);
        for nn=nuser2
            YI=YI+UTT(:,1:NR-1,nuser)'*H(ns,:)*W(:,:,nn);
        end
    else
        YI = zeros(NR,NR);
        for nn=nuser2
            YI=YI+UTT(:,1:NR,nuser)'*H(ns,:)*W(:,:,nn);
        end
    end
    RIP(:,nuser) = sum(abs(YI).^2,2); % 干渉波電力
    if NR > 1
        YS=UTT(:,1:NR,nuser)'*H(ns,:)*W(:,:,nuser);
    else
        YS=UTT(:,1:NR,nuser)'*H(ns,:)*W(:,:,nuser);
    end
    SP(:,nuser) = sum(abs(YS).^2,2); % 所望波電力
    
end
% YI3 = YI
% YS3 = YS
% RIP3 = RIP
% 確認用：現在はコメント
% abs(H*W2)