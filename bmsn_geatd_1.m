function [W,UTT,STT,RIP,SP] = bmsn_geatd_1(NT,NR,NU,H,a)
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
    [D1,IN] = sort(diag(abs(D)).','descend'); %~ = D1（下のセクション）
    EW = EW(:,IN);
    %切り替え条件%
    %第3固有値のSNR%
    SNR2=D1(:,2)./a;
    SNR3=D1(:,3)./a;
    %第1,2固有値の比%
    D12=D1(:,1)/D1(:,2);
    %条件分岐%
    if SNR3<=20 && D12<=5 && SNR2>=1
    EW(:,3) = EW(:,2);
    end
    if SNR3<=20 && D12>5 && D12<=8
    EW(:,3) = EW(:,1);
    end
    if SNR2<=7
     EW(:,3) = EW(:,1);
     EW(:,2) = EW(:,1);
    end

    %%
    Wopt(:,:,nuser) = EW(:,1:NR);
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
        YI = zeros(NR,NR);
        for nn=nuser2
            YI=YI+UTT(:,1:NR,nuser)'*H(ns,:)*W(:,:,nn);
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