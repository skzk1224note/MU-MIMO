% bmsn3_gev_SINR.m
% BMSN(GEV) algorithm for NU-user
% Condition : NT >= NR*NU

function [W,UTT,STT,RIP,SP] = bmsn_ge_atd2(NT,NR,NU,H,a)
Wopt=zeros(NT,NR,NU);W=zeros(NT,NR,NU);
UTT=zeros(NR,NR,NU);STT=zeros(NR,NR,NU);VTT=zeros(NR,NR,NU);
% ne = NR*(NU-1)+1:NT;
D2=zeros(NR,NR,NU);


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
    [D1,IN] = sort(diag(abs(D)).','descend');    %~ = D1（下のセクション）
    D12=D1(:,1)/D1(:,2);
    SNR1=D1(:,1)./a;
    SNR2=D1(:,2)./a;
    EWs = EW(:,IN);
    %固有値のルートの値をD3として格納
    D3=sqrt(D1);
    %D3をEの行列にかけれるように行列として対角化
    D2(:,:,nuser) = diag(D3(:,1:NR));
    %% Eigen SNR2が使えない状態であると判断されるdB以下で第一固有値の道へ変更
    if SNR2<18.36
    %if
       % D1(:,NR) < 0.1
        EWs(:,NR) = EWs(:,NR-1);
    end
    %BFへの切り替え
    %if SNR2<18.36
    Wopt(:,:,nuser) = EWs(:,1:NR);
    %else
    %Wopt(:,:,nuser) = EWs(:,1:NR)*(D2(:,:,nuser));
    %end
    Wopt(:,:,nuser) = Wopt(:,:,nuser)/norm(Wopt(:,:,nuser),'fro')*sqrt(NR);
    % Normalization
    
    
    % Woptをnuerのチャネル行列に乗算 -> Block channel matrix: HTT                                 
    HTT=Hu*Wopt(:,:,nuser);
    
    % HTT3 = HTT
    % Block channel matrixをSVD
    [UTT(:,:,nuser),STT(:,:,nuser),VTT(:,:,nuser)]=svd(HTT);
    % nuserのウエイト(信号部分空間を利用)
    if NR > 1
        W(:,:,nuser)=Wopt(:,:,nuser)*VTT(:,1:NR,nuser); %lambda1に対応する固有ベクトルのみをふたつ使用
    else
        W(:,:,nuser)=Wopt(:,:,nuser)*VTT(:,1:NR,nuser); %lambda1に対応する固有ベクトルのみをふたつ使用
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
SP = zeros(NR,NU);
RIP = zeros(NR,NU);
for nuser=1:NU
    % nuserにおける受信アンテナ番号
    ns = NR*(nuser-1)+1:NR*nuser; 
        
    nuser2=1:NU;
    nuser2(nuser)=[];
    if NR>1
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