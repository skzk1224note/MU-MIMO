% cap_bd_muser.m
% NU-ユーザでのBD法の計算
% TDMA, Upper boundと比較

clear;%close;


% パラメータ条件 NT >= NR*NU
K_tar = 20;      % ターゲットK[dB]

CDF = 50;         % ABRのターゲットCDF
SNR_tar = 10;     % ターゲットSNR[dB]
Nt = 16;          % 送信素子数
Nr = 2;           % 受信素子数 (アンテナ選択の場合1, そうでない場合は2)
Nu = 8;           % ユーザ数
SIMU = 1000;      % 試行回数（通常 1000）

An = 2; % 指向性関数の係数(cosパターンの場合は2)

I = eye(Nt,Nt);
Nru = Nr*Nu;


dr_min = 0.5;
dr_max = 1.0;
dr_box = (dr_min:0.05:dr_max); % 送信アンテナ間隔（in wavelength)
Ldr=length(dr_box);        % dt_boxの大きさ

d_t = 0.5;        % 受信アンテナ間隔（in wavelength)
derad = pi/180;   % degree -> rad



%T 所望のチャネル行列
T = zeros(Nr,Nr);
for inu = 1:Nu
    T(:,:,inu) = eye(Nr,Nr);
end


Q = zeros(SIMU,5,Nu);
QmC = zeros(Ldr,5);
QmC_cos = zeros(Ldr,5);
% 入力ファイル基幹名
% folder='CSV3/';
% fn1 = 'Eig16x2x8u_BD_SNR';
% fn2 = 'Eig16x2x8u_BMSN1(GEV)_SNR';
% fn3 = 'Eig16x2x8u_BMSN2a_SNR';  % alpha = 1e-2
% fn4 = 'Eig16x2x8u_BMSN2b_SNR';  % alpha = 1e-6
% fn5 = 'Eig16x2x8u_BMSN3_SNR';   % alpha = sigma2Nt equal to MMSE
% fn6 = 'Eig16x2x8u_ZF_SNR';
% fn7 = 'Eig16x2x8u_MMSE_SNR';

% 入力ファイル基幹名
% folder='CSV4/';
% fn1 = 'Eig16x2x8u_ZF_SNR';
% fn2 = 'Eig16x2x8u_BD_SNR';
% %fn3 = 'Eig16x2x8u_BMSN2_SNR'; % 2: alpha = 1e-2, 3: alpha=sigma2xNT
% %fn4 = 'Eig16x2x8u_BDAS_SNR';  % antenna selection
% %fn5 = 'Eig16x2x8u_BMSNAS2_SNR';  % alpha = 1e-2 and antenna selection
% fn3 = 'Eig16x2x8u_BMSN3_SNR';   % alpha = sigma2Nt equal to MMSE
% fn4 = 'Eig16x2x8u_BMSN3(GEV)_SNR';   % alpha = sigma2Nt equal to MMSE
% fn5 = 'Eig16x2x8u_MMSE_SNR';

% 出力ファイル名
%folder= 'CSV/';
%cdfn1 = 'ABRCDFvsSNR_16x2x8u_BD_MMSE_BMSNs_USERave_CDF';
%cdfn2 = 'Eig_CDF_16x2x8u_BMSN_SNR';
%%
if CDF < 10
   target_CDF=strcat('CDF=',num2str(CDF,'%01d'),'%');
else
   target_CDF=strcat('CDF=',num2str(CDF,'%02d'),'%');
end

if SNR_tar < 10
   target_SNR=strcat('SNR=',num2str(SNR_tar,'%01d'),'dB');
else
   target_SNR=strcat('SNR=',num2str(SNR_tar,'%02d'),'dB');
end

if K_tar < 10
   target_K=strcat('{\it{K}}=',num2str(K_tar,'%01d'),'dB');
else
   target_K=strcat('{\it{K}}=',num2str(K_tar,'%02d'),'dB');
end

if d_t < 10
   target_d_t=strcat('d_t=',num2str(d_t,'%.3g'),'\lambda');
else
   target_d_t=strcat('d_t=',num2str(d_t,'%02d'),'\lambda');
end
% 出力ファイル名 with SNR in dB
%cdffile1 = strcat(folder,cdfn1,num2str(CDF,'%02d'),'_1000itr.csv');
%cdffile2 = strcat(folder,cdfn2,num2str(SN_tar,'%02d'),'dB_1000itr.csv');
%%
H_los = zeros(Nu*Nr,Nt); % 伝搬チャネル行列の直接波成分(LOS チャネル)
sigma2 = 1/(10^(SNR_tar/10)); % noise power
a = sigma2*Nt; % 擬似雑音

St = zeros(Nr,Nu);
MSt_BMSN_BF = zeros(SIMU, Nr*Nu);
MSt_BMSN_GE = zeros(SIMU, Nr*Nu);
MSt_MMSE = zeros(SIMU, Nr*Nu);
MSt_GMI1 = zeros(SIMU, Nr*Nu);
MSt_GMI2 = zeros(SIMU, Nr*Nu);
% MSt_BD = zeros(SIMU, Nr*Nu);
% MSt_ZF = zeros(SIMU, Nr*Nu);

for Directivity_switch = 0:1 % 送受信素子の指向性考慮の有無 0:無,1:有
    for idr = 1:Ldr 
      dr_tar = dr_box(idr);
      % 出力ファイル名 with alpha and SNR in dB
      % evfile1 = strcat(folder,fn1,num2str(snr(isnr),'%02d'),'dB_1000itr.csv'); 
        for isimu = 1:SIMU      
          % 伝搬チャネル行列のマルチパス成分 (i.i.d.Rayleigh)
          H_iid = (randn(Nr*Nu,Nt)+1j*randn(Nr*Nu,Nt))/sqrt(2);
          % LOS チャネル      
            if Directivity_switch == 1 % 送受信素子の指向性考慮有り
             Theta_t = (rand(1,Nu)-0.5)*180; % ユーザ毎の送信角 指向性:(-90deg - 90deg)
             Theta_r = (rand(1,Nu)-0.5)*180; % ユーザ毎の受信角 指向性:(-90deg - 90deg)
             for n = 1 : Nu
                a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad))*An*cos(Theta_t(1,n)*derad); % ユーザ毎の送信モードベクトル
                a_r = exp(-1j*2*pi*dr_tar*(0:Nr-1).'*sin(Theta_r(1,n)*derad))*An*cos(Theta_r(1,n)*derad); % ユーザ毎の受信モードベクトル
                H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t';                 % ユーザ毎のLOSチャネル行列
             end
            elseif Directivity_switch == 0 % 送受信素子の指向性考慮無し
              Theta_t = (rand(1,Nu)-0.5)*360; % ユーザ毎の送信角 等方性:(-180deg - 180deg) 
              Theta_r = (rand(1,Nu)-0.5)*360; % ユーザ毎の受信角 等方性:(-180deg - 180deg) 
             for n = 1 : Nu
             a_t = exp(-1j*2*pi*d_t*(0:Nt-1).'*sin(Theta_t(1,n)*derad)); % ユーザ毎の送信モードベクトル
             a_r = exp(-1j*2*pi*dr_tar*(0:Nr-1).'*sin(Theta_r(1,n)*derad)); % ユーザ毎の受信モードベクトル
             H_los((n-1)*Nr+1:(n-1)*Nr+Nr,:) = a_r*a_t';                 % ユーザ毎のLOSチャネル行列
             end
            end
            % Kを真値にする
            K = 10^(K_tar/10);
            
            % H=[sqrt(K/(K+1))*(LOS チャネル)]+[sqrt(1/(K+1))*(NLOS チャネル)]
            H = sqrt(K/(K+1))*H_los + sqrt(1/(K+1))*H_iid;
            
            % BMSN-BF
            [~,~,STT,RIP,~] = bmsn_bf(Nt,Nr,Nu,H,a,T);
            for inu=1:Nu
                St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
            end
            %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
            MSt_BMSN_BF(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
            
            % BMSN-GE
            [~,~,STT,RIP,~] = bmsn_gev(Nt,Nr,Nu,H,a);
            for inu=1:Nu
                St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
            end
            %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
            MSt_BMSN_GE(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
            
            % MMSE-CI
            [~,~,STT,RIP,~] = mmse(Nt,Nr,Nu,H,a);
            for inu=1:Nu
                St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
            end
            %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
            MSt_MMSE(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
            
            % GMI1
            [~,~,STT,RIP,~] = gmmse_m1(Nt,Nr,Nu,H,a);
            for inu=1:Nu
                St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
            end
            %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
            MSt_GMI1(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
            
            % GMI2
            [~,~,STT,RIP,~] = gmmse_m2(Nt,Nr,Nu,H,a);
            for inu=1:Nu
                St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
            end
            %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
            MSt_GMI2(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
            
            % BD
            %     [~,~,STT,RIP,~] = bd(Nt,Nr,Nu,H);
            %     for inu=1:Nu
            %         St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
            %     end
            %     %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
            %     MSt_BD(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
            
            % ZF-CI
            %     [~,~,STT,RIP,~] = zf(Nt,Nr,Nu,H);
            %     for inu=1:Nu
            %         St(:,inu) = diag(STT(1:Nr,1:Nr,inu));
            %     end
            %     %受信アンテナ（nr)の固有値分布 for B-MSN（干渉波成分を考慮）
            %     MSt_ZF(isimu,:)=((reshape(St,[Nru,1]).').^2)./(reshape(RIP,[Nru,1]).'+Nt*sigma2);
            
            % fprintf('Iteration = %d / %d\n',isimu, SIMU);
            
        end % isimu
        
        % 各アルゴリズムのABR
        for nuser=1:Nu
            ns = Nr*(nuser-1)+1:Nr*nuser;
            Q(:,1,nuser)=sort(sum(log2(1 + MSt_BMSN_BF(:,ns)),2));
            Q(:,2,nuser)=sort(sum(log2(1 + MSt_BMSN_GE(:,ns)),2));
            Q(:,3,nuser)=sort(sum(log2(1 + MSt_MMSE(:,ns)),2));
            Q(:,4,nuser)=sort(sum(log2(1 + MSt_GMI1(:,ns)),2));
            Q(:,5,nuser)=sort(sum(log2(1 + MSt_GMI2(:,ns)),2));
            %         Q(:,6,nuser)=sort(sum(log2(1 + MSt_BD(:,ns)),2));
            %         Q(:,7,nuser)=sort(sum(log2(1 + MSt_ZF(:,ns)),2));
        end % nuser end
        
        
        
        Qm = mean(Q,3); % Qのユーザ回数平均    
        if Directivity_switch == 0 % 指向性無しの場合の行列
            QmC(idr,:)=Qm(round(CDF*SIMU/100),:);% Channel capacity
        elseif Directivity_switch == 1 % 指向性有りの場合の行列
            QmC_cos(idr,:)=Qm(round(CDF*SIMU/100),:);% Channel capacity
        end
        fprintf('dr = %.3g lambda \n',dr_box(idr));
    end % idt end
    %csvwrite(cdffile1,[dr,QmC]);
end
%% グラフ表示 figure ABRvsK
figure;
mycol = [1 0 1;0 1 0;1 0 0;0 0 1;0 1 1;0 0 0]; % 色
set(groot,'defaultAxesColorOrder',mycol)
%plot(K,QmC(:,1),K,QmC(:,2),K,QmC(:,3),K,QmC(:,4),'Linewidth',2);

plot(dr_box,QmC(:,1),'r--',dr_box,QmC(:,2),'b--',dr_box,QmC(:,3),'g--',dr_box,QmC(:,4),'y--',dr_box,QmC(:,5),'k--',...
    dr_box,QmC_cos(:,1),'r-',dr_box,QmC_cos(:,2),'b-',dr_box,QmC_cos(:,3),'g-',dr_box,QmC_cos(:,4),'y-',dr_box,QmC_cos(:,5),'k-','Linewidth',2);
axis([dr_min,dr_max,min(min(QmC))-0.05,max(max(QmC_cos))+2])
set(gca,'Fontsize',18,'Fontname','Times New Roman');
lgd = legend;
lgd.NumColumns = 2; % 凡例の列数を指定
legend('BMSN-BF','BMSN-GE','MMSE-CI','GMI1','GMI2',...
    'BMSN-BFcos','BMSN-GEcos','MMSE-CIcos','GMI1cos','GMI2cos','Location','Northwest');
xlabel('{\it{d_r}} [\lambda]','Fontsize',6,'Fontname','Times New Roman');
ylabel('Channel capacity [bits/s/Hz]','Fontsize',6,'Fontname','Times New Roman');
set(gca,'Fontsize',18,'Fontname','Times New Roman');
title(strcat(target_CDF,',',target_SNR,',',target_K,',',target_d_t));
grid on;
hold on;
% End