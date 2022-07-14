% Directional pattern
clear;

derad = pi/180;
theta = -90:0.1:90; % -90degから90degまで0.1deg刻みで
theta_iso = -180:0.1:180; % -180degから180degまで0.1deg刻みで
Ltheta = length(theta);
Ltheta_iso = length(theta_iso);
% cosθのn乗（nは整数値）
cos_n = 0;
An = sqrt(2 * (cos_n + 1));

E0 = zeros(1,Ltheta_iso);
E01 = zeros(1,Ltheta_iso);
E02 = zeros(1,Ltheta_iso);

E1 = zeros(1,Ltheta);
E11 = zeros(1,Ltheta);
E2 = zeros(1,Ltheta);
E3 = zeros(1,Ltheta);
E4 = zeros(1,Ltheta);
E5 = zeros(1,Ltheta);

Nt = 16;
Nr = 2;
d_t = 0.5;
d_r = 0.5;
%%
for ithe_iso = 1:Ltheta_iso
    u_t  = 2*pi*d_t*sin(theta_iso(ithe_iso)*derad); % u=kdcosθ
    u_r  = 2*pi*d_r*sin(theta_iso(ithe_iso)*derad); % u=kdcosθ
    E0(ithe_iso) = 1;    
    E01(ithe_iso) = abs(sin(Nt*u_t/2)/sin(u_t/2));    % 無指向性送信アレー指向性
    E02(ithe_iso) = abs(sin(Nr*u_r/2)/sin(u_r/2));    % 無指向性受信アレー指向性
end
for ithe = 1:Ltheta
    u_t = 2*pi*d_t*sin(theta(ithe)*derad); % u=kdcosθ
    u_r = 2*pi*d_r*sin(theta(ithe)*derad); % u=kdcosθ
    E1(ithe) = sqrt(2 * (0 + 1)) * cos(theta(ithe)*derad)^(0/2);
    E2(ithe) = sqrt(2 * (1 + 1)) * cos(theta(ithe)*derad)^(1/2);
    E3(ithe) = sqrt(2 * (2 + 1)) * cos(theta(ithe)*derad)^(2/2);
    E4(ithe) = sqrt(2 * (3 + 1)) * cos(theta(ithe)*derad)^(3/2);
    E5(ithe) = sqrt(2 * (4 + 1)) * cos(theta(ithe)*derad)^(4/2);
end
%% 正規化

Pst = E2.^2;    % Power pattern
Pst = Pst/max(Pst);    % Power pattern
Psr = E3.^2;    % Power pattern
Psr = Psr/max(Psr);    % Power pattern

E01 = E01/max(E01);
E02 = E02/max(E02);

PtdB = 10*log10(Pst);  % dB値
PrdB = 10*log10(Psr);  % dB値
%%
% figure % 直角座標表示
% plot(theta, PtdB,'-b','LineWidth',1)
% axis([-90 90 -100 0]);
% set(gca, 'XTick',-90:30:90,'FontSize',10,'FontName','Arial')
% xlabel('theta angle [degree]','FontSize',12,'FontName','Arial');
% ylabel('relative amplitude [dB]','FontSize',12,'FontName','Arial');
% title('Directional Pattern','FontSize',14,'FontName','Arial');
% grid on
% hold on
% plot(theta, PrdB,'-m','LineWidth',1)
% axis([-90 90 -100 0]);
% set(gca, 'XTick',-90:30:90,'FontSize',10,'FontName','Arial')
% xlabel('theta angle [degree]','FontSize',12,'FontName','Arial');
% ylabel('relative amplitude [dB]','FontSize',12,'FontName','Arial');
% title('Directional Pattern','FontSize',14,'FontName','Arial');
% grid on
% hold on
%%
figure; % 極座標表示
h=polarplot(theta_iso*derad,E0,'k-');
set(h,'LineWidth',2)
rlim([0 4])
grid on
hold on

h=polarplot(theta*derad,E1,'-m');
set(h,'LineWidth',2)
grid on
hold on

h=polarplot(theta*derad,E2,'-b');
set(h,'LineWidth',2)
grid on
hold on

h=polarplot(theta*derad,E3,'-g');
set(h,'LineWidth',2)
grid on
hold on

h=polarplot(theta*derad,E4,'-c');
set(h,'LineWidth',2)
grid on
hold on

h=polarplot(theta*derad,E5,'-r');
set(h,'LineWidth',2)
grid on
hold on

lgd = legend('isotropic','cos^0 \theta','cos^1 \theta','cos^2 \theta','cos^3 \theta','cos^4 \theta','Location','northeast');
lgd.FontSize = 14;
% End