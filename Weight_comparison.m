clear;

lam = 1e-4:1e-4:100;
lamx = 10*log10(lam);

y_gev = lam./lam;
y_ge4 = log2(1+lam);
y_ge5 = sqrt(lam);
y_ge6 = (log2(1+lam)).^2;

figure
plot(lamx,y_gev,'r-',lamx,y_ge4,'b-',lamx,y_ge5,'g-',lamx,y_ge6,'k-')
set(gca,'Fontsize',14,'Fontname','Times New Roman');
legend('1','log2(1+λ)','√λ','{log2(1+λ)}^2');
xlabel('Eigenvalue (λ) [dB]','Fontsize',14,'Fontname','Times New Roman');
ylabel('Weight','Fontsize',16,'Fontname','Times New Roman');
set(gca,'Fontsize',16,'Fontname','Times New Roman');
grid on;
hold on;