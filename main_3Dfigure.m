clear;
derad = pi/180;
theta = (-90:1:90).*derad; % -90degから90degまで0.1deg刻みで
phi = (-90:1:90).*derad;

n = 0;
r = 2*(n+1);
Ex = r^(1/2) * cos(theta).^(n/2) .* cos(phi);
Ey = r^(1/2) * cos(theta).^(n/2) .* sin(phi);

[X,Y] = meshgrid(Ex,Ey);
Z = r * cos(theta).^(n/2);
surf(X,Y,abs(Z))