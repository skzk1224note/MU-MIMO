function y = Mapping(x,a)

switch a   
case 0
    y = zeros(size(x));
case 1
    y = x*2-1;
case 2
    y = (rem(x,2)*2-1+1i*(floor(x/2)*2-1))/1.4142; %sqrt(2)
case 3
    y = (rem(x,4)*2-3+1i*(floor(x/4)*2-1))/2.2882; %(sqrt(2)+sqrt(10))/2
case 4
    y = (rem(x,4)*2-3+1i*(floor(x/4)*2-3))/2.9954; %sqrt(2)+sqrt(10)/2
case 5
    y = (rem(x+1+(x>3)+(x>27),6)*2-5+1i*(floor(x/6)*2-5))/4.2302; %(2*sqrt(2)+sqrt(10)+sqrt(26)+sqrt(34))/4
case 6
    y = (rem(x,8)*2-7+1i*(floor(x/8)*2-7))/6.0869; %(13*sqrt(2)+sqrt(10)+sqrt(26)+sqrt(34)+sqrt(58)+sqrt(74))/8
case 7
    y = (rem(x+2+(x>7)*4+(x>15)*2+(x>111)*2+(x>119)*4,12)*2-11+1i*(floor(x/12)*2-11))/7.8896; %(13*sqrt(2)+4*sqrt(10)+10*sqrt(13)+10*sqrt(17)+sqrt(26)+sqrt(34)+sqrt(58)+sqrt(74)+sqrt(82)+sqrt(106)+sqrt(122)+sqrt(130)+sqrt(146))/24
case 8
    y = (rem(x,16)*2-15+1i*(floor(x/16)*2-15))/12.2253;
end