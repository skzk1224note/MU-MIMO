function y = Decode(x,a)

switch a
case 0
    y = x;
case 1
    y = (sign(real(x))+1)/2;
case 2
    x = x*1.4142;
    y = (real(x)>=0)+(imag(x)>=0)*2;
case 3
    x = x*2.2882;
    x1 = real(x);
    x1 = x1+(x1<-3).*(-x1-3)+(x1>3).*(-x1+3);
    x1 = round((x1+3)/2);
    x2 = imag(x)>=0;
    y = x1+x2*4;
case 4
    x = x*2.9954;
    x1 = real(x);
    x2 = imag(x);
    x1 = x1+(x1<-3).*(-x1-3)+(x1>3).*(-x1+3);
    x2 = x2+(x2<-3).*(-x2-3)+(x2>3).*(-x2+3);
    x1 = round((x1+3)/2);
    x2 = round((x2+3)/2);
    y = x1+x2*4;
case 5
    x = x*4.2302;
    x1 = real(x);
    x2 = imag(x);
    if (abs(x1)>=4) && (abs(x2)>=4)
    if abs(x1)>=abs(x2)
        x = (sign(x1)*5+1i*sign(x2)*3);
    else
        x = (sign(x1)*3+1i*sign(x2)*5);
    end
        x1 = real(x);
        x2 = imag(x);
    end
    x1 = x1+(x1<-5).*(-x1-5)+(x1>5).*(-x1+5);
    x2 = x2+(x2<-5).*(-x2-5)+(x2>5).*(-x2+5);
    x1 = round((x1+5)/2);
    x2 = round((x2+5)/2);
    y = x1+x2*6;
    y = y-1-(y>4)-(y>29);
case 6
    x = x*6.0869;
    x1 = real(x);
    x2 = imag(x);
    x1 = x1+(x1<-7).*(-x1-7)+(x1>7).*(-x1+7);
    x2 = x2+(x2<-7).*(-x2-7)+(x2>7).*(-x2+7);
    x1 = round((x1+7)/2);
    x2 = round((x2+7)/2);
    y = x1+x2*8;
case 7
    x = x*7.8896;
    x1 = real(x);
    x2 = imag(x);
    if (abs(x1)>=8) && (abs(x2)>=8)
    if (abs(x1)<10) && (abs(x2)<10)
    if abs(x1)>=abs(x2)
        x = (sign(x1)*9+1i*sign(x2)*7);
    else
        x = (sign(x1)*7+1i*sign(x2)*9);
    end
    else
    if abs(x1)>=abs(x2)
        x = (sign(x1)*11+1i*sign(x2)*7);
    else
        x = (sign(x1)*7+1i*sign(x2)*11);
    end
    end
        x1 = real(x);
        x2 = imag(x);
    end
    x1 = x1+(x1<-5).*(-x1-5)+(x1>5).*(-x1+5);
    x2 = x2+(x2<-5).*(-x2-5)+(x2>5).*(-x2+5);
    x1 = round((x1+5)/2);
    x2 = round((x2+5)/2);
    y = x1+x2*6;
    y = y-1-(y>4)-(y>29);
case 8
    x = x*12.2253;
    x1 = real(x);
    x2 = imag(x);
    x1 = x1+(x1<-15).*(-x1-15)+(x1>15).*(-x1+15);
    x2 = x2+(x2<-15).*(-x2-15)+(x2>15).*(-x2+15);
    x1 = round((x1+15)/2);
    x2 = round((x2+15)/2);
    y = x1+x2*16;
end