function y = peruser(x,a)

b = size(x,1)/a;
y = zeros(b,size(x,2),a);
for i = 1:a
for ii = 1:b
    y(ii,:,i) = x((ii+b*(i-1)),:);
end
end