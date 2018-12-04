function y = hermite(n,x)
if n == 0
    y = 1;
elseif n==1
    y = x;
elseif n==2
    y = (x^2-1)/sqrt(2); 
end
end