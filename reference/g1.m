function g1=g1(x,b,A,B)
beta = 0.004645639*(log(b*1000000/0.26111998)-1);
if x==0
    g1=(A*beta*(x - 2)^2)/4  - (A*beta^2*(16*dilog(x)))/32 ;
else
    logx=log(x);
    log1x=log(1-x);
    x2=x^2;
    x3=x^3;
deltaprime = 0.006480359964569+(0.75-0.541217.*beta).*beta;
g1=  (B*beta^2*(18*x + 18*log1x - 36*x2*logx + 12*x3*logx - 36*x*log1x - 18*x2 + 2*x3 + 27*x2*log1x - 9*x3*log1x))/72 + (A*beta*(x - 2)^2)/4 + A*x^beta*(deltaprime + 1) - (A*beta^2*(10*x + 18*log1x + 16*dilog(x) - 8*x2*logx - 24*x*log1x + 32*x*logx - x2 + 6*x2*log1x + 16*log1x*logx))/32 + (B*beta*x2*(x - 3))/6 + (2*B*beta*x^(beta + 1)*(deltaprime + 1))/(2*beta + 2);
end