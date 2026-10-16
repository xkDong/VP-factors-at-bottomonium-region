function PIL=PIL(s,m2)
alpha=1.D0/137.035999074;
if s<4.0*m2
    a=(4.0*m2/s-1)^0.5;
    PIL=alpha/pi*(8.0/9.0+a*a/3.0-2.0*(0.5+a*a/6.0)*a*acot(a));
    return
else 
    a=(1-4.0*m2/s)^0.5;
    b=(1-a)/(1+a);
    PILR=alpha/pi*(8.0/9.0-a*a/3.0+(0.5-a*a/6.0)*a*log(b));
    PILI=a*alpha/3*(1+2*m2/s);
    PIL=PILR+1i*PILI;
    return
end
end