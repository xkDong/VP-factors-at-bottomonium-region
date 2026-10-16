function sigmu=sigmu(s)
beta=(s/4-0.105658^2).^0.5./(s/4).^0.5;
sigmu=4*pi*1/137.036^2/3./s.*beta.*(1.5-0.5*beta.^2)*389379000;
 %sigmu=4*pi*1/137.036^2/3./s*389379;
end