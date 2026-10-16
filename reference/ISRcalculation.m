function [sigma1,delta]=ISRcalculation(datacombis,n,l,curve)
m=length(datacombis(1,:));
sigma=zeros(l,m,n);
sigma(:,:,1)=Newsigma(datacombis,l,curve);
sigma1=sigma(:,:,1);
for i=1:l
    sigma(i,23:end,1)=smooth(datacombis(1,23:end),sigma(i,23:end,1),0.045,'lowess');
end
delta=iteration(datacombis(1,:),sigma,l,m,n);
end
