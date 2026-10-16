clear;
clc;
load datacombi
datacombi(:,1)=[];
datacombi(1,1)=10.55852;
datacombi(2:6,1)=zeros(5,1);

datacombis(1,:)=datacombi(1,:);
for i=2:3
    datacombis(i,:)=datacombi(i,:).*86.85./datacombi(1,:).^2;
end
datacombis(6,1)=0;
n=30;
l=10000;
m=length(datacombis(1,:));

curve=1;
[sigma1,delta]=ISRcalculation(datacombis,n,l,curve);
onedeltasamp=zeros(l,m);
sigmadresamp=zeros(l,m);
for i=1:m
    onedeltasamp(:,i)=delta(:,i,n);
    sigmadresamp(:,i)=sigma1(:,i)./onedeltasamp(:,i);
end
[sigmadre(1,:),sigmadre(2,:)]=normfit(sigmadresamp);
[onedelta(1,:),onedelta(2,:)]=normfit(onedeltasamp);
onedeltacov=cov(onedeltasamp);
Rbdresamp=zeros(l,m);
for i=1:m
Rbdresamp(:,i)=sigmadresamp(:,i)*datacombi(1,i)^2/86.85;
end
[Rbdre(1,:),Rbdre(2,:)]=normfit(Rbdresamp);
Rbdrecov=cov(Rbdresamp);
save 'ISR_factor_sta_err.mat'
% %iteration convergence
% iterationdelta=ones(n,m);
% for i=1:n
%     for j=1:m
%         iterationdelta(i,j)=delta(1,k*j-j+1,i);
%     end
% end
% figure
% hold on
% for i=1:n
%     plot(energy,iterationdelta(i,:));
% end
% hold off