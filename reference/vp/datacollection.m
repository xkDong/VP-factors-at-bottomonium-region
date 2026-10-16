vp1=real(vp1);
plot(Etot(2:end),vp1(2:end,1))
hold on
plot(Etot(2:end),vp1(2:end,2))
plot(Etot(2:end),vp1(2:end,3))
% plot(Etot(2:5601),vp1(2:5601,4))
% plot(Etot(2:5601),vp1(2:5601,5))
hold off
nEbmin=4250;%mim 4145
nEbmax=5350;%max 5429
EShowIter=Etot(nEbmin:nEbmax);
delvpIter21=(vp1(nEbmin:nEbmax,2)-vp1(nEbmin:nEbmax,1))./vp1(nEbmin:nEbmax,1);
delvpIter32=(vp1(nEbmin:nEbmax,3)-vp1(nEbmin:nEbmax,2))./vp1(nEbmin:nEbmax,2);
figure
plot(EShowIter,vp1(nEbmin:nEbmax,1),EShowIter,vp1(nEbmin:nEbmax,2),EShowIter,vp1(nEbmin:nEbmax,3))
figure
plot(EShowIter,delvpIter21,EShowIter,delvpIter32)
%vpplot
Ebmin=4250;%mim 4145
Ebmax=5350;%max 5429
vp2b1=vp2(:,Ebmin:Ebmax,1);
vp2b1mean=mean(vp2b1(2:end,:));
vp2b1err=std(vp2b1(2:end,:));
vp2b2=vp2(:,Ebmin:Ebmax,2);
vp2b2mean=mean(vp2b2(2:end,:));
vp2b2err=std(vp2b2(2:end,:));
vp2b3=vp2(:,Ebmin:Ebmax,3);
vp2b3mean=mean(vp2b3(2:end,:));
vp2b3err=std(vp2b3(2:end,:));
% figure
% hold on
% EShowIter=Etot(Ebmin:Ebmax);
% plot(Etot(Ebmin:Ebmax),vp2b1mean)
% 
% plot(Etot(Ebmin:Ebmax),vp2b2mean)
% 
% plot(Etot(Ebmin:Ebmax),vp2b3mean)
% hold off
figure
% errorbar(Etot(Ebmin:Ebmax),vp2b1mean,vp2b1err)
% hold on
% errorbar(Etot(Ebmin:Ebmax),vp2b2mean,vp2b2err)
errorbar(Etot(Ebmin:Ebmax),vp2b3mean,vp2b3err);
%hold off

%vpb-belleenergy
nEbmin=4300;
vpBelleSamp=zeros(10000,210);
for i=1:210
    vpBelleSamp(:,i)=vp2(:,nEbmin+5*i-4,3);
end
vpbmean_belle_energy=mean(vpBelleSamp)';
vpberr_belle_energy=std(vpBelleSamp)';
vpbcov_belle_energy=cov(vpBelleSamp);
RbBornSamp=Rbdresamp./vpBelleSamp;
RbBornmean=mean(RbBornSamp)';
RbBornerr=std(RbBornSamp)';
RbBorncov=cov(RbBornSamp);
Rbdremean=mean(Rbdresamp)';
Rbdreerr=std(Rbdresamp)';
Rbdrecov=cov(Rbdresamp);


