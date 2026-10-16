
function PIH=PIH(s,n)
global Sigtot Stot;
alpha  =1.D0/137.035999074D0;
PImass=0.13957018;

RESmass=[9.590888,13.587399,89.491600,100.465741,107.230167];%mass of J/psi psi upsilon^2 in GeV^2
REStotwid=[.0929E-3, .304E-3,.05402E-3,.0319E-3, .02032E-3];%total width in GeV
RESlepwid=[5.55E-6*0.95, 2.35E-6*0.95,1.34E-6*0.93, .612E-6*0.93, .443E-6*0.93];%leptonic width in GeV
resi=RESlepwid./RESmass.^0.5.*(s-RESmass)./((s-RESmass).^2+RESmass.*REStotwid.^2);
RES=-3.0*s/alpha*sum(resi);

s1=1.0D4;
highene=-alpha/3.0/pi*RpQCD(s1)*log(abs(s-s1)/s1);
singu=s*Sigtot(n)/4/pi^2/alpha*log(abs((s1-s)/(4*PImass^2-s)));


integ=(Sigtot-Sigtot(n))./(Stot-Stot(n));
integ(n)=0;
numint=trapz(Stot,integ)*Stot(n)/4/pi/pi/alpha;
%NumInt=integral(integ,4*PImass^2,s1,'AbsTol',10e-12);
% PIH=1i*s/4/pi*sigma(s)/389379;
PIH=RES+highene+(singu+numint+1i*s/4/pi/alpha*Sigtot(n))/389379000;
end