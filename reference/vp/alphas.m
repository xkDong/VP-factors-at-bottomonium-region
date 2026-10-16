function alphas=alphas(s)
Cmass=1.85;
Bmass=5.6;
TAUmass=1.7769;
alphas0=0.327;%at S=mtau^2
% alphas0=0.3;%at S=mtau^2
a0=alphas0/pi;
zeta3=1.202056903;
if s<4.0*Cmass^2
    nf=3.0;
else if s<4*Bmass^2
        nf=4.0;
    else
        nf=5.0;
    end
end

beta0=0.25*(11.0-2.0/3.0*nf);
beta1=1.0/16.0*(102.0-38.0/3.0*nf);
beta2=1.0/64.0*(2857.0/2.0-5033.0/18.0*nf+325.0/54.0*nf^2);
beta3=1.0/256.0*(149753.0/6.0+3564.0*zeta3-(1078361.0/162.0+6508.0/27.0*zeta3)*nf ...
    +(50065.0/162.0+6472.0/81.0*zeta3)*nf^2+1093.0/729.0*nf^3);

K=1+beta0*a0*log(s/TAUmass^2);
L=log(K);
b1=beta1/beta0;
b2=beta2/beta0;
b3=beta3/beta0;

a2mu=a0/(K+b1*a0*L+b1^2*a0^2*(1-K+L)/K);
c2mu=b2*(1-K);
c3mu=b3/2.0*(1-K^2)+b1*b2*K*(K-1-L)+b1^3/2.0*(L^2-(1-K)^2);
as=a2mu*(1+c2mu*a2mu^2+c3mu*a2mu^3);
alphas=as*pi;
end