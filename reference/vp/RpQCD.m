function RpQCD=RpQCD(s)
%PDG Eq.(9.8,9.9)
% alpha=7.297352533E-3;
% PImsss=0.13957018;
Cmass=1.85;
Bmass=5.62;%avoid the last two points of Belle and BaBar.

if s<4.0*Cmass^2
    nf=3.0;
    R=2.0;
    eta=0;
else if s<4.0*Bmass^2
        nf=4.0;
        R=10.0/3.0;
        eta=2.0/15.0;
    else
        nf=5.0;
        R=11.0/3.0;
        eta=1.0/33.0;
    end
end

a=alphas(s)/pi;%G.Rodrigo hep-ph 9707474
%a=0.1736/pi;
c2=1.9857-0.1152*nf;
c3=-6.63694-1.20013*nf-0.00518*nf^2-1.24*eta;
c4=-156.61+18.775*nf-0.7974*nf^2+0.0215*nf^3+(17.828-0.575*nf)*eta;
RpQCD=R*(1+a+c2*a^2+c3*a^3+c4*a^4);

end