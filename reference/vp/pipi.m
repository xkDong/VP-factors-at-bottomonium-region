function pipi=pipi(s)
global mpi alpha;
         r2pi = 0.439D0/0.197327^2;    % !in GeV^-2, Table 3 of EPJC27,497
         c1   = 6.8D0;
         c2   =-0.7D0;
         mpi = 0.13957018D0;
         beta = (1D0-4D0*mpi^2./s).^0.5;
         Fpi  = 1D0 + r2pi*s/6D0 + c1*s.^2 + c2*s.^3 ;         %!Eq.29 of EPJC27,497
         XSPP = Fpi.^2*pi*alpha^2.*beta.^3/3D0./s ; %!Eq.28 of EPJC27,497
         pipi=XSPP./sigmu(s)*380379000;
end