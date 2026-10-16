function vacuum2 = vacuum(s,n)
% alpha=1.D0/137.035999074;
% Emass=0.000510998928;
% MUmass=0.1056583715;
% TAUmass=1.77682;
global PILemutau;

PIs=-PILemutau(n)-PIH(s,n);
% vacuum1=PIs;
vacuum2=1/((1-PIs)*conj((1-PIs)));
%[vacuum1,vacuum2]=[PIs,1/((1-PIs)*conj((1-PIs)))];
 %vacuum=1/((1-PIs)*conj((1-PIs)));
%vacuum=PIL(s,Emass^2)+PIL(s,MUmass^2)+PIL(s,TAUmass^2);
end