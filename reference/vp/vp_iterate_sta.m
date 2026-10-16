clear
clc

load 'inputR_for_MC_sta_err_0320.mat';%Ebb,Ebes,Epdg;Rbb,Rbes,Rpdg
global mpi alpha Sigtot Stot;
mpi=0.13957018D0;
alpha=1.D0/137.035999074;

for i=1:10000
    Rbb(i,22:end)=smooth(Ebb(22:end),Rbb(i,22:end),0.045,'lowess');
end
Epipi=2*mpi:0.001:0.36;
Rpipi=pipi(Epipi.^2);

EPDG=Newenergy(Epdg,5);
EBES=Newenergy(Ebes,5);
EBB=Newenergy(Ebb,5);

EpQCD1=2.02:0.02:3.68;
RpQCD1=zeros(1,length(EpQCD1));
for i=1:length(EpQCD1)
    RpQCD1(i)=RpQCD(EpQCD1(i)^2);
end

EpQCD2=4.65:0.05:9.0;
RpQCD2=zeros(1,length(EpQCD2));
for i=1:length(EpQCD2)
    RpQCD2(i)=RpQCD(EpQCD2(i)^2);
end

EpQCD4=9.01:0.01:10.55;
RpQCD4=zeros(1,length(EpQCD4));
for i=1:length(EpQCD4)
    RpQCD4(i)=RpQCD(EpQCD4(i)^2);
end

EpQCD3=[11.25:0.05:15,15.5:0.5:100];
RpQCD3=zeros(1,length(EpQCD3));
for i=1:length(EpQCD3)
    RpQCD3(i)=RpQCD(EpQCD3(i)^2);
end

RpQCDb=zeros(1,length(Ebb));
for i=1:length(Ebb)
    RpQCDb(i)=RpQCD(Ebb(i)^2);
end

Etot=[Epipi,EPDG,EpQCD1,EBES,EpQCD2,EpQCD4,EBB,EpQCD3];
Stot=Etot.^2;

Emass2=0.000510998928^2;
MUmass2=0.1056583715^2;
TAUmass2=1.77682^2;

global PILemutau;
PILemutau=zeros(1,length(Etot));
for i=1:length(Etot)
PILemutau(i)=PIL(Stot(i),Emass2)+PIL(Stot(i),MUmass2)+PIL(Stot(i),TAUmass2);
end


a=length([Epipi,EPDG,EpQCD1,EBES,EpQCD2,EpQCD4]);
b=a+length(EBB);
sigmamu=sigmu(Etot.^2);
% PIS1=zeros(length(Etot),3);
% PIS2=zeros(10000,length(Etot),3);
vp1=zeros(length(Etot),3);
vp2=zeros(10000,length(Etot),3);
RpQCDB=Newenergy(RpQCDb,5);
RBB=Newenergy(Rbb(1,:),5);
RBtot=RpQCDB+RBB;
RBES=Newenergy(Rbes(1,:),5);
RPDG=Newenergy(Rpdg(1,:),5);
Rtot=[Rpipi,RPDG,RpQCD1,RBES,RpQCD2,RpQCD4,RBtot,RpQCD3];
Sigtot=Rtot.*sigmamu;
sigmuBB=sigmu(EBB.^2);
SigpQCDB=RpQCDB.*sigmuBB;
SigBB=RBB.*sigmuBB;
for j=1:5
%     if j==1
%         c=2;
%         d=length(Etot);
%     else
%         c=a+1;
%         d=b;
%     end
    for k=2:length(Etot)-1
        vp1(k,j)=vacuum(Stot(k),k);
    end
    Sigtot(a+1:b)=SigpQCDB+SigBB./vp1(a+1:b,j)';
end

for i=1:10000
    RBB=Newenergy(Rbb(i,:),5);
    RBtot=RpQCDB+RBB;
    RBES=Newenergy(Rbes(i,:),5);
    RPDG=Newenergy(Rpdg(i,:),5);
    Rtot=[Rpipi,RPDG,RpQCD1,RBES,RpQCD2,RpQCD4,RBtot,RpQCD3];
    Sigtot=Rtot.*sigmamu;
    SigBB=RBB.*sigmuBB;
    for j=1:3
        for k=4145:5419
        vp2(i,k,j)=vacuum(Stot(k),k);
        end
        Sigtot(a+1:b)=SigpQCDB+SigBB./(vp2(i,a+1:b,j));
    end
end
save 'vp_smo_sta_err_0320.mat';