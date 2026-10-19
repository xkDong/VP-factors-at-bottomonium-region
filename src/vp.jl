using LinearAlgebra
include("vp_pqcd.jl")
export VPKernel, vp_replay, vp_literal, vp_grid, vp_smooth
const VP_ALPHA=1/137.035999074
const VP_MPI=.13957018
function vp_mu(s)
    beta=sqrt(s/4-.105658^2)/sqrt(s/4)
    4pi/137.036^2/3/s*beta*(1.5-.5beta^2)*389379000
end
function vp_pion(s)
    beta=sqrt(max(0.,1-4VP_MPI^2/s))
    f=1+(.439/.197327^2)*s/6+6.8s^2-.7s^3
    f^2*pi*VP_ALPHA^2*beta^3/3/s/vp_mu(s)*380379000
end
function vp_lepton(s,m2)
    if s<4m2
        a=sqrt(4m2/s-1)
        return complex(VP_ALPHA/pi*(8/9+a*a/3-2*(.5+a*a/6)*a*atan(1/a)))
    end
    a=sqrt(1-4m2/s);b=(1-a)/(1+a)
    complex(VP_ALPHA/pi*(8/9-a*a/3+(.5-a*a/6)*a*log(b)),a*VP_ALPHA/3*(1+2m2/s))
end
function vp_fixed(s)
    m2=(9.590888,13.587399,89.491600,100.465741,107.230167)
    widths=(.0929e-3,.304e-3,.05402e-3,.0319e-3,.02032e-3)
    ee=(5.55e-6*.95,2.35e-6*.95,1.34e-6*.93,.612e-6*.93,.443e-6*.93)
    res=-3s/VP_ALPHA*sum(ee[j]/sqrt(m2[j])*(s-m2[j])/((s-m2[j])^2+m2[j]*widths[j]^2) for j in 1:5)
    res-VP_ALPHA/(3pi)*r_pqcd(1e4)*log(abs(s-1e4)/1e4)
end
function vp_refine(x,k=5)
    v=zeros((length(x)-1)*k+1)
    for i in 1:length(x)-1
        dl=(x[i+1]-x[i])/k
        for j in 0:k-1;v[(i-1)*k+j+1]=x[i]+dl*j;end
    end
    v[end]=x[end];v
end
function refine_matrix(n,k=5)
    A=zeros((n-1)*k+1,n)
    for i in 1:n-1, j in 0:k-1
        A[(i-1)*k+j+1,i]=1-j/k
        A[(i-1)*k+j+1,i+1]=j/k
    end
    A[end,end]=1;A
end
vp_smooth(energy,draws)=smoothing_matrix(energy;first=22,span=.045)*draws
function vp_grid(epdg,ebes,ebb)
    pp=collect(2VP_MPI:.001:.36)
    pdg=vp_refine(epdg);bes=vp_refine(ebes);bb=vp_refine(ebb)
    q1=collect(2.02:.02:3.68);q2=collect(round(last(ebes)+.05;digits=8):.05:9.0)
    q4=collect(9.01:.01:10.55);q3=vcat(collect(11.25:.05:15),collect(15.5:.5:100))
    e=vcat(pp,pdg,q1,bes,q2,q4,bb,q3)
    n=length(e);background=zeros(n)
    ipdg=length(pp).+(1:length(pdg))
    ibes=(last(ipdg)+length(q1)).+(1:length(bes))
    ibb=(last(ibes)+length(q2)+length(q4)).+(1:length(bb))
    background[1:length(pp)]=vp_pion.(pp.^2)
    background[last(ipdg)+1:first(ibes)-1]=r_pqcd.(q1.^2)
    background[last(ibes)+1:first(ibb)-1]=r_pqcd.(vcat(q2,q4).^2)
    background[ibb]=vp_refine(r_pqcd.(ebb.^2))
    background[last(ibb)+1:end]=r_pqcd.(q3.^2)
    (;energy=e,ipdg,ibes,ibb,background)
end
struct VPKernel
    grid::Vector{Float64}
    ibb::UnitRange{Int}
    mu::Vector{Float64}
    constant::Vector{ComplexF64}
    low::Matrix{Float64}
    charm::Matrix{Float64}
    bottom::Matrix{Float64}
    local_imag::Vector{Float64}
    smooth::Matrix{Float64}
    expand::Matrix{Float64}
    background::Vector{Float64}
    ipdg::UnitRange{Int}
    ibes::UnitRange{Int}
end
function VPKernel(epdg,ebes,ebb)
    g=vp_grid(epdg,ebes,ebb);e=g.energy;s=e.^2;n=length(s);b=g.ibb
    all(diff(e).>0) || error("Invalid VP grid")
    w=vcat((s[2]-s[1])/2,(s[3:end].-s[1:end-2])/2,(s[end]-s[end-1])/2)
    K=zeros(length(b),n)
    for (row,i) in enumerate(b)
        for j in 1:n
            j==i && continue
            K[row,j]=w[j]/(s[j]-s[i])
        end
        K[row,i]=log(abs((1e4-s[i])/(4VP_MPI^2-s[i])))-sum(K[row,:])
        K[row,:].*=s[i]/(4pi^2*VP_ALPHA*389379000)
    end
    mu=vp_mu.(s)
    K .*= transpose(mu)
    leptons=[sum(vp_lepton(s[i],m^2) for m in (.000510998928,.1056583715,1.77682)) for i in b]
    imK=s[b]./(4pi*389379000).*mu[b]
    constant=1 .+leptons.+vp_fixed.(s[b]).+K*g.background.+im.*imK.*g.background[b]
    VPKernel(e,b,mu,constant,K[:,g.ipdg]*refine_matrix(length(epdg)),
        K[:,g.ibes]*refine_matrix(length(ebes)),K[:,b],imK,
        smoothing_matrix(ebb;first=22),refine_matrix(length(ebb)),g.background,g.ipdg,g.ibes)
end
"""Three VP iterations for historical toys, with only bottom R updated."""
function vp_replay(k::VPKernel,pdg,bes,bb;iterations=3,chunk=128)
    iterations>=1 && chunk>=1 || throw(ArgumentError("Positive iterations and chunk required"))
    n=size(bb,2);size(pdg,2)==size(bes,2)==n || throw(DimensionMismatch())
    result=zeros(size(bb));smooth=k.smooth*bb
    for first in 1:chunk:n
        ix=first:min(first+chunk-1,n)
        original=k.expand*smooth[:,ix]; born=copy(original)
        base=k.constant .+ k.low*pdg[:,ix] .+ k.charm*bes[:,ix]
        re=similar(original);f=similar(original)
        for iteration in 1:iterations
            mul!(re,k.bottom,born)
            @. f=inv(abs2(base+re+im*k.local_imag*born))
            @. born=original/f
        end
        result[:,ix]=f[1:5:end,:]
    end
    result
end
"""Literal PIH.m trapezoid sum; diagnostic reference for selected samples."""
function vp_literal(k,pdg,bes,bb;iterations=3)
    original=k.expand*(k.smooth*bb)
    r=copy(k.background)
    r[k.ipdg].+=vp_refine(pdg);r[k.ibes].+=vp_refine(bes)
    r[k.ibb].+=original
    s=k.grid.^2;n=length(s);history=zeros(length(k.ibb),iterations)
    for t in 1:iterations
        sigma=r.*k.mu
        for (j,i) in enumerate(k.ibb)
            integ=(sigma.-sigma[i])./(s.-s[i]);integ[i]=0
            integral=sum((integ[1:end-1].+integ[2:end]).*diff(s))/2
            sing=s[i]*sigma[i]/(4pi^2*VP_ALPHA)*log(abs((1e4-s[i])/(4VP_MPI^2-s[i])))
            pih=vp_fixed(s[i])+(sing+integral*s[i]/(4pi^2*VP_ALPHA)+im*s[i]/(4pi)*sigma[i])/389379000
            pil=sum(vp_lepton(s[i],m^2) for m in (.000510998928,.1056583715,1.77682))
            history[j,t]=inv(abs2(1+pil+pih))
        end
        r[k.ibb]=k.background[k.ibb]+original./history[:,t]
    end
    history[1:5:end,:]
end
