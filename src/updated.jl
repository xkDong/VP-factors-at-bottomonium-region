# Updated numerical model. Historical functions above remain unchanged.
export updated_isr, updated_vp_kernel, updated_vp, updated_samples, covariance_summary
const UPDATE_ALPHA=VP_ALPHA
const UPDATE_ME=.000510998928
const UPDATE_CONV=389379000.0

updated_mu(s)=4pi*UPDATE_ALPHA^2/(3s)*sqrt(1-4*.105658^2/s)*(1+2*.105658^2/s)*UPDATE_CONV
function gauss_rule(n)
    # Golub-Welsch, standard library only.
    eig=eigen(SymTridiagonal(zeros(n),[i/sqrt(4i^2-1) for i in 1:n-1]))
    eig.values,2 .*eig.vectors[1,:].^2
end
function updated_radiator(s)
    beta=2UPDATE_ALPHA/pi*(log(s/UPDATE_ME^2)-1)
    delta=UPDATE_ALPHA/pi*(pi^2/3-.5)+.75beta+beta^2*(9/32-pi^2/12)
    beta,delta
end
function regular_radiator(x,b)
    -b*(1-x/2)+b^2/8*(4*(2-x)*(-log(x))-(1+3*(1-x)^2)*log1p(-x)/x-6+x)
end
"""Piecewise-linear sigma(s), soft term integrated analytically, regular term Gaussian."""
function updated_isr_matrix(e;order=48)
    s=e.^2;n=length(s);W=zeros(n,n);nodes,weights=gauss_rule(order)
    for i in 2:n
        beta,delta=updated_radiator(s[i])
        for j in 1:i-1
            xa=1-s[j+1]/s[i];xb=1-s[j]/s[i]
            # Shape functions linear in radiation fraction x.
            for (col,A,B) in ((j,(s[j+1]-s[i])/(s[j+1]-s[j]),s[i]/(s[j+1]-s[j])),
                             (j+1,(s[i]-s[j])/(s[j+1]-s[j]),-s[i]/(s[j+1]-s[j])))
                v=(1+delta)*(A*(xb^beta-xa^beta)+B*beta/(beta+1)*(xb^(beta+1)-xa^(beta+1)))
                for q in eachindex(nodes)
                    u=(nodes[q]+1)/2
                    x=xa==0 ? xb*u^4 : (xa+xb)/2+(xb-xa)/2*nodes[q]
                    jac=xa==0 ? 2xb*u^3 : (xb-xa)/2
                    v+=weights[q]*jac*regular_radiator(x,beta)*(A+B*x)
                end
                W[i,col]+=v
            end
        end
    end
    W
end
function energy_smoothing(e;start=10.6659,span=.045)
    first=searchsortedfirst(e,start)
    first<=length(e)-2 || throw(ArgumentError("Too few points above smoothing start"))
    smoothing_matrix(e;first,span)
end
"""Returns a signed linear estimator: never clip or remove negative unfolded draws."""
function updated_isr(e,legacy_sigma;start=10.6659,span=.045,order=48)
    # Stored samples encode R via 86.85/s. Reconstruct R before using a consistent mu cross section.
    rawR=legacy_sigma.*(e.^2 ./86.85)
    mu=updated_mu.(e.^2)./1000
    observed=rawR.*mu
    S=energy_smoothing(e;start,span)
    y=S*observed;W=updated_isr_matrix(e;order)
    d=zeros(size(y));d[2:end,:]=LowerTriangular(W[2:end,2:end])\y[2:end,:]
    inverse=ones(size(y));inverse[2:end,:]=d[2:end,:]./y[2:end,:]
    # Consistent single-smoothed estimator, not inverse factor times unsmoothed data.
    R=d./mu
    residual=maximum(abs.(W*d.-y))
    (;dressed=R,inverse,observed_R=rawR,residual,
      negative_count=count(<(0),R[2:end,:]),W,smoothed_sigma=y,mu)
end
function updated_lepton(s,m)
    a=sqrt(1-4m^2/s)
    l=log(4m^2/s)-2log1p(a) # stable for electron mass at bottomonium energies
    complex(UPDATE_ALPHA/pi*(8/9-a*a/3+(.5-a*a/6)*a*l),a*UPDATE_ALPHA/3*(1+2m^2/s))
end
"""Weights for analytic PV integral of a continuous piecewise-linear function in s."""
function analytic_pv_weights(x,s)
    x[1]<s<x[end] || throw(DomainError(s))
    w=zeros(length(x));cf=log((x[end]-s)/(s-x[1]))
    for j in 1:length(x)-1
        a=x[j];b=x[j+1];h=b-a
        w[j]-=1;w[j+1]+=1
        if !(a<=s<=b)
            L=log1p(h/(a-s))
            w[j]+=(1-(s-a)/h)*L;w[j+1]+=(s-a)/h*L
            cf-=L
        end
    end
    i=clamp(searchsortedlast(x,s),1,length(x)-1);f=(s-x[i])/(x[i+1]-x[i])
    w[i]+=cf*(1-f);w[i+1]+=cf*f
    w
end
function updated_grid(ep,ec,eb;refinement=20,upper=200.,tail_step=.125,background_refinement=2)
    pp=collect(2VP_MPI:(.001/background_refinement):.36)
    pdg=vp_refine(ep,refinement);bes=vp_refine(ec,refinement);bb=vp_refine(eb,refinement)
    q1=collect(2.02:(.02/background_refinement):3.68);q2=collect(4.65:(.05/background_refinement):9.)
    q4=collect(9.01:(.01/background_refinement):10.55)
    q3=vcat(collect(11.25:(.05/background_refinement):15),collect(15+tail_step:tail_step:upper))
    e=vcat(pp,pdg,q1,bes,q2,q4,bb,q3)
    ip=length(pp).+(1:length(pdg));ic=(last(ip)+length(q1)).+(1:length(bes))
    ib=(last(ic)+length(q2)+length(q4)).+(1:length(bb))
    bg=zeros(length(e))
    for i in eachindex(pp)
        s=pp[i]^2;f=1+(.439/.197327^2)*s/6+6.8s^2-.7s^3
        sig=f^2*pi*UPDATE_ALPHA^2*max(0.,1-4VP_MPI^2/s)^1.5/(3s)*UPDATE_CONV
        bg[i]=sig/updated_mu(s)
    end
    bg[last(ip)+1:first(ic)-1]=r_pqcd.(q1.^2)
    bg[last(ic)+1:first(ib)-1]=r_pqcd.(vcat(q2,q4).^2)
    bg[ib]=vp_refine(r_pqcd.(eb.^2),refinement)
    bg[last(ib)+1:end]=r_pqcd.(q3.^2)
    (;energy=e,ip,ic,ib,bg)
end
function updated_vp_kernel(ep,ec,eb;refinement=20,upper=200.,tail_step=.125,background_refinement=2,
        second_smoothing=false,start=10.6659,span=.045)
    g=updated_grid(ep,ec,eb;refinement,upper,tail_step,background_refinement)
    e=g.energy;s=e.^2;b=g.ib
    K=zeros(length(b),length(e))
    for (j,i) in enumerate(b)
        K[j,:]=analytic_pv_weights(s,s[i]).*(s[i]/(4pi^2*UPDATE_ALPHA*UPDATE_CONV))
    end
    mu=updated_mu.(s);K .*=transpose(mu)
    imK=s[b]./(4pi*UPDATE_ALPHA*UPDATE_CONV).*mu[b]
    leptons=[sum(updated_lepton(s[i],m) for m in (UPDATE_ME,.1056583715,1.77682)) for i in b]
    # Same narrow-resonance convention as baseline; high-energy constant-R tail at chosen cutoff.
    fixed=[vp_fixed(v)+UPDATE_ALPHA/(3pi)*r_pqcd(1e4)*log(abs(v-1e4)/1e4)-
        UPDATE_ALPHA/(3pi)*r_pqcd(upper^2)*log1p(-v/upper^2) for v in s[b]]
    constant=1 .+leptons.+fixed.+K*g.bg.+im.*imK.*g.bg[b]
    S=second_smoothing ? energy_smoothing(eb;start,span) : Matrix{Float64}(I,length(eb),length(eb))
    VPKernel(e,b,mu,constant,K[:,g.ip]*refine_matrix(length(ep),refinement),
        K[:,g.ic]*refine_matrix(length(ec),refinement),K[:,b],imK,
        S,refine_matrix(length(eb),refinement),g.bg,g.ip,g.ic)
end
"""Converge VP to a tolerance; return factor and self-consistency residual for every batch."""
function updated_vp(k,pdg,bes,bb;rtol=1e-11,maxiter=100,chunk=128)
    maxiter>0 && rtol>0 || throw(ArgumentError("Invalid solver tolerances"))
    n=size(bb,2);out=zeros(size(bb));worst=0.;used=0
    step=(length(k.ibb)-1)÷(size(bb,1)-1)
    for first in 1:chunk:n
        ix=first:min(first+chunk-1,n)
        original=k.expand*(k.smooth*bb[:,ix]);born=copy(original)
        base=k.constant .+k.low*pdg[:,ix].+k.charm*bes[:,ix]
        re=similar(original);f=ones(size(original));prev=similar(original);converged=false
        for iteration in 1:maxiter
            prev.=f;mul!(re,k.bottom,born)
            @. f=inv(abs2(base+re+im*k.local_imag*born))
            @. born=original/f
            err=maximum(abs.(f.-prev)./max.(abs.(f),1.))
            if err<rtol
                converged=true;used=max(used,iteration);break
            end
        end
        converged || error("VP failed to converge within $maxiter iterations")
        mul!(re,k.bottom,born)
        residual=maximum(abs.(inv.(abs2.(base.+re.+im.*k.local_imag.*born)).-f))
        worst=max(worst,residual)
        out[:,ix]=f[1:step:end,:]
    end
    (;factor=out,residual=worst,iterations=used)
end
function covariance_summary(samples;exclude_nominal=true)
    x=exclude_nominal ? samples[:,2:end] : samples
    mu,sd=ensemble_summary(x);d=x.-mu
    (;mean=mu,sd,covariance=d*transpose(d)/(size(x,2)-1),mean_se=sd/sqrt(size(x,2)),n=size(x,2))
end
"""Gaussian measurement toys; optional conditional truncation preserves shared nuisances."""
function updated_samples(input;toys=10000,seed=2020,positivity=:none)
    positivity in (:none,:conditional) || throw(ArgumentError("Unknown positivity model"))
    rng=Xoshiro(seed);n=size(input,1);a=zeros(n,toys)
    for t in 1:toys
        common=randn(rng)*input[:,5]+randn(rng)*input[:,6]
        for j in 1:n
            v=input[j,2]+common[j]+randn(rng)*input[j,3]+randn(rng)*input[j,4]
            tries=0
            while positivity==:conditional && v<0
                tries+=1;tries<100000 || error("Positivity conditioning infeasible")
                v=input[j,2]+common[j]+randn(rng)*input[j,3]+randn(rng)*input[j,4]
            end
            a[j,t]=v
        end
    end
    a[:,1]=input[:,2];a
end
