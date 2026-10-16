module PaperISR
using Random
export read_csv, write_csv, load_draws, sample_inputs, Kernel, unfold, unfold_literal, replay, ensemble_summary

function read_csv(path)
    lines=readlines(path)
    reduce(vcat,permutedims.(map(l->parse.(Float64,split(l,',')),lines[2:end])))
end
function write_csv(path,header,a)
    mkpath(dirname(path))
    open(path,"w") do io
        println(io,join(header,","))
        for r in eachrow(a); println(io,join(r,",")); end
    end
end

# g1.m uses MATLAB dilog(x) = Li₂(1-x). Here 0 ≤ x < 0.2.
function li2_small(x)
    z=x; total=0.0
    for k in 1:1000
        term=z/(k*k); next=total+term
        next==total && return total
        total=next; z*=x
    end
    error("dilogarithm series did not converge")
end
spence(x)=x==0 ? pi^2/6 : pi^2/6-log(x)*log1p(-x)-li2_small(x)

function primitive(x,s,A,B)
    beta=.004645639*(log(s*1000000/.26111998)-1)
    x==0 && return A*beta-A*beta^2*spence(0)/2
    lx=log(x); l1=log1p(-x); x2=x*x; x3=x2*x
    dp=.006480359964569+(.75-.541217*beta)*beta
    (B*beta^2*(18*x+18*l1-36*x2*lx+12*x3*lx-36*x*l1-18*x2+2*x3+27*x2*l1-9*x3*l1))/72+
    A*beta*(x-2)^2/4+A*x^beta*(dp+1)-
    (A*beta^2*(10*x+18*l1+16*spence(x)-8*x2*lx-24*x*l1+32*x*lx-x2+6*x2*l1+16*l1*lx))/32+
    B*beta*x2*(x-3)/6+B*beta*x^(beta+1)*(dp+1)/(beta+1)
end

"""Ordinary local-linear LOWESS coefficients, matching smooth(..., .045, 'lowess')."""
function smoothing_matrix(e;first=23,span=.045)
    n=length(e); S=zeros(n,n)
    for i in 1:first-1; S[i,i]=1; end
    x=e[first:end]; k=clamp(ceil(Int,span*length(x)),3,length(x))
    for i in eachindex(x)
        d=abs.(x.-x[i]); h=partialsort(copy(d),k)
        h>0 || error("Degenerate LOWESS neighbourhood")
        dx=x.-x[i]; w=(1 .-min.(d./h,1).^3).^3
        a=sum(w); b=sum(w.*dx); c=sum(w.*dx.^2); det=a*c-b*b
        det>0 || error("Singular LOWESS fit")
        S[first+i-1,first:end]=w.*(c.-b.*dx)./det
    end
    S
end

struct Kernel
    energy::Vector{Float64}
    smoothing::Matrix{Float64}
    left::Matrix{Float64}
    right::Matrix{Float64}
end
function segment(s,i,j,left,right)
    slope=(right-left)/(s[j+1]-s[j])
    A=left+(s[i]-s[j])*slope; B=-s[i]*slope
    primitive(1-s[j]/s[i],s[i],A,B)-primitive(1-s[j+1]/s[i],s[i],A,B)
end
function Kernel(energy;first=23,span=.045)
    e=Float64.(energy); n=length(e)
    all(diff(e).>0) || error("Energy grid must be strictly increasing")
    # Deliberately limited to the verified bottomonium domain.
    0<1-(e[1]/e[end])^2<.2 || error("Outside verified radiator domain")
    L=zeros(n,n); R=zeros(n,n); s=e.^2
    for i in 2:n, j in 1:i-1
        L[i,j]=segment(s,i,j,1.,0.); R[i,j]=segment(s,i,j,0.,1.)
    end
    Kernel(e,smoothing_matrix(e;first,span),L,R)
end

function smoothed(k,y)
    length(y)==length(k.energy) || throw(DimensionMismatch())
    all(isfinite,y) && all(y.>=0) || error("Input must be finite and nonnegative")
    y[1]==0 || error("Threshold anchor must be zero")
    v=k.smoothing*y
    all(v[2:end].>0) || error("Nonpositive smoothed cross section")
    v
end
"""Historical sequential energy iteration: 30 states = 29 updates."""
function unfold(k::Kernel,observed;states=30,literal=false)
    states>=2 || throw(ArgumentError("states must be at least 2"))
    y=smoothed(k,observed); n=length(y)
    history=zeros(n,states); history[:,1]=y
    factor=ones(n); s=k.energy.^2
    for i in 2:n
        older=0.0
        # Same descending segment accumulation as iteration.m.
        for j in i-2:-1:1
            older+=literal ? segment(s,i,j,history[j,states],history[j+1,states]) :
                k.left[i,j]*history[j,states]+k.right[i,j]*history[j+1,states]
        end
        for t in 2:states
            prev=history[i,t-1]; left=history[i-1,t-1]
            adjacent=literal ? segment(s,i,i-1,left,prev) :
                k.left[i,i-1]*left+k.right[i,i-1]*prev
            delta=(older+adjacent)/prev
            isfinite(delta) && delta>0 || error("Invalid ISR iteration at point $i")
            history[i,t]=y[i]/delta
            factor[i]=1/delta
        end
    end
    factor
end
unfold_literal(k,y;states=30)=unfold(k,y;states,literal=true)

function load_draws(path,n)
    bytes=filesize(path)
    bytes%(8*n)==0 || error("Invalid binary length")
    v=Vector{UInt64}(undef,bytes÷8)
    read!(path,v)
    reshape(reinterpret(Float64,ltoh.(v)),n,:)
end
"""Newsigma.m prescription, including independent redraws for negative entries."""
function sample_inputs(input;toys=10000,seed=2020)
    toys>=1 || throw(ArgumentError("toys must be positive"))
    size(input,2)==6 || throw(DimensionMismatch())
    rng=Xoshiro(seed); n=size(input,1); draws=zeros(n,toys)
    for t in 1:toys
        a=randn(rng); b=randn(rng)
        q=randn(rng,n); w=randn(rng,n)
        for j in 1:n
            v=input[j,2]+q[j]*input[j,3]+w[j]*input[j,4]+a*input[j,5]+b*input[j,6]
            while v<0
                v=input[j,2]+sum(randn(rng)*input[j,c] for c in 3:6)
            end
            draws[j,t]=v
        end
    end
    draws[:,1]=input[:,2]
    draws
end
function replay(k,draws;states=30)
    ans=similar(draws)
    Threads.@threads for t in axes(draws,2)
        ans[:,t]=unfold(k,view(draws,:,t);states)
    end
    ans
end
function ensemble_summary(a)
    n=size(a,2); n>1 || error("At least two toys needed")
    mu=vec(sum(a;dims=2))/n
    sd=sqrt.(vec(sum((a.-mu).^2;dims=2))/(n-1))
    mu,sd
end
include("vp.jl")
include("updated.jl")
end
