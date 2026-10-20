# Article-method revision: robust LOWESS, iterative ISR, Gaussian likelihood summaries.
export robust_lowess, paper_isr, gaussian_summary
function robust_lowess(e,y;start=10.6659,span=.045,robust_steps=3)
 first=searchsortedfirst(e,start);x=e[first:end];n=length(x);nnear=clamp(ceil(Int,span*n),3,n)
 neighbourhoods=Vector{Vector{Int}}(undef,n);weights=Vector{Vector{Float64}}(undef,n)
 for i in 1:n
  dist=abs.(x.-x[i]);h=partialsort(copy(dist),nnear);ids=findall(dist.<h)
  neighbourhoods[i]=ids;weights[i]=(1 .-(dist[ids]./h).^3).^3
 end
 out=copy(y)
 for col in axes(y,2)
  z=view(y,first:length(e),col);rw=ones(n);fit=zeros(n)
  for step in 0:robust_steps
   for i in 1:n
    ids=neighbourhoods[i];w=weights[i].*rw[ids];dx=x[ids].-x[i]
    a=sum(w);b=dot(w,dx);c=dot(w,dx.^2);det=a*c-b*b
    if det>eps(Float64)*max(a*c,1e-30)
     fit[i]=dot(w.*(c.-b.*dx),z[ids])/det
    else
     # Degenerate robust support: local ordinary LOWESS fallback, explicitly counted by tests.
     w=weights[i];a=sum(w);b=dot(w,dx);c=dot(w,dx.^2)
     fit[i]=dot(w.*(c.-b.*dx),z[ids])/(a*c-b*b)
    end
   end
   step==robust_steps && break
   residual=abs.(z.-fit);sr=sort(residual);scale=(sr[(n+1)÷2]+sr[(n+2)÷2])/2
   scale<=eps(Float64)*max(maximum(abs.(z)),1.) && break
   @. rw=(1-min(residual/(6scale),1)^2)^2
  end
  out[first:end,col]=fit
 end
 out
end
function paper_isr(e,raw;rtol=1e-4,maxiter=500,order=48)
 rawR=raw.*(e.^2 ./86.85);mu=updated_mu.(e.^2)./1000
 observed=rawR.*mu;y=robust_lowess(e,observed);W=updated_isr_matrix(e;order)
 # Signed Jacobi fixed-point iteration for the same piecewise-linear convolution.
 d=copy(y);d[1,:].=0.;f=ones(size(y));iterations=zeros(Int,size(y,2));errfinal=zeros(size(y,2));active=trues(size(y,2))
 for it in 1:maxiter
  ids=findall(active);isempty(ids) && break
  old=d[:,ids];next=old.+(y[:,ids].-W*old)./reshape(vcat(1.,diag(W)[2:end]),:,1);next[1,:].=0.
  fn=ones(size(next));fn[2:end,:]=next[2:end,:]./y[2:end,ids]
  fo=f[:,ids];err=vec(maximum(abs.(fn.-fo)./max.(abs.(fn),1e-12);dims=1))
  d[:,ids]=next;f[:,ids]=fn
  for (j,id) in enumerate(ids)
   if err[j]<rtol;active[id]=false;iterations[id]=it;errfinal[id]=err[j];end
  end
 end
 any(active) && error("ISR convergence failed")
 (;dressed=d./mu,inverse=f,iterations,relative_change=maximum(errfinal),residual=maximum(abs.(W*d.-y)),W,smoothed_sigma=y,mu)
end
function gaussian_summary(a)
 # Closed-form unbinned Gaussian maximum-likelihood fit, without assuming normality was established.
 n=size(a,2);m=vec(sum(a;dims=2))/n;z=a.-m;v=vec(sum(z.^2;dims=2))/n;sd=sqrt.(v)
 skew=vec(sum(z.^3;dims=2))/n./max.(sd.^3,eps(Float64)^3)
 kurt=vec(sum(z.^4;dims=2))/n./max.(sd.^4,eps(Float64)^4).-3
 cov=z*transpose(z)/n
 (;mean=m,sd,covariance=cov,skew,kurt,n)
end
