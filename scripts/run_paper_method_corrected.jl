using PaperISR,LinearAlgebra,TOML,SHA,Random
BLAS.set_num_threads(4)
function main()
 root=dirname(@__DIR__);out=joinpath(root,"results/paper_method_corrected");N=10000
 report=Dict{String,Any}("completed"=>false,"rtol"=>1e-4,"samples"=>N,"seed"=>20261020,"low_background_policy"=>"Separate unresolved low-energy combined input uncertainty; never label it statistical or systematic", "central_estimator"=>"Unbinned Gaussian MLE; normality diagnostics reported separately")
 save()=open(joinpath(out,"validation.toml"),"w") do io;TOML.print(io,report);end
 save()
 input=read_csv(joinpath(root,"data/total_input.csv"));e=input[:,1];low=read_csv(joinpath(out,"low_energy_input.csv"));charm=read_csv(joinpath(root,"data/besii/input.csv"));ep=low[:,1];ec=charm[:,1]
 @assert ep≈vec(read_csv(joinpath(root,"data/vp/Epdg.csv")))
 rng=Xoshiro(20261020);ns=randn(rng,length(e),N);nu=randn(rng,length(e),N);c1=randn(rng,1,N);c2=randn(rng,1,N)
 ds=input[:,3].*ns;dy=input[:,4].*nu.+input[:,5]*c1.+input[:,6]*c2
 lowdelta=low[:,3]*randn(rng,1,N)
 cs=charm[:,3]*randn(rng,1,N);cy=charm[:,4]*randn(rng,1,N)
 # Distinct independent nuisances for stat/sys; each background component fully correlated over energy, per article convention.
 println("Building corrected kernel");flush(stdout);k=updated_vp_kernel(ep,ec,e)
 nominal=paper_isr(e,reshape(input[:,2],:,1));vnom=updated_vp(k,reshape(low[:,2],:,1),reshape(charm[:,2],:,1),nominal.dressed;rtol=1e-4)
 summaries=Dict();totalraw=zeros(0,0);totala=nothing;totalv=nothing
 for mode in ("total","stat","sys","low_background")
  println("Starting ",mode);flush(stdout)
  raw=repeat(input[:,2],1,N);pdg=repeat(low[:,2],1,N);bes=repeat(charm[:,2],1,N)
  if mode in ("total","stat");raw.+=ds;bes.+=cs;end
  if mode in ("total","sys");raw.+=dy;bes.+=cy;end
  if mode in ("total","low_background");pdg.+=lowdelta;end
  a=paper_isr(e,raw);v=updated_vp(k,pdg,bes,a.dressed;rtol=1e-4,full_grid=(mode=="total"))
  fac=mode=="total" ? v.factor[1:20:end,:] : v.factor
  dre=raw.*a.inverse;born=dre./fac
  vals=Dict("isr"=>a.inverse,"vp"=>fac,"sigma_dressed"=>dre,"sigma_born"=>born)
  stats=Dict(name=>gaussian_summary(val) for (name,val) in vals);summaries[mode]=stats
  for (name,st) in stats
   @assert all(isfinite,st.mean) && all(isfinite,st.sd)
   @assert eigmin(Symmetric(st.covariance[2:end,2:end])) > -1e-10
   write_csv(joinpath(out,mode*"_"*name*"_covariance.csv"),string.(e[2:end]),st.covariance[2:end,2:end])
   write_csv(joinpath(out,mode*"_"*name*"_gaussian_fit.csv"),["energy_GeV","mu_MLE","sigma_MLE","skewness","excess_kurtosis"],hcat(e[2:end],st.mean[2:end],st.sd[2:end],st.skew[2:end],st.kurt[2:end]))
  end
  report[mode*"_ISR_max_iterations"]=maximum(a.iterations);report[mode*"_ISR_max_relative_change"]=a.relative_change
  report[mode*"_VP_iterations"]=v.iterations;report[mode*"_VP_fixed_point_residual"]=v.residual
  if mode=="total"
   totala=a;totalv=v;totalraw=raw
   # Use every fourth refined point: five subintervals per experimental interval.
   ids=findall(t->10.4-1e-10<=t<e[1],k.grid);s=k.grid.^2
   K=reduce(vcat,[transpose(PaperISR.analytic_pv_weights(s,s[i]).*(s[i]/(4pi^2*PaperISR.UPDATE_ALPHA*PaperISR.UPDATE_CONV)).*k.mu) for i in ids])
   fixed=[PaperISR.vp_fixed(z)+PaperISR.UPDATE_ALPHA/(3pi)*PaperISR.r_pqcd(1e4)*log(abs(z-1e4)/1e4)-PaperISR.UPDATE_ALPHA/(3pi)*PaperISR.r_pqcd(200.0^2)*log1p(-z/200.0^2) for z in s[ids]]
   imK=s[ids]./(4pi*PaperISR.UPDATE_ALPHA*PaperISR.UPDATE_CONV).*k.mu[ids]
   constant=1 .+[sum(PaperISR.updated_lepton(z,m) for m in (PaperISR.UPDATE_ME,.1056583715,1.77682)) for z in s[ids]].+fixed.+K*k.background.+im.*imK.*k.background[ids]
   L=K[:,k.ipdg]*PaperISR.refine_matrix(length(ep),20);C=K[:,k.ibes]*PaperISR.refine_matrix(length(ec),20);B=K[:,k.ibb];below=zeros(length(ids),N)
   for j in 1:128:N
    ii=j:min(j+127,N);br=(k.expand*a.dressed[:,ii])./v.factor[:,ii]
    below[:,ii]=inv.(abs2.(constant.+L*pdg[:,ii].+C*bes[:,ii].+B*br))
   end
   dense=vcat(below,v.factor[1:4:end,:]);dm,dsd=ensemble_summary(dense);dsd.*=sqrt((N-1)/N)
   write_csv(joinpath(out,"vp_dense.csv"),["energy_GeV","VP_factor","VP_total_error"],hcat(vcat(k.grid[ids],k.grid[k.ibb][1:4:end]),dm,dsd))
   # Check the same random inputs against tighter ISR/VP tolerances.
   sel=collect(1:100:N);at=paper_isr(e,raw[:,sel];rtol=1e-6);vt=updated_vp(k,pdg[:,sel],bes[:,sel],at.dressed;rtol=1e-11)
   report["tight_ISR_max_abs_difference"]=maximum(abs.(at.inverse.-a.inverse[:,sel]))
   report["tight_VP_max_abs_difference"]=maximum(abs.(vt.factor.-fac[:,sel]))
   report["tight_ISR_max_difference_in_SD"]=maximum(abs.(at.inverse[2:end,:].-a.inverse[2:end,sel])./max.(stats["isr"].sd[2:end],1e-9))
   @assert report["tight_VP_max_abs_difference"]<1e-5
  end
  save();println("Finished ",mode);flush(stdout)
 end
 st=summaries["total"];ii=st["isr"];vv=st["vp"];dd=st["sigma_dressed"];bb=st["sigma_born"]
 write_csv(joinpath(out,"results.csv"),["energy_GeV","ISR_correction","ISR_total_error","VP_factor","VP_total_error","sigma_dressed_nb","sigma_dressed_error_nb","sigma_Born_nb","sigma_Born_error_nb"],hcat(e[2:end],ii.mean[2:end],ii.sd[2:end],vv.mean[2:end],vv.sd[2:end],dd.mean[2:end],dd.sd[2:end],bb.mean[2:end],bb.sd[2:end]))
 cols=[e[2:end]];names=["energy_GeV"]
 for q in ("isr","vp","sigma_dressed","sigma_born"),mode in ("stat","sys","low_background","total")
  push!(names,q*"_"*mode*"_error");push!(cols,summaries[mode][q].sd[2:end])
 end
 write_csv(joinpath(out,"error_components.csv"),names,hcat(cols...))
 report["nominal_ISR_iterations"]=maximum(nominal.iterations);report["nominal_VP_iterations"]=vnom.iterations
 report["negative_observed_draws"]=count(<(0),totalraw);report["negative_unfolded_draws"]=count(<(0),totala.dressed)
 report["robust_steps"]=3;report["span"]=.045;report["smoothing_start_GeV"]=10.6659
 report["hashes"]=Dict(f=>bytes2hex(sha256(read(joinpath(root,f)))) for f in ["src/paper_method.jl","src/updated.jl","scripts/run_paper_method_corrected.jl","data/total_input.csv","data/besii/input.csv","results/paper_method_corrected/low_energy_input.csv"])
 report["completed"]=true;save();println(report)
end
main()
