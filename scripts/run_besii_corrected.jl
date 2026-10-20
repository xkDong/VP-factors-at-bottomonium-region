using PaperISR, LinearAlgebra, TOML, SHA
BLAS.set_num_threads(4)
function reportwrite(out,r)
 open(joinpath(out,"validation.toml"),"w") do io;TOML.print(io,r);end
end
function main()
 root=dirname(@__DIR__);out=joinpath(root,"results/besii_corrected")
 r=Dict{String,Any}("completed"=>false,"scope"=>"numerical-updates kernels + published BES input; total input-error ensemble means, not nominal values", "numerical_base"=>"805455f", "central_definition"=>"Means of 9999 paired random samples; first fixed bottom sample excluded. Fixed-input checks use background ensemble mean and published BES central values.")
 reportwrite(out,r)
 e=vec(read_csv(joinpath(root,"data/vp/Ebb.csv")));ep=vec(read_csv(joinpath(root,"data/vp/Epdg.csv")));ec=vec(read_csv(joinpath(root,"data/vp/Ebes.csv")))
 raw=load_draws(joinpath(root,"data/total_draws.f64"),length(e))[:,2:end]
 pdg=load_draws(joinpath(root,"data/vp/Rpdg.f64"),length(ep))[:,2:end]
 bes=load_draws(joinpath(root,"data/vp/Rbes.f64"),length(ec))[:,2:end]
 input=read_csv(joinpath(root,"data/total_input.csv"));p0=reshape(vec(sum(pdg;dims=2))/size(pdg,2),:,1);c0=reshape(read_csv(joinpath(root,"data/besii/input.csv"))[:,2],:,1)
 a=updated_isr(e,raw);a0=updated_isr(e,reshape(input[:,2],:,1))
 r["ISR_forward_residual_nb"]=a.residual;@assert a.residual<1e-10
 r["ISR_order48_96_max_weight_difference"]=maximum(abs.(a.W.-PaperISR.updated_isr_matrix(e;order=96)));@assert r["ISR_order48_96_max_weight_difference"]<1e-10
 println("ISR done; building corrected VP kernel");flush(stdout)
 k=updated_vp_kernel(ep,ec,e);v0=updated_vp(k,p0,c0,a0.dressed)
 println("Propagating 9999 paired samples on refined grid");flush(stdout)
 v=updated_vp(k,pdg,bes,a.dressed;full_grid=true);step=20
 factor=v.factor[1:step:end,:];dressed=raw.*a.inverse;born=dressed./factor
 r["VP_residual"]=v.residual;r["VP_iterations"]=v.iterations;@assert v.residual<1e-10
 sumry(x)=covariance_summary(x;exclude_nominal=false)
 ii=sumry(a.inverse);vv=sumry(factor);dd=sumry(dressed);bb=sumry(born)
 write_csv(joinpath(out,"results.csv"),["energy_GeV","ISR_correction","ISR_total_error","VP_factor","VP_total_error","sigma_dressed_nb","sigma_dressed_error_nb","sigma_Born_nb","sigma_Born_error_nb"],hcat(e[2:end],ii.mean[2:end],ii.sd[2:end],vv.mean[2:end],vv.sd[2:end],dd.mean[2:end],dd.sd[2:end],bb.mean[2:end],bb.sd[2:end]))
 write_csv(joinpath(out,"fixed_inputs.csv"),["energy_GeV","ISR_fixed_input","VP_fixed_input"],hcat(e[2:end],a0.inverse[2:end,1],v0.factor[2:end,1]))
 for (name,z) in (("isr",ii),("vp",vv),("sigma_dressed",dd),("sigma_born",bb))
  @assert maximum(abs.(diag(z.covariance).-z.sd.^2))<1e-11
  @assert eigmin(Symmetric(z.covariance[2:end,2:end])) > -1e-10
  write_csv(joinpath(out,name*"_covariance.csv"),string.(e[2:end]),z.covariance[2:end,2:end])
 end
 println("Evaluating VP below bottom threshold");flush(stdout)
 ix=findall(t->10.4-1e-10<=t<e[1],k.grid);s=k.grid.^2
 K=reduce(vcat,[transpose(PaperISR.analytic_pv_weights(s,s[i]).*(s[i]/(4pi^2*PaperISR.UPDATE_ALPHA*PaperISR.UPDATE_CONV)).*k.mu) for i in ix])
 fixed=[PaperISR.vp_fixed(z)+PaperISR.UPDATE_ALPHA/(3pi)*PaperISR.r_pqcd(1e4)*log(abs(z-1e4)/1e4)-PaperISR.UPDATE_ALPHA/(3pi)*PaperISR.r_pqcd(200.0^2)*log1p(-z/200.0^2) for z in s[ix]]
 imK=s[ix]./(4pi*PaperISR.UPDATE_ALPHA*PaperISR.UPDATE_CONV).*k.mu[ix]
 constant=1 .+[sum(PaperISR.updated_lepton(z,m) for m in (PaperISR.UPDATE_ME,.1056583715,1.77682)) for z in s[ix]].+fixed.+K*k.background.+im.*imK.*k.background[ix]
 L=K[:,k.ipdg]*PaperISR.refine_matrix(length(ep),20);C=K[:,k.ibes]*PaperISR.refine_matrix(length(ec),20);B=K[:,k.ibb]
 below=zeros(length(ix),size(raw,2))
 for j in 1:128:size(raw,2)
  ids=j:min(j+127,size(raw,2));br=(k.expand*a.dressed[:,ids])./v.factor[:,ids]
  below[:,ids]=inv.(abs2.(constant.+L*pdg[:,ids].+C*bes[:,ids].+B*br))
 end
 # Same plotted sampling as the original: five subintervals, evaluated from the converged refinement-20 solution.
 dense=vcat(below,v.factor[1:4:end,:]);de=vcat(k.grid[ix],k.grid[k.ibb][1:4:end]);dm,ds=ensemble_summary(dense)
 write_csv(joinpath(out,"vp_dense.csv"),["energy_GeV","VP_factor","VP_total_error"],hcat(de,dm,ds))
 @assert maximum(abs.(dm[length(ix)+6:5:end].-vv.mean[2:end]))<1e-12
 r["dense_points"]=length(de);r["samples"]=size(raw,2);r["pQCD_start_GeV"]=round(last(ec)+.05;digits=8)
 # Compare means consistently using the SAME 9999 input samples in the historical algorithm.
 println("Computing paired historical baseline");flush(stdout)
 ecold=vec(read_csv(joinpath(root,"reference/legacy_charm/Ebes.csv")));cold=load_draws(joinpath(root,"reference/legacy_charm/Rbes.f64"),length(ecold))[:,2:end]
 oi=replay(Kernel(e),raw);od=raw.*oi;ov=vp_replay(VPKernel(ep,ecold,e),pdg,cold,od.*(e.^2 ./86.85));ob=od./ov
 cols=[e[2:end]];names=["energy_GeV"]
 for (name,new,old) in (("ISR",a.inverse,oi),("VP",factor,ov),("dressed",dressed,od),("Born",born,ob))
  nm,ns=ensemble_summary(new);om,os=ensemble_summary(old);delta,deltasd=ensemble_summary(new.-old)
  change=100 .* (nm[2:end]./om[2:end].-1)
  append!(names,[name*"_old_mean",name*"_new_mean",name*"_change_percent",name*"_paired_shift_MC_SE"])
  append!(cols,[om[2:end],nm[2:end],change,deltasd[2:end]./sqrt(size(raw,2))])
  r[name*"_max_abs_change_percent"]=maximum(abs.(change))
 end
 write_csv(joinpath(out,"comparison.csv"),names,hcat(cols...))
 println("Checking resolution variations");flush(stdout)
 for (tag,kw) in (("grid30",(;refinement=30)),("background4",(;background_refinement=4)),("tail400",(;upper=400.)),("tailfine",(;tail_step=.0625)))
  kk=updated_vp_kernel(ep,ec,e;kw...);z=updated_vp(kk,p0,c0,a0.dressed)
  r[tag*"_max_VP_abs"]=maximum(abs.(z.factor.-v0.factor));kk=nothing;GC.gc()
  println(tag," ",r[tag*"_max_VP_abs"]);flush(stdout)
 end
 r["resolution_budget_abs"]=sum(r[t*"_max_VP_abs"] for t in ("grid30","background4","tail400","tailfine"));@assert r["resolution_budget_abs"]<1e-6
 r["tolerance_change_abs"]=maximum(abs.(updated_vp(k,p0,c0,a0.dressed;rtol=1e-13).factor.-v0.factor));@assert r["tolerance_change_abs"]<1e-10
 r["negative_estimates_retained"]=a.negative_count;r["uncertainty_scope"]="Propagated historical input-error model; excludes smoothing/model and missing parameter uncertainties; no claim of a complete statistical/systematic split"
 r["hashes"]=Dict(f=>bytes2hex(sha256(read(joinpath(root,f)))) for f in ["src/updated.jl","scripts/run_besii_corrected.jl","data/total_draws.f64","data/vp/Rpdg.f64","data/vp/Rbes.f64","data/besii/input.csv"])
 r["completed"]=true;reportwrite(out,r);println(r)
end
main()
