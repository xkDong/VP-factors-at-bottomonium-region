using PaperISR,LinearAlgebra,TOML
BLAS.set_num_threads(2)
root=dirname(@__DIR__);out=joinpath(root,"results/paper_method_corrected")
function main(root,out)
 inp=read_csv(joinpath(root,"data/total_input.csv"));e=inp[:,1];low=read_csv(joinpath(out,"low_energy_input.csv"));c=read_csv(joinpath(root,"data/besii/input.csv"))
 ep=low[:,1];ec=c[:,1];p=reshape(low[:,2],:,1);b=reshape(c[:,2],:,1)
 a=paper_isr(e,reshape(inp[:,2],:,1));at=paper_isr(e,reshape(inp[:,2],:,1);rtol=1e-8)
 k=updated_vp_kernel(ep,ec,e);v=updated_vp(k,p,b,a.dressed;rtol=1e-4);vt=updated_vp(k,p,b,at.dressed;rtol=1e-11)
 r=Dict{String,Any}("completed"=>false,"nominal_ISR_1e4_iterations"=>maximum(a.iterations),"nominal_ISR_tight_iterations"=>maximum(at.iterations),"nominal_ISR_tolerance_abs"=>maximum(abs.(at.inverse.-a.inverse)),"nominal_VP_tolerance_abs"=>maximum(abs.(vt.factor.-v.factor)))
 for (tag,kw) in (("grid30",(;refinement=30)),("background4",(;background_refinement=4)),("tail400",(;upper=400.)),("tailfine",(;tail_step=.0625)))
  kk=updated_vp_kernel(ep,ec,e;kw...);vv=updated_vp(kk,p,b,at.dressed;rtol=1e-11)
  r[tag]=maximum(abs.(vv.factor.-vt.factor));kk=nothing;GC.gc()
 end
 r["resolution_sum"]=sum(r[t] for t in ("grid30","background4","tail400","tailfine"));@assert r["resolution_sum"]<1e-6
 @assert r["nominal_VP_tolerance_abs"]<1e-5
 r["completed"]=true
 open(joinpath(out,"resolution.toml"),"w") do io;TOML.print(io,r);end
 println(r)
end
main(root,out)
