using PaperISR, LinearAlgebra, TOML, SHA
BLAS.set_num_threads(4)
function main()
 root=dirname(@__DIR__);out=joinpath(root,"docs/erratum_besii/support")
 ep=vec(read_csv(joinpath(root,"data/vp/Epdg.csv")));ec=vec(read_csv(joinpath(root,"data/vp/Ebes.csv")));eb=vec(read_csv(joinpath(root,"data/vp/Ebb.csv")))
 pdg=load_draws(joinpath(root,"data/vp/Rpdg.f64"),length(ep))
 bes=load_draws(joinpath(root,"data/vp/Rbes.f64"),length(ec));bb=load_draws(joinpath(root,"data/vp/Rbb.f64"),length(eb))
 k=VPKernel(ep,ec,eb);s=k.grid.^2;n=length(s)
 lowix=findall(i->10.4-1e-10<=k.grid[i]<eb[1],eachindex(k.grid));targets=vcat(lowix,collect(k.ibb))
 w=vcat((s[2]-s[1])/2,(s[3:end].-s[1:end-2])/2,(s[end]-s[end-1])/2)
 K=zeros(length(lowix),n)
 for (row,i) in enumerate(lowix)
  for j in 1:n;j==i && continue;K[row,j]=w[j]/(s[j]-s[i]);end
  K[row,i]=log(abs((1e4-s[i])/(4PaperISR.VP_MPI^2-s[i])))-sum(K[row,:])
  K[row,:].*=s[i]/(4pi^2*PaperISR.VP_ALPHA*389379000)
 end
 K .*=transpose(k.mu)
 ik=s[lowix]./(4pi*389379000).*k.mu[lowix]
 c=1 .+[sum(PaperISR.vp_lepton(s[i],m^2) for m in (.000510998928,.1056583715,1.77682))+PaperISR.vp_fixed(s[i]) for i in lowix]+K*k.background+im.*ik.*k.background[lowix]
 L=K[:,k.ipdg]*PaperISR.refine_matrix(length(ep));C=K[:,k.ibes]*PaperISR.refine_matrix(length(ec));B=K[:,k.ibb]
 result=zeros(length(targets),size(bb,2));smooth=k.smooth*bb
 for first in 1:128:size(bb,2)
  ix=first:min(first+127,size(bb,2));original=k.expand*smooth[:,ix];born=copy(original)
  base=k.constant.+k.low*pdg[:,ix].+k.charm*bes[:,ix]
  re=similar(original);f=similar(original)
  for iteration in 1:3
   if iteration==3;result[1:length(lowix),ix]=inv.(abs2.(c.+L*pdg[:,ix].+C*bes[:,ix].+B*born));end
   mul!(re,k.bottom,born)
   @. f=inv(abs2(base+re+im*k.local_imag*born))
   @. born=original/f
  end
  result[length(lowix)+1:end,ix]=f
 end
 mu,sd=ensemble_summary(result);old=read_csv(joinpath(out,"retained_input_only_results.csv"))
 sel=length(lowix).+collect(6:5:length(k.ibb))
 dm=maximum(abs.(mu[sel].-old[:,4]));ds=maximum(abs.(sd[sel].-old[:,5]));@assert max(dm,ds)<1e-12
 check=vp_literal(k,pdg[:,2],bes[:,2],bb[:,2]);@assert maximum(abs.(result[length(lowix)+1:5:end,2].-check[:,3]))<1e-12
 write_csv(joinpath(out,"vp_dense.csv"),["energy_GeV","VP_factor","VP_total_error"],hcat(k.grid[targets],mu,sd))
 report=Dict("passed"=>true,"scope"=>"Historical BES input-only conventions retained; five subdivisions per measured interval; original background grid below threshold", "samples"=>size(bb,2),"grid_points"=>length(targets),"measured_mean_max_abs"=>dm,"measured_sd_max_abs"=>ds,"plot_limits_GeV"=>[10.4,11.2],"input_hashes"=>Dict(f=>bytes2hex(sha256(read(joinpath(root,f)))) for f in ["data/vp/Rpdg.f64","data/vp/Rbes.f64","data/vp/Rbb.f64","src/vp.jl"]))
 open(joinpath(out,"vp_dense_validation.toml"),"w") do io;TOML.print(io,report);end
 println(report)
end
main()
