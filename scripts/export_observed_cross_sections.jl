# Apply paired ISR and VP factors to the unsmoothed experimental cross sections.
using PaperISR, LinearAlgebra, TOML, SHA
BLAS.set_num_threads(4)
function main()
    root=dirname(@__DIR__);out=joinpath(root,"results/observed_cross_sections");mkpath(out)
    e=vec(read_csv(joinpath(root,"data/vp/Ebb.csv")))
    ep=vec(read_csv(joinpath(root,"data/vp/Epdg.csv")));ec=vec(read_csv(joinpath(root,"data/vp/Ebes.csv")))
    raw=load_draws(joinpath(root,"data/total_draws.f64"),length(e))
    pdg=load_draws(joinpath(root,"data/vp/Rpdg.f64"),length(ep))
    bes=load_draws(joinpath(root,"data/vp/Rbes.f64"),length(ec))
    println("Computing ISR and VP factors from the declared smoothed solver spectrum");flush(stdout)
    a=updated_isr(e,raw)
    k=updated_vp_kernel(ep,ec,e)
    v=updated_vp(k,pdg,bes,a.dressed)
    # raw is the original experimental cross section in nb; do not replace it
    # with the smoothed input or the consistently-renormalized solver input.
    dressed_sigma=raw.*a.inverse
    born_sigma=dressed_sigma./v.factor
    mu=PaperISR.updated_mu.(e.^2)./1000
    dressed=dressed_sigma./mu;born=born_sigma./mu
    ss=covariance_summary(dressed_sigma);bs=covariance_summary(born_sigma)
    dr=covariance_summary(dressed);br=covariance_summary(born)
    ii=covariance_summary(a.inverse);vv=covariance_summary(v.factor)
    old=read_csv(joinpath(root,"results/updates/comparison.csv"))
    @assert maximum(abs.(a.inverse[2:end,1].-old[:,2]))<1e-12
    @assert maximum(abs.(v.factor[2:end,1].-old[:,5]))<1e-12
    @assert maximum(abs.(ss.sd.^2 .- diag(ss.covariance)))<1e-12
    @assert maximum(abs.(bs.sd.^2 .- diag(bs.covariance)))<1e-12
    write_csv(joinpath(out,"comparison.csv"),["energy_GeV","isr_nominal","isr_mean","isr_sd","vp_nominal","vp_mean","vp_sd","vp_mean_se","dressed_nominal","dressed_mean","dressed_sd","born_nominal","born_mean","born_sd","appendix_vp"],
        hcat(e[2:end],a.inverse[2:end,1],ii.mean[2:end],ii.sd[2:end],v.factor[2:end,1],vv.mean[2:end],vv.sd[2:end],vv.mean_se[2:end],dressed[2:end,1],dr.mean[2:end],dr.sd[2:end],born[2:end,1],br.mean[2:end],br.sd[2:end],old[:,15]))
    for (name,s) in (("dressed",dr),("born",br),("sigma_dressed_nb",ss),("sigma_born_nb",bs))
        write_csv(joinpath(out,name*"_covariance.csv"),string.(e[2:end]),s.covariance[2:end,2:end])
    end
    write_csv(joinpath(out,"sigma.csv"),["energy_GeV","observed_nominal_nb","dressed_nominal_nb","dressed_mean_nb","dressed_sd_nb","born_nominal_nb","born_mean_nb","born_sd_nb"],hcat(e[2:end],raw[2:end,1],dressed_sigma[2:end,1],ss.mean[2:end],ss.sd[2:end],born_sigma[2:end,1],bs.mean[2:end],bs.sd[2:end]))
    report=Dict{String,Any}("completed"=>true,"definition"=>"Original unsmoothed experimental sigma_nb times paired ISR inverse, divided by paired VP for Born","toys"=>9999,"vp_residual"=>v.residual,"isr_residual"=>a.residual,"vp_iterations"=>v.iterations,
        "input_hashes"=>Dict(f=>bytes2hex(sha256(read(joinpath(root,f)))) for f in ["data/total_draws.f64","data/vp/Rpdg.f64","data/vp/Rbes.f64","src/updated.jl","scripts/export_observed_cross_sections.jl"]),
        "output_hashes"=>Dict(f=>bytes2hex(sha256(read(joinpath(out,f)))) for f in ["comparison.csv","sigma.csv","dressed_covariance.csv","born_covariance.csv"]))
    open(joinpath(out,"validation.toml"),"w") do io;TOML.print(io,report);end
    println("Completed original-data cross sections and paired covariance");flush(stdout)
end
main()
