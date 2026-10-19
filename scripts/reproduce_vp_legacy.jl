using PaperISR, LinearAlgebra, TOML, SHA
BLAS.set_num_threads(4)
function main()
    root=dirname(@__DIR__)
    mode=isempty(ARGS) ? "isr" : ARGS[1]
    mode in ("isr","saved") || error("Usage: reproduce_vp.jl [isr|saved]")
    out=joinpath(root,"results",mode=="isr" ? "vp" : "vp_saved");mkpath(out)
    report=Dict{String,Any}("completed"=>false,"input_mode"=>mode,"julia_version"=>string(VERSION))
    open(joinpath(out,"validation.toml"),"w") do io;TOML.print(io,report);end
    ep=vec(read_csv(joinpath(root,"data/vp/Epdg.csv")))
    ec=vec(read_csv(joinpath(root,"reference/legacy_charm/Ebes.csv")))
    eb=vec(read_csv(joinpath(root,"data/vp/Ebb.csv")))
    pdg=load_draws(joinpath(root,"data/vp/Rpdg.f64"),length(ep))
    bes=load_draws(joinpath(root,"reference/legacy_charm/Rbes.f64"),length(ec))
    bb=load_draws(joinpath(root,"data/vp/Rbb.f64"),length(eb))
    if mode=="isr"
        println("Rebuilding dressed R from historical ISR inputs");flush(stdout)
        isrin=read_csv(joinpath(root,"data/total_input.csv"))
        @assert maximum(abs.(isrin[:,1].-eb))<1e-12
        draws=load_draws(joinpath(root,"data/total_draws.f64"),length(eb))
        factors=replay(Kernel(eb),draws)
        rebuilt=draws.*factors.*(eb.^2 ./86.85)
        report["isr_to_vp_input_max_abs"]=maximum(abs.(rebuilt.-bb))
        @assert report["isr_to_vp_input_max_abs"]<1e-10
        bb=rebuilt
    end
    k=VPKernel(ep,ec,eb)
    grid=read_csv(joinpath(root,"reference/vp/grid.csv"))
    report["toys"]=size(bb,2)
    report["grid_max_abs"]=maximum(abs.(k.grid.-grid[:,1]))
    report["mu_max_relative"]=maximum(abs.(k.mu./grid[:,2].-1))
    smooth=vp_smooth(eb,bb)
    report["smooth_max_abs"]=maximum(abs.(smooth.-load_draws(joinpath(root,"reference/vp/Rbb_smoothed.f64"),length(eb))))
    @assert report["grid_max_abs"]<1e-12
    @assert report["mu_max_relative"]<1e-12
    @assert report["smooth_max_abs"]<1e-10
    println(report);flush(stdout)
    literal=vp_literal(k,pdg[:,1],bes[:,1],bb[:,1];iterations=5)
    nominal=vp_replay(k,pdg[:,1:1],bes[:,1:1],bb[:,1:1];iterations=5)
    report["literal_max_abs"]=maximum(abs.(nominal[:,1].-literal[:,5]))
    report["nominal_saved_max_abs"]=maximum(abs.(literal.-grid[k.ibb[1:5:end],5:9]))
    @assert report["nominal_saved_max_abs"]<1e-8
    @assert report["literal_max_abs"]<1e-12
    write_csv(joinpath(out,"nominal.csv"),["energy_GeV","vp1","vp2","vp3","vp4","vp5"],hcat(eb,literal))
    println("Replaying all VP samples");flush(stdout)
    result=vp_replay(k,pdg,bes,bb)
    saved=load_draws(joinpath(root,"reference/vp/factors.f64"),length(eb))
    report["samples_max_abs"]=maximum(abs.(result.-saved))
    report["samples_max_relative"]=maximum(abs.(result./saved.-1))
    @assert report["samples_max_abs"]<1e-8
    mu,sd=ensemble_summary(result)
    ref=read_csv(joinpath(root,"reference/cpc.csv"))
    @assert maximum(abs.(eb[2:end].-ref[:,1]))<1e-12
    dm=mu[2:end].-ref[:,10];ds=sd[2:end].-ref[:,11]
    report["appendix_mean_max_abs"]=maximum(abs.(dm))
    report["appendix_error_max_abs"]=maximum(abs.(ds))
    @assert maximum(abs.(dm))<1e-8
    @assert maximum(abs.(ds))<1e-8
    report["nominal_3_to_5_max_abs"]=maximum(abs.(literal[:,3].-literal[:,5]))
    write_csv(joinpath(out,"appendix_comparison.csv"),
        ["energy_GeV","vp_mean","vp_sd","appendix_vp","appendix_sd","mean_difference","sd_difference"],
        hcat(eb[2:end],mu[2:end],sd[2:end],ref[:,10:11],dm,ds))
    report["completed"]=true
    report["hashes"]=Dict(f=>bytes2hex(sha256(read(joinpath(root,f)))) for f in
        ["src/PaperISR.jl","src/vp.jl","src/vp_pqcd.jl","scripts/reproduce_vp.jl","data/vp/Rpdg.f64","reference/legacy_charm/Rbes.f64","data/vp/Rbb.f64","data/total_draws.f64"])
    open(joinpath(out,"validation.toml"),"w") do io;TOML.print(io,report);end
    println(report)
end
main()
