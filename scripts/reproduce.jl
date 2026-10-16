using PaperISR, SHA, TOML, LinearAlgebra
BLAS.set_num_threads(1)
function main()
    root=dirname(@__DIR__)
    mode=isempty(ARGS) ? "historical" : ARGS[1]
    mode in ("historical","fresh") || error("Usage: reproduce.jl [historical|fresh] [toys=10000] [seed=2020]")
    toys=length(ARGS)>1 ? parse(Int,ARGS[2]) : 10000
    seed=length(ARGS)>2 ? parse(Int,ARGS[3]) : 2020
    out=joinpath(root,"results",mode); mkpath(out)
    report=Dict{String,Any}("mode"=>mode,"states"=>30,"julia_version"=>string(VERSION),
        "threads"=>Threads.nthreads(),"seed"=>seed,"hashes"=>Dict{String,String}())
    report["completed"]=false
    open(joinpath(out,"validation.toml"),"w") do io; TOML.print(io,report); end
    for f in ("src/PaperISR.jl", "scripts/reproduce.jl", "reference/cpc.csv")
        report["hashes"][f]=bytes2hex(sha256(read(joinpath(root,f))))
    end
    factors=Dict{String,Matrix{Float64}}()
    energy=Float64[]
    for kind in ("stat","total")
        inp=joinpath(root,"data",kind*"_input.csv")
        isfile(inp) || error("Missing imported input $inp. Restore bundled data or rerun the import described in README.")
        input=read_csv(inp); energy=input[:,1]; k=Kernel(energy)
        report["hashes"][relpath(inp,root)]=bytes2hex(sha256(read(inp)))
        if mode=="historical"
            file=joinpath(root,"data",kind*"_draws.f64")
            draws=load_draws(file,length(energy))
            report["hashes"][relpath(file,root)]=bytes2hex(sha256(read(file)))
        else
            draws=sample_inputs(input;toys,seed=seed+(kind=="total"))
        end
        println("Computing ",kind,": ",size(draws,2)," toys"); flush(stdout)
        result=replay(k,draws)
        if mode=="historical"
            saved=load_draws(joinpath(root,"reference",kind*"_forward_factors.f64"),length(energy))
            err=maximum(abs.(result.-1 ./saved))
            report[kind*"_full_replay_max_abs"]=err
            @assert err<1e-8 "Historical sample replay failed"
        end
        # Verbatim arithmetic reference checks representative samples of this input domain.
        ids=unique([1,2,min(101,size(draws,2)),size(draws,2)])
        err=maximum(maximum(abs.(unfold_literal(k,draws[:,t]).-result[:,t])) for t in ids)
        report[kind*"_literal_max_abs"]=err
        @assert err<1e-8
        factors[kind]=result
        report[kind*"_toys"]=size(draws,2)
        nominal=unfold(k,input[:,2])
        write_csv(joinpath(out,kind*"_nominal.csv"),["energy_GeV","inverse"],hcat(energy,nominal))
        # Historical iteration count is retained, convergence changes are reported separately.
        f60=unfold(k,input[:,2];states=60)
        report[kind*"_30_to_60_max_abs"]=maximum(abs.(f60.-nominal))
    end
    mu,_=ensemble_summary(factors["stat"]); _,sd=ensemble_summary(factors["total"])
    ref=read_csv(joinpath(root,"reference","cpc.csv"))
    @assert maximum(abs.(ref[:,1].-energy[2:end]))<1e-12
    dm=mu[2:end].-ref[:,8]; ds=sd[2:end].-ref[:,9]
    report["appendix_mean_max_abs"]=maximum(abs.(dm))
    report["appendix_error_max_abs"]=maximum(abs.(ds))
    if mode=="historical"
        @assert maximum(abs.(dm))<1e-8
        @assert maximum(abs.(ds))<1e-8
    end
    write_csv(joinpath(out,"appendix_comparison.csv"),
        ["energy_GeV","inverse_stat_mean","inverse_total_sd","appendix_inverse","appendix_sd","mean_difference","sd_difference"],
        hcat(energy[2:end],mu[2:end],sd[2:end],ref[:,8:9],dm,ds))
    report["completed"]=true
    open(joinpath(out,"validation.toml"),"w") do io; TOML.print(io,report); end
    println(report)
end
main()
