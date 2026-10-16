using PaperISR, LinearAlgebra, TOML, SHA
BLAS.set_num_threads(4)
function main()
    root=dirname(@__DIR__);out=joinpath(root,"results/updates");mkpath(out)
    e=vec(read_csv(joinpath(root,"data/vp/Ebb.csv")))
    open(joinpath(out,"validation.toml"),"w") do io;TOML.print(io,Dict("completed"=>false));end
    ep=vec(read_csv(joinpath(root,"data/vp/Epdg.csv")));ec=vec(read_csv(joinpath(root,"data/vp/Ebes.csv")))
    raw=load_draws(joinpath(root,"data/total_draws.f64"),length(e))
    pdg=load_draws(joinpath(root,"data/vp/Rpdg.f64"),length(ep))
    bes=load_draws(joinpath(root,"data/vp/Rbes.f64"),length(ec))
    a=updated_isr(e,raw)
    report=Dict{String,Any}("isr_residual_nb"=>a.residual,"negative_dressed_values"=>a.negative_count,
        "minimum_dressed_R"=>minimum(a.dressed[2:end,:]),"n_toys"=>size(raw,2))
    println(report);flush(stdout)
    println("Building production VP kernel");flush(stdout)
    k=updated_vp_kernel(ep,ec,e)
    v=updated_vp(k,pdg,bes,a.dressed)
    report["vp_residual"]=v.residual;report["vp_iterations"]=v.iterations
    println("Total uncertainty propagated");flush(stdout)
    inv=a.inverse;born=a.dressed./v.factor
    stats=covariance_summary(v.factor);rb=covariance_summary(born);ii=covariance_summary(inv)
    dr=covariance_summary(a.dressed)
    ref=read_csv(joinpath(root,"reference/cpc.csv"))
    # Nominal and mean are separate estimators; paired toy propagation throughout.
    write_csv(joinpath(out,"comparison.csv"),["energy_GeV","isr_nominal","isr_mean","isr_sd","vp_nominal","vp_mean","vp_sd","vp_mean_se","dressed_nominal","dressed_mean","dressed_sd","born_nominal","born_mean","born_sd","appendix_vp"],
        hcat(e[2:end],inv[2:end,1],ii.mean[2:end],ii.sd[2:end],v.factor[2:end,1],stats.mean[2:end],stats.sd[2:end],
        stats.mean_se[2:end],a.dressed[2:end,1],dr.mean[2:end],dr.sd[2:end],born[2:end,1],rb.mean[2:end],rb.sd[2:end],ref[:,10]))
    for (name,s) in (("vp",stats),("isr",ii),("dressed",dr),("born",rb))
        write_csv(joinpath(out,name*"_covariance.csv"),string.(e[2:end]),s.covariance[2:end,2:end])
    end
    # Quantify a common incorrect independent-error approximation.
    indep=sqrt.((dr.sd./stats.mean).^2 .+(dr.mean.*stats.sd./stats.mean.^2).^2)
    write_csv(joinpath(out,"propagation.csv"),["energy_GeV","paired_born_sd","independent_delta_method_sd","difference_percent"],
        hcat(e[2:end],rb.sd[2:end],indep[2:end],100 .* (indep[2:end]./rb.sd[2:end].-1)))
    rows=Vector{Float64}[]
    nominal=a.dressed[:,1:1]
    for refn in (5,10,20,30)
        println("Grid check ",refn);flush(stdout)
        kk=updated_vp_kernel(ep,ec,e;refinement=refn)
        vv=updated_vp(kk,pdg[:,1:1],bes[:,1:1],nominal)
        push!(rows,[refn,length(kk.grid),vv.iterations,vv.residual,maximum(abs.(vv.factor.-v.factor[:,1:1]))])
        write_csv(joinpath(out,"vp_grid_$refn.csv"),["energy_GeV","vp"],hcat(e,vv.factor[:,1]))
    end
    write_csv(joinpath(out,"grid_convergence.csv"),["refinement","grid_points","iterations","residual","max_diff_from_ref20"],permutedims(hcat(rows...)))
    report["isr_quadrature48_96_max_abs"]=maximum(abs.(a.W.-PaperISR.updated_isr_matrix(e;order=96)))
    # Settings changes are modelling choices; do not fold into statistical errors.
    rows=Vector{Float64}[]
    for (start,span,second) in ((10.6613,.045,false),(10.6659,.035,false),(10.6659,.055,false),(10.6659,.045,true))
        aa=updated_isr(e,raw[:,1:1];start,span)
        kk=updated_vp_kernel(ep,ec,e;second_smoothing=second,start,span)
        vv=updated_vp(kk,pdg[:,1:1],bes[:,1:1],aa.dressed)
        push!(rows,[start,span,Float64(second),maximum(abs.(vv.factor./v.factor[:,1:1].-1)),maximum(abs.(aa.dressed.-nominal))])
    end
    write_csv(joinpath(out,"smoothing_sensitivity.csv"),["start_GeV","span","second_smoothing","max_vp_relative_change","max_dressed_R_change"],permutedims(hcat(rows...)))
    for (tag,upper,step) in (("tail100",100.,.125),("tail400",400.,.125),("tailcoarse",200.,.25),("tailfine",200.,.0625))
        kk=updated_vp_kernel(ep,ec,e;upper,tail_step=step)
        vv=updated_vp(kk,pdg[:,1:1],bes[:,1:1],nominal)
        report[tag*"_max_vp_abs"]=maximum(abs.(vv.factor.-v.factor[:,1:1]))
    end
    report["vp_tolerance_check"]=maximum(abs.(updated_vp(k,pdg[:,1:1],bes[:,1:1],nominal;rtol=1e-13).factor.-v.factor[:,1:1]))
    for bg in (1,4)
        kk=updated_vp_kernel(ep,ec,e;background_refinement=bg)
        vv=updated_vp(kk,pdg[:,1:1],bes[:,1:1],nominal)
        report["background_refinement_$(bg)_max_abs"]=maximum(abs.(vv.factor.-v.factor[:,1:1]))
    end
    report["isr_condition_number"]=cond(a.W[2:end,2:end])
    report["nominal_min_dressed_R"]=minimum(a.dressed[2:end,1])
    neg=vec(sum(a.dressed[:,2:end].<0;dims=2))
    write_csv(joinpath(out,"negative_estimates.csv"),["energy_GeV","negative_toy_count","fraction"],
        hcat(e[2:end],neg[2:end],neg[2:end]./(size(raw,2)-1)))
    # Statistical ensemble from the matching historical stat-only inputs.
    println("Propagating statistical ensemble");flush(stdout)
    statraw=load_draws(joinpath(root,"data/stat_draws.f64"),length(e))
    # Historical stat-labelled background arrays are identical to total-labelled arrays.
    # For a bottom-stat-only component, freeze BOTH backgrounds at nominal values.
    statpdg=repeat(pdg[:,1:1],1,size(raw,2))
    statbes=repeat(bes[:,1:1],1,size(raw,2))
    sa=updated_isr(e,statraw);sv=updated_vp(k,statpdg,statbes,sa.dressed)
    sb=covariance_summary(sa.dressed./sv.factor);ss=covariance_summary(sv.factor)
    write_csv(joinpath(out,"stat_total_uncertainty.csv"),["energy_GeV","vp_bottom_stat_only_sd","vp_all_input_sd","born_bottom_stat_only_sd","born_all_input_sd"],
        hcat(e[2:end],ss.sd[2:end],stats.sd[2:end],sb.sd[2:end],rb.sd[2:end]))
    report["vp_bottom_stat_variance_above_total_count"]=count(ss.sd[2:end].>stats.sd[2:end])
    # Preserve stat covariance rather than forcing a PSD systematic covariance by subtraction.
    write_csv(joinpath(out,"born_stat_covariance.csv"),string.(e[2:end]),sb.covariance[2:end,2:end])
    println("Propagating background-only component");flush(stdout)
    fixedborn=repeat(a.dressed[:,1:1],1,size(raw,2))
    bgv=updated_vp(k,pdg,bes,fixedborn)
    bgs=covariance_summary(bgv.factor);bgb=covariance_summary(fixedborn./bgv.factor)
    write_csv(joinpath(out,"error_components.csv"),["energy_GeV","vp_bottom_stat_only_sd","vp_background_only_sd","vp_all_input_sd","born_bottom_stat_only_sd","born_background_only_sd","born_all_input_sd"],
        hcat(e[2:end],ss.sd[2:end],bgs.sd[2:end],stats.sd[2:end],sb.sd[2:end],bgb.sd[2:end],rb.sd[2:end]))
    report["background_error_model"]="Stored low/charm combined-error ensembles; no justified statistical/systematic split"
    # MC stability from independent disjoint batches of saved toys, not fixed-output equality.
    batchsd=zeros(length(e),20)
    for b in 1:20
        ix=(2+(b-1)*499):(1+b*499)
        batchsd[:,b]=ensemble_summary(v.factor[:,ix])[2]
    end
    _,sdspread=ensemble_summary(batchsd)
    write_csv(joinpath(out,"mc_precision.csv"),["energy_GeV","vp_sd","vp_mean_se","sd_batch_se"],
        hcat(e[2:end],stats.sd[2:end],stats.mean_se[2:end],sdspread[2:end]./sqrt(20)))
    # Validate corrected negative-redraw sampler with independent seeds.
    toyrows=Vector{Float64}[]
    inp=read_csv(joinpath(root,"data/total_input.csv"))
    for seed in (20261016,20261017)
        newsamples=updated_samples(inp;toys=2000,seed)
        na=updated_isr(e,newsamples)
        nv=updated_vp(k,pdg[:,1:2000],bes[:,1:2000],na.dressed)
        ns=covariance_summary(nv.factor)
        push!(toyrows,[seed,na.residual,nv.residual,na.negative_count,maximum(abs.(ns.mean.-stats.mean)),
            maximum(abs.(ns.mean.-stats.mean)./sqrt.(ns.mean_se.^2 .+ stats.mean_se.^2))])
    end
    write_csv(joinpath(out,"new_seed_check.csv"),["seed","isr_residual","vp_residual","negative_values","max_mean_difference","max_mean_z"],
        permutedims(hcat(toyrows...)))
    report["completed"]=true
    report["scientific_status"]="provisional until resolution and negative-draw audits are assessed"
    report["hashes"]=Dict(f=>bytes2hex(sha256(read(joinpath(root,f)))) for f in ["src/updated.jl","scripts/audit_updates.jl"])
    open(joinpath(out,"validation.toml"),"w") do io;TOML.print(io,report);end
    println(report)
end
main()
