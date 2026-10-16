using PaperISR, LinearAlgebra, TOML, SHA
BLAS.set_num_threads(2)
function main()
    root=dirname(@__DIR__);out=joinpath(root,"results/updates")
    r=TOML.parsefile(joinpath(out,"validation.toml"))
    independent=TOML.parsefile(joinpath(out,"independent_methods.toml"))
    @assert r["completed"] && independent["completed"]
    for (file,expected) in r["hashes"]
        @assert bytes2hex(sha256(read(joinpath(root,file))))==expected
    end
    grid=read_csv(joinpath(out,"grid_convergence.csv"))
    # Combined nominal discretization sensitivity budget, conservatively added.
    bound=grid[end,5]+r["background_refinement_4_max_abs"]+r["tailfine_max_vp_abs"]+r["tail400_max_vp_abs"]
    @assert bound<1e-6
    @assert r["isr_residual_nb"]<1e-11 && r["vp_residual"]<1e-10
    @assert independent["vp_independent_Pi_max_abs"]<1e-10
    @assert independent["isr_independent_forward_abs_nb"]<1e-9
    eigmins=Dict{String,Float64}()
    for name in ("vp","isr","dressed","born","born_stat")
        c=read_csv(joinpath(out,name*"_covariance.csv"))
        eigmins[name]=eigmin(Symmetric(c))
        @assert eigmins[name]>-1e-10
    end
    input=read_csv(joinpath(root,"data/total_input.csv"));e=input[:,1];k=Kernel(e)
    n=length(e);W=zeros(n,n)
    for i in 2:n,j in 1:i-1
        W[i,j]+=k.left[i,j];W[i,j+1]+=k.right[i,j]
    end
    y=k.smoothing*input[:,2];exact=zeros(n);exact[2:end]=LowerTriangular(W[2:end,2:end])\y[2:end]
    target=ones(n);target[2:end]=exact[2:end]./y[2:end]
    rows=Vector{Float64}[]
    for states in (30,60,120,240)
        f=unfold(k,input[:,2];states)
        push!(rows,[states,maximum(abs.(f.-target)),maximum(abs.(W*(y.*f).-y))])
    end
    write_csv(joinpath(out,"legacy_iteration_convergence.csv"),["states","max_inverse_error_vs_exact","max_forward_residual_nb"],permutedims(hcat(rows...)))
    comparison=read_csv(joinpath(out,"comparison.csv"))
    stat=read_csv(joinpath(out,"stat_total_uncertainty.csv"))
    z=(stat[:,2].-stat[:,3])./sqrt.((stat[:,2].^2 .+stat[:,3].^2)./(2*9998))
    propagation=read_csv(joinpath(out,"propagation.csv"))
    verdict=Dict{String,Any}("passed"=>true,
        "scope"=>"Given historical input uncertainty model and declared single-smoothing signed linear estimator",
        "nominal_vp_resolution_budget_abs"=>bound,"resolution_target_abs"=>1e-6,
        "covariance_eigenvalue_minima"=>eigmins,
        "negative_draws_policy"=>"Retain signed estimates; no clipping, dropping, or Gaussian-prior positivity claim",
        "nominal_positive"=>r["nominal_min_dressed_R"]>0,
        "vp_bottom_stat_above_all_input_max_gaussian_mc_z"=>maximum(z),
        "max_independent_error_approximation_difference_percent"=>maximum(abs.(propagation[:,4])),
        "max_nominal_vp_change_vs_appendix_percent"=>maximum(abs.(100 .* (comparison[:,5]./comparison[:,15].-1))),
        "not_included"=>"Smoothing model bias, pQCD/model/parameter uncertainty, newly inferred experimental correlations",
        "settings"=>Dict("isr_smoothing_start_GeV"=>10.6659,"span"=>.045,"second_smoothing"=>false,
            "isr_order"=>48,"vp_refinement"=>20,"vp_background_refinement"=>2,"vp_upper_GeV"=>200.,
            "vp_tail_step_GeV"=>.125,"vp_rtol"=>1e-11,"historical_toys"=>10000,"uncertainty_samples"=>9999))
    verdict["input_hashes"]=Dict(f=>bytes2hex(sha256(read(joinpath(root,f)))) for f in
        ["data/stat_draws.f64","data/total_draws.f64","data/vp/Rpdg.f64","data/vp/Rbes.f64","data/vp/stat_provenance.json"])
    verdict["evidence_hashes"]=Dict(f=>bytes2hex(sha256(read(joinpath(out,f)))) for f in
        ["validation.toml","independent_methods.toml","comparison.csv","grid_convergence.csv"])
    open(joinpath(out,"acceptance.toml"),"w") do io;TOML.print(io,verdict);end
    println(verdict)
end
main()
