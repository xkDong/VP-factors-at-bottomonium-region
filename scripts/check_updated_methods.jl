using PaperISR, LinearAlgebra, TOML
BLAS.set_num_threads(2)
function main()
    root=dirname(@__DIR__);out=joinpath(root,"results/updates");mkpath(out)
    e=vec(read_csv(joinpath(root,"data/vp/Ebb.csv")))
    ep=vec(read_csv(joinpath(root,"data/vp/Epdg.csv")));ec=vec(read_csv(joinpath(root,"data/vp/Ebes.csv")))
    input=read_csv(joinpath(root,"data/total_input.csv"));raw=reshape(input[:,2],:,1)
    a=updated_isr(e,raw)
    # Independent convolution: soft contribution is integrated after t=x^beta,
    # instead of using analytic moments of linear sigma(s).
    nodes,weights=PaperISR.gauss_rule(64);s=e.^2;d=a.dressed[:,1].*a.mu
    forward=zeros(length(e))
    for i in 2:length(e)
        beta,delta=PaperISR.updated_radiator(s[i])
        for j in 1:i-1
            xa=1-s[j+1]/s[i];xb=1-s[j]/s[i];ta=xa^beta;tb=xb^beta
            for q in eachindex(nodes)
                t=(ta+tb)/2+(tb-ta)/2*nodes[q];x=t^(1/beta)
                f=(s[i]*(1-x)-s[j])/(s[j+1]-s[j])
                forward[i]+=weights[q]*(tb-ta)/2*(1+delta)*(d[j]*(1-f)+d[j+1]*f)
                u=(nodes[q]+1)/2;x=xa==0 ? xb*u^4 : xa+(xb-xa)*u
                jac=xa==0 ? 2xb*u^3 : (xb-xa)/2
                f=(s[i]*(1-x)-s[j])/(s[j+1]-s[j])
                forward[i]+=weights[q]*jac*PaperISR.regular_radiator(x,beta)*(d[j]*(1-f)+d[j+1]*f)
            end
        end
    end
    report=Dict{String,Any}("isr_independent_forward_abs_nb"=>maximum(abs.(forward.-a.smoothed_sigma[:,1])))
    @assert report["isr_independent_forward_abs_nb"]<1e-9
    k=updated_vp_kernel(ep,ec,e;refinement=20,upper=200.,tail_step=.125)
    pdg=load_draws(joinpath(root,"data/vp/Rpdg.f64"),length(ep))[:,1:1]
    bes=load_draws(joinpath(root,"data/vp/Rbes.f64"),length(ec))[:,1:1]
    v=updated_vp(k,pdg,bes,a.dressed)
    R=copy(k.background)
    R[k.ipdg].+=PaperISR.vp_refine(pdg[:,1],20)
    R[k.ibes].+=PaperISR.vp_refine(bes[:,1],20)
    R[k.ibb].+=k.expand*a.dressed[:,1] # initial spectrum is enough to test the PV functional
    x=k.grid.^2;y=R.*k.mu
    # Direct quadrature of subtracted rational integrand versus analytic PV weights.
    rows=Vector{Float64}[]
    for en in (10.6063,10.6128,10.6458,10.9077,11.2062)
        i=argmin(abs.(k.grid.-en));sn=x[i];ys=y[i]
        ref=ys*log((x[end]-sn)/(sn-x[1]))
        for j in 1:length(x)-1
            h=x[j+1]-x[j]
            if x[j]<=sn<=x[j+1]
                ref+=y[j+1]-y[j]
            else
                for q in eachindex(nodes)
                    t=(nodes[q]+1)/2
                    ref+=weights[q]*h/2*((y[j]-ys)+(y[j+1]-y[j])*t)/((x[j]-sn)+h*t)
                end
            end
        end
        analytical=dot(PaperISR.analytic_pv_weights(x,sn),y)
        pierr=abs(ref-analytical)*sn/(4pi^2*PaperISR.UPDATE_ALPHA*PaperISR.UPDATE_CONV)
        push!(rows,[en,analytical,ref,pierr])
    end
    report["vp_independent_Pi_max_abs"]=maximum(r[4] for r in rows)
    println(rows);println(report);flush(stdout)
    @assert report["vp_independent_Pi_max_abs"]<1e-10
    write_csv(joinpath(out,"independent_pv.csv"),["energy_GeV","analytic_integral","quadrature_integral","Pi_abs_difference"],permutedims(hcat(rows...)))
    # Compare general dense solve and triangular solve independently at the algebraic level.
    dense=a.W[2:end,2:end]\a.smoothed_sigma[2:end,1]
    report["isr_dense_triangular_abs_nb"]=maximum(abs.(dense.-d[2:end]))
    @assert report["isr_dense_triangular_abs_nb"]<1e-10
    report["completed"]=true
    open(joinpath(out,"independent_methods.toml"),"w") do io;TOML.print(io,report);end
    println(report)
end
main()
