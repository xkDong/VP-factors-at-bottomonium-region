# Diagnostic only: isolate the missing 1/alpha in the historical hadronic imaginary term.
using PaperISR, LinearAlgebra, TOML
BLAS.set_num_threads(1)
function main()
    root=dirname(@__DIR__)
    ep=vec(read_csv(joinpath(root,"data/vp/Epdg.csv")))
    ec=vec(read_csv(joinpath(root,"data/vp/Ebes.csv")))
    eb=vec(read_csv(joinpath(root,"data/vp/Ebb.csv")))
    k=VPKernel(ep,ec,eb)
    pdg=load_draws(joinpath(root,"data/vp/Rpdg.f64"),length(ep))[:,1:1]
    bes=load_draws(joinpath(root,"data/vp/Rbes.f64"),length(ec))[:,1:1]
    bb=load_draws(joinpath(root,"data/vp/Rbb.f64"),length(eb))[:,1:1]
    newim=k.local_imag./PaperISR.VP_ALPHA
    constant=k.constant.+im.*(newim.-k.local_imag).*k.background[k.ibb]
    corrected=VPKernel(k.grid,k.ibb,k.mu,constant,k.low,k.charm,k.bottom,newim,
        k.smooth,k.expand,k.background,k.ipdg,k.ibes)
    old=vec(vp_replay(k,pdg,bes,bb))
    new=vec(vp_replay(corrected,pdg,bes,bb))
    d=100 .* (new./old.-1)
    write_csv(joinpath(root,"results/vp/imaginary_convention.csv"),
        ["energy_GeV","historical_vp","paper_imaginary_term_vp","difference_percent"],hcat(eb,old,new,d))
    report=Dict("scope"=>"Nominal input; only imaginary hadronic 1/alpha changed",
        "max_abs_percent"=>maximum(abs.(d)),"energy_GeV"=>eb[argmax(abs.(d))])
    open(joinpath(root,"results/vp/imaginary_convention.toml"),"w") do io;TOML.print(io,report);end
    println(report)
end
main()
