using PaperISR, Test, LinearAlgebra
BLAS.set_num_threads(1)
@testset "Numerical core (synthetic inputs)" begin
    e=collect(range(10.55852,11.2;length=210))
    y=@. .1+exp(-((e-10.87)/.08)^2)
    y[1]=0
    k=Kernel(e)
    f=unfold(k,y); ref=unfold_literal(k,y)
    @test all(isfinite,f)
    @test f[1]==1
    @test maximum(abs.(f.-ref))<1e-8
    @test maximum(abs.(unfold(k,2y).-f))<1e-10
    @test maximum(abs.(k.smoothing*ones(210).-ones(210)))<1e-12
    @test maximum(abs.((k.smoothing*e)[23:end].-e[23:end]))<1e-12
    @test_throws ErrorException Kernel(reverse(e))
    @test_throws ErrorException unfold(k,-y)
    input=hcat(e,y,fill(.001,210),zeros(210,3)); input[1,2:end].=0
    d=sample_inputs(input;toys=30,seed=1)
    @test d==sample_inputs(input;toys=30,seed=1)
    @test d[:,1]==y
    @test all(d.>=0)
    @test replay(k,d)[:,1]≈f
    m,s=ensemble_summary(hcat(f,f))
    @test m==f
    @test all(iszero,s)
    mktemp() do p,io
        close(io)
        vals=reshape(collect(1.:12.),3,4)
        write(p,htol.(reinterpret(UInt64,vec(vals))))
        @test load_draws(p,3)==vals
    end
end
root=dirname(@__DIR__)
if isfile(joinpath(root,"data","total_input.csv"))
    @testset "Recovered MATLAB nominal result" begin
        a=read_csv(joinpath(root,"data","total_input.csv"))
        ref=read_csv(joinpath(root,"reference","central_replay.csv"))
        @test a[:,1]≈ref[:,1]
        @test maximum(abs.(unfold(Kernel(a[:,1]),a[:,2]).-ref[:,3]))<1e-8
    end
else
    @warn "Historical-data tests NOT RUN: external-drive input import incomplete."
end
@testset "Historical VP" begin
    ep=vec(read_csv(joinpath(root,"data/vp/Epdg.csv")))
    ec=vec(read_csv(joinpath(root,"data/vp/Ebes.csv")))
    eb=vec(read_csv(joinpath(root,"data/vp/Ebb.csv")))
    k=VPKernel(ep,ec,eb)
    grid=read_csv(joinpath(root,"reference/vp/grid.csv"))
    @test maximum(abs.(k.grid.-grid[:,1]))<1e-12
    @test maximum(abs.(k.mu./grid[:,2].-1))<1e-12
    ids=[1,2,10000]
    pdg=load_draws(joinpath(root,"data/vp/Rpdg.f64"),length(ep))[:,ids]
    bes=load_draws(joinpath(root,"data/vp/Rbes.f64"),length(ec))[:,ids]
    bb=load_draws(joinpath(root,"data/vp/Rbb.f64"),length(eb))[:,ids]
    expected=load_draws(joinpath(root,"reference/vp/factors.f64"),length(eb))[:,ids]
    smooth=load_draws(joinpath(root,"reference/vp/Rbb_smoothed.f64"),length(eb))[:,ids]
    @test maximum(abs.(vp_smooth(eb,bb).-smooth))<1e-12
    got=vp_replay(k,pdg,bes,bb)
    @test maximum(abs.(got.-expected))<1e-12
    for i in 1:3
        @test maximum(abs.(got[:,i].-vp_literal(k,pdg[:,i],bes[:,i],bb[:,i])[:,3]))<1e-12
    end
    @test_throws ArgumentError vp_replay(k,pdg,bes,bb;iterations=0)
end
include("updated.jl")
