@testset "Updated numerical model" begin
    # Exact Cauchy-PV identities on a nonuniform grid.
    x=[.1,.3,1.,2.,4.,10.]
    for s in (.3,1.7,4.)
        w=PaperISR.analytic_pv_weights(x,s)
        L=log((x[end]-s)/(s-x[1]))
        @test dot(w,ones(length(x)))≈L atol=1e-12
        @test dot(w,x)≈(x[end]-x[1])+s*L atol=1e-12
    end
    e=vec(read_csv(joinpath(root,"data/vp/Ebb.csv")))
    input=read_csv(joinpath(root,"data/total_input.csv"))
    a=updated_isr(e,reshape(input[:,2],:,1))
    @test a.residual<1e-12
    @test minimum(a.dressed[2:end,:])>0
    @test maximum(abs.(a.W.-PaperISR.updated_isr_matrix(e;order=96)))<1e-11
    # Inversion must recover a known, unsmoothed spectrum; checks orientation and units.
    trial=@. exp(-((e-10.87)/.1)^2)+.1
    trial[1]=0
    y=a.W*trial
    recovered=LowerTriangular(a.W[2:end,2:end])\y[2:end]
    @test maximum(abs.(recovered.-trial[2:end]))<1e-12
    # One declared smoothing threshold used across implementations.
    @test PaperISR.energy_smoothing(e;start=10.6659)[22,22]==1
    @test PaperISR.energy_smoothing(e;start=10.6659)[23,23]<1
    ep=vec(read_csv(joinpath(root,"data/vp/Epdg.csv")))
    ec=vec(read_csv(joinpath(root,"data/vp/Ebes.csv")))
    k=updated_vp_kernel(ep,ec,e;refinement=3)
    pdg=load_draws(joinpath(root,"data/vp/Rpdg.f64"),length(ep))[:,1:1]
    bes=load_draws(joinpath(root,"data/vp/Rbes.f64"),length(ec))[:,1:1]
    @test k.smooth==Matrix{Float64}(I,length(e),length(e))
    expected=k.grid[k.ibb].^2 ./ (4pi*PaperISR.UPDATE_ALPHA*PaperISR.UPDATE_CONV).*k.mu[k.ibb]
    @test k.local_imag≈expected
    v=updated_vp(k,pdg,bes,a.dressed)
    @test v.residual<1e-10
    @test_throws ErrorException updated_vp(k,pdg,bes,a.dressed;maxiter=1)
    # Joint ratio propagation retains perfect correlation and cancels it.
    z=[1. 2 3 4;2 4 6 8]
    c=covariance_summary(z;exclude_nominal=false)
    @test c.covariance[1,2]≈2c.covariance[1,1]
    @test ensemble_summary(reshape(z[2,:]./z[1,:],1,:))[2][1]==0
    nominal=copy(input);nominal[:,2].=10.;nominal[:,3:4].=0.;nominal[:,5].=.01;nominal[:,6].=0.
    t=updated_samples(nominal;toys=100,seed=9)
    @test maximum(abs.(t[1,:].-t[2,:]))==0
    @test t==updated_samples(nominal;toys=100,seed=9)
end
