using PaperISR,LinearAlgebra,Test
root=dirname(@__DIR__)
@testset "Requested paper methods" begin
 e=collect(range(10.56,11.2;length=210));y=reshape(@.(1+.2*(e-10.56)),:,1)
 f=robust_lowess(e,y);@test maximum(abs.(f.-y))<1e-12
 trial=copy(y);trial[:,1].+=.001 .*sin.(collect(1:210));trial[110,1]+=1
 ordinary=PaperISR.energy_smoothing(e)*trial;robust=robust_lowess(e,trial)
 @test abs(robust[110]-y[110])<abs(ordinary[110]-y[110])
 @test robust[1:22]==trial[1:22]
 a=[1. 2. 3. 4.;2. 4. 6. 8.];g=gaussian_summary(a)
 @test g.mean≈[2.5,5.]
 @test g.sd≈sqrt.([1.25,5.])
 @test diag(g.covariance)≈g.sd.^2
 inp=read_csv(joinpath(root,"data/total_input.csv"));en=inp[:,1];p=paper_isr(en,reshape(inp[:,2],:,1));tight=paper_isr(en,reshape(inp[:,2],:,1);rtol=1e-8)
 @test p.relative_change<1e-4
 @test maximum(abs.(tight.W*(tight.dressed.*tight.mu).-tight.smoothed_sigma))<1e-9
 @test maximum(abs.(p.inverse.-tight.inverse))<.001
 println("iterations: ",maximum(p.iterations)," -> ",maximum(tight.iterations))
end
