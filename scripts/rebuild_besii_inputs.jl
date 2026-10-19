using PaperISR, SHA, TOML
root=dirname(@__DIR__)
a=read_csv(joinpath(root,"data/besii/input.csv"))
z=read_csv(joinpath(root,"data/besii/shared_nuisance.csv"))
r=a[:,2].+a[:,3]*transpose(z[:,1]).+a[:,4]*transpose(z[:,2])
@assert all(r.>0)
# Retain bundled bytes if already equivalent (cross-language last-bit rounding).
p=joinpath(root,"data/vp/Rbes.f64")
if isfile(p)
    @assert maximum(abs.(load_draws(p,size(r,1)).-r))<1e-13
else
    write(p,htol.(reinterpret(UInt64,vec(r))))
end
ep=joinpath(root,"data/vp/Ebes.csv")
if isfile(ep)
    @assert vec(read_csv(ep))==a[:,1]
else
    write_csv(ep,["energy_GeV"],a[:,1:1])
end
println("BES input arrays regenerated/verified from bundled central values and nuisances")
