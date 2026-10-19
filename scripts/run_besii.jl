using PaperISR, LinearAlgebra, TOML, SHA
BLAS.set_num_threads(4)
function save_report(path,report)
    open(path,"w") do io; TOML.print(io,report); end
end
function savecov(out,name,energy,a)
    m,sd=ensemble_summary(a); centered=a.-m
    c=centered*transpose(centered)/(size(a,2)-1)
    @assert maximum(abs.(diag(c).-sd.^2))<1e-10
    write_csv(joinpath(out,name*"_covariance.csv"),string.(energy[2:end]),c[2:end,2:end])
end
function main()
    root=dirname(@__DIR__);out=joinpath(root,"results/besii");mkpath(out)
    report=Dict{String,Any}("completed"=>false,"base_commit"=>"5cabb61","julia_version"=>string(VERSION),"scope"=>"BES input replacement only; historical numerical conventions retained")
    save_report(joinpath(out,"validation.toml"),report)
    ep=vec(read_csv(joinpath(root,"data/vp/Epdg.csv")))
    e=vec(read_csv(joinpath(root,"data/vp/Ebb.csv")))
    ec=vec(read_csv(joinpath(root,"data/vp/Ebes.csv")))
    ec0=vec(read_csv(joinpath(root,"reference/legacy_charm/Ebes.csv")))
    pdg=load_draws(joinpath(root,"data/vp/Rpdg.f64"),length(ep))
    bes=load_draws(joinpath(root,"data/vp/Rbes.f64"),length(ec))
    legacy=load_draws(joinpath(root,"reference/legacy_charm/Rbes.f64"),length(ec0))
    raw=load_draws(joinpath(root,"data/total_draws.f64"),length(e))
    stat=load_draws(joinpath(root,"data/stat_draws.f64"),length(e))
    n=size(raw,2);@assert n==10000
    println("Recomputing both ISR ensembles: ",n," samples each");flush(stdout)
    k=Kernel(e);isr=replay(k,raw);isrstat=replay(k,stat)
    expected=load_draws(joinpath(root,"reference/total_forward_factors.f64"),length(e))
    report["isr_legacy_max_abs"]=maximum(abs.(isr.-1 ./expected));@assert report["isr_legacy_max_abs"]<1e-8
    irmean,_=ensemble_summary(isrstat);_,irsd=ensemble_summary(isr)
    appendix=read_csv(joinpath(root,"reference/cpc.csv"))
    report["isr_appendix_mean_max_abs"]=maximum(abs.(irmean[2:end].-appendix[:,8]))
    report["isr_appendix_sd_max_abs"]=maximum(abs.(irsd[2:end].-appendix[:,9]))
    @assert max(report["isr_appendix_mean_max_abs"],report["isr_appendix_sd_max_abs"])<1e-8
    dressed=raw.*isr;bb=dressed.*(e.^2 ./86.85)
    savedbb=load_draws(joinpath(root,"data/vp/Rbb.f64"),length(e))
    report["dressed_input_max_abs"]=maximum(abs.(bb.-savedbb));@assert report["dressed_input_max_abs"]<1e-10
    println("Recomputing legacy and BES VP with paired backgrounds");flush(stdout)
    k0=VPKernel(ep,ec0,e);kn=VPKernel(ep,ec,e)
    old=vp_replay(k0,pdg,legacy,bb);new=vp_replay(kn,pdg,bes,bb)
    expectedvp=load_draws(joinpath(root,"reference/vp/factors.f64"),length(e))
    report["vp_legacy_samples_max_abs"]=maximum(abs.(old.-expectedvp));@assert report["vp_legacy_samples_max_abs"]<1e-8
    ids=[1,2,101,10000];worst=0.
    for i in ids
        literal=vp_literal(kn,pdg[:,i],bes[:,i],bb[:,i];iterations=6)
        worst=max(worst,maximum(abs.(literal[:,3].-new[:,i])))
        if i==1
            report["vp_3_to_6_fixed_input_max_abs"]=maximum(abs.(literal[:,3].-literal[:,6]))
            write_csv(joinpath(out,"fixed_input_iterations.csv"),["energy_GeV";["vp_iteration_$j" for j in 1:6]],hcat(e,literal))
        end
    end
    report["vp_literal_max_abs"]=worst;@assert worst<1e-12
    # Paired toy cross sections: raw experimental points, NOT smoothed solver output.
    born0=dressed./old;born=dressed./new
    vm0,vs0=ensemble_summary(old);vm,vs=ensemble_summary(new)
    dm,ds=ensemble_summary(dressed);bm0,bs0=ensemble_summary(born0);bm,bs=ensemble_summary(born)
    @assert all(isfinite,new) && all(new.>0)
    write_csv(joinpath(out,"comparison.csv"),["energy_GeV","isr_stat_mean","isr_total_sd","vp_legacy_mean","vp_legacy_sd","vp_besii_mean","vp_besii_sd","sigma_dressed_mean_nb","sigma_dressed_sd_nb","sigma_born_legacy_mean_nb","sigma_born_legacy_sd_nb","sigma_born_besii_mean_nb","sigma_born_besii_sd_nb","vp_change_percent","born_change_percent"],hcat(e[2:end],irmean[2:end],irsd[2:end],vm0[2:end],vs0[2:end],vm[2:end],vs[2:end],dm[2:end],ds[2:end],bm0[2:end],bs0[2:end],bm[2:end],bs[2:end],100 .* (vm[2:end]./vm0[2:end].-1),100 .* (bm[2:end]./bm0[2:end].-1)))
    rdm,rds=ensemble_summary(bb);rbm,rbs=ensemble_summary(born.*(e.^2 ./86.85))
    write_csv(joinpath(out,"results.csv"),["energy_GeV","ISR_correction","ISR_total_error","VP_factor","VP_total_error","sigma_dressed_nb","sigma_dressed_error_nb","sigma_Born_nb","sigma_Born_error_nb","R_dressed","R_dressed_error","R_Born","R_Born_error"],hcat(e[2:end],irmean[2:end],irsd[2:end],vm[2:end],vs[2:end],dm[2:end],ds[2:end],bm[2:end],bs[2:end],rdm[2:end],rds[2:end],rbm[2:end],rbs[2:end]))
    for (name,v) in (("isr",isr),("vp",new),("sigma_dressed",dressed),("sigma_born",born));savecov(out,name,e,v);end
    delta=new.-old;delm,dels=ensemble_summary(delta)
    write_csv(joinpath(out,"paired_vp_shift.csv"),["energy_GeV","mean_shift","MC_standard_error"],hcat(e[2:end],delm[2:end],dels[2:end]./sqrt(n)))
    # Validate the historical two-common-nuisance input model.
    central=read_csv(joinpath(root,"data/besii/input.csv"));z=read_csv(joinpath(root,"data/besii/shared_nuisance.csv"))
    reconstructed=central[:,2].+central[:,3]*transpose(z[:,1]).+central[:,4]*transpose(z[:,2])
    report["bes_input_reconstruction_max_abs"]=maximum(abs.(bes.-reconstructed));@assert report["bes_input_reconstruction_max_abs"]<1e-13
    report["toys"]=n;report["charm_points"]=length(ec);report["pQCD_restart_GeV"]=5.05
    report["max_vp_change_percent"]=maximum(abs.(100 .* (vm[2:end]./vm0[2:end].-1)))
    report["max_born_change_percent"]=maximum(abs.(100 .* (bm[2:end]./bm0[2:end].-1)))
    report["vp_sd_change_percent_range"]=[minimum(100 .* (vs[2:end]./vs0[2:end].-1)),maximum(100 .* (vs[2:end]./vs0[2:end].-1))]
    report["max_paired_shift_MC_SE"]=maximum(dels[2:end])/sqrt(n)
    report["historical_ISR_30_to_60_max_abs"]=maximum(abs.(unfold(k,raw[:,1];states=60).-isr[:,1]))
    report["uncertainty_scope"]="Historical 10000-sample convention retained, including fixed first bottom sample; two correlated charm nuisances; sample mean is reported, first background sample is not nominal. No model/systematic uncertainty added for numerical conventions."
    report["hashes"]=Dict(f=>bytes2hex(sha256(read(joinpath(root,f)))) for f in ["src/vp.jl","src/PaperISR.jl","scripts/run_besii.jl","data/besii/input.csv","data/besii/shared_nuisance.csv","data/vp/Rbes.f64","data/vp/Ebes.csv","data/vp/Rpdg.f64","data/total_draws.f64"])
    report["completed"]=true;save_report(joinpath(out,"validation.toml"),report)
    println("COMPLETED ",report);flush(stdout)
end
main()
