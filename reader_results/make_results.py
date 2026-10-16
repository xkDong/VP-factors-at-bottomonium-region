"""Standalone table/figure builder: python make_results.py (NumPy, Matplotlib)."""
from pathlib import Path
import hashlib, json, re, shutil, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
for sub in ('data','figures','verification'):
    (ROOT/sub).mkdir(exist_ok=True)
source = json.loads((ROOT/'source/provenance.json').read_text())
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
for name, meta in source.items():
    assert sha(ROOT/'source'/name) == meta['sha256'], name
for name, key in [('old_isr_validation.toml','completed'),('old_vp_validation.toml','completed'),('new_acceptance.toml','passed')]:
    assert re.search(r'^'+key+r' = true$',(ROOT/'source'/name).read_text(),re.M)
accept = (ROOT/'source/new_acceptance.toml').read_text()
assert source['new_results.csv']['sha256'] == re.search(r'^"comparison.csv" = "([a-f0-9]{64})"$',accept,re.M).group(1)
def read(name):
    return np.genfromtxt(ROOT/'source'/name,delimiter=',',names=True)
oi,on,ov,vn,new = (read(n) for n in ['old_isr.csv','old_isr_nominal.csv','old_vp.csv','old_vp_nominal.csv','new_results.csv'])
applied=read('new_applied_results.csv')
applied_validation=(ROOT/'source/new_applied_validation.toml').read_text()
assert re.search(r'^completed = true$',applied_validation,re.M)
for name in ('comparison.csv','dressed_covariance.csv','born_covariance.csv'):
    source_name={'comparison.csv':'new_applied_results.csv','dressed_covariance.csv':'new_dressed_covariance.csv','born_covariance.csv':'new_born_covariance.csv'}[name]
    assert source[source_name]['sha256']==re.search(r'^"'+re.escape(name)+r'" = "([a-f0-9]{64})"$',applied_validation,re.M).group(1)
for name in ('energy_GeV','isr_nominal','isr_mean','isr_sd','vp_nominal','vp_mean','vp_sd'):
    assert np.allclose(new[name],applied[name],rtol=1e-11,atol=1e-12),name
new=applied
on,vn=on[1:],vn[1:]
e = new['energy_GeV']
assert len(e)==209 and np.all(np.diff(e)>0)
for tab in (oi,on,ov,vn):
    assert np.allclose(e,tab['energy_GeV'],rtol=0,atol=1e-12)
def write(name, headers, cols):
    a=np.column_stack(cols)
    assert np.isfinite(a).all()
    np.savetxt(ROOT/'data'/name,a,delimiter=',',header=','.join(headers),comments='',fmt='%.17g')

old={'isr':(on['inverse'],oi['inverse_stat_mean'],oi['inverse_total_sd']),
     'vp':(vn['vp3'],ov['vp_mean'],ov['vp_sd'])}
updated={k:tuple(new[f'{k}_{n}'] for n in ('nominal','mean','sd')) for k in ('isr','vp')}
for version, data in [('paper_reproduction',old),('numerical_updates',updated)]:
    for k,(nom,mean,sd) in data.items():
        assert np.all(nom>0) and np.all(mean>0) and np.all(sd>=0)
        label='stat_toy_mean' if version=='paper_reproduction' and k=='isr' else 'total_toy_mean'
        write(f'{version}_{k}.csv',['energy_GeV','factor_nominal',f'factor_{label}','factor_input_sd'],[e,nom,mean,sd])
write('numerical_updates_cross_sections.csv',['energy_GeV','R_dressed_nominal','R_dressed_total_toy_mean','R_dressed_input_sd','R_born_nominal','R_born_total_toy_mean','R_born_input_sd'],[e]+[new[f'{q}_{n}'] for q in ('dressed','born') for n in ('nominal','mean','sd')])
old_r = read('old_dressed_R_summary.csv')
assert np.allclose(old_r['energy_GeV'],e,rtol=0,atol=1e-12)
# Exact normalization of each branch; new implementation returns pb, converted to nb.
mu_old = 86.85/e**2
alpha = 1/137.035999074
mu_new = 4*np.pi*alpha**2/(3*e**2)*np.sqrt(1-4*.105658**2/e**2)*(1+2*.105658**2/e**2)*389379000/1000
old_sigma = [old_r[n]*mu_old for n in ('R_nominal','R_total_toy_mean','R_input_sd')]
new_sigma = [new['dressed_'+n]*mu_new for n in ('nominal','mean','sd')]
sigma_direct=read('new_applied_sigma.csv')
for col,name in zip(new_sigma,('dressed_nominal_nb','dressed_mean_nb','dressed_sd_nb')):
    assert np.allclose(col,sigma_direct[name],rtol=1e-12,atol=1e-13)
for version,cols,mu in [('paper_reproduction',old_sigma,mu_old),('numerical_updates',new_sigma,mu_new)]:
    write(version+'_sigma_dressed_nb.csv',['energy_GeV','sigma_dressed_nominal_nb','sigma_dressed_total_toy_mean_nb','sigma_dressed_input_sd_nb','sigma_mumu_reference_nb'],[e]+cols+[mu])
write('comparison_sigma_dressed_nb.csv',['energy_GeV','old_nominal_nb','new_nominal_nb','old_total_toy_mean_nb','new_total_toy_mean_nb','old_input_sd_nb','new_input_sd_nb'],[e,old_sigma[0],new_sigma[0],old_sigma[1],new_sigma[1],old_sigma[2],new_sigma[2]])
rd_cov=np.loadtxt(ROOT/'source/new_dressed_covariance.csv',delimiter=',',skiprows=1)
sigma_cov=rd_cov*mu_new[:,None]*mu_new[None,:]
assert np.allclose(np.diag(sigma_cov),new_sigma[2]**2,rtol=1e-10,atol=1e-24)
np.savetxt(ROOT/'data/numerical_updates_sigma_dressed_covariance_nb2.csv',np.column_stack([e,sigma_cov]),delimiter=',',header='energy_GeV,'+','.join('%.17g'%v for v in e),comments='',fmt='%.17g')
old_born_r=read('old_born_R_summary.csv')
assert np.allclose(old_born_r['energy_GeV'],e,rtol=0,atol=1e-12)
old_born_sigma=[old_born_r[n]*mu_old for n in ('R_nominal','R_total_toy_mean','R_input_sd')]
new_born_sigma=[new['born_'+n]*mu_new for n in ('nominal','mean','sd')]
assert np.allclose(old_born_sigma[0],old_sigma[0]/old['vp'][0],rtol=1e-12,atol=1e-13)
assert np.allclose(new_born_sigma[0],new_sigma[0]/updated['vp'][0],rtol=1e-12,atol=1e-13)
for version,cols in [('paper_reproduction',old_born_sigma),('numerical_updates',new_born_sigma)]:
    write(version+'_sigma_born_nb.csv',['energy_GeV','sigma_born_nominal_nb','sigma_born_total_toy_mean_nb','sigma_born_input_sd_nb'],[e]+cols)
write('comparison_sigma_born_nb.csv',['energy_GeV','old_nominal_nb','new_nominal_nb','old_total_toy_mean_nb','new_total_toy_mean_nb','old_input_sd_nb','new_input_sd_nb'],[e,old_born_sigma[0],new_born_sigma[0],old_born_sigma[1],new_born_sigma[1],old_born_sigma[2],new_born_sigma[2]])
born_cov=np.loadtxt(ROOT/'source/new_born_covariance.csv',delimiter=',',skiprows=1)*mu_new[:,None]*mu_new[None,:]
assert np.allclose(np.diag(born_cov),new_born_sigma[2]**2,rtol=1e-10,atol=1e-24)
np.savetxt(ROOT/'data/numerical_updates_sigma_born_covariance_nb2.csv',np.column_stack([e,born_cov]),delimiter=',',header='energy_GeV,'+','.join('%.17g'%v for v in e),comments='',fmt='%.17g')
for k in ('isr','vp','dressed','born'):
    matrix=np.loadtxt(ROOT/'source'/f'new_{k}_covariance.csv',delimiter=',',skiprows=1)
    energies=np.loadtxt(ROOT/'source'/f'new_{k}_covariance.csv',delimiter=',',max_rows=1)
    assert matrix.shape==(209,209) and np.allclose(energies,e,rtol=0,atol=1e-12)
    assert np.allclose(matrix,matrix.T,atol=1e-16,rtol=1e-12)
    assert np.allclose(np.diag(matrix),new[f'{k}_sd']**2,rtol=1e-10,atol=1e-24)
    # Add row energies: first cell is zero, remaining header cells are column energies.
    np.savetxt(ROOT/'data'/f'numerical_updates_{k}_covariance.csv',np.column_stack([e,matrix]),delimiter=',',header='energy_GeV,'+','.join('%.17g'%v for v in e),comments='',fmt='%.17g')
components=read('new_error_components.csv')
component_names=[n for n in components.dtype.names if n=='energy_GeV' or n.startswith('vp_')]
write('numerical_updates_error_components.csv',component_names,[components[n] for n in component_names])
plt.rcParams.update({'font.size':11,'axes.titlesize':11,'axes.labelsize':11,'legend.fontsize':10,'xtick.labelsize':9,'ytick.labelsize':9})
blue,orange='#0072B2','#E69F00'
figures=[]
quantities=[]
summary={}
def scatter(a,y,label,color,marker='o'):
    a.scatter(e,y,s=10,label=label,c=color,marker=marker,linewidths=.7)
def finish(fig,axes,name):
    for a in axes.flat:
        a.set_xlabel(r'$\sqrt{s}$ [GeV]');a.grid(alpha=.15)
    fig.savefig(ROOT/'figures'/f'{name}.png',dpi=180)
    fig.savefig(ROOT/'figures'/f'{name}.pdf')
    fig.canvas.draw()
    for i,a in enumerate(axes.flat):
        bbox=a.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(ROOT/'verification'/f'{name}_panel{i+1}.png',dpi=120,bbox_inches=bbox.expanded(1.03,1.04))
    figures.append((name,fig))
for k,label in [('isr',r'$C_{\mathrm{ISR}}=1/(1+\delta)$'),('vp',r'$F_{\mathrm{VP}}=|1-\Pi|^{-2}$')]:
    a,b=old[k],updated[k]
    delta=100*(b[0]/a[0]-1); ds=b[2]-a[2]
    write(f'comparison_{k}.csv',['energy_GeV','old_nominal','new_nominal','nominal_change_percent','old_input_sd','new_input_sd','sd_change_absolute'],[e,a[0],b[0],delta,a[2],b[2],ds])
    idx=int(np.argmax(abs(delta)))
    summary[k]={'max_abs_nominal_change_percent':float(abs(delta[idx])),'energy_GeV':float(e[idx]),'max_abs_sd_change':float(max(abs(ds)))}
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for vals,col,tag,marker in [(a,blue,'Paper reproduction','o'),(b,orange,'Numerical updates','x')]:
        scatter(axes[0,0],vals[0],tag,col,marker)
        scatter(axes[1,0],vals[2],tag,col,marker)
    axes[0,0].set(title='a   Nominal factors',ylabel=label);axes[0,0].legend()
    scatter(axes[0,1],delta,'Nominal change',orange,'x')
    axes[0,1].set(title='b   Same-definition factor comparison',ylabel='Updated / reproduction − 1 (%)')
    axes[1,0].set(title='c   Propagated input uncertainties',ylabel='Factor standard deviation');axes[1,0].legend()
    scatter(axes[1,1],ds,'SD difference',orange,'x')
    axes[1,1].set(title='d   Absolute uncertainty difference',ylabel='Updated SD − reproduction SD')
    for axis in (axes[0,1],axes[1,1]):axis.axhline(0,c='0.6',lw=.7)
    finish(fig,axes,f'comparison_{k}')
    quantities += [f'{k}_old_nominal',f'{k}_new_nominal',f'{k}_old_sd',f'{k}_new_sd',f'{k}_nominal_change',f'{k}_sd_change']
for version,data,color in [('paper_reproduction',old,blue),('numerical_updates',updated,orange)]:
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for j,k in enumerate(('isr','vp')):
        vals=data[k]
        # Historical panels use published center convention; updated panels use nominal center.
        center=vals[1] if version=='paper_reproduction' else vals[0]
        centerlabel=('Stat-toy mean' if k=='isr' else 'Total-toy mean') if version=='paper_reproduction' else 'Nominal'
        axes[0,j].errorbar(e,center,yerr=vals[2],fmt='.',ms=3,elinewidth=.6,c=color,label=centerlabel+' ± input SD')
        axes[0,j].set(title=f'{k.upper()} factor',ylabel='ISR correction factor' if k=='isr' else 'VP factor');axes[0,j].legend()
        scatter(axes[1,j],vals[2],'Input SD',color)
        axes[1,j].set(title=f'{k.upper()} input uncertainty',ylabel='Factor standard deviation')
        quantities.append(f'{version}_{k}_display_center')
    finish(fig,axes,version)
# Main GitHub figures: same nominal-center convention for both branches.
# A complete immutable commit identifies code, defaults and dependencies;
# source hashes and this transform's hash also enter every evaluator fingerprint.
fingerprints={}
series_rows=[]
main_specs=[('isr',old['isr'][0],updated['isr'][0],old['isr'][2],updated['isr'][2],r'$C_{\mathrm{ISR}}=1/(1+\delta)$'),
            ('sigma_dressed',old_sigma[0],new_sigma[0],old_sigma[2],new_sigma[2],r'$\sigma_{\mathrm{dressed}}$ [nb]'),
            ('sigma_born',old_born_sigma[0],new_born_sigma[0],old_born_sigma[2],new_born_sigma[2],r'$\sigma_{\mathrm{Born}}$ [nb]'),
            ('vp',old['vp'][0],updated['vp'][0],old['vp'][2],updated['vp'][2],r'$F_{\mathrm{VP}}=|1-\Pi|^{-2}$')]
for key,old_y,new_y,old_sd,new_sd,label in main_specs:
    for version,y,sd in [('paper_reproduction',old_y,old_sd),('numerical_updates',new_y,new_sd)]:
        sid=key+'_'+version
        metadata={'version':version,'quantity':key,'sources':source,'transform_sha256':sha(Path(__file__)),'center':'nominal','uncertainty':'sample SD of historical input ensemble; old 10000, new 9999','connection':'straight segments at measured energies, no new evaluator or fit'}
        fp=hashlib.sha256(json.dumps(metadata,sort_keys=True).encode()).hexdigest()
        fingerprints[sid]={'fingerprint':fp,'metadata':metadata}
        series_rows += [(sid,fp,float(x),float(v),float(u)) for x,v,u in zip(e,y,sd)]
with (ROOT/'data/github_series.csv').open('w',newline='') as f:
    w=csv.writer(f,lineterminator="\n");w.writerow(['series_id','evaluator_fingerprint','energy_GeV','nominal','input_sd']);w.writerows(series_rows)
for sid,info in fingerprints.items():
    assert len({r[1] for r in series_rows if r[0]==sid})==1
(ROOT/'verification/github_fingerprints.json').write_text(json.dumps(fingerprints,indent=2)+'\n')
for key,old_y,new_y,old_sd,new_sd,label in main_specs:
    fig,axes=plt.subplots(2,1,figsize=(10,6.6),sharex=True,gridspec_kw={'height_ratios':[3,1]},layout='constrained')
    for y,sd,color,tag,style in [(old_y,old_sd,blue,'Paper reproduction','--'),(new_y,new_sd,orange,'Numerical updates','-')]:
        if key=='vp':
            assert np.all(np.isfinite(y-sd)) and np.all(sd>=0)
            axes[0].fill_between(e,y-sd,y+sd,color=color,alpha=.20,linewidth=0)
        axes[0].plot(e,y,linestyle=style,color=color,lw=1.3,marker='.' if style=='-' else None,ms=3,label=tag)
    axes[0].set_ylabel(label)
    axes[0].legend(loc='best',framealpha=.95)
    if key=='vp':axes[0].set_title('Nominal factors with ±1 input-SD bands',loc='left')
    else:axes[0].set_title('Nominal results at 209 measured energies',loc='left')
    delta=100*(new_y/old_y-1)
    axes[1].scatter(e,delta,c=orange,s=12,marker='x',linewidths=.8)
    axes[1].axhline(0,color='0.5',lw=.7)
    axes[1].set_ylabel('Change [%]')
    finish(fig,axes,'github_'+key)
    quantities += [key+'_github_old_nominal',key+'_github_new_nominal',key+'_github_relative_change']
    if key=='vp':quantities += ['vp_github_old_band_limits','vp_github_new_band_limits']
(ROOT/'verification/summary.json').write_text(json.dumps(summary,indent=2)+'\n')
(ROOT/'verification/quantity_verdict.json').write_text(json.dumps({'schema_id':'quantity_verdict_v1','schema_version':1,'quantities':quantities,'verdict':'pass'},indent=2)+'\n')
# Separate overview binding for each figure, with the complete checked observable list.
for name,fig in figures:
    manifest={'figure':f'../figures/{name}.pdf','data':[f'../data/{p.name}' for p in sorted((ROOT/'data').glob('*.csv'))],'script':'../make_results.py','command':'python make_results.py','style':'rcParams in make_results.py; points only; blue/orange','source':'../source/provenance.json','display_acceptance':{'plotted_quantities':quantities,'verdict_bindings':[{'quantity':q,'verdict_path':'quantity_verdict.json','verdict_sha256':sha(ROOT/'verification/quantity_verdict.json')} for q in quantities],'overview_figure':{'path':f'../figures/{name}.png','archived':True}},'sha256':{f'../figures/{name}.png':sha(ROOT/'figures'/f'{name}.png'),'../make_results.py':sha(Path(__file__))}}
    if name.startswith('github_'):
        manifest['style']='rcParams in script; blue dashed old, orange solid updated; VP bands nominal +/- input SD'
        manifest['series_provenance']={'table':'../data/github_series.csv','series_column':'series_id','fingerprint_column':'evaluator_fingerprint','metadata':'github_fingerprints.json'}
    (ROOT/'verification'/f'{name}.provenance.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(summary,indent=2))
