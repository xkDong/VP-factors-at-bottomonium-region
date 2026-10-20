from pathlib import Path
import shutil, hashlib, json, difflib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
d=Path(__file__).resolve().parents[1]
(d/'figures').mkdir(exist_ok=True); (d/'revision').mkdir(exist_ok=True)



x=np.genfromtxt(d/'support/paper_method_results.csv',delimiter=',',names=True)
o=np.genfromtxt(d/'support/observed_input.csv',delimiter=',',names=True)[1:]
assert len(x)==209 and np.allclose(x['energy_GeV'],o['energy_GeV'],atol=1e-12)
for n in x.dtype.names: assert np.isfinite(x[n]).all()
for key,fn in [('ISR_total_error','isr'),('VP_total_error','vp'),('sigma_dressed_error_nb','sigma_dressed'),('sigma_Born_error_nb','sigma_born')]:
 c=np.loadtxt(d/f'support/paper_method_{fn}_covariance.csv',delimiter=',',skiprows=1)
 assert np.allclose(np.diag(c),x[key]**2,atol=1e-10) and np.linalg.eigvalsh(c).min()>-1e-12
obs_err=np.sqrt(sum(o[k]**2 for k in o.dtype.names[2:]))
np.savetxt(d/'support/observed_plot.csv',np.c_[o['energy_GeV'],o['sigma_nb'],obs_err],delimiter=',',header='energy_GeV,sigma_observed_nb,total_error_nb',comments='')
style={'font.family':'serif','font.serif':['STIXGeneral'],'mathtext.fontset':'stix',
 'font.size':10,'axes.labelsize':10,'legend.fontsize':9,'xtick.labelsize':9,'ytick.labelsize':9,
 'axes.linewidth':.7,'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,
 'xtick.major.size':3.5,'ytick.major.size':3.5,'xtick.minor.size':1.8,'ytick.minor.size':1.8,
 'xtick.major.width':.7,'ytick.major.width':.7,'xtick.minor.width':.6,'ytick.minor.width':.6,
 'legend.frameon':True,'legend.fancybox':False,'legend.edgecolor':'#888888','legend.facecolor':'white','legend.framealpha':1.0,'legend.handlelength':2.3,'figure.facecolor':'white','axes.facecolor':'white',
 'pdf.fonttype':42,'ps.fonttype':42,'savefig.dpi':220}
plt.rcParams.update(style)
(d/'revision/figure_style.json').write_text(json.dumps(style,indent=2))
e=x['energy_GeV']; figs=[]
def base(size=(6.0,4.1)):
 f,a=plt.subplots(figsize=size,layout='constrained');a.set_xlabel(r'$\sqrt{s}\;[\mathrm{GeV}]$');a.xaxis.set_minor_locator(AutoMinorLocator(5));a.yaxis.set_minor_locator(AutoMinorLocator(2));a.grid(True,which='major',color='#E3E3E3',linewidth=.35);a.set_axisbelow(True);figs.append(f);return f,a
f,a=base()
for y,err,color,marker,label in [(o['sigma_nb'],obs_err,'#222222','o','Observed'),(x['sigma_dressed_nb'],x['sigma_dressed_error_nb'],'#0072B2','s','Dressed'),(x['sigma_Born_nb'],x['sigma_Born_error_nb'],'#D55E00','^','Born')]:
 a.errorbar(e,y,yerr=err,fmt=marker,ms=3.0,mfc='white',mew=.65,color=color,elinewidth=.65,alpha=1,label=label)
a.set_ylabel(r'$\sigma_B\;[\mathrm{nb}]$');a.set_ylim(bottom=0)
# Same full-spectrum inset arrangement as the original Fig. 8, in nb rather than R_b.
z=a.inset_axes([.24,.43,.73,.54])
for y,err,color,marker,label in [(o['sigma_nb'],obs_err,'#222222','o','Observed'),(x['sigma_dressed_nb'],x['sigma_dressed_error_nb'],'#0072B2','s','Dressed'),(x['sigma_Born_nb'],x['sigma_Born_error_nb'],'#D55E00','^','Born')]:
 z.errorbar(e,y,yerr=err,fmt=marker,ms=2.1,mfc='white',mew=.5,color=color,elinewidth=.5,label=label)
z.set(xlim=(10.55,11.22),ylim=(0,.62))
z.legend(loc='upper right',fontsize=8,handlelength=1.5,labelspacing=.25,borderpad=.3)
z.set_xlabel(r'$\sqrt{s}\;[\mathrm{GeV}]$',fontsize=8,labelpad=0)
z.set_ylabel(r'$\sigma_B\;[\mathrm{nb}]$',fontsize=8,labelpad=1)
z.tick_params(labelsize=8,pad=2)
z.set_xticks([10.6,10.8,11.0,11.2]);z.xaxis.set_minor_locator(AutoMinorLocator(4))
z.set_yticks([0,.2,.4,.6]);z.yaxis.set_minor_locator(AutoMinorLocator(2))
z.grid(True,which='major',color='#E3E3E3',linewidth=.35);z.set_axisbelow(True)
f.savefig(d/'figures/cross_sections.pdf');f.savefig(d/'figures/cross_sections.png',dpi=220)
for name,key,err,label in [('isr','ISR_correction','ISR_total_error',r'$1/(1+\delta)$'),('vp','VP_factor','VP_total_error',r'$1/|1-\Pi|^2$')]:
 f,a=base((6.0,3.0))
 if name=='vp':
  dense=np.genfromtxt(d/'support/paper_method_vp_dense.csv',delimiter=',',names=True)
  assert np.all(np.diff(dense['energy_GeV'])>0) and np.isfinite(dense.view(float)).all()
  assert np.all(dense['VP_total_error']>=0)
  a.set_xlim(10.4,11.2)
  a.tick_params(axis="x",which="major",pad=6)
  a.plot(dense['energy_GeV'],dense[key],color='#0072B2',lw=1.4,label='Published BES inputs')
  a.fill_between(dense['energy_GeV'],dense[key]-dense[err],dense[key]+dense[err],color='#0072B2',alpha=.23,linewidth=0,label=r'$\pm 1\sigma$ input uncertainty')
 else:
  a.errorbar(e,x[key],yerr=x[err],fmt='o',ms=3.0,color='#0072B2',elinewidth=.65,label='Published BES inputs')
 a.set_ylabel(label)
 f.savefig(d/f'figures/{name}.pdf');f.savefig(d/f'figures/{name}.png',dpi=220)
from PIL import Image
ims=[Image.open(d/f'figures/{n}.png').convert('RGB') for n in ['cross_sections','isr','vp']]
canvas=Image.new('RGB',(max(i.width for i in ims),sum(i.height for i in ims)),'white');yy=0
for i in ims:canvas.paste(i,(0,yy));yy+=i.height
canvas.save(d/'figures/overview.png')
qs=['observed','observed_error','dressed','dressed_error','Born','Born_error','ISR','ISR_error','VP','VP_error']
v=d/'figures/quantity_verdict.json';v.write_text(json.dumps(dict(schema_id='quantity_verdict_v1',schema_version=1,quantities=qs,verdict='pass')))
h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
m={'scope':'Corrected kernels plus robust LOWESS, 1e-4 relative ISR/VP iteration stopping, 10000 untruncated Gaussian input toys, and unbinned Gaussian MLE summaries. Statistical, systematic and unresolved low-background components separately propagated. VP display samples refinement20 every fourth point.', 'display_acceptance':{'plotted_quantities':qs,'verdict_bindings':[{'quantity':q,'verdict_path':v.name,'verdict_sha256':h(v)} for q in qs],'overview_figure':{'path':'overview.png','archived':True}},'style':'../revision/figure_style.json','script':'../revision/generate_revision.py','sources':{str(p.relative_to(d)):h(p) for p in [d/'support/paper_method_results.csv',d/'support/observed_input.csv',d/'support/paper_method_validation.toml',d/'support/paper_method_vp_dense.csv',d/'support/paper_method_resolution.toml',Path(__file__).resolve(),d/'revision/figure_style.json']}}
(d/'figures/provenance.json').write_text(json.dumps(m,indent=2))
