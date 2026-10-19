from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/besii';F=OUT/'figures';F.mkdir(exist_ok=True)
a=np.genfromtxt(OUT/'comparison.csv',delimiter=',',names=True)
b=np.genfromtxt(ROOT/'data/besii/input.csv',delimiter=',',names=True)
u=np.genfromtxt(ROOT/'reference/legacy_charm/central.csv',delimiter=',',names=True)
assert len(a)==209 and len(b)==77 and np.all(np.diff(a['energy_GeV'])>0)
assert np.allclose(a['vp_change_percent'],100*(a['vp_besii_mean']/a['vp_legacy_mean']-1),atol=1e-11)
assert np.allclose(a['born_change_percent'],100*(a['sigma_born_besii_mean_nb']/a['sigma_born_legacy_mean_nb']-1),atol=1e-11)
plt.rcParams.update({'font.size':10,'axes.labelsize':10,'legend.fontsize':9,'xtick.labelsize':8,'ytick.labelsize':8})
blue='#0072B2';orange='#D55E00';e=a['energy_GeV']
fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained',sharex=True)
def band(ax,m,sd,color,label,style='-'):
 ax.plot(e,m,color=color,lw=1,ls=style,label=label)
 ax.fill_between(e,m-sd,m+sd,color=color,alpha=.17,lw=0)
for ax in axes.flat:ax.grid(alpha=.15);ax.set_xlabel('Energy (GeV)')
band(axes[0,0],a['isr_stat_mean'],a['isr_total_sd'],blue,'Both inputs (identical)')
axes[0,0].set_ylabel('ISR correction factor')
band(axes[0,1],a['vp_legacy_mean'],a['vp_legacy_sd'],blue,'Historical inputs')
band(axes[0,1],a['vp_besii_mean'],a['vp_besii_sd'],orange,'Published BES inputs','--')
axes[0,1].set_ylabel('VP factor')
band(axes[1,0],a['sigma_dressed_mean_nb'],a['sigma_dressed_sd_nb'],blue,'Both inputs (identical)')
axes[1,0].set_ylabel('Dressed cross section (nb)')
band(axes[1,1],a['sigma_born_legacy_mean_nb'],a['sigma_born_legacy_sd_nb'],blue,'Historical inputs')
band(axes[1,1],a['sigma_born_besii_mean_nb'],a['sigma_born_besii_sd_nb'],orange,'Published BES inputs','--')
axes[1,1].set_ylabel('Born cross section (nb)')
for ax in axes.flat:ax.legend(loc='best')
fig.savefig(F/'overview.png',dpi=180);fig.savefig(F/'overview.pdf')
for i,ax in enumerate(axes.flat):fig.savefig(F/f'overview_panel_{i+1}.png',dpi=130,bbox_inches=ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted()).expanded(1.02,1.05))
fig2,ax=plt.subplots(figsize=(10,4),layout='constrained')
ax.errorbar(u['energy_GeV'],u['R'],yerr=np.hypot(u['error1'],u['error2']),fmt='.',ms=4,color=blue,elinewidth=.6,alpha=.65,label='Historical inputs')
ax.errorbar(b['energy_GeV'],b['R'],yerr=np.hypot(b['stat'],b['sys']),fmt='o',ms=3.5,mfc='white',color=orange,elinewidth=.8,capsize=2,label='BES 2002 + BES 2000 at 5.0 GeV')
ax.set(xlabel='Energy (GeV)',ylabel='Inclusive R');ax.legend();ax.grid(alpha=.15)
fig2.savefig(F/'input_comparison.png',dpi=180);fig2.savefig(F/'input_comparison.pdf')
fig3,axs=plt.subplots(2,1,figsize=(9,6),layout='constrained',sharex=True)
axs[0].plot(e,a['vp_change_percent'],color=orange,label='VP factor')
axs[0].plot(e,a['born_change_percent'],color=blue,ls='--',label='Born cross section')
axs[0].set_ylabel('Change in mean (%)');axs[0].legend()
axs[1].plot(e,100*(a['vp_besii_sd']/a['vp_legacy_sd']-1),color=orange,label='VP uncertainty')
axs[1].set_ylabel('Change in standard deviation (%)');axs[1].set_xlabel('Energy (GeV)');axs[1].legend()
for ax in axs:ax.grid(alpha=.15)
fig3.savefig(F/'relative_changes.png',dpi=180);fig3.savefig(F/'relative_changes.pdf')
from PIL import Image
panels=[Image.open(F/n).convert('RGB') for n in ['input_comparison.png','overview.png','relative_changes.png']]
w=max(im.width for im in panels)
scaled=[im.resize((w,round(im.height*w/im.width))) for im in panels]
canvas=Image.new('RGB',(w,sum(im.height for im in scaled)), 'white')
y=0
for im in scaled:canvas.paste(im,(0,y));y+=im.height
canvas.save(F/'review_overview.png')
qs=['input_R','input_errors','ISR','ISR_error','VP','VP_error','dressed','dressed_error','Born','Born_error','relative_changes']
assert all(np.isfinite(a[n]).all() for n in a.dtype.names)
assert all(np.all(a[n]>=0) for n in a.dtype.names if n.endswith('_sd') or n.endswith('_sd_nb'))
v=F/'quantity_verdict.json';v.write_text(json.dumps({'schema_id':'quantity_verdict_v1','schema_version':1,'quantities':qs,'verdict':'pass'}))
h=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
manifest={'sources':{str(p.relative_to(ROOT)):h(p) for p in [OUT/'comparison.csv',OUT/'validation.toml',ROOT/'data/besii/input.csv',ROOT/'reference/legacy_charm/central.csv',Path(__file__)]},'caption':'209 energies; bands are propagated historical-model total SD, not MC error or numerical accuracy. ISR mean uses statistical toys. Input errors combine stat/sys only for display. All comparison entries are paired.','display_acceptance':{'plotted_quantities':qs,'verdict_bindings':[{'quantity':q,'verdict_path':v.name,'verdict_sha256':h(v)} for q in qs],'overview_figure':{'path':'review_overview.png','archived':True}},'series_fingerprints':{'historical':'5cabb61 kernels + old inputs + same 10000 samples','bes':'same kernels + data/besii/provenance.json + same 10000 nuisances'},'outputs':{p.name:h(p) for p in F.glob('*.png')}}
(F/'provenance.json').write_text(json.dumps(manifest,indent=2))
