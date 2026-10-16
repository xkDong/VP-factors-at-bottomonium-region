"""Compare committed VP results; run from any directory with numpy/matplotlib."""
from pathlib import Path
import hashlib, io, json, subprocess, sys, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/branch_comparison'
OUT.mkdir(exist_ok=True)
OLD, NEW = '5cabb61', '27d9bd1'
inputs = {}
def read(ref, path):
    raw = subprocess.check_output(['git', 'show', f'{ref}:{path}'], cwd=ROOT)
    inputs[f'{ref}:{path}'] = hashlib.sha256(raw).hexdigest()
    return raw
def table(ref, path):
    return np.genfromtxt(io.BytesIO(read(ref, path)), delimiter=',', names=True)
old = table(OLD, 'results/vp/appendix_comparison.csv')
nom = table(OLD, 'results/vp/nominal.csv')[1:]
new = table(NEW, 'results/updates/comparison.csv')
validation = read(OLD, 'results/vp/validation.toml').decode()
accept = read(NEW, 'results/updates/acceptance.toml').decode()
assert re.search(r'^completed = true$', validation, re.M)
assert re.search(r'^passed = true$', accept, re.M)
expected_hash = re.search(r'^"comparison.csv" = "([a-f0-9]{64})"$', accept, re.M).group(1)
scope = re.search(r'^scope = "([^"]+)"$', accept, re.M).group(1)
assert inputs[f'{NEW}:results/updates/comparison.csv'] == expected_hash
e = new['energy_GeV']
assert len(e) == 209 and np.allclose(e, old['energy_GeV'], atol=1e-12, rtol=0)
assert np.allclose(e, nom['energy_GeV'], atol=1e-12, rtol=0)
om, nm = old['vp_mean'], new['vp_mean']
on, nn = nom['vp3'], new['vp_nominal']
os, ns = old['vp_sd'], new['vp_sd']
assert all(np.all(np.isfinite(x)) and np.all(x > 0) for x in (om,nm,on,nn,os,ns))
dm, dn, ds = 100*(nm/om-1), 100*(nn/on-1), 100*(ns/os-1)
columns = ['energy_GeV','old_mean','new_mean','old_nominal','new_nominal','old_sd','new_sd','mean_change_percent','nominal_change_percent','sd_change_percent']
values = np.column_stack([e,om,nm,on,nn,os,ns,dm,dn,ds])
np.savetxt(OUT/'comparison.csv',values,delimiter=',',header=','.join(columns),comments='',fmt='%.16g')
summary = {}
for name, delta in [('mean',dm),('nominal',dn),('sd',ds)]:
    k = np.argmax(abs(delta))
    summary[name] = dict(max_abs_percent=float(abs(delta[k])),signed_percent=float(delta[k]),energy_GeV=float(e[k]),median_abs_percent=float(np.median(abs(delta))),min_percent=float(min(delta)),max_percent=float(max(delta)))
summary['old_sd_range'] = [float(min(os)),float(max(os))]
summary['new_sd_range'] = [float(min(ns)),float(max(ns))]
summary['old_relative_sd_percent_range'] = [float(min(100*os/om)),float(max(100*os/om))]
summary['new_relative_sd_percent_range'] = [float(min(100*ns/nm)),float(max(100*ns/nm))]
summary['largest_sd_change_row'] = dict(zip(columns,values[np.argmax(abs(ds))].tolist()))
summary['largest_mean_change_row'] = dict(zip(columns,values[np.argmax(abs(dm))].tolist()))
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
plt.rcParams.update({'font.size':11,'axes.titlesize':11,'axes.labelsize':11,'legend.fontsize':10,'xtick.labelsize':9,'ytick.labelsize':9})
fig, ax = plt.subplots(2,2,figsize=(12,8),layout='constrained')
blue, orange = '#0072B2','#E69F00'
def dots(a,y,label,color,marker='o'):
    a.scatter(e,y,s=11,label=label,color=color,marker=marker,linewidths=.6)
dots(ax[0,0],om,'Reproduction: MC mean',blue)
dots(ax[0,0],nm,'Updated: MC mean',orange,'x')
ax[0,0].set(ylabel=r'$F_{\mathrm{VP}}$',title='a   VP factors: mean versus mean')
ax[0,0].legend(loc='lower right')
dots(ax[0,1],dm,'Mean versus mean',orange,'x')
dots(ax[0,1],dn,'Nominal versus nominal',blue)
ax[0,1].axhline(0,color='0.6',lw=.7)
ax[0,1].set(ylabel='Updated / reproduction − 1 (%)',title='b   Factor changes at matching center definitions')
ax[0,1].legend(loc='best')
dots(ax[1,0],1e3*os,'Reproduction: 10,000 toys',blue)
dots(ax[1,0],1e3*ns,'Updated: 9,999 random toys',orange,'x')
ax[1,0].set(ylabel=r'$\sigma(F_{\mathrm{VP}})\;[10^{-3}]$',title='c   Propagated input uncertainties (sample SD)')
ax[1,0].legend(loc='best')
dots(ax[1,1],ds,'SD change',orange,'x')
ax[1,1].axhline(0,color='0.6',lw=.7)
ax[1,1].set(ylabel='Updated SD / reproduction SD − 1 (%)',title='d   Relative uncertainty changes')
for a in ax.flat:
    a.set_xlabel(r'$\sqrt{s}$ [GeV]')
    a.grid(alpha=.16)
fig.savefig(OUT/'vp_branches.png',dpi=180)
fig.savefig(OUT/'vp_branches.pdf')
sys.path.insert(0,'/Users/dongxiangkun/.codex/skills/figure-hygiene/scripts/bin')
from figure_qa import check_figure
findings = check_figure(fig)
(OUT/'geometry_qa.json').write_text(json.dumps(findings,indent=2,default=str))
assert not findings, findings
for k,a in enumerate(ax.flat):
    bbox = a.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(OUT/f'panel_{k+1}.png',bbox_inches=bbox.expanded(1.03,1.04),dpi=160)
quantities = columns[1:]
verdict = {'schema_id':'quantity_verdict_v1','schema_version':1,'quantities':quantities,'verdict':'pass'}
(OUT/'quantity_verdict.json').write_text(json.dumps(verdict,indent=2)+'\n')
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
manifest = {'figure':'vp_branches.pdf','data':['comparison.csv'],'script':'../../scripts/compare_vp_branches.py','command':'python scripts/compare_vp_branches.py','style':'matplotlib rcParams in script; blue/orange; scatter only','source_commits':{'paper-reproduction':OLD,'numerical-updates':NEW},'source_hashes':inputs,'scope':scope,'display_acceptance':{'plotted_quantities':quantities,'verdict_bindings':[{'quantity':q,'verdict_path':'quantity_verdict.json','verdict_sha256':sha(OUT/'quantity_verdict.json')} for q in quantities],'overview_figure':{'path':'vp_branches.png','archived':True}},'sha256':{p.name:sha(p) for p in [OUT/'comparison.csv',OUT/'vp_branches.pdf',OUT/'vp_branches.png']}}
manifest['sha256']['../../scripts/compare_vp_branches.py']=sha(Path(__file__))
(OUT/'provenance.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(summary,indent=2))
