"""Freeze reader inputs from the two audited commits (requires NumPy)."""
from pathlib import Path
import hashlib, json, subprocess
import io
import numpy as np

root = Path(__file__).resolve().parents[1]
out = root / 'reader_results'
sources = {
    'old_isr.csv': ('5cabb61', 'results/historical/appendix_comparison.csv'),
    'old_isr_nominal.csv': ('5cabb61', 'results/historical/total_nominal.csv'),
    'old_vp.csv': ('5cabb61', 'results/vp/appendix_comparison.csv'),
    'old_vp_nominal.csv': ('5cabb61', 'results/vp/nominal.csv'),
    'new_results.csv': ('27d9bd1', 'results/updates/comparison.csv'),
    'old_isr_validation.toml': ('5cabb61', 'results/historical/validation.toml'),
    'old_vp_validation.toml': ('5cabb61', 'results/vp/validation.toml'),
    'new_acceptance.toml': ('27d9bd1', 'results/updates/acceptance.toml'),
}
for name in ('isr', 'vp', 'dressed', 'born'):
    sources[f'new_{name}_covariance.csv'] = ('27d9bd1', f'results/updates/{name}_covariance.csv')
sources['new_error_components.csv'] = ('27d9bd1', 'results/updates/error_components.csv')
for name in ('test_summary.json','independent_methods.toml','smoothing_sensitivity.csv','legacy_iteration_convergence.csv'):
    sources['new_'+name] = ('27d9bd1', 'results/updates/'+name)
manifest = {}
for name, (ref, path) in sources.items():
    raw = subprocess.check_output(['git', 'show', f'{ref}:{path}'], cwd=root)
    (out/'source'/name).write_bytes(raw)
    commit = subprocess.check_output(['git', 'rev-parse', ref], cwd=root, text=True).strip()
    manifest[name] = {'commit':commit,'path':path,'sha256':hashlib.sha256(raw).hexdigest()}
raw = subprocess.check_output(['git','show','5cabb61:data/vp/Rbb.f64'],cwd=root)
energy_raw = subprocess.check_output(['git','show','5cabb61:data/vp/Ebb.csv'],cwd=root)
energy = np.loadtxt(io.BytesIO(energy_raw),delimiter=',',skiprows=1)
samples = np.frombuffer(raw,dtype='<f8').reshape(10000,210)
assert hashlib.sha256(raw).hexdigest() == '1bd10251507ce2e19f26b5cfdef4805e188cfa17a0e1b657b21ece572e8621db'
name='old_dressed_R_summary.csv'
np.savetxt(out/'source'/name,np.column_stack([energy[1:],samples[0,1:],samples[:,1:].mean(axis=0),samples[:,1:].std(axis=0,ddof=1)]),delimiter=',',header='energy_GeV,R_nominal,R_total_toy_mean,R_input_sd',comments='',fmt='%.17g')
manifest[name]={'commit':subprocess.check_output(['git','rev-parse','5cabb61'],cwd=root,text=True).strip(),'path':'data/vp/Rbb.f64','input_sha256':hashlib.sha256(raw).hexdigest(),'energy_sha256':hashlib.sha256(energy_raw).hexdigest(),'transformation':'Shape 10000 x 210; remove threshold; first nominal, mean and sample SD of all 10000 rows; before second VP smoothing','sha256':hashlib.sha256((out/'source'/name).read_bytes()).hexdigest()}
vp_raw = subprocess.check_output(['git','show','5cabb61:reference/vp/factors.f64'],cwd=root)
vp = np.frombuffer(vp_raw,dtype='<f8').reshape(10000,210)
vp_summary=np.genfromtxt(out/'source/old_vp.csv',delimiter=',',names=True)
assert np.all(vp>0) and np.allclose(vp[:,1:].mean(axis=0),vp_summary['vp_mean'],rtol=0,atol=1e-13)
assert np.allclose(vp[:,1:].std(axis=0,ddof=1),vp_summary['vp_sd'],rtol=0,atol=1e-14)
born = samples/vp
name='old_born_R_summary.csv'
np.savetxt(out/'source'/name,np.column_stack([energy[1:],born[0,1:],born[:,1:].mean(axis=0),born[:,1:].std(axis=0,ddof=1)]),delimiter=',',header='energy_GeV,R_nominal,R_total_toy_mean,R_input_sd',comments='',fmt='%.17g')
manifest[name]={'commit':manifest['old_dressed_R_summary.csv']['commit'],'path':['data/vp/Rbb.f64','reference/vp/factors.f64'],'input_sha256':[hashlib.sha256(raw).hexdigest(),hashlib.sha256(vp_raw).hexdigest()],'transformation':'Derived Born R = paired unsmoothed dressed R / saved VP factor; nominal first row, mean and sample SD over 10000 rows; not independently matched to a published Born table','sha256':hashlib.sha256((out/'source'/name).read_bytes()).hexdigest()}
for name,path in {'new_applied_results.csv':'comparison.csv','new_applied_sigma.csv':'sigma.csv','new_applied_validation.toml':'validation.toml','new_dressed_covariance.csv':'dressed_covariance.csv','new_born_covariance.csv':'born_covariance.csv'}.items():
    location=root/'results/observed_cross_sections'/path
    payload=location.read_bytes()
    (out/'source'/name).write_bytes(payload)
    manifest[name]={'path':str(location.relative_to(root)),'definition':'Paired corrections applied to original unsmoothed experimental sigma_nb','script_sha256':hashlib.sha256((root/'scripts/export_observed_cross_sections.jl').read_bytes()).hexdigest(),'sha256':hashlib.sha256(payload).hexdigest()}
(out/'source/provenance.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Exported',len(sources),'committed inputs')
