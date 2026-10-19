from pathlib import Path
import hashlib,json,tomllib
import numpy as np
r=Path(__file__).resolve().parents[1];o=r/'results/besii'
v=tomllib.loads((o/'validation.toml').read_text());assert v['completed']
for f,h in v['hashes'].items():assert hashlib.sha256((r/f).read_bytes()).hexdigest()==h,f
x=np.genfromtxt(o/'results.csv',delimiter=',',names=True);c=np.genfromtxt(o/'comparison.csv',delimiter=',',names=True)
assert len(x)==209 and np.array_equal(x['energy_GeV'],c['energy_GeV'])
checks={}
for name,key in [('isr','ISR_total_error'),('vp','VP_total_error'),('sigma_dressed','sigma_dressed_error_nb'),('sigma_born','sigma_Born_error_nb')]:
 a=np.loadtxt(o/(name+'_covariance.csv'),delimiter=',',skiprows=1)
 assert a.shape==(209,209) and np.all(np.isfinite(a))
 assert np.max(abs(a-a.T))<1e-12
 assert np.max(abs(np.diag(a)-x[key]**2))<1e-10
 eig=np.linalg.eigvalsh(a);assert eig[0]>-1e-12
 checks[name]={'min_eigenvalue':float(eig[0]),'diagonal_sd_max_abs':float(np.max(abs(np.sqrt(np.diag(a))-x[key])))}
mu=86.85/x['energy_GeV']**2
assert np.allclose(x['R_Born']*mu,x['sigma_Born_nb'],rtol=1e-12)
assert np.allclose(x['R_dressed']*mu,x['sigma_dressed_nb'],rtol=1e-12)
result={'passed':True,'runtime_hashes_verified':len(v['hashes']),'points':len(x),'covariance_checks':checks,'limitations':'Historical algorithm retained; no general numerical convergence claim'}
(o/'output_checks.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
