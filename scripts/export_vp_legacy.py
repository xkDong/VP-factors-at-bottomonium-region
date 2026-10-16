from pathlib import Path
import hashlib,json,shutil
import numpy as np
from scipy.io import loadmat
root=Path(__file__).resolve().parents[1]
src=Path('/Volumes/Seagate_xkDong/Academic/works/VP_B/matlab/vp_iterate')
a=loadmat(src/'inputR_for_MC_tot_err_0320.mat',variable_names=['Rpdg','Rbes','Rbb','Epdg','Ebes','Ebb','Rbbdreunsmo'])
keys=['Rbb','Ebb','Etot','PILemutau','Rpipi','RpQCD1','RpQCD2','RpQCD3','RpQCD4','RpQCDb','sigmamu','vp1','vpBelleSamp','vpbmean_belle_energy','vpberr_belle_energy','a','b']
b=loadmat(src/'vp_smo_tot_err_0320.mat',variable_names=keys)
manifest={'sources':{},'files':{},'a':int(b['a'].item()),'b':int(b['b'].item())}
def binary(path,v):
 v=np.asarray(v,dtype='<f8');v.tofile(root/path)
 manifest['files'][path]={'shape':list(v.shape),'sha256':hashlib.sha256((root/path).read_bytes()).hexdigest()}
for k in ['Rpdg','Rbes','Rbb']:
 binary('data/vp/'+k+'.f64',a[k])
for k in ['Epdg','Ebes','Ebb']:
 np.savetxt(root/'data/vp'/f'{k}.csv',a[k].ravel(),delimiter=',',header='energy_GeV',comments='')
binary('reference/vp/Rbb_smoothed.f64',b['Rbb'])
binary('reference/vp/factors.f64',b['vpBelleSamp'])
np.savetxt(root/'reference/vp/grid.csv',np.column_stack([b['Etot'].ravel(),b['sigmamu'].ravel(),b['PILemutau'].real.ravel(),b['PILemutau'].imag.ravel(),b['vp1']]),delimiter=',',header='energy_GeV,sigma_mu_pb,lepton_real,lepton_imag,vp1,vp2,vp3,vp4,vp5',comments='')
for name in ['vp_iterate_sta.m','vp_iterate.m','PIH.m','PIL.m','vacuum.m','Newenergy.m','sigmu.m','pipi.m','RpQCD.m','alphas.m','datacollection.m']:
 shutil.copyfile(src/name,root/'reference/vp'/name)
for name in ['inputR_for_MC_tot_err_0320.mat','vp_smo_tot_err_0320.mat']:
 manifest['sources'][name]=hashlib.sha256((src/name).read_bytes()).hexdigest()
ref=np.genfromtxt(root/'reference/cpc.csv',delimiter=',',names=True)
print('appendix mean',max(abs(b['vpbmean_belle_energy'].ravel()[1:]-ref['vp'])))
print('appendix error',max(abs(b['vpberr_belle_energy'].ravel()[1:]-ref['vp_error'])))
print('a,b',manifest['a'],manifest['b'])
print('Ebb first',a['Ebb'].ravel()[:3], 'smooth start',a['Ebb'].ravel()[21])
# Establish relationship with the existing ISR module, independently of VP outputs.
isr=loadmat(src.parent.parent/'R_b/matlab/codes1/ISR_factor_total_err.mat',variable_names=['Rbdresamp'])
print('ISR input link max',np.max(abs(a['Rbb']-isr['Rbdresamp'])))
manifest['isr_dressed_input_max_abs']=float(np.max(abs(a['Rbb']-isr['Rbdresamp'])))
(root/'data/vp/provenance.json').write_text(json.dumps(manifest,indent=2)+'\n')
