"""One-time import; runtime reproduction needs only Julia and bundled data."""
from pathlib import Path
import hashlib,json,shutil
import numpy as np
from scipy.io import loadmat
root=Path(__file__).resolve().parents[1]
source=Path('/Volumes/Seagate_xkDong/Academic/works/VP_B/R_b/matlab/codes1')
manifest={'source_directory':str(source),'arrays':{},'source_files':{}}
for mode,name in [('stat','ISR_factor_sta_err.mat'),('total','ISR_factor_total_err.mat')]:
    p=source/name
    d=loadmat(p,variable_names=['datacombis','sigma1','onedeltasamp'])
    inp=d['datacombis'].T
    np.savetxt(root/'data'/f'{mode}_input.csv',inp,delimiter=',',header='energy_GeV,sigma_nb,stat_nb,sys_uncorrelated_nb,sys_common1_nb,sys_common2_nb',comments='')
    for variable,label,folder in [('sigma1','draws','data'),('onedeltasamp','forward_factors','reference')]:
        arr=np.asarray(d[variable],dtype='<f8')
        dest=root/folder/f'{mode}_{label}.f64'
        arr.tofile(dest)
        manifest['arrays'][str(dest.relative_to(root))]={'shape_in_file':[int(v) for v in arr.shape],'encoding':'little-endian Float64, toy-major: all energies of toy 1, then toy 2', 'sha256':hashlib.sha256(dest.read_bytes()).hexdigest()}
    manifest['source_files'][name]=hashlib.sha256(p.read_bytes()).hexdigest()
for name in ['ISR_staerr.m','ISRcalculation.m','iteration.m','g1.m','Newsigma.m']:
    dest=root/'reference'/name
    shutil.copyfile(source/name,dest)
    manifest['source_files'][name]=hashlib.sha256(dest.read_bytes()).hexdigest()
shutil.copyfile(root.parent/'BottomoniumReproduction.jl/data/reference/cpc.csv',root/'reference/cpc.csv')
(root/'data/provenance.json').write_text(json.dumps(manifest,indent=2)+'\n')
print({k:v['shape_in_file'] for k,v in manifest['arrays'].items()})
