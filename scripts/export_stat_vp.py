"""Audit the stat-labelled VP backgrounds against total-labelled inputs."""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.io import loadmat
root=Path(__file__).resolve().parents[1]
source=Path('/Volumes/Seagate_xkDong/Academic/works/VP_B/matlab/vp_iterate/inputR_for_MC_sta_err_0320.mat')
data=loadmat(source,variable_names=['Rpdg','Rbes'])
meta={'source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'arrays':{}}
for key in ('Rpdg','Rbes'):
    arr=data[key];dest=root/'data/vp'/f'{key}_stat.f64'
    # Audit only: both archived backgrounds are identical, so do not duplicate binary inputs.
    current=np.fromfile(root/'data/vp'/f'{key}.f64',dtype='<f8').reshape(arr.shape)
    assert np.array_equal(arr,current)
    meta['arrays'][key]={'identical_to_total':True,'shape':list(arr.shape),'sha256':hashlib.sha256((root/'data/vp'/f'{key}.f64').read_bytes()).hexdigest()}
(root/'data/vp/stat_provenance.json').write_text(json.dumps(meta,indent=2)+'\n')
