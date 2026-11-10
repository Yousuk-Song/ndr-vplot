"""Stage each result in a sibling temporary directory, validate, then publish."""
import argparse, json, os, subprocess, sys, tempfile, uuid
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser()
for key in ['engine','sample','cram','reference','bed','output']:
    p.add_argument('--'+key,required=True)
for key,default in [('threads',4),('half-window',1000),('x-bin',4),('min-length',100),('max-length',200),('mapq',30)]:
    p.add_argument('--'+key,type=int,default=default)
a=p.parse_args(); out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True)
stage=Path(tempfile.mkdtemp(prefix='.'+a.sample+'.running.',dir=out.parent))
cmd=[sys.executable,a.engine,'--sample',a.sample,'--cram',a.cram,'--reference',a.reference,'--ndr','STAD='+a.bed,'--output',str(stage)]
for key in ['threads','half_window','x_bin','min_length','max_length','mapq']:
    cmd.extend(['--'+key.replace('_','-'),str(getattr(a,key))])
try:
    subprocess.run(cmd,check=True)
    with np.load(stage/'vplot.npz',allow_pickle=False) as z:
        raw=z['raw_counts']; norm=z['minmax']
        assert raw.shape==(a.max_length-a.min_length,2*a.half_window//a.x_bin,1)
        assert raw.sum()>0 and np.isfinite(norm).all()
        expected=(raw-raw.min())/(raw.max()-raw.min()) if raw.max()>raw.min() else np.zeros_like(norm)
        assert np.allclose(norm,expected)
        assert json.loads((stage/'COMPLETE.json').read_text())['assignments']==int(raw.sum())
    meta=json.loads((stage/'run.json').read_text())
    meta['parameters']['output']=str(out)
    meta['workflow']='Snakemake; staged output validation'
    (stage/'run.json').write_text(json.dumps(meta,indent=2))
    if out.exists():
        backup=out.parent/'.previous'/f'{a.sample}.{uuid.uuid4().hex}'
        backup.parent.mkdir(exist_ok=True)
        out.rename(backup)
        print('Previous remaining files:',backup,flush=True)
    stage.rename(out)
    print('COMPLETE:',out,flush=True)
except BaseException:
    print('Failed staging directory retained:',stage,file=sys.stderr)
    raise
