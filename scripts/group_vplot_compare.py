#!/usr/bin/env python3
"""Descriptive group-average and difference V-plots from existing sample NPZ files."""
import argparse,csv,itertools,json,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--manifest',required=True,help='TSV columns sample, group, npz; optional run_json')
 p.add_argument('--output',required=True,help='New or empty output directory')
 p.add_argument('--channel',default='STAD',help='NPZ channel name (default: STAD)')
 a=p.parse_args();mf=Path(a.manifest).resolve();out=Path(a.output).resolve()
 if out.exists() and any(out.iterdir()): raise SystemExit('Output directory must be new or empty')
 with mf.open(newline='') as f: rows=list(csv.DictReader(f,delimiter='\t'))
 if not rows or not {'sample','group','npz'}.issubset(rows[0]):raise SystemExit('Manifest requires sample, group, npz columns')
 if any(not r['sample'].strip() or not r['group'].strip() or not r['npz'].strip() for r in rows):raise SystemExit('Missing sample, group or npz')
 if len({r['sample'] for r in rows})!=len(rows):raise SystemExit('Sample IDs must be unique')
 groups=list(dict.fromkeys(r['group'] for r in rows))
 if not 2<=len(groups)<=6:raise SystemExit('Provide 2–6 groups')
 matrices=[];edges=None;hashes=set();metadata_missing=[]; audit=[]
 for r in rows:
  path=Path(r['npz']).expanduser();path=path if path.is_absolute() else mf.parent/path
  with np.load(path,allow_pickle=False) as z:
   raw=z['raw_counts'];channels=z['channels'].astype(str).tolist();x=z['x_edges'];y=z['length_edges']
   if a.channel not in channels:raise SystemExit(f'{r["sample"]}: channel {a.channel} missing')
   if raw.ndim!=3 or raw.shape!=(len(y)-1,len(x)-1,len(channels)) or not np.isfinite(raw).all() or (raw<0).any():raise SystemExit(f'{r["sample"]}: invalid matrix')
   if edges is None:edges=(x.copy(),y.copy())
   elif not np.array_equal(edges[0],x) or not np.array_equal(edges[1],y):raise SystemExit('Inconsistent bin edges')
   v=raw.astype(float);lo=v.min();hi=v.max()
   if hi<=lo:raise SystemExit(f'{r["sample"]}: min–max normalization undefined')
   mat=((v-lo)/(hi-lo))[:,:,channels.index(a.channel)];matrices.append(mat)
  meta=Path(r.get('run_json','') or path.parent/'run.json');meta=meta if meta.is_absolute() else mf.parent/meta
  bedhash=''
  if meta.is_file():
   info=json.loads(meta.read_text()); bedhash=next((v for k,v in info.get('bed_sha256',{}).items() if k.split('=')[0]==a.channel),'')
  if bedhash:hashes.add(bedhash)
  else:metadata_missing.append(r['sample'])
  audit.append({'sample':r['sample'],'group':r['group'],'npz':str(path.resolve()),'npz_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'bed_sha256':bedhash})
 if len(hashes)>1:raise SystemExit('Different BED hashes: do not compare incompatible regions')
 out.mkdir(parents=True,exist_ok=True);means=[np.mean([v for r,v in zip(rows,matrices) if r['group']==g],axis=0) for g in groups];pairs=list(itertools.combinations(range(len(groups)),2));diffs=[means[j]-means[i] for i,j in pairs]; maxdiff=max(float(np.abs(d).max()) for d in diffs) or 1.;maxmean=max(float(v.max()) for v in means) or 1.
 cols=max(len(groups),len(pairs));fig,axs=plt.subplots(2,cols,figsize=(4.3*cols,8),squeeze=False,layout='constrained');extent=[edges[0][0],edges[0][-1],edges[1][0],edges[1][-1]]
 for ax in axs.flat:ax.set_visible(False)
 for i,(g,v) in enumerate(zip(groups,means)):
  ax=axs[0,i];ax.set_visible(True);im=ax.imshow(v,origin='lower',extent=extent,aspect='auto',vmin=0,vmax=maxmean,cmap='viridis',interpolation='nearest');ax.axvline(0,color='white',ls=':',lw=.7);ax.set_title(f'{g} (n={sum(r["group"]==g for r in rows)})')
 for i,((lo,hi),v) in enumerate(zip(pairs,diffs)):
  ax=axs[1,i];ax.set_visible(True);dm=ax.imshow(v,origin='lower',extent=extent,aspect='auto',norm=TwoSlopeNorm(vmin=-maxdiff,vcenter=0,vmax=maxdiff),cmap='RdBu_r',interpolation='nearest');ax.axvline(0,color='black',ls=':',lw=.7);ax.set_title(f'{groups[hi]} − {groups[lo]}')
 for ax in axs.flat:
  if ax.get_visible():ax.set_xlabel('Distance from NDR midpoint (bp)');ax.set_ylabel('Fragment length (bp)')
 fig.colorbar(im,ax=[x for x in axs[0] if x.get_visible()],label='Mean sample-wise min–max count',shrink=.8);fig.colorbar(dm,ax=[x for x in axs[1] if x.get_visible()],label='Difference in group means',shrink=.8)
 fig.suptitle('NDR V-plot: group means and descriptive differences',fontsize=15);fig.supxlabel('Descriptive only: differences are not pixel-wise significance tests.',fontsize=9)
 fig.savefig(out/'group_vplot_means_and_differences.png',dpi=220,bbox_inches='tight');fig.savefig(out/'group_vplot_means_and_differences.pdf',bbox_inches='tight');plt.close(fig)
 with (out/'sample_audit.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=audit[0],delimiter='\t');w.writeheader();w.writerows(audit)
 np.savez_compressed(out/'group_matrices.npz',group_names=np.array(groups),means=np.stack(means),difference_names=np.array([groups[j]+' - '+groups[i] for i,j in pairs]),differences=np.stack(diffs),x_edges=edges[0],length_edges=edges[1])
 (out/'comparison_metadata.json').write_text(json.dumps({'groups':{g:sum(r['group']==g for r in rows) for g in groups},'channel':a.channel,'normalization':'sample-wide min-max on raw_counts across all channels','bed_hashes':list(hashes),'bed_hash_unverified_samples':metadata_missing,'pairwise_statistics':'not performed','difference_color_limit':maxdiff,'mean_color_max':maxmean,'smoothing':False},indent=2)+'\n')
 if metadata_missing:print('WARNING: BED provenance unverified for',len(metadata_missing),'samples; inspect comparison_metadata.json')
 print('Saved:',out)
if __name__=='__main__':main()
