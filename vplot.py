#!/usr/bin/env python3
"""Methods-based cfDNA NDR midpoint/length histogram; not author-exact cfWGS preprocessing."""
import argparse,collections,hashlib,json,math,os,sys
from pathlib import Path
import numpy as np
import pysam

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def read_sites(specs,refs,half):
 names=[];sites=collections.defaultdict(list);audit=[]
 for ch,spec in enumerate(specs):
  name,path=spec.split('=',1)
  if name in names:raise ValueError('Duplicate channel: '+name)
  names.append(name);seen=set()
  with open(path) as f:
   for line in f:
    if not line.strip() or line.startswith(('#','track','browser')):continue
    a=line.split();chrom=a[0];start,end=map(int,a[1:3]);key=(chrom,start,end)
    if key in seen:raise ValueError('Duplicate BED interval: '+str(key))
    seen.add(key)
    if chrom not in refs:raise ValueError('BED contig absent in CRAM: '+chrom)
    if not 0<=start<end<=refs[chrom]:raise ValueError('Invalid BED coordinates: '+str(key))
    center=(start+end)/2;valid=center-half>=0 and center+half<=refs[chrom]
    audit.append([name,chrom,start,end,center,'included' if valid else 'excluded_boundary'])
    if valid:sites[chrom].append((center,ch))
  if not seen:raise ValueError('Empty BED: '+path)
 for chrom in sites:sites[chrom].sort()
 if set(ch for v in sites.values() for _,ch in v)!=set(range(len(names))):raise ValueError('Channel has no usable sites')
 return names,sites,audit

def windows(centers,half,pad,limit):
 out=[]
 for c in centers:
  s=max(0,math.floor(c-half-pad));e=min(limit,math.ceil(c+half+pad))
  if out and s<=out[-1][1]:out[-1][1]=max(out[-1][1],e)
  else:out.append([s,e])
 return out

def run(a):
 out=Path(a.output)
 if out.exists() and any(out.iterdir()):raise ValueError('Output must be new or empty: '+str(out))
 if not Path(a.reference+'.fai').exists():raise ValueError('Reference .fai required; input reference is not modified')
 if a.half_window<=0 or a.x_bin<=0 or 2*a.half_window%a.x_bin:raise ValueError('Window must be divisible by x-bin')
 if not 0<a.min_length<a.max_length or a.threads<1:raise ValueError('Invalid length or threads')
 mode='rc' if a.cram.lower().endswith('.cram') else 'rb'
 with pysam.AlignmentFile(a.cram,mode,reference_filename=a.reference,threads=a.threads) as bam:
  if not bam.has_index():raise ValueError('CRAM/BAM index required')
  refs=dict(zip(bam.references,bam.lengths))
  with pysam.FastaFile(a.reference) as fa:
   for chrom,n in refs.items():
    if chrom not in fa.references or fa.get_reference_length(chrom)!=n:raise ValueError('Reference dictionary mismatch: '+chrom)
  names,sites,audit=read_sites(a.ndr,refs,a.half_window)
  raw=np.zeros((a.max_length-a.min_length,2*a.half_window//a.x_bin,len(names)),dtype=np.int64)
  qc=collections.Counter();site_hits=collections.Counter()
  for chrom,ss in sites.items():
   centers=np.array([v[0] for v in ss]);channels=np.array([v[1] for v in ss]);qc['fetch_windows']+=len(windows(centers,a.half_window,a.max_length,refs[chrom]))
   for lo,hi in windows(centers,a.half_window,a.max_length,refs[chrom]):
    pending={}
    for r in bam.fetch(chrom,lo,hi):
     qc['fetched_alignment_records']+=1
     if r.is_unmapped or r.mate_is_unmapped or not r.is_paired or not r.is_proper_pair or r.reference_id!=r.next_reference_id:
      qc['reject_pair_flags']+=1;continue
     if r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail:
      qc['reject_alignment_flags']+=1;continue
     if r.mapping_quality<a.mapq:qc['reject_mapq']+=1;continue
     length=abs(r.template_length)
     if not a.min_length<=length<a.max_length:qc['reject_length']+=1;continue
     key=(r.query_name,r.get_tag('RG') if r.has_tag('RG') else '')
     mate=pending.pop(key,None)
     if mate is None:pending[key]=r;continue
     if mate.is_read1==r.is_read1 or mate.reference_start!=r.next_reference_start or r.reference_start!=mate.next_reference_start or mate.template_length!=-r.template_length:
      qc['reject_inconsistent_pair']+=1;continue
     start=min(r.reference_start,mate.reference_start);end=max(r.reference_end,mate.reference_end)
     if end-start!=length:qc['reject_span_tlen_mismatch']+=1;continue
     mid=(start+end)/2
     # x in [-half, half); include every overlapping NDR, once per fragment/site.
     left=np.searchsorted(centers,mid-a.half_window,side='right');right=np.searchsorted(centers,mid+a.half_window,side='right')
     if right==left:continue
     qc['fragments_assigned']+=1
     for idx in range(left,right):
      col=int(math.floor((mid-centers[idx]+a.half_window)/a.x_bin));ch=int(channels[idx])
      raw[length-a.min_length,col,ch]+=1;site_hits[(chrom,idx)]+=1;qc['fragment_site_assignments']+=1
    qc['unpaired_after_filters_or_fetch']+=len(pending)
   print('Processed',chrom,'assigned',qc['fragments_assigned'],flush=True)
 if not raw.sum():raise ValueError('No usable fragment/site counts; no successful output generated')
 out.mkdir(parents=True,exist_ok=True)
 norm=(raw-raw.min())/float(raw.max()-raw.min()) if raw.max()>raw.min() else np.zeros_like(raw,dtype=float)
 x=np.arange(-a.half_window,a.half_window+a.x_bin,a.x_bin);y=np.arange(a.min_length,a.max_length+1)
 np.savez_compressed(out/'vplot.npz',raw_counts=raw,minmax=norm,x_edges=x,length_edges=y,channels=np.array(names))
 with (out/'regions.tsv').open('w') as f:
  f.write('channel\tchrom\tstart\tend\tcenter\tstatus\n');f.writelines('\t'.join(map(str,r))+'\n' for r in audit)
 qc['regions_total']=len(audit);qc['regions_used']=sum(r[-1]=='included' for r in audit);qc['regions_with_counts']=len(site_hits)
 for ch,name in enumerate(names):qc['assignments_'+name]=int(raw[:,:,ch].sum())
 meta=dict(parameters=vars(a),qc=dict(qc),shape=list(raw.shape),bed_sha256={s:sha(s.split('=',1)[1]) for s in a.ndr},reference_fai_sha256=sha(a.reference+'.fai'),input_stat={'size':os.stat(a.cram).st_size,'mtime_ns':os.stat(a.cram).st_mtime_ns},versions={'pysam':pysam.__version__,'numpy':np.__version__},doi='10.1038/s41467-023-37768-3',limitations=['Midpoint counting is an explicit interpretation; author fragment-halving code unavailable.','Default +/-1000 bp at 4 bp gives 500 columns, not reported 250.','Length interval is [100,200) by default; author exact boundary convention unavailable.','hg38 inputs are an adaptation of the hg19 publication.','No GC correction, smoothing, classifier training, or tissue-specific NDR selection performed.','Sample-wide min-max follows paper prose, not verified author model preprocessing.'])
 (out/'run.json').write_text(json.dumps(meta,indent=2))
 if not a.no_plot:
  import matplotlib;matplotlib.use('Agg')
  import matplotlib.pyplot as plt
  for ch,name in enumerate(names):
   fig,ax=plt.subplots(figsize=(10,4.8),layout='constrained')
   im=ax.imshow(norm[:,:,ch],origin='lower',aspect='auto',extent=[x[0],x[-1],y[0],y[-1]],interpolation='nearest',cmap='viridis',vmin=0,vmax=1)
   ax.set(xlabel='Fragment midpoint distance from NDR midpoint (bp)',ylabel='Fragment length (bp)',title=f'{a.sample} | {name} | {sum(r[0]==name and r[-1]=="included" for r in audit):,} NDRs')
   fig.colorbar(im,ax=ax,label='Sample-wide min–max normalized count')
   fig.savefig(out/f'vplot_channel{ch+1:02d}.png',dpi=300);plt.close(fig)
 (out/'COMPLETE.json').write_text(json.dumps({'status':'complete','assignments':int(raw.sum()),'shape':list(raw.shape)}))

def parser():
 p=argparse.ArgumentParser(description=__doc__)
 for k in ['cram','reference','output','sample']:p.add_argument('--'+k,required=True)
 p.add_argument('--ndr',action='append',required=True,help='NAME=BED; repeat for multiple tissue channels; same assembly as CRAM')
 p.add_argument('--half-window',type=int,default=1000);p.add_argument('--x-bin',type=int,default=4)
 p.add_argument('--min-length',type=int,default=100);p.add_argument('--max-length',type=int,default=200,help='Exclusive upper boundary')
 p.add_argument('--mapq',type=int,default=30);p.add_argument('--threads',type=int,default=4);p.add_argument('--no-plot',action='store_true')
 return p
if __name__=='__main__':run(parser().parse_args())
