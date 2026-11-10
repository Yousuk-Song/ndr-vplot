import tempfile,subprocess,sys,json
from pathlib import Path
import numpy as np,pysam
P=Path(__file__).resolve().parent
with tempfile.TemporaryDirectory() as td:
 d=Path(td);ref=d/'ref.fa';ref.write_text('>chr1\n'+'A'*10000+'\n');pysam.faidx(str(ref))
 header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':10000}]};records=[]
 for name,start,L,q2,dup in [('ok100',2950,100,60,False),('ok199',2900,199,60,False),('exclude200',2900,200,60,False),('lowmate',2950,100,10,False),('duplicate',2950,100,60,True)]:
  for read1 in [True,False]:
   r=pysam.AlignedSegment();r.query_name=name;r.query_sequence='A'*50;r.query_qualities=pysam.qualitystring_to_array('I'*50)
   r.flag=(99 if read1 else 147)+(1024 if dup else 0);r.reference_id=0;r.reference_start=start if read1 else start+L-50;r.mapping_quality=60 if read1 else q2;r.cigarstring='50M';r.next_reference_id=0;r.next_reference_start=start+L-50 if read1 else start;r.template_length=L if read1 else -L;records.append(r)
 cram=d/'test.cram'
 with pysam.AlignmentFile(str(cram),'wc',header=header,reference_filename=str(ref)) as f:
  for r in sorted(records,key=lambda r:r.reference_start):f.write(r)
 pysam.index(str(cram));bed=d/'sites.bed';bed.write_text('chr1\t2990\t3010\nchr1\t3090\t3110\n')
 cmd=[sys.executable,str(P/'vplot.py'),'--cram',str(cram),'--reference',str(ref),'--ndr','test='+str(bed),'--sample','synthetic','--output',str(d/'out')]
 subprocess.run(cmd,check=True)
 z=np.load(d/'out/vplot.npz');a=z['raw_counts'];assert a.shape==(100,500,1);assert a.sum()==4;assert a[0].sum()==2;assert a[99].sum()==2
 assert np.allclose(z['minmax'], (a-a.min())/(a.max()-a.min()))
 assert z['x_edges'][0]==-1000 and z['x_edges'][-1]==1000
 assert z['length_edges'][0]==100 and z['length_edges'][-1]==200
 assert (d/'out/COMPLETE.json').exists()
 assert a[0,250,0]==1 and a[0,225,0]==1
 q=json.load(open(d/'out/run.json'))['qc'];assert q['fragments_assigned']==2
 assert (d/'out/vplot_channel01.png').stat().st_size>1000
 # Existing outputs cannot be silently overwritten.
 assert subprocess.run(cmd,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL).returncode!=0
 print('PASS: indexed CRAM decoding, paired-once counting, mate MAPQ, duplicates, length boundaries, overlapping regions, midpoint bins, plotting, overwrite guard')
