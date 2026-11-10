#!/usr/bin/env python3
"""Reusable Snakemake entry point: external run config + sample manifest."""
import argparse,csv,hashlib,json,re,shutil,subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parent

def write_changed(path,text):
    if not path.exists() or path.read_text()!=text: path.write_text(text)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',required=True,type=Path)
    p.add_argument('--run',action='store_true',help='Default is dry-run')
    a=p.parse_args(); cfgpath=a.config.resolve(); c=json.loads(cfgpath.read_text())
    def path(value):
        v=Path(value).expanduser();return v if v.is_absolute() else cfgpath.parent/v
    manifest=path(c['samples']);reference=path(c['reference']);out=path(c['output_dir']).resolve()
    cores=int(c.get('cores',8));threads=int(c.get('threads_per_sample',4));jobs=int(c.get('jobs',2))
    memory=int(c.get('memory_mb',16000));per_mem=int(c.get('mem_mb_per_sample',8000))
    if min(cores,threads,jobs,memory,per_mem)<1 or per_mem>memory:raise ValueError('Invalid resource settings')
    rows=list(csv.DictReader(manifest.open(),delimiter='\t'));seen=set();ready=[];audit=[]
    for r in rows:
        s=r['sample']; cr=Path(r['cram']).expanduser()
        if not cr.is_absolute():cr=manifest.parent/cr
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',s) or s in seen:raise ValueError('Invalid/duplicate sample '+s)
        seen.add(s)
        ix=next((v for v in [Path(str(cr)+'.crai'),cr.with_suffix('.crai')] if v.is_file() and v.stat().st_size>0),None)
        status='ready' if cr.is_file() and cr.stat().st_size>0 and ix else 'missing_cram_or_index'
        audit.append([s,str(cr),str(ix or ''),status])
        if status=='ready':ready.append([s,str(cr.resolve())])
    if not rows:raise ValueError('Empty manifest')
    state=out/'workflow_inputs';state.mkdir(parents=True,exist_ok=True)
    with (state/'input_audit.tsv').open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['sample','cram','index','status']);w.writerows(audit)
    print(f'Selected={len(rows)} Ready={len(ready)} Missing={len(rows)-len(ready)}',flush=True)
    print('Audit:',state/'input_audit.tsv',flush=True)
    if not ready:raise SystemExit('No runnable samples. No analysis started.')
    if not reference.is_file() or not Path(str(reference)+'.fai').is_file():raise SystemExit('Reference FASTA and .fai required')
    exe=shutil.which(c.get('snakemake','snakemake'))
    if not exe:raise SystemExit('Snakemake not found; activate runtime or configure snakemake path')
    python=c.get('python',sys.executable)
    subprocess.run([python,'-c','import numpy,pysam,matplotlib'],check=True)
    bed=path(c['bed']) if c.get('bed') else BASE/'resources/STAD_NDR_hg38_autosomes.bed'
    if not bed.is_file():raise SystemExit('BED missing')
    selected=state/'samples.tsv'
    write_changed(selected,'sample\tcram\n'+''.join('\t'.join(r)+'\n' for r in ready))
    resolved={k:c[k] for k in ['half_window','x_bin','min_length','max_length','mapq'] if k in c}
    resolved.update(samples=str(selected),reference=str(reference.resolve()),output_dir=str(out/'samples'),python=python,bed=str(bed.resolve()),threads_per_sample=threads,mem_mb_per_sample=per_mem)
    config=state/'config.json';write_changed(config,json.dumps(resolved,indent=2))
    (state/'run_request.json').write_text(json.dumps({'requested_config':c,'pipeline':str(BASE),'version':(BASE/'VERSION').read_text().strip(),'manifest_sha256':hashlib.sha256(manifest.read_bytes()).hexdigest()},indent=2))
    cmd=[exe,'--snakefile',str(BASE/'Snakefile'),'--directory',str(out),'--configfile',str(config),'--cores',str(cores),'--jobs',str(jobs),'--resources',f'mem_mb={memory}','--rerun-incomplete','--keep-going']
    # NDR_SNAKEMAKE_CACHE_COMPAT_V1: optional flag is absent in older releases.
    help_result=subprocess.run([exe,'--help'],capture_output=True,text=True,check=True)
    if '--runtime-source-cache-path' in help_result.stdout + help_result.stderr:
        cmd.extend(['--runtime-source-cache-path',str(out/'.snakemake/source-cache')])
    if not a.run:cmd.append('--dry-run')
    raise SystemExit(subprocess.run(cmd).returncode)
if __name__=='__main__':main()
