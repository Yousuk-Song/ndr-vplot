import csv
import re
from pathlib import Path

BASE = Path(workflow.basedir)
# Supply --configfile explicitly; relative paths are relative to the working directory.
required = ['samples', 'reference', 'output_dir']
for key in required:
    if key not in config:
        raise ValueError('Missing config: ' + key + '; use --configfile config/config.yaml')
with open(config['samples'], newline='') as handle:
    rows = list(csv.DictReader(handle, delimiter='\t'))
SAMPLES = {}
for row in rows:
    if not {'sample', 'cram'} <= row.keys():
        raise ValueError('samples TSV requires sample and cram columns')
    name = row['sample']
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', name) or name in SAMPLES:
        raise ValueError('Invalid or duplicate sample: ' + name)
    SAMPLES[name] = row['cram']
if not SAMPLES:
    raise ValueError('Empty sample manifest')
REF = config['reference']
BED = config.get('bed', str(BASE / 'resources/STAD_NDR_hg38_autosomes.bed'))
OUT = str(Path(config['output_dir']).absolute())
PYTHON = config.get('python', 'python3')

def index(wc):
    cram = Path(SAMPLES[wc.sample])
    candidates = [Path(str(cram)+'.crai'), cram.with_suffix('.crai')]
    return str(next((p for p in candidates if p.exists()), candidates[0]))

rule all:
    input:
        expand(OUT+'/{sample}/'+f, sample=list(SAMPLES))
        for f in ['vplot.npz', 'regions.tsv', 'run.json', 'vplot_channel01.png', 'COMPLETE.json']

rule vplot:
    input:
        cram=lambda wc: SAMPLES[wc.sample],
        crai=index,
        reference=REF,
        fai=REF+'.fai',
        bed=BED,
        engine=str(BASE/'vplot.py'),
        runner=str(BASE/'scripts/run_sample.py')
    output:
        npz=OUT+'/{sample}/vplot.npz',
        regions=OUT+'/{sample}/regions.tsv',
        metadata=OUT+'/{sample}/run.json',
        png=OUT+'/{sample}/vplot_channel01.png',
        complete=OUT+'/{sample}/COMPLETE.json'
    threads: int(config.get('threads_per_sample', 4))
    resources:
        mem_mb=int(config.get('mem_mb_per_sample', 8000))
    params:
        python=PYTHON,
        out=lambda wc: OUT+'/'+wc.sample,
        half=int(config.get('half_window', 1000)),
        bin=int(config.get('x_bin', 4)),
        minimum=int(config.get('min_length', 100)),
        maximum=int(config.get('max_length', 200)),
        mapq=int(config.get('mapq', 30))
    log: OUT+'/logs/{sample}.log'
    benchmark: OUT+'/benchmarks/{sample}.tsv'
    shell:
        '''
        {params.python:q} {input.runner:q} --engine {input.engine:q} \
          --sample {wildcards.sample:q} --cram {input.cram:q} \
          --reference {input.reference:q} --bed {input.bed:q} --output {params.out:q} \
          --threads {threads} --half-window {params.half} --x-bin {params.bin} \
          --min-length {params.minimum} --max-length {params.maximum} --mapq {params.mapq} \
          > {log:q} 2>&1
        '''
