# coding: utf-8
import pandas as pd
from pathlib import Path
import os, shutil
import subprocess
import re
import yaml
raw = Path('download_raw/')
files = list(raw.glob('*/*.fastq.gz'))
samples = [re.search('/([^/]+_S[0-9]+)[^/]+.fastq.gz', str(f)).groups()[0] for f in files]
df=pd.DataFrame({'files': [str(f) for f in files], 'samples': samples})
df['is_R1'] = df['files'].astype(str).str.contains('_R1_')

output = {'samples': {group: files[files.is_R1]['files'].tolist() for group, files in df.groupby('samples')},
          'paired': {group: files[~files.is_R1]['files'].tolist() for group, files in df.groupby('samples')}
          }

with open('config.yaml', 'w') as f:
    yaml.dump(output, f)

# df.groupby('samples')['files'].apply(lambda group: subprocess.run(f'cat {" ".join([str(f) for f in group])} > fastq/{group.name}.fastq.gz', shell=True))
