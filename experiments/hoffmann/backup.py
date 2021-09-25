# coding: utf-8
pd
import pandas as pd
df = pd.read_csv('SraRunTable.txt', sep='\t')
df.head()
df = pd.read_csv('SraRunTable.txt', sep=',')
df.head()
df.iloc[0]
get_ipython().system('ls fastq')
df.iloc[0]
get_ipython().system('ls fastq')
df.iloc[0]
df.iloc[1]
df.iloc[2]
df.iloc[1]
df.iloc[0]
df.iloc[0]
df.iloc[1]
df.iloc[2]
df.iloc[3]
df.iloc[2]
df.iloc[1]
df.iloc[0]
df.iloc[3]
df.iloc[4]
df.iloc[3]
df.iloc[0]
df.iloc[3]
df.iloc[4]
df.iloc[3]
df.iloc[4]
df.shape
df.groupby(['source_name', 'Temperature', 'Treatment', 'Replicate']).count()
df.iloc[4]
dfsars = df[df.Treatment=='SARS-CoV-2']
dfsars.head()
dfsars.shape
dfsars.groupby(['source_name', 'Temperature', 'Treatment', 'Replicate']).count()
dfsars = df[df.Treatment.isin('mock', 'SARS-CoV-2')]
dfsars = df[df.Treatment.isin(['mock', 'SARS-CoV-2'])]
dfsars.groupby(['source_name', 'Temperature', 'Treatment', 'Replicate']).count()
dfsars.head()
get_ipython().set_next_input("dfsars.groupby(['source_name', 'Temperature', 'Treatment', 'Replicate']).apply");get_ipython().run_line_magic('pinfo', 'apply')
get_ipython().set_next_input("dfsars.groupby(['source_name', 'Temperature', 'Treatment', 'Replicate']).apply");get_ipython().run_line_magic('pinfo', 'apply')
dfsars.groupby(['source_name', 'Temperature', 'Treatment', 'Replicate'])
def combine_fastqs(group):
    import subprocess
    filenames = [f'fastq/{g["Run"]}.fastq.gz' for g in group]
    cmd = f'cat {" ".join(filenames)} > final/{group.name}.fastq.gz'
    print(cmd)
    
dfsars.groupby(['source_name', 'Temperature', 'Treatment', 'Replicate']).apply(combine_fastqs)
dfsars.groupby(['source_name', 'Temperature', 'Treatment', 'Replicate']).aggregate(combine_fastqs)
groups=dfsars.groupby(['source_name', 'Temperature', 'Treatment', 'Replicate'])
next(groups)
next(iter(groups))
def combine_fastqs(group):
    import subprocess
    print(group)
    filenames = [f'fastq/{g["Run"]}.fastq.gz' for g in group]
    cmd = f'cat {" ".join(filenames)} > final/{group.name}.fastq.gz'
    print(cmd)
    
dfsars.groupby(['source_name', 'Temperature', 'Treatment', 'Replicate']).aggregate(combine_fastqs)
def combine_fastqs(group):
    import subprocess
    print(group)
    print('lala')
    filenames = [f'fastq/{g["Run"]}.fastq.gz' for g in group]
    cmd = f'cat {" ".join(filenames)} > final/{group.name}.fastq.gz'
    print(cmd)
    
dfsars.groupby(['source_name', 'Temperature', 'Treatment', 'Replicate']).aggregate(combine_fastqs)
get_ipython().run_line_magic('save', 'backup.py 1-55')
l
dfsars.groupby(['source_name', 'Temperature', 'Treatment', 'Replicate']).aggregate(combine_fastqs)
def combine_fastqs(group):
    filenames = [f'fastq/{g["Run"]}.fastq.gz' for i, g in group.iterrows()]
    cmd = f'cat {" ".join(filenames)} > final/{group.name}.fastq.gz'
    print(cmd)
    
dfsars.groupby(['source_name', 'Temperature', 'Treatment', 'Replicate']).aggregate(combine_fastqs)
dfsars.groupby(['source_name', 'Temperature', 'Treatment', 'Replicate']).apply(combine_fastqs)
def combine_fastqs(group):
    filenames = [f'fastq/{g["Run"]}.fastq.gz' for i, g in group.iterrows()]
    cmd = f'cat {" ".join(filenames)} > final/{"hoffmann_" + "_".join(group.name)}.fastq.gz'
    print(cmd)
    
dfsars.groupby(['source_name', 'Temperature', 'Treatment', 'Replicate']).apply(combine_fastqs)
def combine_fastqs(group):
    filenames = [f'fastq/{g["Run"]}.fastq.gz' for i, g in group.iterrows()]
    cmd = f'cat {" ".join(filenames)} > final/{"hoffmann_" + "_".join([str(v) for v in group.name])}.fastq.gz'
    print(cmd)
    
dfsars.groupby(['source_name', 'Temperature', 'Treatment', 'Replicate']).apply(combine_fastqs)
get_ipython().run_line_magic('ls', 'fastq')
get_ipython().run_line_magic('mkdir', 'final')
dfsars.groupby(['source_name', 'Temperature', 'Treatment', 'Replicate']).apply(combine_fastqs)
def combine_fastqs(group):
    filenames = [f'fastq/{g["Run"]}.fastq.gz' for i, g in group.iterrows()]
    cmd = f'cat {" ".join(filenames)} > final/{"hoffmann_" + "_".join([str(v) for v in group.name])}.fastq.gz'
    subprocess.run(cmd, shell=True)
    
    
dfsars.groupby(['source_name', 'Temperature', 'Treatment', 'Replicate']).apply(combine_fastqs)
import subprocess
dfsars.groupby(['source_name', 'Temperature', 'Treatment', 'Replicate']).apply(combine_fastqs)
def combine_fastqs(group):
    filenames = [f'fastq/{g["Run"]}.fastq.gz' for i, g in group.iterrows()]
    cmd = f'cat {" ".join(filenames)} > final/{"hoffmann_" + "_".join([str(v) for v in group.name])}.fastq.gz'
    print(cmd)
    subprocess.run(cmd, shell=True)
    
    
dfsars.groupby(['source_name', 'Temperature', 'Treatment', 'Replicate']).apply(combine_fastqs)
def combine_fastqs(group):
    filenames = [f'fastq/{g["Run"]}.fastq.gz' for i, g in group.iterrows()]
    cmd = f'cat {" ".join(filenames)} > final/{"hoffmann_" + "_".join([str(v).replace(' ', '-') for v in group.name])}.fastq.gz'
    print(cmd)
    subprocess.run(cmd, shell=True)
    
    
def combine_fastqs(group):
    filenames = [f'fastq/{g["Run"]}.fastq.gz' for i, g in group.iterrows()]
    cmd = f'cat {" ".join(filenames)} > final/{"hoffmann_" + "_".join([str(v).replace(" ", "-") for v in group.name])}.fastq.gz'
    print(cmd)
    subprocess.run(cmd, shell=True)
    
    
dfsars.groupby(['source_name', 'Temperature', 'Treatment', 'Replicate']).apply(combine_fastqs)
get_ipython().run_line_magic('save', 'backup.py 1-77')
