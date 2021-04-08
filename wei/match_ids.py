# coding: utf-8
from Bio import SeqIO
import modin.pandas as pd
import re


def parse_gene_name(description):
    try:
        return re.search('gene_symbol:(\w+) ', description).groups()[0]
    except (IndexError, TypeError, AttributeError):
        return re.search('gene:(ENSCSAG[0-9.]+) ', description).groups()[0]
        
x=SeqIO.parse('./Chlorocebus_sabaeus.ChlSab1.1.cdna.abinitio.fa', 'fasta')

import pdb; pdb.set_trace()
seqs = {s.name: str(s.seq) for s in x}
    
merged_counts = pd.read_csv('./sgrna_counts_merged.tsv', sep='\t')
matches = merged_counts['sgRNA'].apply(lambda sg: [gene for gene, seq in seqs.items() if re.search(f'{sg}.GG', seq)])
merged_counts['new_gene'] = matches.apply(lambda row: re.match('[^.]+', row[0]).group() if len(row) == 1 else (None if len(row) == 0 else f'multiKO_{"_".join([re.match("[^.]+", v).group() for v in row])}'))

