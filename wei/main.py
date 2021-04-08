#!/usr/bin/env python3

import subprocess

tests = {
        'wei_vero_v1_d10_lowMOI': ('Cas9-v1 D10 Mock', 'Cas9-v1 D10 Lo-MOI'),
        'wei_vero_v1_d10_highMOI': ('Cas9-v1 D10 Mock', 'Cas9-v1 D10 Hi-MOI'),
        'wei_vero_v2_d10_5e6_highMOI': ('Cas9-v2 D5 Mock', 'Cas9-v2 D10 5e6 Hi-MOI'),
        'wei_vero_v2_d5_5e6_highMOI': ('Cas9-v2 D5 Mock', 'Cas9-v2 D5 5e6 Hi-MOI'),
        'wei_vero_v2_d2_5e6_highMOI': ('Cas9-v2 D5 Mock', 'Cas9-v2 D2 5e6 Hi-MOI'),
        'wei_vero_v2_d5_2.5e6_lowMOI': ('Cas9-v2 D5 Mock', 'Cas9-v2 D5 2.5e6 Lo-MOI'),
        'wei_vero_v2_d5_2.5e6_highMOI': ('Cas9-v2 D5 Mock', 'Cas9-v2 D5 2.5e6 Hi-MOI')}
template = open('template.conf', 'r').read()
for name, samples in tests.items():
    subprocess.run(f'mageck test -k merged_counts_with_replaced_ensembl_mapping.tsv -t "{samples[1]}" -c "{samples[0]}" -n {name} --normcounts-to-file', shell=True)
    with open(f'{name}.vispr.yaml', 'w') as f:
        f.write(template.replace('myexperiment', name))
    print(f'{name} done. {samples}')
