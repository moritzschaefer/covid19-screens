from pathlib import Path
import pandas as pd 


count_files = [str(p) for p in Path('.').glob('*/workflow/results/count/all.count.txt')]
normalized_count_files = [str(p) for p in Path('.').glob('*/workflow/results/count/all.count_normalized.txt')]
# print(count_files)

# make separate columns for vero:
human_count_files = [f for f in count_files if not f.startswith('wei') and not f.startswith('vero')]
normalized_human_count_files = [f for f in normalized_count_files if not f.startswith('wei') and not f.startswith('vero')]


vero_count_files = [f for f in count_files if f.startswith('wei') or f.startswith('vero')]
normalized_vero_count_files = [f for f in normalized_count_files if f.startswith('wei') or f.startswith('vero')]

for files, output in [
    # (human_count_files, 'supp_data_human_read_counts.xlsx'),
    # (normalized_human_count_files, 'supp_data_human_norm_read_counts.xlsx'),
    (vero_count_files, 'supp_data_vero_read_counts.xlsx'),
    # (normalized_vero_count_files, 'supp_data_vero_norm_read_counts.xlsx')
    ]:
    with pd.ExcelWriter(output) as writer:  
        for f in files:
            pd.read_csv(f, sep='\t').to_excel(writer, sheet_name=f[:f.find('/')], index=False)
