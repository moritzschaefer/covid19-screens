awk -F '\t' '{print $29}' E-MTAB-9638.sdrf.txt | xargs -n 1 -P 7 wget -r
