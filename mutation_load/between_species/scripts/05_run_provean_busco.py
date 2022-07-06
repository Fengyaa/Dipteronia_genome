import sys
import os
import pandas as pd

OG_list = sys.argv[1]
provean_path = '/home/xubo/Software/provean-1.1.5/bin/provean.sh'

provean_out_df = pd.DataFrame(columns=['OG_name','HGVS','Provean_score'])
provean_out = open(OG_list+'provean_out',"w",buffering=1) 
with open(OG_list) as in_list:
    for line in in_list:
        OG_name = line.strip()
        provean_cmd =f'{provean_path} -q {OG_name}.anc.pep -v {OG_name}.var --num_threads 2'
        provean_output = os.popen(provean_cmd).readlines()
        try:
            header_index = provean_output.index('# VARIATION\tSCORE\n')
            provean_scores = [x.strip().split() for x in provean_output[header_index+1:]]
            for var_info in provean_scores:
                var_name = var_info[0]
                provean_score = float(var_info[1])
                provean_out.write(f'{OG_name}\t{var_name}\t{provean_score}\n')
                # provean_out_df = provean_out_df.append({
                #     'OG_name':OG_name, 
                #     'HGVS':var_name,
                #     'Provean_score': provean_score,
                #     },ignore_index=True)
            print(f'{OG_name} finished')
        except ValueError:
            print(f'{OG_name} failed')

# provean_out_df_indexed = provean_out_df.set_index(['OG_name','HGVS'])
# provean_out_df.to_csv(input_vep.replace("txt", "_provean_out.txt"),sep="\t")
