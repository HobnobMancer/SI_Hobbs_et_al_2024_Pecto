# run coinfinder to pectobacteriaceae, and the pectobacterium and dickeya data sets

mkdir results/pectobact/coinfinder

coinfinder \
    -i data/pecto_dic/cazomes/pd_fam_genomes_taxs \
    -p data/pecto_dic/tree/pecto_dic_bestTree_taxs \
    -a \
    -o results/pecto_dic/coinfinder/pectodic_coinfinder_taxs_rectangular \
    > results/pecto_dic/coinfinder/coinfinder_pectodic_taxs_run.log
