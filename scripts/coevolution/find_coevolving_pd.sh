# run coinfinder to pectobacteriaceae, and the pectobacterium and dickeya data sets

mkdir results/pectobact/coinfinder

coinfinder \
    -i data/pecto_dic/cazomes/pd_fam_genomes \
    -p data/pecto_dic/tree/pecto_dic_bestTree \
    -a \
    -o results/pecto_dic/coinfinder/pectodic_coinfinder_ \
    > results/pecto_dic/coinfinder/pectodic_pectobact_run.log
