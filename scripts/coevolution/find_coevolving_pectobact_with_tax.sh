# run coinfinder to pectobacteriaceae, and the pectobacterium and dickeya data sets

mkdir results/pectobact/coinfinder

coinfinder \
    -i data/pectobact/cazomes/coinfinder_pecto_fam_genomes \
    -p data/pectobact/tree/ani_tree/pyani_ani_tree.new \
    -a \
    -o results/pectobact/coinfinder/pectobact_coinfinder_taxs_rectangular \
    > results/pectobact/coinfinder/coinfinder_pectobact_taxs_run.log