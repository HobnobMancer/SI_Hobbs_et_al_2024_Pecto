# Supplementary Information: Hobbs _et al._, 2023

## Association between niche adaptation and evolution of carbohydrate active enzymes in _Pectobacteriaceae_

This repository contains supplementary information for analyses reported in Hobbs et al. (2023), describing the evolution of the CAZome.

Run all commands provided in the walkthrough from the root of this directory.

## Online supplementary

Owing to the size of the data sets used, the figures are consequently compressed in the final manuscript. This remote repository contains the original full size, high resolution figures.

Additionally, some analyses are only briefly mentioned in the manuscript. The full method and results of these analyses are stored in this repository.

For the complete analysis and results for each data set, follow the links:
* [_Pectobacteriaceae_ GenBank genomes](https://hobnobmancer.github.io/SI_Hobbs_et_al_2023_Pecto/notebooks/explore_pectobact_cazomes.html)
* [_Pectobacterium_ and _Dickeya_ reference genomes](https://hobnobmancer.github.io/SI_Hobbs_et_al_2023_Pecto/notebooks/explore_pecto_dic_cazomes.html)
* [Intracellular and extracellular CAZyme](https://hobnobmancer.github.io/SI_Hobbs_et_al_2023_Pecto/notebooks/explore_IE_cazomes.html)

Find a full list of the results [here](#results).

## How to use this repository.

You can use this repository like a website, to browse and see how we performed the analysis, or you can download it to inspect, verify, reproduce and build on our analysis.

## Downloading this repository

You can use `git` to _clone_ this repository to your local hard drive:

```bash
git clone git@github.com:HobnobMancer/SI_Hobbs_et_al_2023_Pecto.git
```

Or you can download it as a compressed `.zip` archive from [this link](https://github.com/HobnobMancer/SI_Hobbs_et_al_2023_Pecto/archive/refs/heads/master.zip).

## If you have problems with this repository

Please raise an issue at the corresponding `GitHub` page:

* [Issues for this repository](https://github.com/HobnobMancer/SI_Hobbs_et_al_2023_Pecto/issues)

## Repo structure

The structure of this repository:

```bash
.
├── _config.yml
├── data
│   ├── bac120_taxonomy.tsv
│   ├── genomic_accessions
│   │   ├── genomes_for_coinfinder.txt
│   │   ├── pectobact_accessions
│   │   └── pecto_dic_accessions
│   ├── pectobact
│   │   ├── cazomes
│   │   │   ├── coinfinder_pecto_fam_genomes
│   │   │   ├── coinfinder_pecto_fam_genomes_taxs
│   │   │   ├── pecto_fam_genomes
│   │   │   ├── pecto_fam_genomes_proteins
│   │   │   ├── pecto_fam_genomes_proteins_taxs
│   │   │   └── pecto_fam_genomes_taxs
│   │   ├── genomes
│   │   │   ├── classes.txt
│   │   │   └── labels.txt
│   │   ├── missing_genomes
│   │   ├── proteomes
│   │   └── tree
│   │       └── ani_tree
│   │           ├── anim_matrix.tab
│   │           ├── genomes.tab
│   │           ├── logs
│   │           │   ├── pyani_01_createdb.log
│   │           │   ├── pyani_02_index.log
│   │           │   ├── pyani_03_anim.log
│   │           │   ├── pyani_04_plot.log
│   │           │   ├── pyani_05_reports.log
│   │           │   └── temp.log
│   │           ├── matrix_aln_lengths_4.tab
│   │           ├── matrix_aln_lengths_run4.pdf
│   │           ├── matrix_coverage_4.tab
│   │           ├── matrix_coverage_run4.pdf
│   │           ├── matrix_hadamard_4.tab
│   │           ├── matrix_hadamard_run4.pdf
│   │           ├── matrix_identity_4.tab
│   │           ├── matrix_identity_run4.pdf
│   │           ├── matrix_sim_errors_4.tab
│   │           ├── matrix_sim_errors_run4.pdf
│   │           ├── pyani_ani_tree.new
│   │           ├── pyani_ani_tree_taxs.new
│   │           └── scatter_identity_vs_coverage_run4.pdf
│   ├── pecto_dic
│   │   ├── cazomes
│   │   │   ├── pd_fam_genomes
│   │   │   ├── pd_fam_genomes_proteins
│   │   │   ├── pd_fam_genomes_proteins_taxs
│   │   │   └── pd_fam_genomes_taxs
│   │   └── tree
│   │       └── pecto_dic_bestTree
│   └── README.md
├── LICENSE
├── notebooks
│   ├── explore_pectobact_cazomes.ipynb
│   └── explore_pecto_dic_cazomes.ipynb
├── README.md
├── requirements.txt
├── results
│   ├── pectobact
│   │   ├── cazome_size
│   │   │   └── cazome_sizes.csv
│   │   ├── cazy_classes
│   │   │   └── cazy_class_sizes.csv
│   │   ├── cazy_families
│   │   │   ├── cazy_fam_freqs.csv
│   │   │   ├── fam_freq_clustermap.svg
│   │   │   ├── genus_species_fam_freq_clustermap.svg
│   │   │   ├── paper_fam_freq_clustermap_FILTERED.svg
│   │   │   ├── paper_fam_freq_clustermap.png
│   │   │   ├── paper_fam_freq_clustermap.svg
│   │   │   ├── paper_genus_species_fam_freq_clustermap.png
│   │   │   └── paper_pheno_genus_fam_freq_clustermap.png
│   │   ├── coinfinder
│   │   │   ├── pectobact_coinfinder__components.tsv
│   │   │   ├── pectobact_coinfinder__edges.tsv
│   │   │   ├── pectobact_coinfinder__heatmap0.pdf
│   │   │   ├── pectobact_coinfinder__heatmap_annotated.svg
│   │   │   ├── pectobact_coinfinder__network.gexf
│   │   │   ├── pectobact_coinfinder__nodes_in.csv
│   │   │   ├── pectobact_coinfinder__nodes.tsv
│   │   │   ├── pectobact_coinfinder_output.png
│   │   │   ├── pectobact_coinfinder__pairs.tsv
│   │   │   ├── pectobact_coinfinder_taxs_rectangular__components.tsv
│   │   │   ├── pectobact_coinfinder_taxs_rectangular__edges.tsv
│   │   │   ├── pectobact_coinfinder_taxs_rectangular__heatmap0.pdf
│   │   │   ├── pectobact_coinfinder_taxs_rectangular__network.gexf
│   │   │   ├── pectobact_coinfinder_taxs_rectangular__nodes_in.csv
│   │   │   ├── pectobact_coinfinder_taxs_rectangular__nodes.tsv
│   │   │   ├── pectobact_coinfinder_taxs_rectangular__pairs.tsv
│   │   │   ├── pectobact_coinfinder_taxs_rectangular__tree25.pdf
│   │   │   └── pectobact_coinfinder__tree25.pdf
│   │   ├── cooccurring_families
│   │   │   ├── cooccurring_fams_freqs.csv
│   │   │   ├── fam_corr_M_filled.csv
│   │   │   ├── paper-cooccurring_fams_freqs.csv
│   │   │   ├── paper-pecto-cooccurring-families.svg
│   │   │   └── pecto-cooccurring-families.svg
│   │   ├── core_cazome
│   │   │   ├── core_cazome_freqs.csv
│   │   │   ├── genera_core_cazome.svg
│   │   │   ├── genera_soft_hard_core_cazome.svg
│   │   │   └── unique_grp_fams.tsv
│   │   └── pca
│   │       ├── pc1_pc2_loadings_plot.png
│   │       ├── pca_explained_variance.png
│   │       ├── pca_pc1_vs_pc2_musicola_annotated.png
│   │       ├── pca_pc1_vs_pc2.png
│   │       └── pectobact_pca_pca_scree.png
│   └── pecto_dic
│       ├── cazy_classes
│       │   └── cazy_class_sizes.csv
│       ├── cazy_families
│       │   ├── cazy_fam_freqs.csv
│       │   ├── paper_pd_fam_freq_clustermap.png
│       │   ├── paper_pd_genus_species_fam_freq_clustermap.png
│       │   ├── pd_fam_freq_clustermap.svg
│       │   └── pd_genus_species_fam_freq_clustermap.svg
│       ├── cooccurring_families
│       │   ├── cooccurring_fams_freqs.csv
│       │   ├── fam_corr_M_filled.csv
│       │   └── pd_cooccurring_fams.svg
│       ├── core_cazome
│       │   ├── core_cazome_freqs.csv
│       │   ├── pd_all_core_cazome_freqs.csv
│       │   └── pd_core_cazome_freqs.csv
│       └── pca
│           ├── dickeya
│           │   ├── dic_explained_variance.png
│           │   ├── dic_pc1_pc2.png
│           │   ├── dic_pc1_pc3.png
│           │   ├── dic_pc1_pc4.png
│           │   ├── dic_pc2_pc3.png
│           │   ├── dic_pc2_pc4.png
│           │   ├── dic_pc3_pc4_LABELLED.png
│           │   ├── dic_pc3_pc4.png
│           │   └── dic_scree_plot.png
│           ├── pecto
│           │   ├── pecto_explained_variance.png
│           │   ├── pecto_pc1_pc2_LABELLED.png
│           │   ├── pecto_pc1_pc2.png
│           │   ├── pecto_pc1_pc3.png
│           │   ├── pecto_pc1_pc4.png
│           │   ├── pecto_pc2_pc3.png
│           │   ├── pecto_pc2_pc4.png
│           │   ├── pecto_pc3_pc4_LABELLED.png
│           │   ├── pecto_pc3_pc4.png
│           │   └── pecto_scree_plot.png
│           └── pecto_dic
│               ├── pd_pca_explained_variance.png
│               └── pd_pca_scree.png
└── scripts
    ├── annotate_cazome
    │   ├── get_cazy_cazymes_pd.sh
    │   ├── get_cazy_cazymes_pecto.sh
    │   ├── get_dbcan_cazymes_pectobact.sh
    │   ├── get_dbcan_cazymes_pecto_dic.sh
    │   ├── run_dbcan_pectobact.sh
    │   └── run_dbcan_pecto_dic.sh
    ├── coevolution
    │   ├── find_coevolving.sh
    │   ├── pd_network.R
    │   ├── pectobact_circular_network.R
    │   ├── pectobact_taxs_rectangular_network.R
    │   └── remove_excluded_genomes.py
    ├── download
    │   ├── annotate_genomes-copy.sh
    │   ├── annotate_genomes.sh
    │   ├── build_cazyme_database.sh
    │   ├── download_pd_genomes.sh
    │   ├── download_pecto_genomes.sh
    │   ├── download_same_pecto_genomes.sh
    │   ├── ident_missing_proteomes.py
    │   └── install_pyani_v0-3x.sh
    ├── README.md
    ├── signalp
    │   ├── gather_cazyme_seqs.py
    │   ├── get_ie_cazymes.py
    │   └── run_signalp.sh
    ├── taxs
    │   ├── add_ani_tax.py
    │   ├── add_tax_phylotree.py
    │   └── add_taxs.sh
    └── tree
        ├── analysis
        │   └── build_tanglegrams.R
        ├── ani
        │   ├── build_anim_tree.R
        │   ├── build_anim_tree.sh
        │   ├── parse_anim_tab.py
        │   └── run_anim.sh
        ├── phylo
        │   ├── align_scos.sh
        │   ├── annotate_genomes.sh
        │   ├── backtranslate.sh
        │   ├── concatenate_cds.py
        │   ├── extract_cds.py
        │   ├── find_orthologues.sh
        │   └── raxml_ng_build_tree.sh
        └── README.md
```

## Set up and reproducing analyses (quickstart)

You can use this archive to browse, validate, reproduce, or build on the phylogenomics analysis for the Hobbs et al. (2023) manuscript.

We recommend creating a conda environment specific for this activity, for example using the commands:
```bash
conda create -n pectobacteriaceae python=3.9 -y
conda activate pectobacteriaceae
conda install --file requirements.txt -y -c bioconda -c conda-forge -c predector
```

To use `pyani` in this analysis, version 0.3+ must be installed. At the time of development, `pyani` v0.3+ must be installed from `source`, this can be done by using the bash script `install_pyani_v0-3x.sh` (run from the root of this repository):
```bash
scripts/download/install_pyani_v0-3x.sh
```

### Install dbCAN version 3
The CAZyme classifier `dbCAN` verions >= 3.0.6 can be installed via Bioconda (recommended). The full installation instructions can be found [here](https://github.com/linnabrown/run_dbcan/tree/c9bad4835b2bc1a9685d693237f1d6a9d56ff3a1) -- and must be followed to ensure all additional database files are downloaded and compiled correctly.

### Install dbCAN version 2
The installation instructions for `dbCAN` v==2.0.11 can be found [here](https://github.com/linnabrown/run_dbcan/tree/fde6d7225441ef3d4cb29ea29e39cfdcc41d8b19) and were followed to install dbCAN for the analysis presented in the manuscript.

**Note:** If using `cazomevolve` to run `dbCAN` after installing `dbCAN` version 2, ensure `saintBioutils` and `cazomevolve` are installed in the same virtual environment as `dbCAN` version 2.

### Install signalP
Installation instructions for `signalP6` can be found [here](https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md).

### All scripts
1. Download datasets
    1. `download_pecto_genomes.sh` - Downloads _Pectobacteriaceae_ genomes
    2. `download_same_pecto_genomes.sh` - Download the genomes used in the manuscript
    3. `ident_missing_protomes.py` - Identify genomes were a .faa file was not available
    4. `annotate_genomes.sh` - Predicte proteome using Prodigal
    5. `download_pd_genomes.sh` - Download _Pectobacterium_ and _Dickeya_ reference genomes
    6. `build_cazyme_db.sh` - Build a local CAZyme db
2. Annotate CAZomes
    1. `get_cazy_cazymes_pecto.sh` - Retrieve CAZy family annotations from the local CAZyme db for _Pectobacteriaceae_
    2. `get_cazy_cazymes_pd.sh` - Retrieve CAZy family annotations from the local CAZyme db for _Pecto&Dickeya_ data set
    3. `run_dbcan_dbcan_pectobact.sh` - Run dbCAN version 2 on _Pectobacteriaceae_ proteomes
    4. `run_dbcan_pecto_dic.sh` - Run dbCAN version 3 on _Pectobacterium_ and _Dickeya_ proteomes
    5. `get_dbcan_cazymes_pectobact.sh` - Parse dbCAN output for _Pectobacteriaceae_
    6. `get_dbcan_cazymes_pecto_dic.sh` - Parse dbCAN output for _Pectobacterium_ and _Dickeya_
3. Run ANI analysis and build dendrogram
    1. `run_anim.sh`
    2. `build_anim_tree.sh`
    3. `build_anim_tree.R`
    4. `parse_anim_tab.py`
4. Reconstruct phylogenetic tree
    1. `annotate_genomes.sh`
    2. `find_orthologues.sh`
    3. `align_scos.sh`
    4. `extract_cds.py`
    5. `backtranslates.sh`
    6. `concatenate_cds.py`
    7. `raxml_ng_build_tree.sh`
5. Annotate intracellular and extracellular CAZymes
    1. `gather_cazyme_seqs.py`
    2. `run_signalp.sh`
6. Add taxonomic classifications
    1. `add_taxs.sh`
    2. `add_ani_tax.py`
    3. `add_tax_phylotree.py`
7. Explore CAZome composition
    1. `explore_pectobact_cazome.ipynb`
    2. `explore_pecto_dic_cazome.ipynb`
8. Compare trees
    1. `build_tanglegrams.R`
9. Identify networkds of co-evolving CAZy families
    1. `find_colevolving_pectobact.sh`
    2. `find_colevolving_pd.sh`

# Results and online supplementary
<a id="results"></a>

Owing to the size of the data sets used, the figures are consequently compressed in the final manuscript. This remote repository contains the original full size, high resolution figures.

The original figures are found in the `results` directory, and contained within the `jupyter notebooks` used to run the analyses, which can be found here (the raw notebooks are for downloading and re-running locally, the website version are for viewing the results):
* _Pectobacteriaceae_: 
    * [raw notebook](https://github.com/HobnobMancer/SI_Hobbs_et_al_2023_Pecto/blob/master/notebooks/explore_pectobact_cazomes.ipynb)
    * [website](https://hobnobmancer.github.io/SI_Hobbs_et_al_2023_Pecto/notebooks/explore_pectobact_cazomes.html)
* _Pectobacterium_ and _Dickeya_ reference data set:
    * [raw notebook](https://github.com/HobnobMancer/SI_Hobbs_et_al_2023_Pecto/blob/master/notebooks/explore_pecto_dic_cazomes.ipynb)
    * [website](https://hobnobmancer.github.io/SI_Hobbs_et_al_2023_Pecto/notebooks/explore_pecto_dic_cazomes.html)
* Intracellular and extracellular CAZymes:
    * [raw notebook](https://github.com/HobnobMancer/SI_Hobbs_et_al_2023_Pecto/blob/master/notebooks/explore_IE_cazomes.ipynb)
    * [website](https://hobnobmancer.github.io/SI_Hobbs_et_al_2023_Pecto/notebooks/explore_IE_cazomes.html)


# Method: Reproducing the analyses

## Build a local CAZyme database

Configure using [`cazy_webscraper`](https://hobnobmancer.github.io/cazy_webscraper/) (Hobbs _et al., 2022) to download all data from the CAZy database, and compile the data into a local CAZyme database.

> cazy_webscraper: local compilation and interrogation of comprehensive CAZyme datasets
Emma E. M. Hobbs, Tracey M. Gloster, Leighton Pritchard
bioRxiv 2022.12.02.518825; doi: https://doi.org/10.1101/2022.12.02.518825

```bash
# create a local CAZyme database
scripts/build_cazyme_db.sh <email>
```

This generated the local CAZyme database `data/cazy/cazy_db`.

## _Pectobacteriaceae_

### Download data

`cazomevolve` was used to download all complete _Pectobacteriaceae_ genomic assemblies (in genome sequence and protein sequence FASTA file format) from NCBI Assembly, by querying the NCBI Taxonomy database and retrieving all genomic assemblies linked to the _Pectobacteriaceae_ (NCBI:txid1903410). To repeat this method, run the following command from the root of this directory:

```bash
# download Pectobacteriaceae genomes from GenBank
scripts/download/download_pecto_genomes.sh <email>
```

**Note:** With the continual addition of new genomic assemblies to the NCBI Assembly database, repeating the download of _Pectobacteriaceae_ genomes may generate a different dataset to that presented in Hobbs _et al._. To repeat the analysis presented in the manuscript, run the following command from the root of the directory to configure `ncbi-genome-download` to download the 660 genomic assemblies of the genomes used in the manuscript:

```bash
scripts/download/download_same_pecto_genomes.sh
```

In both cases, the downloaded genomic sequence files were written to the dir `data/pectobact/genomes`, the downloaded protein FASTA files were written to `data/pectobact/proteomes`.

### Predict proteomes

Not all genomic assemblies in NCBI are annotated, i.e. a proteome FASTA file (`.faa` file) is not available for all genomic sequences in NCBI.

To identify those genomes were a proteome FASTA file was not available, and thus was not downloaded, the Python script `ident_missing_protomes.py` was run.

```bash
scripts/download/ident_missing_proteomes.py
```

The script generated a text file listing the genomic accession of each assembly for which a proteome FASTA file (`.faa`) file was not downloaded. The file was written to `data/pectobact/missing_genomes`. 

If using the 660 assemblies presented in the manuscript, proteome FASTA files were not available for 104 assemblies.

The script `annotate_genomes.sh` coordinates running `prodigal` on all genome sequences were a proteome FASTA file could not be retrieved, and copies the predicted proteome FASTA file to the `data/pectobact/proteome` directory.

```bash
scripts/download/annotates_genomes.sh
```

### Annotate CAZomes

#### Get CAZy annotated CAZymes

Configure using `cazomevolve` to identify CAZymes classified in the local CAZyme database, for both the _Pectobacteriaceae_.

```bash
scripts/annotate_cazome/get_cazy_cazymes_pecto.sh
```

Two tab delimited lists were created:
1. Listing the CAZy family accession and genomic accession per line: `data/pectobact/cazomes/pecto_fam_genomes`
2. Listing the CAZy family, genomic accession and protein accession per line: `data/pectobact/cazomes/pecto_fam_genomes_proteins`

Proteins in the download protein FASTA files that were not listed in the local CAZyme database were written to `data/pectobact/cazomes/dbcan_input` for _Pectobacteriaceae_.

#### Get dbCAN predicted CAZymes

To retrieve the most comprehensive CAZome for each genome, protein sequences not found in the local CAZyme database were parsed by the CAZyme classifier [`dbCAN`](https://github.com/linnabrown/run_dbcan) (Zhang _et al._ 2018), configured using `cazomevolve`.

> Han Zhang, Tanner Yohe, Le Huang, Sarah Entwistle, Peizhi Wu, Zhenglu Yang, Peter K Busk, Ying Xu, Yanbin Yin; dbCAN2: a meta server for automated carbohydrate-active enzyme annotation, Nucleic Acids Research, Volume 46, Issue W1, 2 July 2018, Pages W95–W101, https://doi.org/10.1093/nar/gky418

Run the following command from the root of this directory. *Note: depending on the computer resources this may take multiple days to complete*

**Note:** The following commands MUST be run from the same directory containing the `db` directory created when installing `dbCAN` - the following commands presumt the `db` dir is located in the root of this repository.

```bash
# run dbcan for the pectobacteriaceae data set - using dbCAN version >= 2
scripts/run_dbcan_pectobact.sh
```

After running dbCAN, use the following commands to parse the output from dbCAN and add the predicted CAZy family annotations, protein accessions and genomic accessions to the tab delimited lists created above.

The command runs the `cazomevolve` command `cazevolve_get_dbcan` which can be used to parse the output from `dbCAN` version 2 and version 3.

```bash
# parse pectobacteriaceae dbcan output
scripts/get_dbcan_cazymes_pectobact.sh
```

At the end, two plain text files will be generated, containing tab separated data:

The _Pectobacteriaceae_ lists were written to:
1. `data/pectobact/cazomes/pecto_fam_genomes`
2. `data/pectobact/cazomes/pecto_fam_genomes_proteins`

### Run ANI analysis and construct dendrogram

The software package `pyani` [Pritchard et al](https://doi.org/10.1039/C5AY02550H) was used to perform an average nucleotide identify (ANI) comparison between all pairs of _Pectobacteriaceae_ genomes, using the ANIm method.

> Pritchard et al. (2016) "Genomics and taxonomy in diagnostics for food security: soft-rotting enterobacterial plant pathogens" Anal. Methods 8, 12-24

```bash
scripts/tree/ani/run_anim.sh
```

This created a pyani database in `data/pectobact/tree`. Graphical outputs summarising the pyani analysis were written to `results/pectobact/tree/anim`.

**May need to run to double check what output is created, and how to get the tsv file for generating the dendrogram**.

A dendrogram was reconstructed from the ANIm analysis using the bash script `build_anim_tree.sh`, which coordinated the extraction of the calclated ANI values from the local `pyani` database, replacing the `pyani` genome IDs with the NCBI genomic version accessions using the Python script `parse_anim_tab.py`, and coordinated the R script `build_anim_tree.R`, which build a distance matrix and used hierarchical clustering (using the 'single' method) to build a dendorgram that was written in Newick format to `data/pectobact/tree/pyani_ani_tree.new`:
```bash
scripts/tree/ani/build_anim_tree.sh
```


## _Pectobacterium_ and _Dickeya_

### Download data

Configure `ncbi-genome-download` to download genomic assemblies (genome sequence and protein sequence FASTA files) for the assembly accessions listed in `supplementary_files/supplementary_file_1`, which include the genomes for _Pectobacterium_, _Dickeya_, and _Hafnia alvei_ ATCC 13337, and decompresses the downloaded files.

```bash
# download Pectobacterium and dickeya genomes from RefSeq
scripts/download_pd_genomes.sh
```

The downloaded genomic sequence files were written to the dir `data/pecto_dic/genomes`, the downloaded protein FASTA files were written to `data/pecto_dic/proteomes`.

### Annotate CAZomes

#### Get CAZy classified CAZymes

Configure using `cazomevolve` to identify CAZymes classified in the local CAZyme database.

```bash
scripts/annotate_cazome/get_cazy_cazymes_pd.sh
```

Two tab delimited lists were created:
1. Listing the CAZy family accession and genomic accession per line: `data/pecto_dic/cazomes/pd_fam_genomes`
2. Listing the CAZy family, genomic accession and protein accession per line: `data/pecto_dic/cazomes/pd_fam_genomes_proteins`

Proteins in the download protein FASTA files that were not listed in the local CAZyme database were written to `data/pecto_dic/cazomes/dbcan_input`.

#### Get dbCAN predicted CAZymes

To retrieve the most comprehensive CAZome for each genome, protein sequences not found in the local CAZyme database were parsed by the CAZyme classifier [`dbCAN`](https://github.com/linnabrown/run_dbcan) (Zhang _et al._ 2018), configured using `cazomevolve`.

> Han Zhang, Tanner Yohe, Le Huang, Sarah Entwistle, Peizhi Wu, Zhenglu Yang, Peter K Busk, Ying Xu, Yanbin Yin; dbCAN2: a meta server for automated carbohydrate-active enzyme annotation, Nucleic Acids Research, Volume 46, Issue W1, 2 July 2018, Pages W95–W101, https://doi.org/10.1093/nar/gky418

Run the following command from the root of this directory. *Note: depending on the computer resources this may take multiple days to complete*

**Note:** The following commands MUST be run from the same directory containing the `db` directory created when installing `dbCAN` - the following commands presumt the `db` dir is located in the root of this repository.

```bash
# run dbcan for the pectobacterium&dickeya data set - using dbCAN version >= 3
scripts/run_dbcan_pecto_dic.sh
```

After running dbCAN, use the following commands to parse the output from dbCAN and add the predicted CAZy family annotations, protein accessions and genomic accessions to the tab delimited lists created above.

The commands run the `cazomevolve` command `cazevolve_get_dbcan` which can be used to parse the output from `dbCAN` version 2 and version 3.

```bash
# parse pectobacterium and dickeya dbcan output
scripts/get_dbcan_cazymes_pecto_dic.sh
```

At the end, two plain text files will be generated, containing tab separated data:
1. `data/pecto_dic/cazomes/pd_fam_genomes`
2. `data/pecto_dic/cazomes/pd_fam_genomes_proteins`

### Reconstruct phylogenetic tree

To reconstruct the phylogenetic tree of _Pectobacterium_, _Dickeya_ and _Musicola_, the method presented in 
[Hugouvieux-Cotte-Pattat et al.](https://pure.strath.ac.uk/ws/portalfiles/portal/124038859/Hugouvieux_Cotte_Pattat_etal_IJSEM_2021_Proposal_for_the_creation_of_a_new_genus_Musicola_gen_nov_reclassification_of_Dickeya_paradisiaca.pdf) was used.

The specific method is found in the [Hugouvieux-Cotte-Pattat et al. supplementary](https://widdowquinn.github.io/SI_Hugouvieux-Cotte-Pattat_2021/).

#### CDS prediction

To ensure consistency of nomenclature and support back threading the nucleotides sequences onto 
aligned single-copy orthologues, all downloaded RefSeq genomes were reannotated using 
[`prodigal`](https://github.com/hyattpd/Prodigal)

> Hyatt D, Chen GL, Locascio PF, Land ML, Larimer FW, Hauser LJ. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics. 2010 Mar 8;11:119. doi: 10.1186/1471-2105-11-119. PMID: 20211023; PMCID: PMC2848648.

```bash
scripts/tree/phylo/annotate_genomes.sh
```

The annotate features were written to the following directories:  
Proteins: `data/pecto_dic/tree/genomes/proteins`  
CDS: `data/pecto_dic/tree/genomes/cds`  
GBK: `data/pecto_dic/tree/genomes/gbk`  


#### Identify Single Copy Orthologues (SCOs)

Orthologues present in the genomes were identified using [`orthofinder`](https://github.com/davidemms/OrthoFinder).

> Emms, D.M. and Kelly, S. (2019) OrthoFinder: phylogenetic orthology inference for comparative genomics. [Genome Biology 20:238](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y)

```bash
scripts/tree/phylo/find_orthologues.sh
```

The Orthofind output was written to `data/pecto_dic/tree/orthologues`. The SCO sequences are written to the dir `data/pecto_dic/tree/orthologues/Results_<date>/Single_Copy_Orthologue_Sequences`.

#### Multiple Sequence Alignment

Each collection of single-copy orthologous was aligned using [`MAFFT`](https://mafft.cbrc.jp/alignment/software/).

The output from `MAFFT` (the aligned files) are placed in the `data/pecto_dic/tree/sco_proteins_aligned` directory.

```bash
scripts/tree/phylo/align_scos.sh <path to dir containing SCO seqs from orthofinder>
```

#### Collect Single-Copy Orthologues CDS sequences

The CDS sequences corresponding to each set of single-copy orthologues are identified and extracted with the Python script `extract_cds.py`. To reproduce this analysis, ensure the `PROTDIR` constant in the script is directed to the correct output directory for orthofinder. The script can then be run from the current directory with:

```bash
python3 scripts/tree/phylo/extract_cds.py
```

The output is a set of unaligned CDS sequences corresponding to each single-copy orthologue, which are 
placed in the `data/pecto_dic/tree/sco_cds` directory

#### Back-translate Aligned Single-Copy Orthologues

The single-copy orthologue CDS sequences were threaded onto the corresponding aligned protein sequences using [`t-coffee`](http://www.tcoffee.org/Projects/tcoffee/).

> T-Coffee: A novel method for multiple sequence alignments. Notredame, Higgins, Heringa, JMB, 302(205-217)2000

The results can be reproduced by executing the `backtranslate.sh` script from this directory.

```bash
scripts/tree/phylo/backtranslate.sh
```

The backtranslated CDS sequences are placed in the `data/pecto_dic/tree/sco_cds_aligned` directory.

#### Concatenating CDS into a Multigene Alignment

The threaded single-copy orthologue CDS sequences were concatenated into a single sequence per input organism using the Python script `concatenate_cds.py`. To reproduce this, execute the script from this directory with:

```bash
python3 scripts/tree/phylo/concatenate_cds.py
```

Two files are generated, a FASTA file with the concatenated multigene sequences, and a partition file allowing a different set of model parameters to be fit to each gene in phylogenetic reconstruction, written to `data/pecto_dic/tree/concatenated.fasta` and `data/pecto_dic/tree/concatenated.part`.

#### Phylogenetic reconstruction

To reconstruct the phylogenetic tree, the bash script `raxml_ng_build_tree.sh` is used, and is 
run from the root of this repository. This executes a series of [`raxml-ng`](https://github.com/amkozlov/raxml-ng) commands.

All genes were considered as separate partitions in the reconstuction, 
with parameters estimated  for the `GTR+FO+G4m+B` model (as recommended by `raxml-ng check`).

Tree reconstructions are placed in the `tree` directory. The best estimate tree is `03_infer.raxml.bestTree` and the midpoint-rooted, manually-annotated/coloured tree (using [`figtree`](http://tree.bio.ed.ac.uk/software/figtree/)) is `03_infer.raxml.bestTree.annotated`)

> Alexey M. Kozlov, Diego Darriba, Tomáš Flouri, Benoit Morel, and Alexandros Stamatakis (2019) RAxML-NG: A fast, scalable, and user-friendly tool for maximum likelihood phylogenetic inference. Bioinformatics, btz305 [doi:10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305)

```bash
scripts/tree/phylo/raxml_ng_build_tree.sh
```

## Annotate intracellular and extracellular CAZymes

To predict which CAZymes were intracellular and which were extracellular, the protein sequences of all CAZymes were gathered into a single FASTA file using the python script `gather_cazyme_seqs.py.

```bash
python3 scripts/signalp/gather_cazyme_seqs.py
```

`signalP` (version 6) ([Teufel et al.](https://www.nature.com/articles/s41587-021-01156-3)) was used to predict the presence of signal peptides in the protein sequences of the identified CAZymes.

> Teufel, F., Almagro Armenteros, J.J., Johansen, A.R. et al. SignalP 6.0 predicts all five types of signal peptides using protein language models. Nat Biotechnol 40, 1023–1025 (2022). https://doi.org/10.1038/s41587-021-01156-3

```bash
scripts/signalp/run_signalp.sh
```

The output from `signalP` was written to `data/pecto_dic/signalp/pd_signalp_output`.

The Python script `get_ie_cazymes.py` was used to add the intracellular/extracellular annotations to the tab delimited file `data/pecto_dic/cazomes/pd_fam_genomes_proteins`.

```bash
scripts/signalp/get_ie_cazymes.py
```

`data/pecto_dic/cazomes/pd_fam_genomes_proteins` includes the headers 'Fam', 'Genome', and 'Protein'.

## Add taxonomic classifications

Download the GTDB database dump from the [GTDB repository](https://data.gtdb.ecogenomic.org/releases/). Release 202.0 was used in the manuscript Hobbs et al. Save the database dump (TSV file) to `data/gtdb/` directory.

The bash script `add_tax.sh` was used to coordinate running `cazomevolve` to add taxonomic information to each genomic accession, in every tab delimited list of (i) CAZy family and genomic accession, and (ii) CAZy family, genomic accession and protein accession that was generated.

```bash
scripts/taxs/add_tax.sh <use email address> <path to gtdb tsv file>
```

Use Python scripts `add_ani_tax.py` and `add_tax_phylotree.py` to add the taxonomic information to the reconstructed ANI and phylogenetic trees, respectively.

```bash
scripts/taxs/add_ani_tax.py
scripts/tax/add_tax_phylotree.py
```

## Explore CAZome composition

Exploration of the CAZomes for each data set (_Pectobacteriacea_, _Pectobacterium_ & _Dickeya_, and intra- and extra-cellular CAZymes) were preformed within `jupyter notebooks`, which are available in this repository  (the raw notebooks are for downloading and re-running locally, the website version are for viewing the results):
* _Pectobacteriaceae_: 
    * [raw notebook](https://github.com/HobnobMancer/SI_Hobbs_et_al_2023_Pecto/blob/master/notebooks/explore_pectobact_cazomes.ipynb)
    * [website](https://hobnobmancer.github.io/SI_Hobbs_et_al_2023_Pecto/notebooks/explore_pectobact_cazomes.html)
* _Pectobacterium_ and _Dickeya_ reference data set:
    * [raw notebook](https://github.com/HobnobMancer/SI_Hobbs_et_al_2023_Pecto/blob/master/notebooks/explore_pecto_dic_cazomes.ipynb)
    * [website](https://hobnobmancer.github.io/SI_Hobbs_et_al_2023_Pecto/notebooks/explore_pecto_dic_cazomes.html)
* Intracellular and extracellular CAZymes:
    * [raw notebook](https://github.com/HobnobMancer/SI_Hobbs_et_al_2023_Pecto/blob/master/notebooks/explore_IE_cazomes.ipynb)
    * [website](https://hobnobmancer.github.io/SI_Hobbs_et_al_2023_Pecto/notebooks/explore_IE_cazomes.html) 


Specifically, the analyses performed in the notebooks were executed using the module `cazomevolve.cazome.explore`, which contains functions for exploring the CAZome annotated by `cazomevolve`. These are:

```python
from cazomevolve.cazome.explore import (
    cazome_sizes,
    identify_families,
    parse_data,
    pca,
    plot,
    taxonomies,
)

# parse the output from cazomevolve tab delimited lists
parse_data.get_dbcan_fams_data()
parse_data.build_fam_freq_df()
parse_data.index_df()  # index genome, genus and species to be row names

# add taxonomic information for taxonomic context
taxonomies.add_tax_info()
taxonomies.get_gtdb_db_tax_dict()  # in development
taxonomies.get_gtdb_search_tax_dict()
taxonomies.get_ncbi_tax_dict()  # in development
taxonomies.get_group_sample_sizes()  # returns the number of genomes per group (genus or species)

# summarise the size of the cazomes
cazome_sizes.get_cazome_size_df()
cazome_sizes.get_proteome_size()
cazome_sizes.get_cazome_proportion_df()
cazome_sizes.get_num_of_fams_per_group()

# identify the core CAZome, i.e. families that appear in every genome
identify_families.identify_core_cazome()
identify_families.get_core_fam_freqs()

# identify families that are specific to a group (i.e. genus or species)
identify_families.get_group_specific_fams()

# identity families that always appear together 
identify_families.get_cooccurring_fams()  # across all genomes
identify_families.get_grps_cooccurring_fams()  # in a specific group (i.e. genus or species)

# visually summarise the data
plot.get_clustermap()  # clustermap of cazy family freqs - potentially add clustering by cazy class freqs

# perform and visualise PCA
pca.perform_pca()
pca.plot_explained_variance()
pca.plot_spree()
pca.plot_pca()  # project genomes onto PCs
pca.plot_loadings()
```

## Compare trees

??? to include here or in the thesis supplementary ???


## Identify associating CAZy families using `coinfinder`

Use the tool [`coinfinder`](https://github.com/fwhelan/coinfinder) (Whelan et al.) to identify CAZy families that are present in the genome together more often than expected by chance and lineage.

> Fiona J. Whelan, Martin Rusilowicz, & James O. McInerney. "Coinfinder: detecting significant associations and dissociations in pangenomes." doi: [https://doi.org/10.1099/mgen.0.000338](https://doi.org/10.1099/mgen.0.000338)

**Generate circular trees and heatmaps:**

To reproduce the output from `coinfinder` in the same structure as presented in the manuscript (i.e. a circular tree surrounded by a heatmap), overwrite the file `network.R` in `coinfinder` with the respective R script in `scripts/coevolution`, and use the corresponding bash script:

* _Pectobacteriaceae_
    * `network.R`: `scripts/coevolution/pectobact_circular_network.R`
    * bash: `scripts/coevolution/find_coevolving_pectobact.sh`
* _Pectobacterium_ and _Dickeya_
    * `network.R`: `scripts/coevolution/pd_circular_network.R`
    * bash: `scripts/coevolution/find_coevolving_pd.sh`

**Generate linear trees and heatmaps, with taxonomy information:**

These circular heatmaps annotate each leaf of the tree with only the respective genomic version accession. To list the taxonomic infomration as well, on each leaf of the tree, overwrite the contents in the file `network.R` in `coinfinder` with the respective R script in `scripts/coevolution`, and use the respective bash script to configure `coinfinder`:

* _Pectobacteriaceae_
    * `network.R`: `scripts/coevolution/pectobact_taxs_rectangular_network.R`
    * bash: `scripts/coevolution/find_coevolving_pectobact_with_tax.sh`
* _Pectobacterium_ and _Dickeya_
    * `network.R`: `scripts/coevolution/pd_taxs_rectangular_network.R`
    * bash: `scripts/coevolution/find_coevolving_pd_with_tax.sh`
