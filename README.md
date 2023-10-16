# Supplementary Information: Hobbs _et al._, 2023:
## Association between niche adaptation and evolution of carbohydrate active enzymes in _Pectobacteriaceae_

[![DOI](https://zenodo.org/badge/609818373.svg)](https://zenodo.org/badge/latestdoi/609818373)

This repository contains supplementary information for analyses reported in Hobbs et al. (2023), exploring the diversity in the 
Carbohydrate Active enZyme (CAZyme) complement and association with the plant host range of _Pectobacteriaceae_.

## You can find the full report, exploring the CAZomes [here](https://hobnobmancer.github.io/SI_Hobbs_et_al_2023_Pecto/notebooks/explore_pectobact_cazomes.html).

A citation for this work will be added once available. At the present please cite this repository as the source, the DOI:10.5281/zenodo.7699655, and the authors (in order):
Emma E. M. Hobbs<sup>1,2,3</sup>, Tracey, M. Gloster<sup>1</sup>, Leighton Prichard<sup>2</sup>.

> 1. School of Biology and Biomedical Sciences Research Complex, University of St Andrews, North Haugh, St Andrews, Fife, KY16 9ST, UK
> 2. Strathclyde Institute of Pharmacy and Biomedical Sciences, University of Strathclyde, Glasgow G4 ORE, UK
> 3. Cell and Molecular Sciences, James Hutton Institute, Invergowrie, Dundee DD2 5DA, UK

```latex
@misc{Hobbs2023,
author = {Emma E. M. Hobbs and Tracey M. Gloster and Leighton Pritchard},
title = {Association between niche adaptation and evolution of carbohydrate active enzymes in Pectobacteriaceae},
howpublished = {\url{https://hobnobmancer.github.io/SI_Hobbs_et_al_2023_Pecto/}},
year = {2023},
note = {Version 1. DOI:10.5281/zenodo.7699655}
}
```

To repeat analyses, run all commands provided in the walkthrough from the root of this directory.

All raw figure files presented in the complete report in the manuscript can be found in the [`results/` directory](https://github.com/HobnobMancer/SI_Hobbs_et_al_2023_Pecto/tree/master/results).

## Online supplementary

Owing to the size of the data sets used, the figures are consequently compressed in the final manuscript. This remote repository contains the original full size, high resolution figures.

Additionally, some analyses are only briefly mentioned in the manuscript. The full method and results of these analyses are stored in this repository.

For the complete analysis of the CAZyme complements (i.e. the CAZomes) are available in the [`jupyter notebook`](https://hobnobmancer.github.io/SI_Hobbs_et_al_2023_Pecto/notebooks/explore_pectobact_cazomes.html).

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
├── LICENSE
├── README.md
├── _config.yml
├── data
│   ├── README.md
│   ├── cazomes
│   │   ├── coinfinder_pecto_fam_genomes
│   │   ├── coinfinder_pecto_fam_genomes_taxs
│   │   ├── pecto_fam_genomes
│   │   ├── pecto_fam_genomes_proteins
│   │   ├── pecto_fam_genomes_proteins_taxs
│   │   └── pecto_fam_genomes_taxs
│   ├── genomes
│   │   ├── classes.txt
│   │   └── labels.txt
│   ├── genomic_accessions
│   │   ├── genomes_for_coinfinder.txt
│   │   └── pectobact_accessions
│   ├── missing_genomes
│   └── tree
│       └── ani_tree
│           ├── anim_matrix.tab
│           ├── genomes.tab
│           ├── logs
│           ├── matrix_aln_lengths_4.tab
│           ├── matrix_aln_lengths_run4.pdf
│           ├── matrix_coverage_4.tab
│           ├── matrix_coverage_run4.pdf
│           ├── matrix_hadamard_4.tab
│           ├── matrix_hadamard_run4.pdf
│           ├── matrix_identity_4.tab
│           ├── matrix_identity_run4.pdf
│           ├── matrix_sim_errors_4.tab
│           ├── matrix_sim_errors_run4.pdf
│           ├── pyani_ani_tree.new
│           ├── pyani_ani_tree_taxs.new
│           └── scatter_identity_vs_coverage_run4.pdf
├── notebooks
│   ├── explore_pecto_dic_cazomes.html
│   ├── explore_pectobact_cazomes.html
│   └── explore_pectobact_cazomes.ipynb
├── requirements.txt
├── results
│   ├── cazome_size
│   │   └── cazome_sizes.csv
│   ├── cazy_classes
│   │   └── cazy_class_sizes.csv
│   ├── cazy_families
│   │   ├── cazy_fam_freqs.csv
│   │   ├── fam_freq_clustermap.svg
│   │   ├── paper_fam_freq_clustermap.png
│   │   ├── paper_fam_freq_clustermap_FILTERED.svg
│   │   ├── paper_genus_species_fam_freq_clustermap.png
│   │   ├── paper_pheno_genus_fam_freq_clustermap.png
│   │   └── unique_grp_fams.tsv
│   ├── cooccurring_families
│   │   ├── cooccurring_fams_freqs.csv
│   │   ├── fam_corr_M_filled.csv
│   │   ├── paper-cooccurring_fams_freqs.csv
│   │   ├── paper-pecto-cooccurring-families.svg
│   │   └── pecto-cooccurring-families.svg
│   ├── core_cazome
│   │   ├── core_cazome_freqs.csv
│   │   ├── genera_core_cazome.svg
│   │   └── genera_soft_hard_core_cazome.svg
│   └── pca
│       ├── PC1-vs-PC2
│       │   ├── pca_pc1_vs_pc2-genus.png
│       │   ├── pca_pc1_vs_pc2-loadings_plot.png
│       │   └── pca_pc1_vs_pc2-species.png
│       ├── PC1-vs-PC3
│       │   ├── pca_pc1_vs_pc3-genus.png
│       │   ├── pca_pc1_vs_pc3-loadings_plot.png
│       │   └── pca_pc1_vs_pc3-species.png
│       ├── PC1-vs-PC4
│       │   ├── pca_pc1_vs_pc4-genus.png
│       │   ├── pca_pc1_vs_pc4-loadings_plot.png
│       │   └── pca_pc1_vs_pc4-species.png
│       ├── PC2-vs-PC3
│       │   ├── pca_pc2_vs_pc3-genus.png
│       │   ├── pca_pc2_vs_pc3-loadings_plot.png
│       │   └── pca_pc2_vs_pc3-species.png
│       ├── PC2-vs-PC4
│       │   ├── pca_pc2_vs_pc4-genus.png
│       │   ├── pca_pc2_vs_pc4-loadings_plot.png
│       │   └── pca_pc2_vs_pc4-species.png
│       ├── PC3-vs-PC4
│       │   ├── pca_pc3_vs_pc4-genus.png
│       │   ├── pca_pc3_vs_pc4-loadings_plot.png
│       │   └── pca_pc3_vs_pc4-species.png
│       ├── pca_explained_variance.png
│       ├── pca_pc_screen_genus.svg
│       ├── pca_pc_screen_species.svg
│       └── pectobact_pca_scree.png
├── scripts
│   ├── README.md
│   ├── annotate_cazome
│   │   ├── get_cazy_cazymes.sh
│   │   ├── get_dbcan_cazymes.sh
│   │   └── run_dbcan.sh
│   ├── coevolution
│   │   ├── find_coevolving_pectobact.sh
│   │   ├── find_coevolving_pectobact_with_tax.sh
│   │   ├── pectobact_circular_network.R
│   │   └── pectobact_taxs_rectangular_network.R
│   ├── download
│   │   ├── annotate_genomes.sh
│   │   ├── build_cazyme_database.sh
│   │   ├── download_genomes.sh
│   │   ├── download_ms_genomes.sh
│   │   └── ident_missing_proteomes.py
│   ├── taxs
│   │   ├── add_ani_tax.py
│   │   └── add_taxs.sh
│   └── tree
│       ├── README.md
│       └── ani
│           ├── build_anim_tree.R
│           ├── build_anim_tree.sh
│           ├── parse_anim_tab.py
│           └── run_anim.sh
└── structure

28 directories, 94 files
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

### Install dbCAN version 2
The installation instructions for `dbCAN` v==2.0.11 can be found [here](https://github.com/linnabrown/run_dbcan/tree/fde6d7225441ef3d4cb29ea29e39cfdcc41d8b19) and were followed to install dbCAN for the analysis presented in the manuscript.

### All scripts
1. Download datasets
    1. `download_genomes.sh` - Search and download all _Pectobacteriaceae_ genomes in NCBI
    2. `download_ms_genomes.sh` - Download the genomes used in the manuscript
    3. `ident_missing_protomes.py` - Identify genomes were a .faa file was not available
    4. `annotate_genomes.sh` - Predicte proteome using Prodigal
    5. `build_cazyme_db.sh` - Build a local CAZyme db
2. Annotate CAZomes
    1. `get_cazy_cazymes.sh` - Retrieve CAZy family annotations from the local CAZyme db for _Pectobacteriaceae_
    3. `run_dbcan_dbcan.sh` - Run dbCAN version 2 on _Pectobacteriaceae_ proteomes
    5. `get_dbcan_cazymes.sh` - Parse dbCAN output
3. Run ANI analysis and build dendrogram
    1. `run_anim.sh`
    2. `build_anim_tree.sh`
    3. `build_anim_tree.R`
    4. `parse_anim_tab.py
4. Add taxonomic classifications
    1. `add_taxs.sh`
    2. `add_ani_tax.py`
7. Explore CAZome composition
    1. `explore_pectobact_cazome.ipynb`
8. Compare trees
    1. `build_tanglegrams.R`
9. Identify networkds of co-evolving CAZy families
    1. `find_colevolving_pectobact.sh`
    2. `find_colevolving_pectobact_with_tax.sh`


# Results and online supplementary
<a id="results"></a>

Owing to the size of the data sets used, the figures are consequently compressed in the final manuscript. This remote repository contains the original full size, high resolution figures.

The original figures are found in the `results` directory, and contained within the `jupyter notebooks` used to run the analyses, which can be found here (the raw notebooks are for downloading and re-running locally, the website version are for viewing the results):
* [raw notebook](https://github.com/HobnobMancer/SI_Hobbs_et_al_2023_Pecto/blob/master/notebooks/explore_pectobact_cazomes.ipynb)
* [website](https://hobnobmancer.github.io/SI_Hobbs_et_al_2023_Pecto/notebooks/explore_pectobact_cazomes.html)


# Method: Reproducing the analyses

Several of the data files required to repeat the analyses presented in the manuscript are stored (available for use) in the repo. These files are stored in the `data/` directory:

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

## Download genomes

`cazomevolve` was used to download all complete _Pectobacteriaceae_ genomic assemblies (in genome sequence and protein sequence FASTA file format) from NCBI Assembly, by querying the NCBI Taxonomy database and retrieving all genomic assemblies linked to the _Pectobacteriaceae_ (NCBI:txid1903410). To repeat this method, run the following command from the root of this directory:

```bash
# download Pectobacteriaceae genomes from GenBank
scripts/download/download_genomes.sh <email>
```

**Note:** With the continual addition of new genomic assemblies to the NCBI Assembly database, repeating the download of _Pectobacteriaceae_ genomes may generate a different dataset to that presented in Hobbs _et al._. To repeat the analysis presented in the manuscript, run the following command from the root of the directory to configure `ncbi-genome-download` to download the 660 genomic assemblies of the genomes used in the manuscript:

```bash
scripts/download/download_ms_genomes.sh
```

In both cases, the downloaded genomic sequence files were written to the dir `data/genomes`, the downloaded protein FASTA files were written to `data/proteomes`.

## Predict proteomes

Not all genomic assemblies in NCBI are annotated, i.e. a proteome FASTA file (`.faa` file) is not available for all genomic sequences in NCBI.

To identify those genomes were a proteome FASTA file was not available, and thus was not downloaded, the Python script `ident_missing_protomes.py` was run.

```bash
scripts/download/ident_missing_proteomes.py
```

The script generated a text file listing the genomic accession of each assembly for which a proteome FASTA file (`.faa`) file was not downloaded. The file was written to `data/missing_genomes`. 

If using the 717 assemblies presented in the manuscript, proteome FASTA files were not available for 107 assemblies.

The script `annotate_genomes.sh` coordinates running `prodigal` on all genome sequences were a proteome FASTA file could not be retrieved, and copies the predicted proteome FASTA file to the `data/proteome` directory.

```bash
scripts/download/annotates_genomes.sh
```

## Annotate CAZomes

### Get CAZy annotated CAZymes

Configure using `cazomevolve` to identify CAZymes classified in the local CAZyme database, for both the _Pectobacteriaceae_.

```bash
scripts/annotate_cazome/get_cazy_cazymes.sh
```

Two tab delimited lists were created:
1. Listing the CAZy family accession and genomic accession per line: `data/cazomes/pecto_fam_genomes`
2. Listing the CAZy family, genomic accession and protein accession per line: `data/cazomes/pecto_fam_genomes_proteins`

Proteins in the download protein FASTA files that were not listed in the local CAZyme database were written to `data/cazomes/dbcan_input` for _Pectobacteriaceae_.

### Get dbCAN predicted CAZymes

To retrieve the most comprehensive CAZome for each genome, protein sequences not found in the local CAZyme database were parsed by the CAZyme classifier [`dbCAN`](https://github.com/linnabrown/run_dbcan) (Zhang _et al._ 2018), configured using `cazomevolve`.

> Han Zhang, Tanner Yohe, Le Huang, Sarah Entwistle, Peizhi Wu, Zhenglu Yang, Peter K Busk, Ying Xu, Yanbin Yin; dbCAN2: a meta server for automated carbohydrate-active enzyme annotation, Nucleic Acids Research, Volume 46, Issue W1, 2 July 2018, Pages W95–W101, https://doi.org/10.1093/nar/gky418

Run the following command from the root of this directory. *Note: depending on the computer resources this may take multiple days to complete*

**Note:** The following commands MUST be run from the same directory containing the `db` directory created when installing `dbCAN` - the following commands presumt the `db` dir is located in the root of this repository.

```bash
scripts/run_dbcan.sh
```

After running dbCAN, use the following commands to parse the output from dbCAN and add the predicted CAZy family annotations, protein accessions and genomic accessions to the tab delimited lists created above.

The command runs the `cazomevolve` command `cazevolve_get_dbcan` which can be used to parse the output from `dbCAN` version 2 and version 3.

```bash
scripts/get_dbcan_cazymes.sh
```

At the end, two plain text files will be generated, containing tab separated data:

The _Pectobacteriaceae_ lists were written to:
1. `data/cazomes/pecto_fam_genomes`
2. `data/cazomes/pecto_fam_genomes_proteins`

## Run ANI analysis and construct dendrogram

The software package `pyani` [Pritchard et al](https://doi.org/10.1039/C5AY02550H) was used to perform an average nucleotide identify (ANI) comparison between all pairs of _Pectobacteriaceae_ genomes, using the ANIm method.

> Pritchard et al. (2016) "Genomics and taxonomy in diagnostics for food security: soft-rotting enterobacterial plant pathogens" Anal. Methods 8, 12-24

```bash
scripts/tree/ani/run_anim.sh
```

This created a pyani database in `data/tree`. Graphical outputs summarising the pyani analysis were written to `results/tree/anim`.

**May need to run to double check what output is created, and how to get the tsv file for generating the dendrogram**.

A dendrogram was reconstructed from the ANIm analysis using the bash script `build_anim_tree.sh`, which coordinated the extraction of the calclated ANI values from the local `pyani` database, replacing the `pyani` genome IDs with the NCBI genomic version accessions using the Python script `parse_anim_tab.py`, and coordinated the R script `build_anim_tree.R`, which build a distance matrix and used hierarchical clustering (using the 'single' method) to build a dendorgram that was written in Newick format to `data/tree/pyani_ani_tree.new`:
```bash
scripts/tree/ani/build_anim_tree.sh
```

## Add taxonomic classifications

Download the GTDB database dump from the [GTDB repository](https://data.gtdb.ecogenomic.org/releases/). Release 202.0 was used in the manuscript Hobbs et al. Save the database dump (TSV file) to `data/gtdb/` directory.

The bash script `add_tax.sh` was used to coordinate running `cazomevolve` to add taxonomic information to each genomic accession, in every tab delimited list of (i) CAZy family and genomic accession, and (ii) CAZy family, genomic accession and protein accession that was generated.

```bash
scripts/taxs/add_tax.sh <use email address> <path to gtdb tsv file>
```

Use Python script `add_ani_tax.py` to add the taxonomic information to the reconstructed ANI trees.

```bash
scripts/taxs/add_ani_tax.py
```

## Explore CAZome composition

Exploration of the CAZomes in the data set was preformed within a `jupyter notebook`, which is available in this repository  (the raw notebooks  is for downloading and re-running locally, the website version is for viewing the results):
* [raw notebook](https://github.com/HobnobMancer/SI_Hobbs_et_al_2023_Pecto/blob/master/notebooks/explore_pectobact_cazomes.ipynb)
* [website](https://hobnobmancer.github.io/SI_Hobbs_et_al_2023_Pecto/notebooks/explore_pectobact_cazomes.html)

Specifically, the analyses performed in the notebook was executed using the module `cazomevolve.cazome.explore`, which contains functions for exploring the CAZome annotated by `cazomevolve`.

### Build rare faction plots

The R script `build_rarefaction_plots.R` was used to estimate the degree of diversity and completeness of the CAZome annotations in the dataset, specifically using the R package [`Vegan`](https://github.com/vegandevs/vegan) (Dixon _et al._, 2003).

> Dixon, P. (2003), VEGAN, a package of R functions for community ecology. Journal of Vegetation Science, 14: 927-930. [https://doi.org/10.1111/j.1654-1103.2003.tb02228.x](https://doi.org/10.1111/j.1654-1103.2003.tb02228.x)

To repeat the analysis, run the following bash command from the root of the reposistory **after** having the Jupyter Notebook (otherwise the script will be unable to find the necessary input files):

```bash
scripts/rare_factions/build_rarefaction_plots.R
```

## Identify associating CAZy families using `coinfinder`

Use the tool [`coinfinder`](https://github.com/fwhelan/coinfinder) (Whelan et al.) to identify CAZy families that are present in the genome together more often than expected by chance and lineage.

> Fiona J. Whelan, Martin Rusilowicz, & James O. McInerney. "Coinfinder: detecting significant associations and dissociations in pangenomes." doi: [https://doi.org/10.1099/mgen.0.000338](https://doi.org/10.1099/mgen.0.000338)

**Generate circular trees and heatmaps:**

To reproduce the output from `coinfinder` in the same structure as presented in the manuscript (i.e. a circular tree surrounded by a heatmap), overwrite the file `network.R` in `coinfinder` with the respective R script in `scripts/coevolution`, and use the corresponding bash script:

* `network.R`: `scripts/coevolution_circular_network.R`
* bash: `scripts/coevolution/find_coevolving_pectobact.sh`

**Generate linear trees and heatmaps, with taxonomy information:**

The circular heatmap annotates each leaf of the tree with only the respective genomic version accession. To list the taxonomic infomration as well, on each leaf of the tree, overwrite the contents in the file `network.R` in `coinfinder` with the respective R script in `scripts/coevolution`, and use the respective bash script to configure `coinfinder`:

* `network.R`: `scripts/coevolution_taxs_rectangular_network.R`
* bash: `scripts/coevolution/find_coevolving_pectobact_with_tax.sh`
