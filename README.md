# Supplementary Information: Hobbs _et al._, 2023

## Association between niche adaptation and evolution of carbohydrate active enzymes in _Pectobacteriaceae_

This repository contains supplementary information for analyses reported in Hobbs et al. (2023), describing the evolution of the CAZome.

Run all commands provided in the walkthrough from the root of this directory.

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

## Reproducing analyses (quickstart)

You can use this archive to browse, validate, reproduce, or build on the phylogenomics analysis for the Hobbs et al. (2023) manuscript.

We recommend creating a conda environment specific for this activity, for example using the commands:
```bash
conda create -name pectobacteriaceae python=3.9 -y
conda activate pectobacteriaceae
conda install --file requirements.txt -y
```

1. Download datasets
    1. `download_pecto_genomes.sh` - Downloads all genomes from NCBI Assembly when querying using the term 'Pectobacteriaceae'
    2. `download_same_pecto_genomes.sh` - Downloads the genomes used in the analysis presented in the manuscript
    3. `download_pd_genomes.sh` - Configures `ncbi-genome-download` to retrieve genomes for a provided list of Assembly accessions
    4. `build_cazyme_db.sh` - Build a local CAZyme db
2. Annotate CAZomes
    1. `get_cazy_cazymes.sh` - Retrieve CAZy family annotations from the local CAZyme db for proteins in the download genomes for both the _Pectobacteriaceae_, and _Pectobacterium_ and _Dickeya_ datasets.

# Reproducing the analyses

## Download data

### _Pectobacteriaceae_

`cazomevolve` was used to download all complete _Pectobacteriaceae_ genomic assemblies (in genome sequence and protein sequence FASTA file format) from NCBI Assembly, by querying the NCBI Taxonomy database and retrieving all genomic assemblies linked to the _Pectobacteriaceae_ (NCBI:txid1903410). To repeat this method, run the following command from the root of this directory:

```bash
# download Pectobacteriaceae genomes from GenBank
scripts/download/download_pecto_genomes.sh <email>
```

**Note:** With the continual addition of new genomic assemblies to the NCBI Assembly database, repeating the download of _Pectobacteriaceae_ genomes may generate a different dataset to that presented in Hobbs _et al._. To repeat the analysis presented in the manuscript, run the following command from the root of the directory to configure `ncbi-genome-download` to download the genomic assemblies of the genomes used in the manuscript:

```bash
scripts/download/download_same_pecto_genomes.sh
```

<<<<<<< HEAD
In both cases, the downloaded genomic sequence files were written to the dir `data/pectobact/genomes`, the downloaded protein FASTA files were written to `data/pectobact/proteomes`.
=======
In both cases, the downloaded genomic sequence files were written to the dir `data/pectobacteriaceae/genomes`, the downloaded protein FASTA files were written to `data/pectobacteriaceae/proteomes`.
>>>>>>> Combine READMEs

### _Pectobacterium_ and _Dickeya_

Configure `ncbi-genome-download` to download genomic assemblies (genome sequence and protein sequence FASTA files) for the assembly accessions listed in `supplementary_files/supplementary_file_1`, which include the genomes for _Pectobacterium_, _Dickeya_, and _Hafnia alvei_ ATCC 13337, and decompresses the downloaded files.

```bash
# download Pectobacterium and dickeya genomes from RefSeq
# listed in a supplementary_file_1
scripts/download_pd_genomes.sh
```

<<<<<<< HEAD
The downloaded genomic sequence files were written to the dir `data/pecto_dic/genomes`, the downloaded protein FASTA files were written to `data/pecto_dic/proteomes`.
=======
The downloaded genomic sequence files were written to the dir `data/pd/genomes`, the downloaded protein FASTA files were written to `data/pd/proteomes`.
>>>>>>> Combine READMEs

### Build a local CAZyme database

Configure using [`cazy_webscraper`](https://hobnobmancer.github.io/cazy_webscraper/) (Hobbs _et al., 2022) to download all data from the CAZy database, and compile the data into a local CAZyme database.

> cazy_webscraper: local compilation and interrogation of comprehensive CAZyme datasets
Emma E. M. Hobbs, Tracey M. Gloster, Leighton Pritchard
bioRxiv 2022.12.02.518825; doi: https://doi.org/10.1101/2022.12.02.518825

```bash
# create a local CAZyme database
scripts/build_cazyme_db.sh <email>
```

This generated the local CAZyme database `data/cazy/cazy_db`.


## Annotate CAZomes

### Get CAZy classified CAZymes

Configure using `cazomevolve` to identify CAZymes classified in the local CAZyme database, for both the _Pectobacteriaceae_, and _Pectobacterium_ & _Dickeya data sets.

```bash
scipts/get_cazy_cazymes.sh
```

For _Pectobacteriaceae_, and the _Pectobacterium_ & _Dickeya_ datasets two tab delimited lists were created:
1. Listing the CAZy family accession and genomic accession per line
2. Listing the CAZy family, genomic accession and protein accession per line

The _Pectobacteriaceae_ lists were written out respectively to:
1. data/pectobacteriaceae/cazomes/pecto_fam_genomes \
2. data/pectobacteriaceae/cazomes/pecto_fam_genomes_proteins 

The _Pectobacterium_ and _Dickeya_ lists were written out respectively to:
1. data/pd/cazomes/pd_fam_genomes \
2. data/pd/cazomes/pd_fam_genomes_proteins 

Proteins in the download protein FASTA files that were not listed in the local CAZyme database were written to `data/pectobacteriaceae/cazomes/dbcan_input` (for _Pectobacteriaceae_) and `data/pd/cazomes/dbcan_input` (for _Pectobacterium_ and _Dickeya_).

### Get dbCAN predicted CAZymes

To retrieve the most comprehensive CAZome for each genome, protein sequences not found in the local CAZyme database were parsed by the CAZyme classifier [`dbCAN`](https://github.com/linnabrown/run_dbcan) (Zhang _et al._ 2018), configured using `cazomevolve`.

> Han Zhang, Tanner Yohe, Le Huang, Sarah Entwistle, Peizhi Wu, Zhenglu Yang, Peter K Busk, Ying Xu, Yanbin Yin; dbCAN2: a meta server for automated carbohydrate-active enzyme annotation, Nucleic Acids Research, Volume 46, Issue W1, 2 July 2018, Pages W95–W101, https://doi.org/10.1093/nar/gky418

Run the following command from the root of this directory. *Note: depending on the computer resources this may take multiple days to complete*

**Note:** The following commands MUST be run from the same directory containing the `db` directory created when installing `dbCAN`.

```bash
# run dbcan for the pectobacteriaceae data set
scripts/run_dbcan_pectobacteriaceae.sh
```

```bash
# run dbcan for the pectobacterium&dickeya data set
scripts/run_dbcan_pd.sh
```

After running dbCAN, using the following commands to parse the output from dbCAN and add the predicted CAZy family annotations, protein accessions and genomic accessions to the tab delimited lists created above.

The commands run the `cazomevolve` command `cazevolve_get_dbcan` which can be used to parse the output from `dbCAN` version 2 and version 3.

```bash
# parse pectobacteriaceae dbcan output
scripts/get_dbcan_cazymes_pecto.sh
```

```bash
# parse pectobacterium and dickeya dbcan output
scripts/get_dbcan_cazymes_pd.sh
```

## Explore CAZome composition

The module `cazomevolve.cazome.explore` contains functions for exploring the CAZome annotated by `cazomevolve`. These are:

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

## Run ANI analysis

## Reconstruct phylogenetic tree

## Identify networks of co-evolving CAZy families

