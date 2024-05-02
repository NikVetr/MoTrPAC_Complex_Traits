# The impact of exercise on gene regulation in association with complex trait genetics

This is the GitHub repository associated with the publication: 

>Vetr, Nikolai G., Nicole R. Gay, MoTrPAC Study Group, and Stephen B. Montgomery. "The Impact of Exercise on Gene Regulation in Association with Complex Trait Genetics." *Nature Communications* **15**, no. 1 (May 1, 2024). [https://doi.org/10.1038/s41467-024-45966-w](https://doi.org/10.1038/s41467-024-45966-w)

It primarily hosts scripts written in the R programming language, but interfacing with several other languages, that were used to carry out all analyses and generate all of the paper figures. Upon publication, it will also host the LaTeX files used to generate the final draft version of the paper. The names of all figure-generating scripts are given in figure order in `/scripts/figures/fig*_*`. Analysis scripts should be run first to generate larger intermediate files and are contained in `/scripts/analyses/analysis_*`. Software and data dependencies are listed below. Smaller data dependencies are also provided in the `/data/internal/` folder, where possible.

This repository contains paper main and supplementary figures (`/figures/`), MCMC output (`/data/MCMC_output`), and supplemental files (`/supplemental_files/`) generated by these scripts.

In the interests of reproducibility, we've endeavored for these scripts function as end-to-end pipelines, capable of taking users from public data downloads to final figures and paper, at least on the assumption that all dependencies have been properly installed. However, just because everything runs from fresh install on my machine (16" 2019 MBP, macOS 10.15.7) does not mean it will run seamlessly on yours. Please contact me at nikgvetr@stanford.edu if you encounter difficulties and I'll do my best to help.

GitHub is limited in how much storage can be devoted to each repository, but the total filesize of all files used in this paper sums to >0.5TB. To download these files, I've separated `/data/internal/` from `/data/external/`, and provide information below on where to obtain larger files from their original sources.

I personally stored this repository on an external drive, opting for simplicity to hard-code all paths at `/Volumes/2TB_External/MoTrPAC_Complex_Traits/`. As you are unlikely to store these files on my external drive, please make sure to modify the included scripts prior to use. One way to perform this modification leverages `find` and `sed`. In terminal, navigate to the scripts subdirectory in wherever you've cloned this repository to:

`cd /Path/To/Your/Directory/MoTrPAC_Complex_Traits/scripts/`

then, run the following command, taking care to appropriately escape forward-slashes:

`find . -type f -exec sh -c 'LANG=C sed -i "" "s/\/Volumes\/2TB_External\//\/Path\/To\/Your\/Directory\//g" "$0"' {} \;`

This should replace all instances of `/Volumes/2TB_External/` with the appropriate location on your machine.

## Dependencies

### R (4.0.4) packages used in these scripts include:

* fgsea: 1.16.0
* Cairo: 1.5.15
* DESeq2: 1.30.1
* EnsDb.Hsapiens.v79: 2.99.0
* EnsDb.Hsapiens.v86: 2.99.0
* MASS: 7.3.57
* MotrpacRatTraining6mo: 1.6.4
* MotrpacRatTraining6moData: 1.9.1
* arrow: 8.0.0
* biomaRt: 2.46.3
* caret: 6.0.92
* circlize: 0.4.15
* clusterProfiler: 3.18.1
* cmdstanr: 0.3.0.9000
* data.table: 1.14.8
* doParallel: 1.0.17
* dplyr: 1.1.2
* edgeR: 3.32.1
* ensembldb: 2.14.1
* fitdistrplus: 1.1.8
* foreach: 1.5.2
* ggplot2: 3.4.2
* invgamma: 1.1
* jpeg: 0.1.9
* jsonlite: 1.8.4
* ks: 1.13.5
* limma: 3.46.0
* org.Hs.eg.db: 3.12.0
* org.Rn.eg.db: 3.12.0
* parallel: 4.0.4
* plotrix: 3.8.2
* posterior: 1.2.2
* pracma: 2.3.8
* sparklyr: 1.7.7
* sparklyr.nested: 0.0.3
* testit: 0.13
* xlsx: 0.6.5

### External software and command-line tools

* [plink2](https://www.cog-genomics.org/plink/2.0/) (2.00a3)
* [GCTA](https://yanglab.westlake.edu.cn/software/gcta/#Overview) (1.93.2)
* [MESC](https://github.com/douglasyao/mesc)  
* [LDSC](https://github.com/bulik/ldsc) (1.0.1)
* [Stan](https://mc-stan.org/) (2.21.0)
* [revTWMR](https://github.com/eleporcu/revTWMR)
* [TWMR](https://github.com/eleporcu/TWMR)
* [GTEx Pipeline](https://github.com/broadinstitute/gtex-pipeline/)
* [Fusion](https://github.com/gusevlab/fusion_twas/tree/master)
* [MetaXcan](https://github.com/hakyimlab/MetaXcan)
* [GenArchDB](https://github.com/jlbren/GenArchDB)

Please follow the installation instructions at the above links to install these software. Please clone all GitHub repos into `/data/external/`, or else modify these scripts accordingly.

### Datasets  

The following external files are used by these scripts. Please download them and place them in `/data/external/`.

* `data/external/eqtl/`: Barbeira et al. `S-PrediXcan` results
available [here](https://zenodo.org/record/3518299#:~:text=Download-,spredixcan_eqtl.tar.gz,-md5%3Ac0474256186dc58ed41705475455ebee), provided as a 
supplement to their manuscript: <https://pubmed.ncbi.nlm.nih.gov/33499903/>. Untar and unzip before use.   
* `data/external/imputed_gwas_hg38_1.1/`: GWAS imputed and harmonized for colocalization 
analysis with GTEx v8. Provided [here](https://zenodo.org/record/3629742/files/harmonized_imputed_gwas.tar?download=1) by Barbeira et al. as a supplement to their 
manuscript: <https://pubmed.ncbi.nlm.nih.gov/33499903/>. Untar before use.
* `data/external/RGD_ORTHOLOGS_20201001.txt`: Data downloaded from RGD FTP website 
on 10/01/2020 by Pierre Jean (<ftp://ftp.rgd.mcw.edu/pub/data_release/orthologs/>, <gs://mawg-data/external-datasets/rat-id-mapping/RGD_ORTHOLOGS_20201001.txt>).
* `data/external/GTEx_Analysis_v8_eQTL_all_associations`: eQTL mapping file downloaded from GTEx at <gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_eQTL_all_associations/>
* `data/external/dea`: Differential expression analysis results from MoTrPAC PASS-1B downloaded from <gs://mawg-data/pass1b-06/transcript-rna-seq/dea/>
* `data/external/rna_dea_20210114.RData`: DE results in a compatible data format downloaded from <gs://mawg-data/pass1b-06/transcript-rna-seq/dea/>
* `data/external/transcript_rna_seq_20211008.RData`: Specific version of RNA DEA downloaded from <gs://motrpac-data-freeze-pass/pass1b-06/v1.1/analysis/transcriptomics/transcript-rna-seq/dea/transcript_rna_seq_20211008.Rdata>
* `data/external/RSID_POS_MAP.txt`: Mapping file between RSIDs and genome positions for build used <https://drive.google.com/file/d/1COG3UXpdMtfDgF9QQ3LCE4UrgLQQRCNK/view?usp=drive_link>
* `data/external/old_dea_deseq_20201121`: Earlier DE results that explicitly incorporate a sex term  <https://drive.google.com/file/d/1kC84BEUUWzOEN30HoqwCWXweKyRixZ9X/view?usp=drive_link>
* `data/external/GTEx_Analysis_v8_eQTL`: Filtered eQTL results from GTEx <https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar>

As a reminder, scripts for generating particular figures will sometimes rely on files generated previously by `analyses` scripts. Unfortunately, filesizes for these sometimes exceeded maximum storage allowances on open repositories such as Zenodo. I've tried to organize these in a manner where those intermediate files are dynamically generated during execution of the script, but if something is not cooperating, please do not hesitate to reach out to me at nikgvetr@stanford.edu for assistance!
