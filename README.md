# Reverse Vaccinology for Infectious Coryza: Identifying Vaccine Targets against *Avibacterium paragallinarum*

This repository contains the scripts and workflows used in the identification and evaluation of antigen candidates for the development of a veterinary vaccine against *Avibacterium paragallinarum*, the causative agent of infectious coryza in chickens. 

## Requirements

### Software Requirements

- **Operating System:**
  The scripts were developed and tested on Ubuntu 22.04, but they have also been verified to work on Ubuntu 20.04.


- **Programming Languages and Packages/Libraries:**
  - R (at least v. 4.4.1) and these R packages:
     - [tidyverse](https://www.tidyverse.org/) v 2.0.0
     - [BiocManager](https://www.bioconductor.org/install/) v 1.30.25
     - [rlang](https://rlang.r-lib.org/) v. 1.1.4
     - [Rpdb](https://cran.r-project.org/web/packages/Rpdb/index.html) v. 2.3.4
     - [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html) v. 4.2.36
     - [BALCONY](https://cran.r-project.org/src/contrib/Archive/BALCONY/) v. 0.2.10
         
  - Python (at least v. 3.9.12) and these python packages:
    - [pyperclip v. 1.9.0](https://pypi.org/project/pyperclip/)
    - [Selenium WebDriver v. 4.27.1](https://pypi.org/project/selenium/)


- **Bioinformatics Tools:**
  - [BLAST+ v. 2.12](https://www.ncbi.nlm.nih.gov/books/NBK279690/) and [RPS-BLAST executables v. 2.13.0](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
  - [PSORTb v. 3.0.3](https://www.psort.org/psortb/index.html)
  - [SignalP v. 6.0](https://services.healthtech.dtu.dk/services/SignalP-6.0/)
  - [VaxiJen v. 3.0](https://www.ddg-pharmfac.net/vaxijen3/home/)
  - [SPAAN](https://sourceforge.net/projects/adhesin/files/SPAAN/)
  - [SeqKit v. 2.8.2](https://bioinf.shenwei.me/seqkit/#installation)
  - [DeepTMHMM](https://dtu.biolib.com/DeepTMHMM)
  - [EMBOSS v. 6.6.0](https://emboss.sourceforge.net/)
  - [MAFFT v. 7.505](https://mafft.cbrc.jp/alignment/software/)
  - [CD-HIT v. 4.7](https://sites.google.com/view/cd-hit/home?authuser=0)


- **Other Tools**
   - [Chrome browser](https://www.google.com/intl/en_uk/chrome/)
   - [A ChromeDriver binary](https://googlechromelabs.github.io/chrome-for-testing/#stable) (compatible with your Chrome version)

### Data
- Genomic sequences of *Av. paragallinarum* strains from [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=728)
- [DEG](http://origin.tubic.org/deg/public/index.php/download) (Database of Essential Genes) Bacteria v. 10
- [VFDB](https://www.mgc.ac.cn/VFs/download.htm) (Virulence Factor Database) full (set B)
- COG preformatted database from the "little_endian" directory (Cog_LE.tar.gz) and summary information about the CD models (cddid.tbl.gz) from the [CDD FTP-archive](https://ftp.ncbi.nih.gov/pub/mmdb/cdd/)
- Descriptions of COG (cog-20.def.tab) and COG functional categories (fun-20.tab) from the [COG2020 FTP-archive](https://ftp.ncbi.nih.gov/pub/COG/COG2020/)

## Usage
0. Download all the required data. Ensure all external tools are installed and accessible in your system's PATH and install all necessary R and python packages

1. Create the required blastp databases with 0.DB_creation.sh script:
   ```bash
   ./0.DB_creation.sh <path/to/database_directory>
   ```

2. Create a file_of_paths (See [example]()). You only need one for use in all the scripts

3. Analyze the different Av. paragallinarum strains with 1.AvP_RV_Analysis.sh script:
   ```bash
   ./1.AvP_RV_Analysis.sh <path/to/file_of_paths>
   ```

4. Analyze the different experimentally tested antigens with 2.AgProtect_RV_Analysis.sh script
   ```bash
   ./2.AgProtect_RV_Analysis.sh <path/to/file_of_paths>
   ```

5. Merge and polish all the results with 3.Final_RV_Analysis.sh script
   ```bash
   ./3.Final_RV_Analysis.sh <path/to/file_of_paths>
   ```
   
## File Structure

```
.
├── AgProtect_data
│   ├── Final_results-AgProtect.tsv    #RV Analysis of the experimentally tested antigens
│   ├── protein.faa      # Fasta file of the collection of experimentally tested antigens
│   ├── BacteriaTypes.tsv         # Gram stain and genus of each antigen
│   └── HostList.tsv         # Hosts for each antigen
├── 0.DB_creation.sh
├── 1.AvP_RV_Analysis.sh
├── 1.AvP_RV_Analysis_RScripts
│   ├── COG.R
│   ├── DEG.R
│   ├── DeepTMHMM.R
│   ├── Final_polishing.R
│   ├── Info.R
│   ├── SignalP.R
│   ├── VaxiJen.R
│   ├── Conservation.R
│   ├── EMBOSS.R
│   ├── Homology_Analysis.R
│   ├── PSORTb.R
│   ├── SPAAN.R
│   └── VFDB.R
├── 2.AgProtect_RV_Analysis.sh
├── 2.AgProtect_RV_Analysis_RScripts
│   ├── COG_AgProtect.R
│   ├── DEG_AgProtect.R
│   ├── Final_polishing_AgProtect.R
│   ├── Info_AgProtect.R
│   ├── SignalP_AgProtect.R
│   ├── VaxiJen_AgProtect.R
│   ├── EMBOSS_AgProtect.R
│   ├── Homology_Analysis_AgProtect.R
│   ├── DeepTMHMM_AgProtect.R
│   ├── PSORTb_AgProtect.R
│   ├── SPAAN_AgProtect.R
│   └── VFDB_AgProtect.R
├── 3.Final_RV_Analysis.sh
├── 3.Final_RV_Analysis_RScripts
│   └──Final_merge.R
├── VaxiJen.py
├── Final_results        # Polished output files
│   ├── Conservation_results.tsv
│   ├── Final_results-GCA_020892835.1.tsv
│   ├── Final_results-GCF_000221945.2.tsv
│   ├── Final_results-GCF_000348525.1.tsv
│   ├── Final_results-GCF_000442905.1.tsv
│   ├── Final_results-GCF_000969305.1.tsv
│   ├── Final_results-GCF_000969315.2.tsv
│   ├── Final_results-GCF_002921155.1.tsv
│   ├── Final_results-GCF_003204245.1.tsv
│   ├── Final_results-GCF_003204295.1.tsv
│   ├── Final_results-GCF_003864155.1.tsv
│   ├── Final_results-GCF_004212345.1.tsv
│   ├── Final_results-GCF_004212355.1.tsv
│   ├── Final_results-GCF_004212375.1.tsv
│   ├── Final_results-GCF_004212415.1.tsv
│   ├── Final_results-GCF_004931865.1.tsv
│   ├── Final_results-GCF_008633115.1.tsv
│   ├── Final_results-GCF_011765605.1.tsv
│   ├── Final_results-GCF_012849875.1.tsv
│   ├── Final_results-GCF_012850075.1.tsv
│   ├── Final_results-GCF_012850275.1.tsv
│   ├── Final_results-GCF_012850455.1.tsv
│   ├── Final_results-GCF_012850655.1.tsv
│   ├── Final_results-GCF_012850815.1.tsv
│   ├── Final_results-GCF_012850955.1.tsv
│   ├── Final_results-GCF_015354775.1.tsv
│   ├── Final_results-GCF_015354785.1.tsv
│   ├── Final_results-GCF_015354815.1.tsv
│   ├── Final_results-GCF_015354825.1.tsv
│   ├── Final_results-GCF_015354845.1.tsv
│   ├── Final_results-GCF_015354865.1.tsv
│   ├── Final_results-GCF_015354885.1.tsv
│   ├── Final_results-GCF_015354895.1.tsv
│   ├── Final_results-GCF_015354935.1.tsv
│   ├── Final_results-GCF_015354955.1.tsv
│   ├── Final_results-GCF_015354975.1.tsv
│   ├── Final_results-GCF_015354985.1.tsv
│   ├── Final_results-GCF_015354995.1.tsv
│   ├── Final_results-GCF_015355035.1.tsv
│   ├── Final_results-GCF_015355055.1.tsv
│   ├── Final_results-GCF_015355065.1.tsv
│   ├── Final_results-GCF_015355095.1.tsv
│   ├── Final_results-GCF_015355125.1.tsv
│   ├── Final_results-GCF_019822885.1.tsv
│   ├── Final_results-GCF_020523945.1.tsv
│   ├── Final_results-GCF_022964995.1.tsv
│   ├── Final_results-GCF_025884415.1.tsv
│   ├── Final_results-GCF_025884555.1.tsv
│   ├── Final_results-GCF_026802155.1.tsv
│   ├── Final_results-GCF_026811325.1.tsv
│   ├── Final_results-GCF_036348495.1.tsv
│   ├── Final_results-GCF_036348515.1.tsv
│   ├── Final_results-GCF_036348535.1.tsv
│   ├── Final_results-GCF_036348555.1.tsv
│   ├── Final_results-GCF_036351875.1.tsv
│   ├── Final_results-GCF_036400095.1.tsv
│   ├── Final_results-GCF_036400115.1.tsv
│   ├── Final_results-GCF_036400155.1.tsv
│   ├── Final_results-GCF_036400175.1.tsv
│   ├── Final_results-GCF_036400195.1.tsv
│   ├── Final_results-GCF_036400215.1.tsv
│   ├── Final_results-GCF_036400235.1.tsv
│   ├── Final_results-GCF_036400255.1.tsv
│   ├── Final_results-GCF_036400275.1.tsv
│   ├── Final_results-GCF_036400285.1.tsv
│   ├── Final_results-GCF_036400315.1.tsv
│   ├── Final_results-GCF_036400335.1.tsv
│   ├── Final_results-GCF_036400355.1.tsv
│   ├── Final_results-GCF_036400375.1.tsv
│   ├── Final_results-GCF_036400395.1.tsv
│   ├── Final_results-GCF_036400415.1.tsv
│   ├── Final_results-GCF_036400435.1.tsv
│   ├── Final_results-GCF_036400455.1.tsv
│   ├── Final_results-GCF_036400475.1.tsv
│   ├── Final_results-GCF_036400495.1.tsv
│   ├── Final_results-GCF_036400515.1.tsv
│   ├── Final_results-GCF_036400535.1.tsv
│   ├── Final_results-GCF_036400555.1.tsv
│   ├── Final_results-GCF_036400595.1.tsv
│   ├── Final_results-GCF_036400615.1.tsv
│   ├── Final_results-GCF_036400635.1.tsv
│   ├── Final_results-GCF_036400645.1.tsv
│   ├── Final_results-GCF_036400675.1.tsv
│   ├── Final_results-GCF_036400695.1.tsv
│   ├── Final_results-GCF_036400715.1.tsv
│   ├── Final_results-GCF_036400735.1.tsv
│   ├── Final_results-GCF_036400755.1.tsv
│   ├── Final_results-GCF_036400775.1.tsv
│   ├── Final_results-GCF_036400795.1.tsv
│   ├── Final_results-GCF_036400835.1.tsv
│   ├── Final_results-GCF_036400895.1.tsv
│   ├── Final_results-GCF_036400915.1.tsv
│   ├── Final_results-GCF_036400955.1.tsv
│   ├── Final_results-GCF_037080095.1.tsv
│   ├── Final_results-GCF_900445345.1.tsv
│   └── Final_results-GCF_900450705.1.tsv
└── README.md            # Documentation
```

## Acknowledgments
We would like to express our gratitude to CONICET (Consejo Nacional de Investigaciones Científicas y Técnicas) and to the Universidad Nacional de Córdoba (UNC) for their support and funding.
Also, we wish to thank the developers and contributors of the various tools, software, and databases that were used in the development of this project, without which this work would not have been possible.

## Contact
For questions or feedback, please contact Esperanza Felici at esperanza.felici@unc.edu.ar.

