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
  - [BLAST+ v. 2.12](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
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
- Genomic sequences of *A. paragallinarum* strains from [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=728)
- [DEG (Database of Essential Genes) Bacteria v. 10](http://origin.tubic.org/deg/public/index.php/download)
- [VFDB (Virulence Factor Database) full (set B)](https://www.mgc.ac.cn/VFs/download.htm)
- COG preformatted database from the "little_endian" directory (Cog_LE.tar.gz) and summary information about the CD models (cddid.tbl.gz) from the [CDD FTP-archive](https://ftp.ncbi.nih.gov/pub/mmdb/cdd/)
- Descriptions of COG (cog-20.def.tab) and COG functional categories (fun-20.tab) from the [COG2020 FTP-archive](https://ftp.ncbi.nih.gov/pub/COG/COG2020/)


## Acknowledgments
We would like to express our gratitude to CONICET (Consejo Nacional de Investigaciones Científicas y Técnicas) and to the Universidad Nacional de Córdoba (UNC) for their support and funding.
Also, we wish to thank the developers and contributors of the various tools, software, and databases that were used in the development of this project, without which this work would not have been possible.

## Contact
For questions or feedback, please contact Esperanza Felici at esperanza.felici@unc.edu.ar.

