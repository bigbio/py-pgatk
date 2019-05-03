# Python tools for ProteoGenomics Analysis Toolkit
============================

**pypgatk** is a Python framework for part of the ProteoGenomics Analysis Toolkit. It provides different tools for proteogenomics data analysis such as:

In order to execute a task in `pypgatk` the user should use a `COMMAND` that perform the specific task and the specific task arguments/options:

```bash
$: python3.7 pypgatk -h
Usage: pypgatk.py [OPTIONS] COMMAND [ARGS]...

  This is the main tool that give access to all commands and options provided by the pypgatk

Options:
  -h, --help  Show this message and exit.

Commands:
  cbioportal-downloader    Command to download the the cbioportal studies
  cbioportal-to-proteindb  Command to translate cbioportal mutation data into proteindb
  cosmic-downloader        Command to download the cosmic mutation database
  cosmic-to-proteindb      Command to translate Cosmic mutation data into proteindb
  ensembl-downloader       Command to download the ensembl information

```

Data downloader
----------------

The `Data downloader` is a set of `COMMANDs` to download data from different Genomics data providers such as ENSEMBL, COSMIC or cBioPortal.

#### Downloading data ENSEMBL database.

from different genomics sources:
     - ENSEMBL data downloader
       - GTF
       - Protein Sequence (FASTA),
       - CDS (FASTA)
       - Variation (VCF))
     - Cosmic data downloader
       - Cosmic mutation file (CosmicMutantExport)
       - Cosmic all genes (All_COSMIC_Genes)
     - Cbioportal `Study` data downloader
- Custom protein database generation from Genomics Cancer mutation datasets:
     - CBioportal (https://www.cbioportal.org/)
     - COSMIC (https://cancer.sanger.ac.uk/cosmic)


Contributions
-----------------------

- Yafeng Zhu ([yafeng](http://github.com/yafeng))
- Husen M. Umer ([husensofteng](https://github.com/husensofteng))
- Enrique Audain ([enriquea](https://github.com/enriquea))
- Yasset Perez-Riverol ([ypriverol](https://github.com/ypriverol))
