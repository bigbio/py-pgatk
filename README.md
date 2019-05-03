# Python tools for ProteoGenomics Analysis Toolkit


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

#### Downloading ENSEMBL data.

Downloading data from [ENSEMBL](https://www.ensembl.org/info/data/ftp/index.html) can be done using the command `ensembl-downloader`. The current tool enables to download the following files for each taxonomy:

- GTF
- Protein Sequence (FASTA),
- CDS (FASTA)
- Variation (VCF))

By default the command `ensembl-downloader` download all file types for all the ensembl species.

```bash
$: python3.7 pypgatk.py ensembl-downloader -h
Usage: pypgatk.py ensembl-downloader [OPTIONS]

  This tool enables to download from enseml ftp the FASTA and GTF files

Options:
  -c, --config_file TEXT          Configuration file for the ensembl data
                                  downloader pipeline
  -o, --output_directory TEXT     Output directory for the peptide databases
  -fp, --folder_prefix_release TEXT
                                  Output folder prefix to download the data
  -t, --taxonomy TEXT             Taxonomy List (comma separated) that will be
                                  use to download the data from Ensembl
  -sg, --skip_gtf                 Skip the gtf file during the download
  -sp, --skip_protein             Skip the protein fasta file during download
  -sc, --skip_cds                 Skip the CDS file download
  -snr, --skip_ncrna              Skip the ncRNA file download
  -h, --help                      Show this message and exit.

```

Each of the file types can be skip using the corresponding option. For example, if the user do not want to download the the protein sequence fasta file, it can be skip by using the argument `pypgatk.py ensembl-downloader --skip_protein`

#### Downloading COSMIC data.

Downloading mutation data from [COSMIC](https://cancer.sanger.ac.uk/cosmic) is performed using the COMMAND ` cosmic-downloader`. The current COMMAND allows users to download the following files:

- Cosmic mutation file (CosmicMutantExport)
- Cosmic all genes (All_COSMIC_Genes)

```bash
$: python3.7 pypgatk.py cosmic-downloader -h
Usage: pypgatk.py cosmic-downloader [OPTIONS]

Options:
  -c, --config_file TEXT       Configuration file for the ensembl data
                               downloader pipeline
  -o, --output_directory TEXT  Output directory for the peptide databases
  -u, --username TEXT          Username for cosmic database -- please if you
                               don't have one register here
                               (https://cancer.sanger.ac.uk/cosmic/register)
  -p, --password TEXT          Password for cosmic database -- please if you
                               don't have one register here
                               (https://cancer.sanger.ac.uk/cosmic/register)
  -h, --help                   Show this message and exit.

```

__Note: In order to be able to download COSMIC data, the user should provide a user and password. Please first register in COSMIC database (https://cancer.sanger.ac.uk/cosmic/register).__

#### Downloading cBioPortal data.

Downloading mutation data from [cBioPortal](https://www.cbioportal.org/) is performed using the command `cbioportal-downloader`. cBioPortal store multiple studies (https://www.cbioportal.org/datasets) containing mutation data. Currently is not possible to search the studies by PubMedID, only can be search by study_id.

```bash

$: python3.7 pypgatk.py cbioportal-downloader -h
Usage: pypgatk.py cbioportal-downloader [OPTIONS]

Options:
  -c, --config_file TEXT       Configuration file for the ensembl data
                               downloader pipeline
  -o, --output_directory TEXT  Output directory for the peptide databases
  -l, --list_studies           Print the list of all the studies in cBioPortal
                               (https://www.cbioportal.org)
  -d, --download_study TEXT    Download an specific Study from cBioPortal --
                               (all to download all studies)
  -h, --help                   Show this message and exit.

```

The argument `-l` (`--list_studies`) allow the users to list all the studies store in cBioPortal. If the user is interested in only one study, it can use the argument `-d` (`--download_study`).



Contributions
-----------------------

- Yafeng Zhu ([yafeng](http://github.com/yafeng))
- Husen M. Umer ([husensofteng](https://github.com/husensofteng))
- Enrique Audain ([enriquea](https://github.com/enriquea))
- Yasset Perez-Riverol ([ypriverol](https://github.com/ypriverol))
