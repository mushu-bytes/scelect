# Scelect

This repository contains a snakemake pipeline that implements a grid search of single cell integration and normalization methods that have been prototyped in this repository: https://github.com/mushu-bytes/scRNAseq-benchmarking?tab=readme-ov-file. 

### Dependencies

In order to use this snakemake pipeline, please install snakemake through pip or conda. Installation instructions can be found here: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

In addition, this snakemake runs two jobs in parallel -- integration and normalization methods implemented in Seurat, and integration and normalization methods hosted in Scanpy. In order to run the two jobs, docker containers must be built from Dockerfiles at workflows/scripts/python_scripts and workflows/scripts/r_scripts, which are used to containerize the Python and R environments.

### Execution

After snakemake and the Docker images have been built, the snakefile within workflows/ can be used to run the snakemake pipeline via snakemake --cores=1. Make sure to edit the snakemake file to fit your environment, such as adding specific filepaths to the data 

### Usage and Developer Notes

Given that this pipeline is running jobs through both Seurat and Scanpy libraries, interoperatability between datasets is tricky. So far, only .h5ad and .h5 files from 10x genomics are supported when reading running Scanpy normalization and Integration. On the other hand, 10x formated h5 files, .rds files, .h5Seurat files, and .h5ad files are available for usage through our R scripts.

Although our team aimed to synchronize outputs from jobs run in R and jobs run in Python, the outputs are not standardized. As a result, it is not trivial to see an apples-apples comparison from a singular reports, but we aimed to dump as much information relevant to our metrics and evaluation as outputs.



