# Overview

This repository contains code, data-processing pipelines, and documentation for an analysis project in systems vaccinology. 

This project makes use of transcriptomic data (microarrays / bulk RNA-seq) together with immune-response profiling to compare vaccine-induced transcriptional signatures across studies and to methodologically evaluate approaches for discovering robust biomarkers of vaccine response.

It uses publically available data from the Human Immune Project Consortium's ImmuneSignatures2 project. The raw data used for this project can be found at [immunespace.org](immunespace.org) and is not stored on this repository.

The work focuses on developing and applying reproducible, well-documented analysis workflows for:

 - the preprocessing of public vaccine transcriptomic datasets,
 - the comparison of gene-set differential analysis across multiple vaccines
 - the robust identification of gene-set biomarkers of the downstream vaccine response using supervised learning
 - sensitivity and robustness analyses 
