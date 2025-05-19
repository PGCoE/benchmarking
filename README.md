# PGCoE Benchmarking Study

This repository hosts the proposal and associated materials for a research study aimed at benchmarking whole genome sequencing (WGS) workflows for key respiratory pathogens. The goal is to develop standardized performance metrics that enable comparison across different sequencing and bioinformatics methods, specifically for public health surveillance use (not clinical diagnostics).

## Overview

The study will:

* Benchmark WGS performance for SARS-CoV-2, RSV A/B, and Influenza A/B.
* Include two implementation phases: (1) bioinformatic assembly from public data; (2) end-to-end lab + bioinformatics from provided RNA samples.
* Use a defined set of computational and laboratory standards to ensure reproducibility and comparability.

## Study Phases

1. **Planning** – Define lab and computational standards.
2. **Implementation I** – Participants assemble genomes from public SRA datasets.
3. **Implementation II** – Participants sequence provided RNA and analyze data with their own workflows.

## Key Metrics

* **Computational**: Genome completeness, accuracy, percent target reads, coverage depth, read quality.
* **Laboratory**: cDNA/library concentrations, fragment and read lengths, instrument/run metrics.

## Standards

* **Computational**: Publicly available SRA datasets with representative samples.
* **Laboratory**: Inactivated RNA and pathogen mixtures with known concentrations, provided by the Springer Lab at Harvard.

## Participation

Participants will receive protocols and be asked to contribute sequencing and analysis results, which will be compiled and compared. Results will be made public, though this is not a pass/fail evaluation—rather a tool for understanding variability and informing method selection.
