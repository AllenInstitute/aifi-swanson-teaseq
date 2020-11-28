## aifi-swanson-teaseq

#### Code related to Swanson, *et. al*. for scATAC-seq, ICICLE-seq,   
#### and TEA-seq data processing, analysis and visualization.  

![TEA-seq overview](common/2020-11-20_TEA-seq_overview.png)

<a name="con"></a>

### Contents

#### [Manuscript/Citation](#man)  

#### [Protocols](#pro)  

#### [GEO Repository](#geo)  

#### [General Notes](#gen)  

#### [scATAC-seq](#sca)  
- [Requirements/Inputs](#sca-req)  
- [Datasets](#sca-dat)  
- [Preprocessing](#sca-pre)  
- [Analysis](#sca-ana)  

#### [BarCounter](#bar)  
- [Main Repository](#bar-mai)  
- [Requirements/Inputs](#bar-req)  
- [Usage](#bar-usa)  

#### [ICICLE-seq](#ici)  
- [Requirements/Inputs](#ici-req)  
- [Datasets](#ici-dat)  
- [Preprocessing](#ici-pre)  
- [Analysis](#ici-ana)  

#### [TEA-seq](#tea)  
- [Requirements/Inputs](#tea-req)  
- [Datasets](#tea-dat)  
- [Preprocessing](#tea-pre)  
- [Analysis](#tea-ana)  

#### [References](#ref)

#### [Legal](#leg)
- [License](#leg-lic)  
- [Level of Support](#leg-lev)  
- [Contributions](#leg-con)  

------------

<a name="man"></a>

### Manuscript

The TEA-seq manuscript is currently available on bioRxiv at:  
https://www.biorxiv.org/content/10.1101/2020.09.04.283887v2  

Citation information:  
Elliott Swanson, Cara Lord, Julian Reading, Alexander T. Heubeck, Adam K. Savage, Richard Green, Xiao-jun Li, Troy R. Torgerson, Thomas F. Bumol, Lucas T. Graybuck, Peter J. Skene. *TEA-seq: a trimodal assay for integrated single cell measurement of transcripts, epitopes, and chromatin accessibility* (2020). bioRxiv 2020.09.04.283887; doi: https://doi.org/10.1101/2020.09.04.283887

[Return to Contents](#con)

------------

<a name="pro"></a>

### Protocols

A detailed bench protocol for TEA-seq written by Elliott Swanson is [available on protocols.io](https://www.protocols.io/edit/tea-seq-bpp2mmqe).

[Return to Contents](#con)

------------

<a name="geo"></a>

### GEO Repository

Data from Swanson, *et al.* will be available on GEO at [Series  GSE158013](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158013)

[Return to Contents](#con)

------------

<a name="gen"></a>

### General Notes

**Genome Builds**  
All samples used in Swanson, *et al.* were from human donors. We use *GRCHg38/hg38* genome builds throughout our processing and analysis.

**Computing Environment**
Preprocessing and analysis scripts were generated and run on Linux/Unix-like platforms (Debian/Ubuntu) with R >= v3.6.3 and 4.0.2 .  

[Return to Contents](#con)

------------

<a name="sca"></a>

### scATAC-seq

<a name="sca-req"></a>

#### scATAC-seq Requirements/Inputs

**Sequencing Data**  
The scATAC-seq preprocessing and analysis pipelines are built around libraries with the [10x Genomics scATAC-seq read structure](https://support.10xgenomics.com/single-cell-atac/sequencing/doc/specifications-sequencing-requirements-for-single-cell-atac):  
- `I1`: 16 nt, 10x Cell Barcodes  
- `I2`: 8 nt, i7 Well/sample Indexes  
- `R1`: 50 nt, ATAC-seq fragment insertion  
- `R2`: 50 nt, ATAC-seq fragment insertion  

**Software**  
Alignment: 10x Genomics [`cellranger-atac count` >= v1.0.0](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/count)   
Preprocessing: [`bedtools2`](https://github.com/arq5x/bedtools2/releases)  
Paralellization: [`GNU Parallel`](https://www.gnu.org/software/parallel/)  

[Return to Contents](#con)

<a name="sca-dat"></a>

#### scATAC-seq Datasets  

Because our samples all originate form human donors, raw data (FASTQ files) used for Swanson, *et al.* are in the process of submission to dbGAP for controlled access.  

Processed data from scATAC-seq datasets used in Swanson, *et al.* will be available for download from these GEO Samples:  

| Accession | PBMC Type | Perm/Nuc | Prep | Purification | Pur. Method | WellID |
| ---       | ---       | ---      | ---  | ---          | ---         | ---     |
| [GSM4784064](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784064) |  leukapheresis | Perm | 0.01% Dig. |none | none | B003-W7 |
| [GSM4784065](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784065) |  leukapheresis | Perm | 0.01% Dig. | FACS | D/D-Neu. | B003-W8 |
| [GSM4784066](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784066) |  ficoll | Nuc | 1x 10xNIB |none | none | X024-W1 |
| [GSM4784067](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784067) |  ficoll | Nuc | 0.25x 10xNIB |none | none | X024-W2 |
| [GSM4784068](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784068) |  ficoll | Nuc | 0.1x 10xNIB |none | none | X024-W3 |
| [GSM4784069](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784069) |  ficoll | Nuc | 1x ANIB |none | none | X024-W4 |
| [GSM4784070](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784070) |  ficoll | Perm | 0.01% Dig. |none | none | X025-W1 |
| [GSM4784071](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784071) |  ficoll | Perm | 0.05% Dig. |none | none | X025-W2 |
| [GSM4784072](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784072) |  ficoll | Perm | 0.1% Dig. |none | none | X025-W3 |
| [GSM4784073](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784073) |  ficoll | Perm | 0.2% Dig. |none | none | X025-W4 |
| [GSM4784074](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784074) |  ficoll | Perm | 0.01% Dig. |none | none | X027-W1 |
| [GSM4784075](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784075) |  ficoll | Perm | 0.01% Dig. | FACS | D/D-Neu. | X027-W3 |
| [GSM4784076](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784076) |  ficoll | Perm | 0.01% Dig. | Mag. Bead | Negative | X032-W1 |
| [GSM4784077](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784077) |  ficoll | Perm | 0.01% Dig. | Mag. Bead | anti-CD15 |  X032-W2 |
| [GSM4784078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784078) |  leukapheresis | Perm | 0.01% Dig. | Mag. Bead | Negative | X032-W3 |
| [GSM4784079](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784079) |  leukapheresis | Perm | 0.01% Dig. | Mag. Bead | anti-CD15 | X032-W4 |
| [GSM4784080](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784080) |  ficoll | Nuc | 1x 10xNIB |none | none | X041-W1 |
| [GSM4784081](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784081) |  ficoll | Nuc | 1x 10xNIB | FACS | D/D | X041-W2 |
| [GSM4784082](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784082) |  ficoll | Nuc | 1x 10xNIB | FACS | D/D-Neu. | X041-W3 |
| [GSM4784083](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784083) |  ficoll | Perm | 0.01% Dig. |none | none | X041-W4 |
| [GSM4784084](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784084) |  ficoll | Perm | 0.01% Dig. | FACS | D/D | X041-W5 |
| [GSM4784085](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784085) |  ficoll | Perm | 0.01% Dig. | FACS | D/D-Neu. | X041-W6 |

[Return to Contents](#con)

<a name="sca-pre"></a>

#### scATAC-seq Preprocessing

*1. cellrange-atac count*  

*2. Custom Preprocessing*  

*3. Custom QC filtering and reporting*  

[Return to Contents](#con)

<a name="sca-ana"></a>

#### scATAC-seq Analysis

*Figure 1*  

*Figure 1 - Figure Supplement 2*  

*Figure 1 - Figure supplement 5*  

*Figure 2*  

[Return to Contents](#con)

------------

<a name="bar"></a>

### BarCounter  

Barcode quantification was performed using the C program BarCounter developed by Elliott Swanson.  

We utilize BarCounter for the processing of ADT counts for both the ICICLE-seq and TEA-seq analyses, below.  

<a name="bar-mai"></a>

#### Main Repository  

For convenience, a version of BarCounter is linked to this repository, and will be cloned along with other repository code.  

The main repository for Barcounter releases is at [https://github.com/AllenInstitute/BarCounter](https://github.com/AllenInstitute/BarCounter).

[Return to Contents](#con)

<a name="bar-req"></a>

#### BarCounter Requirements/Inputs

In its current implementation, BarCounter expects reads to conform to the cell and UMI barcode positions used by 10x Genomics for 3' scRNA-seq and the CellRanger Multiome RNA + ATAC kits. In addition, standard Illumina-like FASTQ naming conventions are required.  

**Sequencing Data**  
- `R1`: 28 nt - First 16 nt are 10x Cell Barcodes; next 12 nt are UMIs  
- `R2`: 15 nt (or longer) - first 15 nt are Antibody/ADI barcodes. Reads may be longer if co-sequenced with RNA/ATAC data  

[Return to Contents](#con)

<a name="bar-usa"></a>

#### BarCounter Usage

*TODO*  

[Return to Contents](#con)

------------

<a name="ici"></a>

### ICICLE-seq

*TODO*  

<a name="ici-req"></a>

#### ICICLE-seq Requirements

*TODO*  

[Return to Contents](#con)

#### ICICLE-seq Datasets

Because our samples all originate form human donors, raw data (FASTQ files) used for Swanson, *et al.* are in the process of submission to dbGAP for controlled access.  

Processed data from ICICLE-seq datasets used in Swanson, *et al.* will be available for download from these GEO Samples:  

| Accession | PBMC Type | Perm/Nuc | Prep | Purification | Pur. Method | WellID |
| ---       | ---       | ---      | ---  | ---          | ---         | ---     |
| [GSM4784086](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784086) | leukapheresis | Perm | 0.01% Dig. | FACS | D/D-Neu. | X044-W1 |
| [GSM4784087](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784087) | leukapheresis | Perm | 0.01% Dig. | FACS | D/D-Neu. | X044-W2 |
| [GSM4784088](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4784088) | leukapheresis | Perm | 0.01% Dig. | FACS | D/D-Neu. | X044-W3 |

[Return to Contents](#con)

<a name="ici-pre"></a>

#### ICICLE-seq Preprocessing

*TODO*  

[Return to Contents](#con)

<a name="ici-ana"></a>

#### ICICLE-seq Analysis

*Figure 3*  

*Figure 3 - Figure Supplement 2*  

[Return to Contents](#con)

------------

<a name="tea"></a>

### TEA-seq

*TODO*  

<a name="tea-req"></a>

#### TEA-seq Requirements

*TODO*  

[Return to Contents](#con)

<a name="tea-dat"></a>

#### TEA-seq Datasets

*TODO*  

[Return to Contents](#con)

<a name="tea-pre"></a>

#### TEA-seq Preprocessing

*TODO*  

[Return to Contents](#con)

<a name="tea-ana"></a>

#### TEA-seq Analysis

*Figure 4*  

*Figure 4 - Figure Supplement 1*  

*Figure 4 - Figure Supplement 2*  

[Return to Contents](#con)

<a name="ref"></a>

#### References

ArchR:  

BedTools:  

GNU Parallel:  

SAMTools:  

Seurat:  

Signac:  

scrattch.vis:  

[Return to Contents](#con)

------------

<a name="leg"></a>

### Legal Information

<a name="leg-lic"></a>

#### License

The license for this package is available on Github in the file [LICENSE in this repository](https://github.com/AllenInstitute/aifi-swanson-teaseq/blob/master/LICENSE)

[Return to Contents](#con)

<a name="leg-lev"></a>

#### Level of Support

We are not currently supporting this code, but simply releasing it to the community AS IS but are not able to provide any guarantees of support. The community is welcome to submit issues, but you should not expect an active response.

[Return to Contents](#con)

<a name="leg-con"></a>

#### Contribution Agreement

If you contribute code to this repository through pull requests or other mechanisms, you are subject to the Allen Institute Contribution Agreement, which is available in the file [CONTRIBUTING in this repository](https://github.com/AllenInstitute/aifi-swanson-teaseq/blob/master/CONTRIBUTING)

[Return to Contents](#con)

