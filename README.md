# ROSMAP_MOFA_public_release

Public release of code and analyses related to **Multi-omics Factor Analysis (MOFA)** using data from the **Religious Orders Study (ROS)** and **Rush Memory and Aging Project (MAP)** cohorts.  

This repository accompanies the preprint:

## Citation

> Scheidemantel LP, de Paiva Lopes K, Gaiteri C, Menon V, De Jager PL, Schneider JA, Buchman AS, Wang Y, Tasaki S, Raittz RT, Bennett DA, Vialle RA.\
 **Integration of aged brain multi-omics reveals cross-system mechanisms underlying Alzheimer’s disease heterogeneity.**\
 bioRxiv (2025). doi: [10.1101/2025.09.23.678110](https://www.biorxiv.org/content/10.1101/2025.09.23.678110v1)

---

## Overview

<p align="center">
  <img src="MOFA_Graphical_Abstract.png" alt="Image" width="800"/>
</p>
<p align="right">
*Figure created in  https://BioRender.com*
</p>

Key findings include:

- Identification of **cross-omic, cross-system biological factors** linked to AD phenotypes.
- Discovery of **immune-related factors** with both detrimental and protective associations.
- **Unsupervised clustering** of participants into **11 molecular subtypes**, including three AD-associated clusters with distinct molecular and phenotypic signatures.
- Insights into **neuroinflammatory processes**, **energetic metabolism**, and **cytoskeletal dynamics** as central mechanisms in AD.

---

## Repository Structure

- **`1.MOFA_analysis/`** – Core MOFA implementation and factor extraction.  
- **`2.GSEA_analysis/`** – Gene Set Enrichment Analysis of MOFA factors.  
- **`3.Loadings_analysis/`** – Examination of omics loadings contributing to factors.  
- **`4.Phenotypes_analysis/`** – Associations between MOFA factors and AD phenotypes.  
- **`5.Clustering_analysis/`** – Unsupervised clustering of participants into molecular subtypes.  
- **`6.Manuscript_code/`** – Misc. scripts used to generate other figures and results.  

---


---

Distributed under terms of the [GNU GENERAL PUBLIC LICENSE](/LICENSE).
