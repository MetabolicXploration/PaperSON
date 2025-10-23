# PaperSON

## Overview

**PaperSON** (Paper Serialization Object Notation) is a community-driven effort to build a structured, text-based, and version-controlled repository of scientific knowledge extracted directly from published papers.  
The project’s long-term goal is to **create a general RDF knowledge graph** that captures experimental data, metadata, and context from the scientific literature — beginning with metabolism, but extensible to any domain.

Unlike automated text-mining systems or black-box AI extraction, PaperSON emphasizes **human curation** by domain experts (e.g., PhD students), ensuring interpretability and fidelity to original sources.

---

## Motivation and Objectives

Scientific data in the literature are often scattered across figures, supplementary tables, and textual descriptions.  
This project aims to:

1. **Aggregate experimental data** from publications into structured JSON files.
2. **Preserve contextual metadata** (organism, strain, conditions, protocols, etc.) that are essential for data reuse and integration.
3. **Provide a transparent, traceable, and version-controlled infrastructure** to reproduce, inspect, and interlink datasets.
4. **Enable evaluation of computational models** (e.g., culture simulations, flux balance models) against well-documented empirical states.

Ultimately, the repository could function as a **test battery** for model validation and a **living archive** of experimentally verified biological states.

---

## Design Principles

| Principle | Description |
|------------|--------------|
| **Simplicity** | Use only human-readable **text-based formats** (JSON, Markdown). No binary formats or databases. |
| **Transparency** | Every dataset corresponds directly to published data; no inferred or AI-generated entries. |
| **Human-Centric Input** | Data entry is mainly done manually by researchers translating data "as is" from papers. |
| **Version Control** | All data are stored in a Git repository for collaborative curation and provenance tracking. |
| **Reproducibility** | Each paper has its own directory containing raw data, visual reproductions, and contextual README files. |
| **Extensibility** | Designed to start with metabolism, but easily applicable to other biological or physical systems. |

---

## Repository Structure

> NOTE/WIP: This might not be the case, we are still figuring out the structure.

```

PaperSON/
│
├── papers/
│   ├── 2023_baldazzi_resource_allocation/
│   │   ├── data.json
│   │   ├── README.md
│   │   └── figures/
│   ├── 2019_long_flux_responses/
│   │   ├── data.json
│   │   └── README.md
│   └── ...
│
├── schema/
│   ├── paper_schema.json      # JSON structure for data entry
│   ├── metadata_schema.json   # Schema for contextual data
│
└── docs/
├── CONTRIBUTING.md
├── DESIGN.md
└── PaperSON.md            # Design reference document

````

Each paper’s folder includes:
- A **`raw.json`** file: the entry point containing/linking to all stored data.
- A **`README.md`** summarizing key results, reproducing tables/figures, and linking to sources.
- Optional **vectorized figures** or **external URLs** referencing supplementary material.

---

## Example Use Case

Although PaperSON is general-purpose, the current focus is **metabolic network modeling**.  
For example:
- Quantitative growth phenotyping on various carbon sources for *E. coli* transcription factor deletions.
- Data aggregation from studies such as [@baldazziResourceAllocationAccounts2023] or [@longMetabolicFluxResponses2019].
- Cross-study comparison of growth yields, flux states, and environmental conditions.

These serve as initial templates to define data structure and demonstrate the approach for other fields.

---

## Documentation and Standards

> NOTE/WIP: This might not be the case, we are still figuring out the structure.

1. **Data Schema:**  
   Each JSON file follows a minimal, human-readable schema describing:
   ```json
   {
     "paper_id": "2023_baldazzi_resource_allocation",
     "organism": "E. coli K-12 NCM3722",
     "conditions": {
       "medium": "M9 minimal salts",
       "carbon_sources": ["glucose", "acetate"]
     },
     "measurements": {
       "growth_rate": "...",
       "fluxes": { "acetate": "...", "biomass": "..." }
     },
     "references": ["doi:10.1038/..."]
   }
````

2. **Metadata Schema:**
   Contextual information includes experimental setup, genetic background, and any transformations applied to data.

3. **Validation Tools:**
   Simple JSON schema validators ensure structural consistency without complex dependencies.

4. **Visualization:**
   Each paper’s `README.md` may reproduce relevant plots or tables for visual reference.

---

## Current Status and Next Steps

* Several metabolism-related papers have been **registered** or **included**.
* Structure and naming conventions are being stabilized.
* Schema design and validation scripts are under refinement.
* Planned addition: automatic graph export to **RDF/OWL** format for semantic web applications.

---

## Tags and Cross-References

`#Project/PaperSON` · `#Metabolism` · `#KnowledgeGraph` · `#DataCuration` · `#Reproducibility`

---

*PaperSON is an open, transparent, and minimalist framework for transforming the scattered knowledge of science into structured, reusable data.*


## Incuded papers

- alexeevaQuantitativeAssessmentOxygen2002
- alterProteomeRegulationPatterns2021
- baldazziResourceAllocationAccounts2023
- bauerMaximalExponentialGrowth1974
- begIntracellularCrowdingDefines2007
- brLargescale13CfluxAnalysis2011
- covertTranscriptionalRegulationConstraintsbased2002
- folsomPhysiologicalBiomassElemental2015
- folsomPhysiologicalProteomicAnalysis2014
- hermsenGrowthRateComposition2015
- joyStudyGrowthEscherichia2010
- kayserMetabolicFluxAnalysis2005
- lewisOmicDataEvolved2010
- longMetabolicFluxResponses2019
- maserAvoidingAminoAcid2019
- mccloskeyEvolutionGeneKnockout2018
- meeSyntrophicExchangeSynthetic2014
- monkGenomescaleMetabolicNetwork2022
- nanchenNonlinearDependencyIntracellular2006
- okanoRegulationUnderlyingHierarchical2019
- perrenoudImpactGlobalTranscriptional2005
- sennGrowthEscherichiaColi1994
- sezonovEscherichiaColiPhysiology2007
- vanheerdenContinuousBatchCultures2013
- wintermuteEmergentCooperationMicrobial2010
- wooMachineLearningIdentifies2024

## TODOs

- Add documentation 
- formalize a quality testing protocole
    - reproduce simple figures or tables.
    - create automatic test suite.
- Add a CONTRIBUTING.md example and specify how new papers should be added
- Document the minimal schema in detail
- Provide a command-line or notebook-based validation tool example
