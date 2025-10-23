Here is a README in Markdown format summarizing our discussion about the `pnas.1405641111.sd03` file:

```markdown
# README for pnas.1405641111.sd03

This document provides clarification and context for the data contained in the `pnas.1405641111.sd03` file, which is part of the supporting information for the paper "Syntrophic exchange in synthetic microbial communities" by Mee et al..

## 1. File Purpose and Content

This file contains **fold growth values for 364 three-member synthetic consortia** of *Escherichia coli*. Crucially, it also includes the **fold growth values for the corresponding two-member subsets** that constitute these three-member groups, derived from pairwise interaction experiments.

The table's columns are:
*   `Co-culture ID`: A unique identifier for each consortium.
*   `Strain 1`, `Strain 2`, `Strain 3`: These are single letters representing the three unique amino acids that define the consortium.
*   `"Strain 1+2+3 fold growth (T84)"`: The total observed fold growth of the three-member consortium after 84 hours of coculture.
*   `"Strain 1+2 fold growth (T84)"`: The observed fold growth of the two-member subset defined by "Strain 1" and "Strain 2" after 84 hours.
*   `"Strain 1+3 fold growth (T84)"`: The observed fold growth of the two-member subset defined by "Strain 1" and "Strain 3" after 84 hours.
*   `"Strain 2+3 fold growth (T84)"`: The observed fold growth of the two-member subset defined by "Strain 2" and "Strain 3" after 84 hours.

## 2. Interpretation of Strain IDs

The single letters (e.g., R, C, G) under `Strain 1`, `Strain 2`, and `Strain 3` in the table refer to **amino acids**, not individual strains directly [conversation history].

*   These letters are shorthand for the 14 single-amino acid auxotrophies initially generated (e.g., 'R' for Arginine ΔargA, 'C' for Cysteine ΔcysE, 'G' for Glycine ΔglyA).
*   For the **three-member consortia** represented by this table, each individual strain in the consortium is actually a **double-amino acid auxotroph**.
*   Therefore, if a row lists `R`, `C`, and `G` as `Strain 1`, `Strain 2`, and `Strain 3` respectively, it signifies a consortium composed of **three double-auxotroph strains**:
    *   A strain requiring **R and C** (RC auxotroph).
    *   A strain requiring **R and G** (RG auxotroph).
    *   A strain requiring **C and G** (CG auxotroph) [conversation history].
*   This design ensures a "strictly syntrophic" community where each member relies on the other two partners for its required amino acids.

## 3. Definition of "Fold Growth (T84)"

"Fold growth" is a measure of the **growth yield of a community** over a specified period.
*   It is calculated as the **final cell density divided by the initial cell density** of the community.
*   `T84` indicates that the growth was measured after **84 hours** of coculture.

## 4. Experimental Context and Data Acquisition

The data in this file was generated through a series of experiments:

### A. Two-Member Pairwise Syntrophic Interactions (for "Strain 1+2 fold growth", etc.)
*   **Strains**: Composed of **14 single-amino acid auxotrophs**.
*   **Measurement**: All 91 possible pairwise cocultures were probed in M9-glucose minimal media. Total coculture fold growth and relative abundance of each member were measured after 84 hours.
*   **Significance**: These experiments provided the "cooperativity coefficients" (c) used in modeling and the baseline "pair based fold growth" values for comparison with three-member consortia.

### B. Three-Member Synthetic Consortia (for "Strain 1+2+3 fold growth")
*   **Strains**: Composed of **91 double-amino acid auxotrophic derivatives**. Each strain was auxotrophic for two amino acids (e.g., MF, MK, FK strains for an M, F, K consortium).
*   **Measurement**: All 364 possible three-member consortia were systematically measured for their total community fold growth after 84 hours.
*   **Control**: Individual double-auxotroph strains showed no growth when only one of their two required amino acids was supplemented. Also, "growth was not observed when only two of the three members are cocultured" (referring to double auxotrophs in a three-member configuration). This confirms the strictly syntrophic nature.
*   **Comparison**: The observed growth of these three-member consortia was directly compared with the fold growth of their corresponding two-member (monoauxotrophic) subsets to assess epistasis.

## 5. Initial Culture Conditions

The initial conditions for these coculture experiments were carefully controlled:

*   **Individual Strain Preparation**: Before coculture, all individual strain concentrations (from pre-cultures) were adjusted to **10^7 cells per mL** in M9 media.
*   **Total Community Seeding Density**: Cocultures were initiated such that the **total initial cell density of the *entire community* was 10^7 cells per mL** [32, 48, conversation history].
*   **Initial Density of Each Strain**: In an experimental community with `N` strains, the initial density of *each individual strain* was therefore **10^7 cells per mL / N** [conversation history].
    *   For example, in a three-member coculture, each of the three strains started at approximately 3.33 x 10^6 cells/mL.
*   **Inoculation Volume**:
    *   For two-member and three-member cocultures, the total volume was **200 μL**.
    *   The strains were inoculated by adding an equal volume from each 10^7 cells/mL pre-culture stock. For a 200 μL three-member culture, approximately 66.7 μL from each strain's 10^7 cells/mL stock was added [conversation history].

## 6. Information Not Explicitly Reported in this Table

*   **Individual Strain Abundance in Three-Member Consortia**: While the authors determined relative abundance for pairwise cocultures and for the 13- and 14-member communities, this specific table (`pnas.1405641111.sd03`) for the 364 three-member consortia **does not explicitly report the individual abundance or relative proportion of each double-auxotroph strain *within* these three-member cultures** [conversation history]. The focus here is on total community yield and its comparison to constituent pairs.
*   **Single-Strain Growth**: Single auxotrophic strains showed **no growth** in M9-glucose minimal media without amino acid supplementation, so no meaningful "single strain fold growth" can be recovered from this data.
```