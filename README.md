# Convergent-Evolution-in-Cancer
A computational model simulating tumour evolution to compare the Oncogene Addiction and Convergent Evolution paradigms, demonstrating how genetic diversity is maintained under selection pressure.

## Overview

This repository contains the code and documentation for a Master's research project that investigates the evolutionary dynamics of cancer. Traditional cancer treatment often fails due to the emergence of drug resistance. The "oncogene addiction" model suggests targeting a single driver gene is sufficient, while the "convergent evolution" model proposes that tumors can achieve resistance through multiple, diverse genetic pathways.

This project uses a **discrete-time Moran process** to simulate a population of tumor cells evolving under treatment pressure. By tracking **genotypic** and **phenotypic diversity** using **Shannon entropy**, the model demonstrates distinct evolutionary patterns that support the convergent evolution framework.

## Key Features

-   **Stochastic Simulation:** Implements a Moran process for realistic birth-death evolutionary dynamics.
-   **Flexible Genotype Representation:** Simulates both small (3-bit) and large (12-bit) genomes.
-   **Fitness Landscapes:** Models both additive fitness (Oncogene Addiction) and peak fitness (Convergent Evolution) scenarios.
-   **Diversity Metrics:** Quantifies genotypic and phenotypic diversity using Shannon entropy.
-   **Comparative Analysis:** Directly compares the predictions of two competing cancer evolution paradigms.
-   **Visualization:** Generates plots for diversity trends over time and statistical summaries.

## Usage

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/your-username/BioCancer-Convergent-Evolution.git
    cd BioCancer-Convergent-Evolution
    ```

2.  **Run the main simulation script from the `scripts` directory:**
    ```bash
    cd scripts
    python Convergent_Evolution.py
    ```

**The script will:**
- Run simulations for the Oncogene Addiction model (3-bit and 12-bit genomes).
- Run a simulation for the Convergent Evolution model (12-bit genome with peak fitness).
- Plot the results comparing genotypic and phenotypic diversity over time.
- Print a statistical summary including Pearson correlation coefficients between diversity metrics.

## Key Results

**The simulations demonstrate:**
-   **Oncogene Addiction Model:** Genotypic (Hg) and phenotypic (Hp) diversity are strongly coupled and decline together under selection.
-   **Convergent Evolution Model:** Genotypic diversity remains high while phenotypic diversity declines, showing a decoupling effect. This indicates multiple genotypes are converging on the same high-fitness phenotype.

## Implications

This work provides computational evidence that tumors can exploit numerous genetic pathways to survive treatment, challenging the "one gene, one drug" paradigm. The findings support the need for evolution-guided, adaptive combination therapies that target multiple pathways simultaneously to prevent treatment resistance.

## Documentation

-   `docs/Master_Project_Thesis.pdf`: The complete thesis paper detailing the background, methodology, results, and discussion.
-   `docs/Research_Poster.pdf`: A conference-style poster summarizing the project.

## References

Key literature influencing this model includes work by Moran (1958), Shannon (1948), Chakravarty (2015), and McGranahan & Swanton (2017). Full citations are available in the thesis document.

