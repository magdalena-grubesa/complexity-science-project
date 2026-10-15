# Multi-node alleviating interventions in a depression–anxiety network

This project simulates network-based interventions using NIRA algorithm (Lunansky et al., 2022) in an Ising model of
depression and anxiety symptoms, based on the network reported by Li et al. (2024).

## Research question
Is it more effective to target one symptom or two symptoms with the same total
treatment dose, and how does the effect change when the dose is spread over more symptoms?

## Model
- **Network:** 16 symptoms. Edge weights approximate the strong,
  medium and bridge connections reported by Li et al. (2024); all other pairs have a weak link (0.1)
- **Thresholds:** not reported by Li et al., so symptoms are grouped by how common they were
  in their sample (Table 1, item means). More common symptoms are easier to activate
- **Outcome:** mean number of active symptoms (sum score) of all possible symptom patterns
- **Intervention:** alleviating: thresholds of targeted symptoms are lowered by a
  number of standard deviations (SD) of the thresholds

## Required packages
```r
install.packages(c("IsingSampler", "qgraph", "ggplot2"))
```

## How to run
1. Download `simulation.R` from this repository.
2. Open it in R or RStudio.
3. Set the working directory to the folder containing the script
4. Run Steps 1–2 first 
5. Run the remaining sections in order


## References:
Li, P., Huang, N., Yang, X., Fang, Y., & Chen, Z. (2024). A simulation-based network analysis of intervention targets for adolescent depressive and anxiety symptoms. Asian Journal of Psychiatry, 99, 104152. https://doi.org/10.1016/j.ajp.2024.104152 

Lunansky, G., Naberman, J., van Borkulo, C. D., Chen, C., Wang, L., & Borsboom, D. (2022). Intervening on psychopathology networks: Evaluating intervention targets through simulations. Methods, 204, 29–37. https://doi.org/10.1016/j.ymeth.2021.11.006 
