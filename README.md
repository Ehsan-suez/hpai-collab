# hpai-collab
HPAI Modelling Challenge Repository

This repository contains modelling code, analyses, and outputs for the HPAI Modelling Challenge.

The objective of this project is to simulate, forecast, and analyze the spread of highly pathogenic avian influenza (HPAI) at the farm level using the datasets provided by the challenge organizers.

This repository will include submissions and analyses for:

• Phase 1 (current)
• Phase 2 (to be added)
• Phase 3 (to be added)

Modelling framework

The epidemic simulations in this repository are based on a stochastic, spatially explicit, farm-level SEIDR (Susceptible–Exposed–Infectious–Detected–Removed) transmission model.

The core transmission model framework was originally developed by Sarita Bugalia (bugalias@arizona.edu) as part of the HPAI Modelling Challenge.

The model represents disease transmission between farms through multiple pathways:

• Distance-dependent spatial transmission
• Animal movement–mediated transmission
• Environmental transmission associated with proximity to wetlands and water bodies
• Increased transmission risk within High Risk Zones (HRZ)
• Detection and removal processes
• Preventive culling events
• Farm-specific activity windows

This repository contains an implementation and extension of the model by:

Ehsan Suez

Additional components implemented in this repository include:

• Parameter inference using particle filtering (pomp mif2)
• Ensemble stochastic simulation
• Forecast generation and uncertainty quantification
• Spatial risk estimation
• Transmission attribution analysis
• Counterfactual intervention simulations

Phase 1

Phase 1 focuses on reconstructing and forecasting the epidemic using observed data through the end of Phase 1.

The following objectives were addressed:

Q1: Descriptive epidemiology
Q2: Forecast temporal and spatial epidemic evolution over the next 4 weeks
Q3: Estimate relative contribution of galliform vs palmiped farms

The Phase 1 implementation includes:

• Model calibration using observed case data
• Ensemble forecasting
• Estimation of spatial risk
• Transmission attribution analysis
• Knockout simulations to evaluate group-specific contributions

Contact

Ehsan Suez
PhD Candidate
Institute of Bioinformatics
University of Georgia

Email: ehsan.suez@uga.edu
