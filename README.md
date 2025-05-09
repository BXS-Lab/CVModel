# Cardiovascular Model

**Version:** 3.1 
**Release Date:** May 9th, 2025  
**Lab:** BXS Lab, UC Davis  
**Authors:** R.S. Whittle, A.J. Kondoor, H.S. Vellore  
**Contact:**  
Dr. Rich Whittle  
Department of Mechanical and Aerospace Engineering  
University of California, Davis  
One Shields Avenue, Davis, CA 95616  
📧 rswhittle@ucdavis.edu

---

## Overview

This repository contains a simulation of the human cardiovascular system developed using [ModelingToolkit.jl](https://mtk.sciml.ai/stable/) in Julia. The model includes:

- A four-chamber heart
- Arterial and venous systems
- Microcirculation
- Arterial baroreflex (ABR) and cardiopulmonary reflex (CPR)
- Hydrostatic effects and interstitial fluid dynamics
- External tissue pressure modeling
- Pulmonary circulation

## Simulation Scenarios

The model supports simulation of:

- Tilt-angle protocols
- Altered-gravity environments
- Lower body negative pressure (LBNP) protocols

## Model Basis

The mathematical formulation is based on prior work by:

- Heldt (2004)  
- Zamanian (2007)
- Mynard (2012)  
- Diaz Artiles (2015)
- Albanese (2016)  
- Whittle (2023)

---

## License

This project is licensed under the [MIT License](LICENSE).

## DOI

[![DOI](https://zenodo.org/badge/894082810.svg)](https://doi.org/10.5281/zenodo.15338311)

