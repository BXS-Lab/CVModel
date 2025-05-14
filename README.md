# Cardiovascular Model

**Version:** 3.2.0<br>
**Release Date:** May 13th, 2025<br>
**Lab:** BXS Lab, UC Davis<br>
**Authors:** R.S. Whittle, A.J. Kondoor, H.S. Vellore<br>
**Contact:**<br>
Dr. Rich Whittle<br>
Department of Mechanical and Aerospace Engineering<br>
University of California, Davis<br>
One Shields Avenue, Davis, CA 95616<br>
📧 rswhittle@ucdavis.edu<br>

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

![Reflex Plot](Images/reflex.png)

*Figure: Arterial Baroreflex and Cardiopulmonary Reflex action on arteriole resistance, heart rate, venous tone, and ventricular contractility in a 90° stand test.*

## Simulation Scenarios

The model supports simulation of:

- Tilt-angle protocols
- Altered-gravity environments
- Lower body negative pressure (LBNP) protocols

![Reflex Plot](Images/hemodynamics.png)

*Figure: Interstitial volume (90° stand test), left heart pressures, cardiac elastances, and ventricular outflow.*

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



