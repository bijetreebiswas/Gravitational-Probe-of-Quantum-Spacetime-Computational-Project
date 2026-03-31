# Gravitational Probe of Quantum Spacetime – Computational Projects

This repository contains Python code that implements three computational projects based on the paper:

> N. Herceg, T. Juric, A. Samsarov, I. Smolic, K. S. Gupta  
> *Gravitational probe of quantum spacetime*  
> [arXiv:xxxx.xxxxx] / Physics Letters B

The code computes noncommutative (quantum) corrections to the Regge‑Wheeler potential for a Schwarzschild black hole and analyses quasinormal mode (QNM) frequencies using the WKB approximation.

## Table of Contents

- [Overview](#overview)
- [Projects Included](#projects-included)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Project Details](#project-details)
- [Output Examples](#output-examples)
- [Customisation](#customisation)
- [References](#references)
- [License](#license)

## Overview

In the paper, the authors introduce a noncommutative (NC) deformation of spacetime via a Drinfeld twist and derive a modified Regge‑Wheeler potential that depends on the azimuthal number \(m\) (a “Zeeman‑like” effect). They then compute QNM frequencies using the third‑order WKB method and show that the Schwarzschild black hole remains stable under axial perturbations.

This repository translates the theoretical results into three self‑contained computational projects:

1. **Plotting the quantum Regge‑Wheeler potential** – visualisation of the NC correction.
2. **Time evolution of a wave packet** – a simple finite‑difference simulation of gravitational wave scattering.
3. **QNM frequencies via 3rd‑order WKB** – numerical calculation of the complex quasinormal mode frequencies.

All code is written in Python and uses only standard scientific libraries.

## Projects Included

### 1. Plotting the Potential (`plot_potential`)
- Plots \(V(r)\) and \(V(r_s)\) for a range of \(qm\) values.
- Reproduces the Zeeman‑like splitting shown in Figures 1 and 2 of the paper.

### 2. Time Evolution (`time_evolution`)
- Solves the 1+1 dimensional wave equation \(\partial_t^2 \psi = \partial_{r_s}^2 \psi - V(r_s) \psi\).
- Uses a leapfrog finite‑difference scheme with absorbing boundaries.
- Visualises the ringdown signal at the peak of the potential.

### 3. QNM Frequencies (`compute_qnm_table`)
- Locates the maximum of the potential \(V(r_s)\).
- Computes numerical derivatives (2nd, 3rd, 4th) at the peak.
- Applies the **3rd‑order WKB formula** (Iyer & Will 1987) to obtain \(\omega\).
- Prints a table of \(\operatorname{Re}(\omega)\) and \(-\operatorname{Im}(\omega)\) for different \(qm\).

## Requirements

- Python 3.7 or later
- Required libraries:
  - `numpy`
  - `scipy`
  - `matplotlib`


