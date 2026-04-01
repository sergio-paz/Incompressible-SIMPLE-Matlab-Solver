# 2D Incompressible Flow Solver (SIMPLE Algorithm)

## Overview
This repository contains a 2D Computational Fluid Dynamics (CFD) solver written in MATLAB. It is designed to solve incompressible flows using the **Finite Volume Method (FVM)** and the **SIMPLE** (Semi-Implicit Method for Pressure Linked Equations) algorithm for pressure-velocity coupling. 

This project was originally inspired by numerical modelling coursework at **Cranfield University**.

## Contact
**Sergio Paz Garcia** MSc Computational Fluid Dynamics Student, Cranfield University  
[Connect on LinkedIn](https://www.linkedin.com/in/sergio-paz-garcia-1872a0243/) | [sergiopazgarcia01@gmail.com]

## Features
* **Language:** MATLAB
* **Algorithm:** SIMPLE 
* **Matrix Solvers:** TDMA
* **Mesh:** Staggered grid arrangement

## Validation Cases
Currently under development. The solver will be validated against classic benchmark problems:
1. Lid-driven cavity flow (Ghia et al., 1982)

## Project Structure
* **/src**: Contains the core numerical engine, including the TDMA solver and SIMPLE algorithm functions.
* **/cases**: Practical implementations and setup files for specific CFD benchmarks (e.g., Lid-Driven Cavity).
* **/results**: Numerical data, validation plots, and comparative studies against established literature.
* **/docs**: Technical documentation and brief theoretical background.

## Acknowledgements and Use of AI
* **Academic Foundation**: This project is based on the principles taught at Cranfield University during the Numerical Modelling of Incompressible Flows module.
* **Generative AI**: I acknowledge the use of Google Gemini as a secondary tool for code optimisation, syntax refinement, and assistance in drafting technical documentation. The core numerical logic, algorithm selection, and physical modelling remain my original work.