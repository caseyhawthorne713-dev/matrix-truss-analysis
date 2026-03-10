# matrix-truss-analysis
![MATLAB](https://img.shields.io/badge/MATLAB-Code-orange)
![Course](https://img.shields.io/badge/COE-321K-black)

## Description:
Below are descriptions of each file within this project folder
### [main.m](main.m)
- Program used to model a truss and gather data on displacement, element stress, and reaction forces, as well as to display the deformation due to applied loads.
### [writeup.pdf](writeup.pdf)
- Short writeup detailing the process and conceptual background of the program and its goals. Includes, diagrams, code snippets, and equations.

## Installations
There are no required dependencies for this program

## Usage
The program is ran via the MATLAB (software or online variant) command terminal

- Truss is built by editing out the following matrices
  - Coords → Node locations (x,y)
  - Bounds → Degrees of freedome (1 = fixed, 0 = free)
  - F_Applied → Applied forces at node i (where i = node row)
  - Connectivity → Places element members connecting each node
  - U_prescribed → Initial displacement of each node

### main.m
Expected outputs will be a matrix of the displacement due to applied loads, vectors of the reaction forces at supports and the stresses in each truss member, and a graph comparing the undeformed truss to the deformation caused by the applied loads
