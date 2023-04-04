# Finite difference solvers for PDEs

PDE solvers in 2D for personal research usage, including wave equation, acoustic wave equation, Maxwell's equations, and convection-diffusion equation.
Features:
- Use Eigen3 for dynamic vector and matrix container. 
- Written in C++ with object-oriented programming style.
- OpenMp is used for parallel computing to increase computational efficiency.

Please compile it with openmp: 
```
g++ --std=c++17 -fopenmp .\src\*.cpp -o Solver.exe
```

Next:
- Improve computing efficiency by improving OpenMP code.
- Dedicated vector and matrix container.
- Build Python interfaces for PDE-constrained optimization work.
- Parallel with Cuda (maybe for another project).

## Demo:
### Convection-diffusion equation
![](data\heat.png)
### Maxwell's equations
![](data\wave.png)

Check the Gifs in the data folder for fun.

