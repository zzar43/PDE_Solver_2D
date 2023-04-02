Please compile with openmp: 
g++ --std=c++17 -fopenmp .\src\*.cpp -o WaveEq2D.exe

Task:
- udpate the acoustic wave equation with PML.
- update the Maxwell's equations.

Next version:
- udpate the Model container, such that all equaiton can use the same model.
- put the save function into the model container.