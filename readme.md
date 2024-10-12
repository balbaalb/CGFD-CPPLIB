The codes in this repository are results of a self-educational workshop for my learning on 
C++, computational geometry and finite volume method. So, it provides a library for: 
-  Triangulation of a domain using different Delaunay triangulation methods. Various quality measures for such triangulations are provided as well. The triangulations are stored in Quad-Edge structure.
-  Solving conduction-convection, or full Navier-Stokes equations on such triangulations employing finite volume method. Both vertex-based and cell-based schemes are available. 

To see references and explanation please refer to 

https://publications.waset.org/10013433/delaunay-triangulations-efficiency-for-conduction-convection-problems

Albaalbaki, B. and Khayat, R.E., 2023. Delaunay Triangulations Efficiency for Conduction-Convection Problems. International Journal of Mechanical and Mechatronics Engineering, 17(12), pp.448-454.

This repository was used exclusively in the numerical calculations of the above paper.

## To compile and test
    cmake -S . Bbuild
    cd build
    cmake --build .
    ctest -C <necessary config options depending on the C++ compiler>

