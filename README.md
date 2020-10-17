#  Order2ForIrregularDomains

Depuis "A  Second-Order-Accurate  Symmetric  Discretization  of  the  Poisson  Equation  on  Irregular  Domains"  Frederic  Gibou.

Ce code de projet est composé de plusieurs parties sources :

 - La bibliothèque partagée O2fid
 - Et des exécutables d'exemples : 1 à 10 (10 qui ne fonctionne pas vraiment), et 12.

## Pour compiler le projet

À la racine du répertoire,  vous devez faire  un `mkdir build` suivi d'un `cd build`, et dans ce répertoire un    `cmake ../`	ou `cmake-gui ../`. La différence entre les deux se situe dans l'interface graphique et donc la facilité avec laquelle la création de chacun des éxecutables peut être effectuée (la compilation de chacun d'entre eux est soumise à un booléen).
Une fois que cela s'est bien déroulé, il vous suffira de faire un `make` et le tour est joué !

## La structure du répertoire

 ```
.
├── README.md
|
├── build
|
├── exec
|
├── docs
|   ├── Doxyfile.in
│   └── doc_doxygen
|
├── libs
│   ├── EIGEN
│   └── lib.cmake
|
└── src
    ├── Example_01
    │   ├── CMakeLists.txt
    │   └── main.cpp
    |
    ├── Example_02
    │   ├── CMakeLists.txt
    │   └── main.cpp
    |
    ├── Example_03
    │   ├── CMakeLists.txt
    │   └── main.cpp
    |
    ├── Example_04
    │   ├── CMakeLists.txt
    │   └── main.cpp
    |
    ├── Example_05
    │   ├── CMakeLists.txt
    │   └── main.cpp
    |
    ├── Example_06
    │   ├── CMakeLists.txt
    │   └── main.cpp
    |
    ├── Example_07
    │   ├── builder.cpp
    │   ├── CMakeLists.txt
    │   ├── headers.h
    │   └── main.cpp
    |
    ├── Example_08
    │   ├── builder.cpp
    │   ├── CMakeLists.txt
    │   ├── headers.h
    │   └── main.cpp
    |
    ├── Example_09
    │   ├── builder.cpp
    │   ├── CMakeLists.txt
    │   ├── headers.h
    │   └── main.cpp
    |
    ├── Example_10
    │   ├── builder.cpp
    │   ├── CMakeLists.txt
    │   ├── headers.h
    │   └── main.cpp
    |
    ├── Example_12
    │   ├── CMakeLists.txt
    │   └── main.cpp
    |
    ├── Example_TVD_WENO
    │   ├── CMakeLists.txt
    │   └── main.cpp
    |
    └── O2FID
        ├── CMakeLists.txt
        |
        ├── O2FID.h
        |
        ├── Data
        │   ├── cell.cpp
        │   ├── cell.h
        │   ├── data.h
        │   ├── datatypedefinitions.h
        │   ├── point.cpp
        │   └── point.h
        |
        ├── Mesh
        │   ├── mesh.cpp
        │   └── mesh.h
        |
        ├── Outputs
        │   ├── writer.cpp
        │   └── writer.h
        |
        ├── Stefan
        │   ├── field.cpp
        │   ├── field.h
        │   ├── interface.cpp
        │   ├── interface.h
        │   ├── stefan.h
        │   ├── wfield.cpp
        │   └── wfield.h
        |
        ├── Toolbox
        │   ├── differencefinite.cpp
        │   ├── differencefinite.h
        │   ├── toolbox.cpp
        │   └── toolbox.h
        |
        └── Tools
            ├── border.cpp
            ├── border.h
            ├── errors.cpp
            ├── errors.h
            ├── funtovec.cpp
            ├── funtovec.h
            ├── impose.cpp
            ├── impose.h
            ├── matrix.cpp
            ├── matrix.h
            ├── solver.cpp
            ├── solver.h
            └── tools.h
```



