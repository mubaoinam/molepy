Molepy
======
Molepy is a python based *ab initio* quantum chemistry program. Currently, the program is in an experimental state and so, it can perform calculations only with Restricted Hartree-Fock method using limited basis sets (cc-pVDZ,cc-pVTZ, def2-TZVP, etc.). The performance is not as efficient as other quantum chemistry programs but will be improved over time (hopefully). The main advantage is that it is a python based program and so it is very easy to understand or modify it. Most of the computationally demanding modules are written in C++ but a python counterpart can be found for most of it. 

## Installation
You might want to clone with --recursive option. Molepy contains submodule molepy/ccmints/src/pybind11

* Requirements
    - CMake 3.0 or higher
    - GCC 4.8 or newer
    - Python 2.7
    - Numpy 1.15.1
* Compile C++ modules
         
         cd molepy/ccmints
         mkdir build
         cd build
         cmake ../src
         make
* Include the top level directory in environmental variable "PYTHONPATH".

         export PYTHONPATH=/path/to/pymole:PYTHONPATH

## Usage
Molepy is quite flexible to use, see molepy/example. Also, read molepy/basis/readme.txt for basis sets. If you choose to use pure python codes only, set 'hfparam.ccmints = False' in molepy/lib/mole_param.py.

   
