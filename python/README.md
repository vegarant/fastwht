# Fast Walsh Hadamard transform

This doirectory contains a script that will create python 2 wrapper code to
interface the C++ code in the above directory. The script require that SWIG is
installed. 

## Installation

1. To create the python wrapper function run the `setup.sh` file. This file use
the command *python2* to run python 2. If your system run python 2 using the
*python* command, you will need to change one line in the `setup.sh` file.

2. Run the file `test.py` to verify your implementation. 

3. In order to make the python code callable from other directories add the
directory to your python path.
```export PYTHONPATH=$PYTHONPATH:/path/to/directory/code/hadam/python```


