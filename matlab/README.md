# Install

To compile the `fastwht` function open MATLAB and move to this folder, using, e.g., the `cd` command. Type `ls` to list all files in your current directory, and check that it contains the file `fastwht.cpp`. In the MATLAB shell type
```
mex CFLAGS='-O2' fastwht.cpp ../hadamard.cpp
mex CFLAGS='-O2' had_mat_idx.cpp ../hadamard.cpp 
```
to compile the source. Then add this directory to your [Matlab path](https://se.mathworks.com/help/matlab/matlab_env/add-folders-to-matlab-search-path-at-startup.html)

