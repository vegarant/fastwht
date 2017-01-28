# Install

To compile the fastwht function, open Matlab and move to this folder, using
i.g., the `cd` command. Type `ls` to list all files in your current directory,
and check that it contains the file 'fastwht.cpp'. In the matlab shell type
```
mex CFLAGS='-O3' fastwht.cpp ../hadamard.cpp
```
to  compile the source. Then add this directory to your 
[Matlab path](https://se.mathworks.com/help/matlab/matlab_env/add-folders-to-matlab-search-path-at-startup.html)

