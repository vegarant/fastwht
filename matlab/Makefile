# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2016 Vegard Antun
#


# Note that on many linux distributions the `mex` command will give you a
# version of `pdflatex`. If this is the case you find the corrcet compiler at 
# matlabroot/bin/mex
CC = mex 
#CC = /opt/app-sync/matlab/bin/mex

UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
CPPFLAGS = '-Ofast -fPIC'
LDFLAGS=""
file_ending = mexmaci64
endif

ifeq ($(UNAME), Linux)
CPPFLAGS = '-O3 -fPIC'
#LDFLAGS=''

# Compile with openMP instead
#CPPFLAGS = '-O3 -fPIC -fopenmp'
#LDFLAGS='-fopenmp'
file_ending = mexa64
endif

# Change the path of the mex compiler for this to work on most computers

main: fastwht.$(file_ending) wal.$(file_ending) had_mat_idx.$(file_ending) #verification
	
fastwht.$(file_ending): ../hadamard.cpp ../hadamard.h fastwht.cpp  
	$(CC) CFLAGS=$(CPPFLAGS) LDFLAGS=$(LDFLAGS)  fastwht.cpp ../hadamard.cpp
had_mat_idx.$(file_ending): ../hadamard.cpp ../hadamard.h wal_comp_core.cpp had_mat_idx.cpp 
	$(CC) CFLAGS=$(CPPFLAGS)  had_mat_idx.cpp ../hadamard.cpp

.PHONY: clean
clean:
	rm -f fastwht.$(file_ending) had_mat_idx.$(file_ending) 
