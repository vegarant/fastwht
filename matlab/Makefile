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

CC = mex

UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
CPPFLAGS += -Ofast
file_ending = mexa64
endif

ifeq ($(UNAME), Linux)
CPPFLAGS += -O3
file_ending = mexa64
endif

# Change the path of the mex compiler for this to work on most computers

main: fastwht.$(file_ending) wal.$(file_ending) had_mat_idx.$(file_ending) #verification
	
fastwht.$(file_ending): ../hadamard.cpp ../hadamard.h fastwht.cpp  
	$(CC) CFLAGS=$(CPPFLAG)  fastwht.cpp ../hadamard.cpp
wal.$(file_ending): ../hadamard.cpp ../hadamard.h wal_comp_core.cpp wal.cpp  
	$(CC) CFLAGS=$(CPPFLAG)  wal.cpp ../hadamard.cpp
had_mat_idx.$(file_ending): ../hadamard.cpp ../hadamard.h wal_comp_core.cpp had_mat_idx.cpp 
	$(CC) CFLAGS=$(CPPFLAG)  had_mat_idx.cpp ../hadamard.cpp
#verification: verification.cpp multidimIdx.h multidimIdx.cpp Makefile
#	g++ -ggdb -std=c++11 multidimIdx.cpp  verification.cpp -o verification

.PHONY: clean
clean:
	rm -f fastwht.$(file_ending) had_mat_idx.$(file_ending) wal.$(file_ending)
