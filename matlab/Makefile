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




# Change the path of the mex compiler for this to work on most computers

CC = mex
CPPFLAG='-O3'
main: fastwht.mexa64 wal.mexa64 verification
	
fastwht.mexa64: ../hadamard.cpp ../hadamard.h fastwht.cpp  
	$(CC) CFLAGS=$(CPPFLAG)  fastwht.cpp ../hadamard.cpp
wal.mexa64: ../hadamard.cpp ../hadamard.h wal.cpp  
	$(CC) CFLAGS=$(CPPFLAG)  wal.cpp ../hadamard.cpp
verification: verification.cpp multidimIdx.h multidimIdx.cpp Makefile
	g++ -ggdb -std=c++11 multidimIdx.cpp  verification.cpp -o verification

.PHONY: clean
clean:
	rm -f fastwht.mexa64 fastwhtn.mexa64
