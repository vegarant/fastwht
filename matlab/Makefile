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

CC = /opt/matlab/opt/app-sync/matlab/bin/mex
CFLAGS = -ggdb 

OBJS   = hadamard.o 
HFILES = ../hadamard.h 


all: fastwht.mexa64 

fastwht.mexa64: $(HFILES) $(OBJS) fastwht.cpp  Makefile
	$(CC)   $(OBJS) fastwht.cpp 
hadamard.o : ../hadamard.cpp  
	$(CC) CFLAGS='-ggdb' -c ../hadamard.cpp 

