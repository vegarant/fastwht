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


import os
import multiprocessing

cc = "g++";
cppFlags = ["-std=c++11", '-O3'];

SetOption('num_jobs', multiprocessing.cpu_count())
SetOption('implicit_cache', 1);

src_files = ["hadamard.cpp", "timer.cpp"];
program1 = "main";
program2 = "verification";

if (cc == "g++"):
    cppFlags.append("-fmax-errors=5");
elif (cc == "clang++"):
    cppFlags.append("-ferror-limit=5");
    


env = Environment(CXX = cc);
#env.Append(CPPFLAGS = ['-std=c++11', '-O3', '-fmax-errors=5']); # gcc flags
env.Append(CPPFLAGS = cppFlags); # clang flags
#env.Append(CPPPATH="fxt/");
#env.Append(LINKFLAGS=["-pthread"]);

env.Program(target = program1, source = [program1+".cpp"] + src_files);
env.Program(target = program2, source = [program2+".cpp"] + src_files); 

# vim: set ft=python:











