import os
import multiprocessing

#Decider('timestamp-match');
#Decider('MD5-timestamp');

SetOption('num_jobs', multiprocessing.cpu_count())
SetOption('implicit_cache', 1);

#lib_files = ["../hadam/hadamard.cpp"];
#library_name = "Hadamard";

src_files = ["hadamard.cpp"];
program1 = "main";
program2 = "verification";

#'-pedantic', '-Wall', '-Wextra', '-Wno-unused-parameter'
#env = Environment(LIBS=library_name, LIBPATH='.');
env = Environment();
env.Append(CXXFLAGS = ['-std=c++11', '-ggdb', '-fmax-errors=5']); #])

#if int(ARGUMENTS.get('release', 0)) == 1:
#    env.Append(CXXFLAGS = ['-O3'])
#else:
#    env.Append(CXXFLAGS = ['-ggdb', '-O0'])

#env.Library(target = library_name, source=lib_files);
env.Program(target = program1, source = [program1+".cpp"] + src_files);
env.Program(target = program2, source = [program2+".cpp"] + src_files);


# vim: set ft=python:

