# Degridder - CPU Implementation

**Todo: description about project, what degridding is, and where it fits in the imaging pipeline.**

---
##### Instructions for installation of this software (includes profiling, linting, building, and unit testing):
1. Install [Valgrind](http://valgrind.org/) (profiling, memory checks, memory leaks etc.)
   ```bash
   $ sudo apt install valgrind
   ```
2. Install [Cmake](https://cmake.org/)/[Makefile](https://www.gnu.org/software/make/) (build tools)
   ```bash
   $ sudo apt install cmake
   ```
3. Install [Google Test](https://github.com/google/googletest) (unit testing) - See [this tutorial](https://www.eriksmistad.no/getting-started-with-google-test-on-ubuntu/) for tutorial on using Google Test library
   ```bash
   $ sudo apt install libgtest-dev
   $ cd /usr/src/gtest
   $ sudo cmake CMakeLists.txt
   $ sudo make
   $ sudo cp *.a /usr/lib
   ```
4. Install [Cppcheck](http://cppcheck.sourceforge.net/) (linting)
   ```bash
   $ sudo apt install cppcheck
   ```
5. Configure the code for usage (**modify degridder.c config**)
6. Build the degridder project (from project folder)
   ```bash
   $ mkdir build && cd build
   $ cmake .. -DCMAKE_BUILD_TYPE=Debug && make
   ```
---
##### Instructions for usage of this software (includes executing, testing, linting, and profiling):
To perform memory checking, memory leak analysis, and profiling using [Valgrind](http://valgrind.org/docs/manual/quick-start.html), execute the following (assumes you are in the appropriate *build* folder (see step 5 above):
```bash
$ valgrind --leak-check=yes ./degridder
$ valgrind --leak-check=yes ./tests
```
To execute linting, execute the following commands (assumes you are in the appropriate source code folder):
```bash
$ cppcheck --enable=all main.cpp
$ cppcheck --enable=all degridder.c
$ cppcheck --enable=all unit_testing.cpp
```
To execute unit testing, execute the following (also assumes appropriate *build* folder):
```bash
$ ./tests
````
To execute the direct fourier transform (once configured and built), execute the following command (also assumes appropriate *build* folder):
```bash
$ ./degridder
```
