#Help for building Triplexator on Mac OSX using the gcc 4.5 compiler

# Introduction #

Mac OSX uses Clang as default compiler. To enable/use the gnu compiler it is helpful to install at least gcc 4.3.

# Details #

To install gcc 4.5 on a mac using mac port, type:

```
sudo port install gcc45 gcc_select
```

To activate the version:

```
sudo port select gcc mp-gcc45
```

You may also need to rehash gcc when you use bash
```
hash gcc
hash g++
```

Finally set the CC compiler for before running CMAKE
```
export CC=/usr/bin/gcc
```

For BOOST support use mac ports as well
```
sudo port install boost
```

and cmake if you don't have it yet
```
sudo port install cmake
```

From here the usual path should work:

obtain triplexator:
```
git clone https://code.google.com/p/triplexator/ triplexator
```
change directory:
```
cd triplexator
```
create directory and change into it:
```
mkdir -p build/Release && cd build/Release
```
run cmake and make:
```
cmake ../.. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++ -G "Unix Makefiles" && make
```
change directory:
```
cd ../..
```
try the binary:
```
./bin/triplexator --help
```
run the smoketest:
```
./demos/smoketest_triplexator.sh ./bin/triplexator
```

To install triplexator into /usr/local/bin call:
```
cd build/Release
make install
```