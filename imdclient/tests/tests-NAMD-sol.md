# Testing of various NAMD builds alongside imdcleint on ASU's SOL

Please find detailed build and IMD protocol related information for NAMD can be found [here](https://github.com/amruthesht/namd-3.0/blob/IMDv3-staging/IMDv3-dev.md). A brief summary of build instructions for NAMD and test script execution procedures are provided below.

## Building NAMD for the various architectures

The basic procedure involves building `charm` followed by NAMD itself. The following dependencies are common for all architectures:

#### Installing header-only dependencies

Download and install TCL and FFTW libraries:
  (cd to `namd-3.0` if you're not already there)
```
    wget http://www.ks.uiuc.edu/Research/namd/libraries/fftw-linux-x86_64.tar.gz
    tar xzf fftw-linux-x86_64.tar.gz
    mv linux-x86_64 fftw
    wget http://www.ks.uiuc.edu/Research/namd/libraries/tcl8.6.13-linux-x86_64.tar.gz
    wget http://www.ks.uiuc.edu/Research/namd/libraries/tcl8.6.13-linux-x86_64-threaded.tar.gz
    tar xzf tcl8.6.13-linux-x86_64.tar.gz
    tar xzf tcl8.6.13-linux-x86_64-threaded.tar.gz
    mv tcl8.6.13-linux-x86_64 tcl
    mv tcl8.6.13-linux-x86_64-threaded tcl-threaded
```

### Compiling and Building NAMD
#### Single-node multicore version
Build and test the Charm++/Converse library (testing is optional)

When using single-node multicore version:
```
    cd charm-8.0.0
    ./build charm++ multicore-linux-x86_64 --with-production
    cd multicore-linux-x86_64/tests/charm++/megatest
    make
    ./megatest +p4   (multicore does not support multiple nodes)
    cd ../../../../..
```

Set up build directory and compile for NAMD:

multicore version:
```
  ./config Linux-x86_64-g++ --charm-arch multicore-linux-x86_64
```

Finally, `make` the namd source files to create an executable

```
  cd Linux-x86_64-g++
  make  
```