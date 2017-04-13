# FFITools

The FFITools are a set of scripts to process calibrated, full-frame images
into light curves.  FFITools covers the photmetry and trending portions of
the TESS pipeline.

## Pre-Requisites (software versions are what used for testing in development)

### Development environment

-   gcc compiler etc. (standard on Linux, need Xcode developer tools
    on MacOS X)

-   gfortran backend (gcc-gfortran subpackage on most Linux distros,
    need separate package on MacOS X:
    https://gcc.gnu.org/wiki/GFortranBinariesMacOS)

-   pkg-config (not a standard part of Xcode developer tools on MacOS
    X, need to install MacPorts: macports.org and then run ```sudo
    port install pkgconfig```)

### From tessgit
    tsig
    ###########

### From PyPI
    numpy==1.11.1 
    scipy==0.18.1
    multiprocessing==0.70a1 
    matplotlib==1.5.3
    pyfits==3.4
    lmfit 
    latex (textlive on centos, dvipng is needed) 
    astroquery==0.3.2
    astropy==1.1.2
    ###########

### must be compiled from source (for the time being, more details below)
    fitsh 0.9.2
    astrometry.net 0.6.6
    vartools 1.33
    

###  potential future dependency
    lfov-gaia.sh (developed by Andras)
    lfov-usno.sh (developed by Andras)

## Installation

Eventually we will be able to set this up install using the latest from the repository, using the following:

    pip install git+https://tessgit.mit.edu/tess/FFITools.git

but for the time being use ```virtualenv```

* Start with creating a new virtual environment to leave the environment use ```deactivate``):

  ```
  virtualenv env
  source env/bin/activate
  ```

Note that throughout we use ```<ENV>``` to refer to the full path name to the 'env' directory created in the previous step, e.g. ```/home/username/tessgit.mit.edu/FFITools/```

* clone the repository:

  ```
  git clone https://tessgit.mit.edu/tess/FFITools.git
  ```

* install all the PyPI packages mentioned above via ```pip```, e.g.

  ```
  pip install numpy==1.11.1
  ```

* install ```fitsh```: get 0.9.2 from https://fitsh.net/

 ```
 ./configure --prefix=<ENV>
 make
 make install
 ```

* install ```astrometry.net```:


 * first install ```cfitsio```, download tarball from http://heasarc.gsfc.nasa.gov/fitsio/

  ```
  ./configure --prefix=<ENV>
  make
  make install
  ```

 *  then get astrometry.net tarball from: http://astrometry.net/use.html

 ```
 PKG_CONFIG_PATH=<ENV>/lib/pkgconfig/ make 
 INSTALL_DIR=<ENV> make install
 ```

* vartools: get 1.34 from http://www.astro.princeton.edu/~jhartman/vartools.html

see https://arxiv.org/pdf/1605.06811.pdf for reference.

  ```
  ./configure --prefix=<ENV> CPPFLAGS="-I<ENV>/include" LDFLAGS="-L<ENV>/lib"
  ```

  (```CPPFLAGS``` and ```LDFLAGS``` needed to enable vartools's configure script to locate ```cfitsio```)

* dvipng: see detailed instructions here: https://mikewest.org/2007/04/installing-libgd-from-source-on-os-x

## Examples

These are examples of how to use FFITools.

### Example 1: detrend a light curve with a single core

    lctools-detrend -c test.cfg -n 1

### Example 2: run bls on a single light curve with a single core

    lctools-bls -c test.cfg -n 1


### Example3: run dft on a single light curve with a single core

    lctools-dft -c test.cfg -n 1
    
### Example4: generate one page summary on a single star with a single core
    
    lctools-candsum -c test.cfg -n 1


### Example5: demonstrate how patools fistar (source extractor) works

    patools-sampler fistar -c example.cfg 
    
### Example6: demonstrate how patools grmatch (astrometry) works
    
    patools-sampler astrometry -c example.cfg
    
### Example7: demonstrate how patools fiphot (photometry) works

    patools-sampler phot -c example.cfg

## Development

The code for TSIG is at https://tessgit.mit.edu/tess/FFITools.git

## Copyright and License

FFITools is Copyright(c) Chelsea Huang, all rights reserved

FFITools is distributed under the terms of GPLv3