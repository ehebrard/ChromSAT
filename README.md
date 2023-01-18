

# ChromSAT

## First go to branch "cpaior2019"

git checkout cpaior2019

## Installation

* meson 

pip3 install --user meson

* sparsehash

git clone https://github.com/sparsehash/sparsehash.git

install it

* ninja

brew install ninja

* tclap

(download and make symbolic link to the tclap/ dir in chromsat)

* minicsp

cd somewhere

git clone https://github.com/gkatsi/minicsp.git

cd somewhere/minicsp/minicsp/core

make

(then from the chromsat dir)

ln -s somewhere/minicsp/minicsp ./minicsp

* build

meson setup builddir

meson compile -C builddir

## Command line

* Default Algorithm 

/builddir/gc <path to data file>

* Algorithm for very large graphs (see CPAIOR2019 paper)

./builddir/gc --strategy 7 --sdsaturiter 20 --randwalkiter 0
--dsatlimit 0  --switchdescent --idsaturlimit -1 --core 3  --lsextra
-1 --rw 1 --dynrandpath --rpfactor 2 --rpdiv 1 --dynamiclimit 2
--dynfactor 2 --dyndiv 1 --lsiter 250000 --maxclique --domegatime 600
--seed <some number> <path to data file>

