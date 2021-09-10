#!/bin/sh

curl -LJO https://github.com/cea-hpc/modules/releases/download/v4.7.1/modules-4$
tar xfz modules-4.7.1.tar.gz
cd modules-4.7.1
./configure
make
make install
cd ~/
export LD_LIBRARY_PATH=/usr/lib64/openmpi-1.10/lib:$LD_LIBRARY_PATH
module load mpi/openmpi-x86_64

Rscript -e 'install.packages("Rmpi", configure.args = paste("--with-Rmpi-includ$
Rscript -e 'install.packages("idr")'

echo "Succesful Installation of RMPI and Modules"