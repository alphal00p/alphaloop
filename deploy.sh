cd libraries/fjcore
make
cd ..
cd QGRAF
make
cd ..
git clone https://github.com/embotech/ecos.git
cd ecos
sed -i -e 's/CTRLC=1/CTRLC=0/g' ecos.mk
make
cd ..
git clone https://github.com/cvxgrp/scs.git
cd scs
sed -i -e 's/CTRLC = 1/CTRLC = 0/g' scs.mk
make
cd ..
wget http://www.feynarts.de/cuba/Cuba-4.2.tar.gz
tar -xvf Cuba-4.2.tar.gz
cd Cuba-4.2
./configure
make