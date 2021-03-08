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
wget http://www.feynarts.de/cuba/Cuba-4.2.1.tar.gz
tar -xvf Cuba-4.2.1.tar.gz
cd Cuba-4.2.1
./configure
make
cd ..
rm Cuba-4.2.1.tar.gz
git clone https://github.com/vermaseren/form.git
cd form
autoreconf -i
./configure
make -j4