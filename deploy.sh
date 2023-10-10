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
cd ..
wget https://github.com/bluescarni/mppp/archive/refs/tags/v0.26.tar.gz
tar -xvf v0.26.tar.gz
rm v0.26.tar.gz
cd mppp-0.26
cmake -DMPPP_WITH_QUADMATH=y -DMPPP_WITH_MPFR=y -DMPPP_WITH_MPC=y .
make -j4
sudo make install
