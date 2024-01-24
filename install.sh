#! /bin/bash -f

#### .bashrc and .bash_aliases must have been scp prior to executing
# this script
#scp ../install.sh /home/tarik/.bash* /home/tarik/bin/colors.sh  <instance>:
black="\u1b[39m"
red="\u1b[31m"
green="\u1b[32m"


printf $green"installing python locally\n"$black
### if python
# wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# chmod +x Miniconda3-latest-Linux-x86_64.sh

pip3 install IPython

# sudo apt-get install python3
# sudo apt-get install python3-dev
# printf $green"installing pip\n"$black
# wget https://bootstrap.pypa.io/get-pip.py
# python3 get-pip.py --user
# pip3 install IPython
###### local python installation
# wget https://www.python.org/ftp/python/3.8.5/Python-3.8.5.tgz
# tar zxfv Python-3.8.5.tgz
# find ~/Python-3.8.5 -type d | xargs chmod 0755
# cd ~/Python-3.8.5
# ./configure --prefix=$HOME/Python-3.8.5
# make
# # printf 'alias python3="~/Python-3.8.5/python"' >> ~/.bash_aliases
# # srcAlias




printf $green"Installing and testing fftw\n"$black
wget http://www.fftw.org/fftw-3.3.9.tar.gz
tar -xvzf fftw-3.3.9.tar.gz
if [ -d ~//Documents/fftw ]; then mkdir -p ~/Documents/fftw;fi
mv fftw-3.3.9 Documents/fftw
cd ~/Documents/fftw/fftw-3.3.9
./configure --enable-threads --enable-float
make
sudo make install
# rm fftw-3.3.9.tar.gz
# cd ..
# make
# ./test
cd ~/


printf $green"Installing Felix\n"$black
cd $ED/..
git clone https://www.github.com/ronandrevon/Felix.git
# git branch continuousED; git checkout
git checkout continousED
#git pull origin continousED
sudo apt-get install mpi libopenmpi-dev libblas-dev liblapack-dev
cd src
make
cd ~/bin
ln -s $pdir/src/felix.OPT64NGNU.d
printf "export LD_LIBRARY_PATH=/usr/local/lib/" >>~/.bash_aliases

printf $green"installing python multislice library\n"$black
cd $ED
pip3 install -e .


printf $green"Installing and testing temsim\n"$black
mkdir -p ~/Documents/git/ccp4/src/
cd ~/Documents/git/ccp4/src/
git clone https://www.github.com/ccp4/electron-diffraction.git
cd $multislice/temsim
mkdir -p ../bin/obj
make
cd ~


printf $green"Testing multislice\n"$black
cd $multislice/examples
mkdir -p ../dat/test
python3 test_base.py
