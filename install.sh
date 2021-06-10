#! /bin/bash -f

#### .bashrc and .bash_aliases must have been scp
#### scp ~/.bash* cluster:
# fftw folder must have been scp to ~/Documents
#### scp -r fftw cluster:


#### get the electron-diffraction library
mkdir -p ~/Documents/git/ccp4/src/
cd ~/Documents/git/ccp4/src/
git clone https://www.github.com/ccp4/electron-diffraction.git

#### install and test fftw
cd ~/Documents/fftw/fftw-3.3.8
./configure --enable-threads --enable-float
cd ../
make
./test

#### install temsim
cd $multislice/temsim
mkdir -p ../bin/obj
make


#### install python locally
cd ~
wget https://www.python.org/ftp/python/3.8.5/Python-3.8.5.tgz
tar zxfv Python-3.8.5.tgz
find ~/Python-3.8.5 -type d | xargs chmod 0755
cd ~/Python-3.8.5
./configure --prefix=$HOME/Python-3.8.5
make
printf 'alias python3="~/Python-3.8.5/python"' >> ~/.bash_aliases
srcAlias

#### pip
wget https://bootstrap.pypa.io/get-pip.py
python3 get-pip.py --user


#### test the library
cdtests
python3 -m tests/
