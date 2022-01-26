#/usr/bin/env bash

DIR=`pwd`
cd ~/Downloads
wget https://webfiles.york.ac.uk/INFODATA/44145f0a-5d82-4604-9494-7cf71190bd82/Bond_et_al_2020_data.tar.gz
tar -xf Bond_et_al_2020_data.tar.gz
tar -xf data/full_reduced_set.tar.gz
mv reduced_full $DIR/data/mr
