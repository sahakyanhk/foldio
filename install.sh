#!/bin/bash

mkdir -p bin
wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz; 
tar xvzf foldseek-linux-avx2.tar.gz; 
mv foldseek/bin/foldseek bin/

rm -r foldseek*

wget https://github.com/rcedgar/reseek/releases/download/v3.01/reseek-v3.01-linux-x86 -O bin/reseek

chmod +x bin/*


export PATH=$(pwd)/bin/:$PATH

