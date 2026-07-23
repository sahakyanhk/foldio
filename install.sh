#!/bin/bash

mkdir -p bin
wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz; 
tar xvzf foldseek-linux-avx2.tar.gz; 
mv foldseek/bin/foldseek bin/
rm -r foldseek*

chmod +x bin/foldseek

