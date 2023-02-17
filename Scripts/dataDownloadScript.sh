#!/usr/bin

##clone datasets if needed. Use the right links and remember to change the zenodo link in the echo command
DIR="../Data"

if [ ! -d $DIR ]; then
    echo "==============================================================" 
    echo "Copying datasets from https://zenodo.org/record/6952995"
    echo "=============================================================="
    
    mkdir -p ../Data
    curl https://zenodo.org/record/6952995/files/clover.tar.gz?download=1 -o ../Data/Clover_Data.tar.gz
    tar -zxvf ../Data/Clover_Data.tar.gz -C ../Data/
    curl https://zenodo.org/record/6952995/files/singlecell.tar.gz?download=1 -o ../Data/scrna_Data.tar.gz
    tar -zxvf ../Data/scrna_Data.tar.gz -C ../Data/
    rm -f ../Data/*.tar.gz
else
    echo "========================================================"
    echo "Datasets folder already exists, no need to download it"
    echo "========================================================"
fi