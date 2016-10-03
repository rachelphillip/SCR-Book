#!/bin/sh
cd /home/ben/GitHub/SCR-Book/analysis/code/
DATE=`date +%Y-%m-%d`
OUTFILE="/home/ben/GitHub/SCR-Book/analysis/test-results/out-$DATE.txt"
touch $OUTFILE
echo "Installing packages..." &> $OUTFILE
R -q -e 'source("packages.R")' &> $OUTFILE
echo DONE
for i in *.R; do
    echo item: $i &> $OUTFILE
done
##DATE=`date +%Y-%m-%d`
##R -q -e 'library(testthat); test_file("tests.R")' > ~/GitHub/SCR-Book/analysis/test-results/out-$DATE.txt
