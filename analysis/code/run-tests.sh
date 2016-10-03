#!/bin/bash
DATE=`date +%Y-%m-%d`
OUTFILE="/home/ben/GitHub/SCR-Book/analysis/test-results/out-$DATE.log"
echo -e "Deleting .RData objects..."
cd /home/ben/GitHub/SCR-Book/analysis/objects
rm -rfv *
echo -e "DONE\n"
cd /home/ben/GitHub/SCR-Book/analysis/code/
echo -e "Installing packages..." 
R -q -e 'source("packages.R")' 1>/dev/null
echo -e "DONE\n"
for i in *.R; do
    if [ $i != "packages.R" -a $i != "tests.R" ]; then
	echo -e "Running $i..."
	R -q --no-save < $i 1>/dev/null
	echo -e "DONE\n"
    fi 
done
echo -e "Running tests..."
R -q -e 'library(testthat); test_file("tests.R")'
