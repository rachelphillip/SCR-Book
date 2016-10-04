#!/bin/bash

## Flags:
## -t to run tests.
## -u to update R packages.

## Sorting out arguments.
doupdate=false
dotest=false
while getopts ":tu" opt; do
    case $opt in
	t)
	    dotest=true
	    ;;
	u)
	    doupdate=true
	    ;;
	\?)
	    echo "Invalid option: -${OPTARG}." >&2
	    ;;
    esac
done
## Printing date.
echo -e "$(date)\n"
## Deleting all .RData objects in analysis/objects/.
echo -e "Deleting .RData objects..."
cd /home/ben/GitHub/SCR-Book/analysis/objects
rm -rfv *
echo -e "DONE\n"
cd /home/ben/GitHub/SCR-Book/analysis/code/
## Updating packages.
if [ "$doupdate" = true ]; then
    echo -e "Updating packages..."
    R -q --no-save < packages.R 1>/dev/null --args -u
    echo -e "DONE\n"
fi
## Installing packages.
echo -e "Installing packages..." 
R -q --no-save < packages.R 1>/dev/null
echo -e "DONE\n"
for i in *.R; do
    if [ $i != "packages.R" -a $i != "tests.R" ]; then
	echo -e "Running $i..."
	R -q --no-save < $i 1>/dev/null
	echo -e "DONE\n"
    fi 
done
## Running tests.
if [ "$dotest" = true ]; then
    echo -e "Running tests..."
    R -q -e 'library(testthat); test_file("tests.R")'
fi
