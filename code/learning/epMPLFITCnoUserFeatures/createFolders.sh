#!/bin/bash

rm -rf simulation*/*
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
do
	echo Creating subfolder $i
	mkdir simulation$i
	cp -R code/* simulation$i
	echo $i > simulation$i/simulationNumber.txt
	cd simulation$i
	nice -18 /home/jmh233/R-2.14.1/bin/R --no-save < experiment.R
	cd ..
done
