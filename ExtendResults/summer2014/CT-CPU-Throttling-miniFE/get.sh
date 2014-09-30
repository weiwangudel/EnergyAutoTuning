#!/bin/bash

gnuplot 3d.gnpl
cp output.pdf time.pdf
gnuplot 3d-energy.gnpl
cp output.pdf energy.pdf 
pdftk energy.pdf time.pdf output et.pdf
evince et.pdf 
