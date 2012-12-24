#!/bin/bash

pdflatex pref_learning_2012.tex
bibtex pref_learning_2012
pdflatex pref_learning_2012.tex
bibtex pref_learning_2012
pdflatex pref_learning_2012.tex
evince pref_learning_2012.pdf
