#!/bin/bash

pdflatex poster.tex
bibtex poster
pdflatex poster.tex
bibtex poster
pdflatex poster.tex
evince poster.pdf
