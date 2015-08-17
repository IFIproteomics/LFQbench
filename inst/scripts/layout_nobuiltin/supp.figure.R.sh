#!/bin/bash
cd "`dirname "$0"`"

pdflatex supp.figure.R

rm *.log *aux