#!/usr/bin/env bash
#
# run_ANIm.sh
#
# Run ANIm analysis (using pyani v0.3) on downloaded genomes

# Create database
pyani createdb -l data/pectobact/tree/logs/pyani_01_createdb.log

# Index genomes
pyani index genomes -l data/pectobact/tree/logs/pyani_02_index.log

# Run ANIm analysis
pyani anim -l data/pectobact/tree/logs/pyani_03_anim.log \
  data/pectobact/genomes \
  data/pectobact/tree/anim_output \
  --name "pecto_dic_ANIm"

# Generate graphical anim output
pyani plot -l data/pectobact/tree/logs/pyani_04_plot.log \
  --formats png,pdf \
  --method seaborn \
  data/pectobact/tree/anim_output 1
  