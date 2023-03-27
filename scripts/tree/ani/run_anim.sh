#!/usr/bin/env bash
#
# run_ANIm.sh
#
# Run ANIm analysis (using pyani v0.3) on downloaded genomes

OUTPUT_DIR=data/pectobact/tree/ani_tree
GENOME_DIR=data/pectobact/genomes

# make output dir
mkdir -p $OUTPUT_DIR/logs

# Create database
pyani createdb -l $OUTPUT_DIR/logs/pyani_01_createdb.log

# Index genomes
pyani index \
  -i $GENOME_DIR \
  -l $OUTPUT_DIR/logs/pyani_02_index.log

# Run ANIm analysis
pyani anim \
  -i $GENOME_DIR\
  -o $OUTPUT_DIR/anim_output \
  -l $OUTPUT_DIR/logs/pyani_03_anim.log \
  --name "pecto_dic_ANIm"

# Generate graphical anim output
pyani plot
  -l $OUTPUT_DIR/logs/pyani_04_plot.log \
  --formats png,pdf \
  --method seaborn \
  -o $OUTPUT_DIR/anim_output
  --run_id 1
  