# CEC 2024 bound constrained single objective optimization comparison

This repository contains results of CEC 2024 competition and code for comparison.

The BC-SOPs folder contains the raw results. Use read_and_reformat.m to generate *.txt files.

Next, there are two options: python code in ranking_and_graph.py or jupyter notebook CEC24-FinalAnalysis.ipynb. The code in these two is identical.

# Ranking results:

U-scores:
BlockEA   28188.0
IEACOP   17877.0
RDE   57437.0
mLSHADE_LR   46214.5
L_SRTDE   69302.5
jSOa   43481.0

Friedman ranking:
BlockEA   123.36
IEACOP   139.08
RDE   76.36
mLSHADE_LR   94.24
L_SRTDE   57.08
jSOa   97.88

Mann-Whitney tests:
BlockEA vs L_SRTDE (23+/1=/5-)
IEACOP vs L_SRTDE (27+/1=/1-)
RDE vs L_SRTDE (18+/5=/6-)
mLSHADE_LR vs L_SRTDE (21+/1=/7-)
jSOa vs L_SRTDE (23+/4=/2-)
