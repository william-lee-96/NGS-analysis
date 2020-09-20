#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 12:44:57 2020

@author: williamlee
"""

import os
path = "/Users/williamlee"
os.chdir(path)

import pandas as pd 

ribosomal_rna = pd.read_csv('Downloads/ribosomal_rna.txt', sep='\t')
rRNA_genes=pd.Series(ribosomal_rna['Approved symbol'].unique())
rRNA_genes.to_csv("ribosomal RNA genes", sep='\t', index=False)

translation_init_go = pd.read_csv('Downloads/translation_init_go.txt', sep='\t', header=1)
trans_init_genes_go=pd.Series(translation_init_go['> Translational initiation'].unique())

miRNA = pd.read_csv('Downloads/miRNA.txt', sep='\t')
miRNA_genes=pd.Series(miRNA['Approved symbol'].unique())
miRNA_genes.to_csv("miRNA genes", sep='\t', index=False)

proteasome = pd.read_csv('Downloads/proteasome.txt', sep='\t')
proteasome_genes=pd.Series(proteasome['Approved symbol'].unique())
proteasome_genes.to_csv("proteasome genes", sep='\t', index=False)

ubiquitin_e1 = pd.read_csv('Downloads/ubiquitin_e1.txt', sep='\t')
ubiquitin_e2 = pd.read_csv('Downloads/ubiquitin_e2.txt', sep='\t')
ubiquitin_e3_ubox = pd.read_csv('Downloads/ubiquitin_e3_ubox.txt', sep='\t')
ubiquitin_e3_ring = pd.read_csv('Downloads/ubiquitin_e3_ring.txt', sep='\t')
ubiquitin_e3_peli = pd.read_csv('Downloads/ubiquitin_e3_peli.txt', sep='\t')
ubiquitin = pd.concat([ubiquitin_e1, ubiquitin_e2, ubiquitin_e3_ubox,
                       ubiquitin_e3_ring, ubiquitin_e3_peli])
ubiquitin_genes=pd.Series(ubiquitin['Approved symbol'])

ubiquitin_kegg=pd.read_csv('Downloads/ubiquitin_kegg.txt', sep='\t', header=1)
ubiquitin_genes_kegg=pd.Series(ubiquitin_kegg['> Ubiquitin mediated proteolysis'])
all_ubiquitin_genes=pd.Series(pd.concat([ubiquitin_genes, ubiquitin_genes_kegg]).unique())
all_ubiquitin_genes.to_csv("ubiquitin genes", sep='\t', index=False)

ubiquitin_partial = all_ubiquitin_genes[~all_ubiquitin_genes.isin(ubiquitin_genes_kegg)]
ubiquitin_partial.to_csv("ubiquitin genes (partly confirmed)", sep='\t', index=False)

autophagy_kegg=pd.read_csv('Downloads/autophagy.txt', sep='\t', header=1)
autophagy_genes_kegg=pd.Series(autophagy_kegg['> Regulation of autophagy'].unique())

tRNA = pd.read_csv('Downloads/tRNA.txt', sep='\t')
tRNA_genes=pd.Series(tRNA['Approved symbol'].unique())
tRNA_genes.to_csv("tRNA genes", sep='\t', index=False)

nuclear_pore_go=pd.read_csv('Downloads/nuclear_pore_go.txt', sep='\t', header=1)
nuclear_pore_genes=pd.Series(nuclear_pore_go['> Nuclear pore'].unique())