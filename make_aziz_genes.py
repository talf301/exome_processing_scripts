#!/usr/bin/env python
import os
import sys
import logging

import numpy as np
from argparse import ArgumentParser
from collections import defaultdict

def script(vcf_path, score_path, ages_path, out_folder, score_field, trans, pool_spidex=None):
    logging.basicConfig(level='DEBUG')
    # So we don't have to check if chr is there on every line
    chr_state = 0
    # Populate our dictionary of scores
    score_file = open(score_path)
    scores = {}
    logging.info('Loading score file...')
    for line in score_file:
        if line.startswith('#'): continue
        tokens = line.strip().split('\t')
        var_info = tokens[:4]
        if var_info[0].startswith('chr'): var_info[0] = var_info[0][3:]
        score = float(tokens[score_field])
        if trans:
            scores[''.join(var_info)] = max(1-np.exp(-0.5*score), 0)
        else:
            scores[''.join(var_info)] = score
    score_file.close()
    logging.info('Done loading score file.')

    if pool_spidex:
        score_file = open(pool_spidex)
        logging.info('Loading spidex score file...')
        for line in score_file:
            if line.startswith('#'): continue
            tokens = line.strip().split('\t')
            var_info = tokens[:4]
            if var_info[0].startswith('chr'): var_info[0] = var_info[0][3:]
            if trans:
                score = np.exp(-2/(np.abs(float(tokens[score_field]))+0.001))
            else:
                score = float(tokens[score_field])
            scores[''.join(var_info)] = max(scores.get(''.join(var_info), 0), score)
        score_file.close()
        logging.info('Done loading spidex score file.')
        
    ages = []
    with open(ages_path) as ages_file:
        for line in ages_file:
            ages.append(float(line.strip()))

    inds = [i for i in range(len(ages)) if ages[i] < 46 or ages[i] >= 49]
    # Get the rest of our stuff
    # List of variant info lines
    variant_lines = []
    # List of variant weight lines based on whatever scoring file we are using
    weight_lines = []
    # Genotype matrix, basically
    matrix_lines = []
    # Dictionary of genes to variants in each gene
    genes = defaultdict(list)

    cols = []
    vcf_file = open(vcf_path)
   
    vcf_line_count = 0
    for line in vcf_file:
        if line.startswith('#CHR'):
            # cols is a list of the names of the patients
            cols = line.strip().split('\t')[9:]
        if line.startswith('#'): continue
        tokens = line.strip().split('\t')
        # Strip off the chr part if it's there
        # chr_state = 0 means hasn't been checked yet, 1 means present, 2 means not
        if chr_state == 0:
            if tokens[0].startswith('chr'):
                chrom = tokens[0][3:]
                chr_state = 1
            else:
                chrom = tokens[0]
                chr_state = 2
        elif chr_state == 1:
            chrom = tokens[0][3:]
        elif chr_state == 2:
            chrom = tokens[0]


        pos = tokens[1]
        ref = tokens[3]
        alt = tokens[4]
        genos = [x.split(':')[0] for x in tokens[9:]]
        info = tokens[7].split(';')
        gene = []
        # Assume that refgene is 28th info field, otherwise look for it manually
        if len(info) > 28 and info[28].startswith('Gene.refGene'):
            gene = info[28].split('=')[1].split(',')
        else:
            for field in info:
                if field.startswith('Gene.refGene'):
                    gene = field.split('=')[1].split(',')
        # Assume that func.refgene is 30th info field, otherwise look for it manually
        if len(info) > 30 and info[30].startswith('ExonicFunc.refGene'):
            func = info[30].split('=')[1]
        else:
            for field in info:
                if field.startswith('ExonicFunc.refGene'):
                    func = field.split('=')[1]
        # Do actual insertions
        if ',' in alt:
            alts = alt.split(',')
            for i, allele in enumerate(alts):
                insert_variant(chrom, pos, ref, allele, func, genos, gene, variant_lines, matrix_lines, weight_lines, genes, scores, which_alt=i+1)
        else:
                insert_variant(chrom, pos, ref, alt, func, genos, gene, variant_lines, matrix_lines, weight_lines, genes, scores)

        vcf_line_count += 1
        if vcf_line_count % 100000 == 0:
            logging.info('Finished parsing %d variants' % vcf_line_count)
            logging.info('Currently %d genes in dictionary' % len(genes.keys()))
    

    assert len(ages) == len(cols)
    logging.info('Writing to files')

    # Create out folder if necessary
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
    
    # Dump stuff to file
    variant_file = open(out_folder + '/variants.txt', 'w')
    gene_file = open(out_folder + '/genes.txt', 'w')
    patient_file = open(out_folder + '/patients.txt', 'w')
    exome_file = open(out_folder + '/exome.txt', 'w')

    # Write variants
    for gene in sorted(genes.keys()):
        for i, var_ind in enumerate(genes[gene]):
            info = '_'.join(variant_lines[var_ind][:4])
            annotation = variant_lines[var_ind][4]
            score = weight_lines[var_ind]
            maf = float(sum([int(x) for x in matrix_lines[i] if x == '1' or x == '2'])) / (2.0*len(cols))
            gene_index = str(i+1)
            variant_file.write('\t'.join([gene, info, annotation, str(score), str(maf), gene_index]) + '\n')

    # Write genes
    for i, gene in enumerate(sorted(genes.keys())):
        gene_file.write('\t'.join([gene, str(i+1), str(len(genes[gene]))]) + '\n')

    # Write patients
    for i in inds:
        patient_id = cols[i]
        patient_index = i+1
        pheno = int(ages[i] < 46)
        patient_file.write('\t'.join([patient_id, str(patient_index), str(pheno)]) + '\n')

    # Write exomes
    for i, var in enumerate(variant_lines):
        info = '_'.join(var[:4])
        for ind in inds:
            if matrix_lines[i][ind] == '1' or matrix_lines[i][ind] == '2':
                exome_file.write('\t'.join([info, cols[ind], matrix_lines[i][ind]]) + '\n')
                
    variant_file.close()
    gene_file.close()
    patient_file.close()
    exome_file.close()
    
    logging.info('Done writing!')

def insert_variant(chrom, pos, ref, alt, func, genos, gene, variant_lines, matrix_lines, weight_lines, genes, scores, which_alt=1):
        if not ''.join([chrom, pos, ref, alt]) in scores: return
        weight_lines.append(scores[''.join([chrom, pos, ref, alt])])
        #logging.debug(' '.join([chrom, pos, ref, alt]))
        #logging.debug(genos)
        variant_lines.append([chrom, pos, ref, alt, func])
        matrix_lines.append([str(sum(int(y) == which_alt for y in x.split('/'))) if not (x == '.' or x == './.') else '.' for x in genos])
        # zero indexed counting
        for g in gene:
            genes[g].append(len(variant_lines)-1)

def parse_args(args):
    parser = ArgumentParser(description='Use scoring file and vcf with genotypes to generate gene and phenotype matrices')
    parser.add_argument('vcf_path', metavar='VCF')
    parser.add_argument('score_path', metavar='SCORE')
    parser.add_argument('ages_path', metavar='AGES')
    parser.add_argument('out_folder', metavar='OUT')
    parser.add_argument('--score_field', metavar='FIELD', help='Integer (0-indexed) indicating which field in the tab separated file contains the score we are'
            'interested in. Default is 4 for cadd. Use 5 for spidex z-scores', type=int, default=4)
    parser.add_argument('--trans', help='Do transformation of all scores (assumes main score is cadd, pooled is spidex', action='store_true')
    parser.add_argument('--pool_spidex', help='Spidex score file to optionally use for pooling with cadd scores')
    return parser.parse_args(args)

def main(args = sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
