#!/usr/bin/env python
import os
import sys
import logging

import numpy as np
from argparse import ArgumentParser
from collections import defaultdict

def script(vcf_path, score_path, out_folder, score_field, trans, pool_spidex=None):
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
        

    # Get the rest of our stuff
    # List of variant info lines
    variant_lines = []
    # List of variant weight lines based on whatever scoring file we are using
    weight_lines = []
    # List of genotypes lists, each entry is the list of genotypes of variant for each patient
    matrix_lines = []
    # Dictionary of genes to variants in each gene
    genes = defaultdict(list)
    vcf_file = open(vcf_path)

    vcf_line_count = 0
    for line in vcf_file:
        if line.startswith('#CHR'):
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
        # Do actual insertions
        if ',' in alt:
            alts = alt.split(',')
            for i, allele in enumerate(alts):
                insert_variant(chrom, pos, ref, allele, genos, gene, variant_lines, weight_lines, matrix_lines, genes, scores, which_alt=i+1)
        else:
                insert_variant(chrom, pos, ref, alt, genos, gene, variant_lines, weight_lines, matrix_lines, genes, scores)

        vcf_line_count += 1
        if vcf_line_count % 100000 == 0:
            logging.info('Finished parsing %d variants' % vcf_line_count)
            logging.info('Currently %d genes in dictionary' % len(genes.keys()))

    assert len(variant_lines) == len(weight_lines)
    assert len(variant_lines) == len(matrix_lines)

    logging.info('Writing to files')

    # Create out folder if necessary
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
    
    # Dump stuff to file
    variant_file = open(out_folder + '/variants.txt', 'w')
    weight_file = open(out_folder + '/weights.txt', 'w')
    matrix_file = open(out_folder + '/genotype_matrix.txt', 'w')
    for i in range(len(variant_lines)):
        variant_file.write('\t'.join(variant_lines[i]) + '\n')
        weight_file.write(str(weight_lines[i]) + '\n')
        matrix_file.write(','.join(matrix_lines[i]) + '\n')
    
    variant_file.close()
    weight_file.close()
    matrix_file.close()
    # Create folder for genes
    if not os.path.exists(out_folder + '/genes'):
        os.makedirs(out_folder + '/genes')
    for gene, vars in genes.items():
        gene_file = open(out_folder + '/genes/' + gene + '.txt', 'w')
        for var in vars:
            gene_file.write(str(var) + '\n')
        gene_file.close()

    logging.info('Done writing!')

def insert_variant(chrom, pos, ref, alt, genos, gene, variant_lines, weight_lines, matrix_lines, genes, scores, which_alt=1):
        if not ''.join([chrom, pos, ref, alt]) in scores: return
        weight_lines.append(scores[''.join([chrom, pos, ref, alt])])
        #logging.debug(' '.join([chrom, pos, ref, alt]))
        #logging.debug(genos)
        variant_lines.append([chrom, pos, ref, alt])
        matrix_lines.append([str(sum(int(y) == which_alt for y in x.split('/'))) if not (x == '.' or x == './.') else '.' for x in genos])
        # One indexed counting
        for g in gene:
            genes[g].append(len(variant_lines))

def parse_args(args):
    parser = ArgumentParser(description='Use scoring file and vcf with genotypes to generate gene and phenotype matrices')
    parser.add_argument('vcf_path', metavar='VCF')
    parser.add_argument('score_path', metavar='SCORE')
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
