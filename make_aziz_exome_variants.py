#!/usr/bin/env python
import os
import sys
import logging

from argparse import ArgumentParser
from collections import defaultdict

def script(vcf_path, out_folder):
    logging.basicConfig(level='DEBUG')
    # So we don't have to check if chr is there on every line
    chr_state = 0
    vcf_file = open(vcf_path)

    variant_lines = []
    exome_lines = []
    genes = defaultdict(list)
    weights = []
    cols = []

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
        # Assume that refgene is 4th info field, otherwise look for it manually
        if len(info) > 3 and info[3].startswith('Gene.refGene'):
            gene = info[3].split('=')[1].split(',')
        else:
            for field in info:
                if field.startswith('Gene.refGene'):
                    gene = field.split('=')[1].split(',')
        # Assume that func.refgene is 3rd info field, otherwise look for it manually
        if len(info) > 2 and info[2].startswith('ExonicFunc.refGene'):
            func = info[2].split('=')[1]
        else:
            for field in info:
                if field.startswith('ExonicFunc.refGene'):
                    func = field.split('=')[1]
        
        # Assume that polyphen2 hdiv score is the 10th info field, otherwise look for it manually
        if len(info) > 9 and info[9].startswith('Polyphen2_HDIV_score'):
            score = info[9].split('=')[1]
        else:
            for field in info:
                if field.startswith('Polyphen2_HDIV_score'):
                    score = field.split('=')[1]
        
         # Do actual insertions
        if ',' in alt:
            alts = alt.split(',')
            for i, allele in enumerate(alts):
                insert_variant(chrom, pos, ref, allele, func, genos, gene, variant_lines, exome_lines, weights, score, genes, which_alt=i+1)
        else:
                insert_variant(chrom, pos, ref, alt, func, genos, gene, variant_lines, exome_lines, weights, score, genes)


        vcf_line_count += 1
        if vcf_line_count % 100000 == 0:
            logging.info('Finished parsing %d variants' % vcf_line_count)
            logging.info('Currently %d genes in dictionary' % len(genes.keys()))
        
    # Do writing
    # Create out folder if necessary
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    # Make files
    variant_file = open(out_folder + '/variants.txt', 'w')
    exome_file = open(out_folder + '/exome.txt', 'w')

    # Write variants
    for gene in sorted(genes.keys()):
        for i, var_ind in enumerate(genes[gene]):
            info = '_'.join(variant_lines[var_ind][:4])
            annotation = variant_lines[var_ind][4]
            score = weights[var_ind]
            maf = float(sum([int(x) for x in exome_lines[i] if x == '1' or x == '2'])) / (2.0*len(cols))
            variant_file.write('\t'.join([info, gene, annotation, str(score), str(maf)]) + '\n')

    # Write exomes
    for i, var in enumerate(variant_lines):
        info = '_'.join(var[:4])
        for j, col in enumerate(cols):
            if exome_lines[i][j] == '1' or exome_lines[i][j] == '2':
                exome_file.write('\t'.join([info, col, exome_lines[i][j]]) + '\n')

    # Close files
    variant_file.close()
    exome_file.close()

    logging.info('Done writing!')


def insert_variant(chrom, pos, ref, alt, func, genos, gene, variant_lines, exome_lines, weights, score, genes, which_alt=1):
        weights.append(score)
        #logging.debug(' '.join([chrom, pos, ref, alt]))
        #logging.debug(genos)
        variant_lines.append([chrom, pos, ref, alt, func])
        exome_lines.append([str(sum(int(y) == which_alt for y in x.split('/'))) if not (x == '.' or x == './.') else '.' for x in genos])
        # zero indexed counting
        for g in gene:
            genes[g].append(len(variant_lines)-1)

   
def parse_args(args):
    parser = ArgumentParser(description='Use annovar file to create exome and variant files.')
    parser.add_argument('vcf_path', metavar = 'ANNO', help='The vcf file outputted by running annovar')
    parser.add_argument('out_folder', metavar='OUT')
    return parser.parse_args(args)

def main(args = sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
