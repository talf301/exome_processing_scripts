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

    cols = []
    
    # Do writing
    # Create out folder if necessary
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    # Make files
    variant_file = open(out_folder + '/variants.txt', 'w')
    exome_file = open(out_folder + '/exome.txt', 'w')

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
        
        # Assume that polyphen2 hdiv score is the 12th info field, otherwise look for it manually
        if len(info) > 11 and info[11].startswith('Polyphen2_HVAR_score'):
            score = info[11].split('=')[1]
        else:
            for field in info:
                if field.startswith('Polyphen2_HVAR_score'):
                    score = field.split('=')[1]
        
         # Do actual insertions
        if ',' in alt:
            alts = alt.split(',')
            for i, allele in enumerate(alts):
                insert_variant(chrom, pos, ref, allele, func, genos, gene, score, variant_file, exome_file, cols, which_alt=i+1)
        else:
                insert_variant(chrom, pos, ref, alt, func, genos, gene, score, variant_file, exome_file, cols)


        vcf_line_count += 1
        if vcf_line_count % 100000 == 0:
            logging.info('Finished parsing %d variants' % vcf_line_count)
        
    # Close files
    variant_file.close()
    exome_file.close()

    logging.info('Done writing!')


def insert_variant(chrom, pos, ref, alt, func, genos, gene, variant_file, exome_file, score, cols, which_alt=1):
        if pos == '.':
            return

        # Make patient exome info a little easier to read
        exome_line = [str(sum(int(y) == which_alt for y in x.split('/'))) if not (x == '.' or x == './.') else '.' for x in genos]

        # Write variant info
        info = '_'.join([chrom, pos, ref, alt])
        maf = float(sum([int(x) for x in exome_lines[i] if x == '1' or x == '2'])) / (2.0*len(cols))
        for g in gene:
            variant_file.write('\t'.join([info, g, func, str(score), str(maf)]) + '\n')

        # Write exome info
        for i, col in enumerate(cols):
            if exome_line[i] == '1' or exome_line[i] == '2':
                exome_file.write('\t'.join([info, col, exome_line[i]]) + '\n')
        
   
def parse_args(args):
    parser = ArgumentParser(description='Use annovar file to create exome and variant files.')
    parser.add_argument('vcf_path', metavar = 'ANNO', help='The vcf file outputted by running annovar with the -vcf flag')
    parser.add_argument('out_folder', metavar='OUT', help='The directory in which to put output files')
    return parser.parse_args(args)

def main(args = sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
