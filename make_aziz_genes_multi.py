#!/usr/bin/env python
import os
import sys
import logging

import numpy as np
from argparse import ArgumentParser
from collections import defaultdict

def script(vcf_path, score_path, gender_path, cc_path,  out_folder, score_field, pool_spidex=None):
    logging.basicConfig(level='DEBUG')
    # So we don't have to check if chr is there on every line
    chr_state = 0
    # Populate our dictionary of scores
    # First is for CADD, second is spidex
    score_file = open(score_path)
    scores_gen_gen = {}
    scores_gen_con = {}
    scores_con_gen = {}
    scores_con_con = {}
    logging.info('Loading score file...')
    for line in score_file:
        if line.startswith('#'): continue
        tokens = line.strip().split('\t')
        var_info = tokens[:4]
        if var_info[0].startswith('chr'): var_info[0] = var_info[0][3:]
        score = float(tokens[score_field])
        scores_gen_gen[''.join(var_info)] = max(1-np.exp(-0.5*score), 0)
        scores_gen_con[''.join(var_info)] = max(1-np.exp(-0.5*score), 0)
        scores_con_gen[''.join(var_info)] = max(1-np.exp(-0.21*score), 0)
        scores_con_con[''.join(var_info)] = max(1-np.exp(-0.21*score), 0)
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
            con_score = np.exp(-1/(np.abs(float(tokens[5]))+0.001))
            gen_score = 1-np.exp(-0.5*np.abs(float(tokens[5])))
            scores_gen_gen[''.join(var_info)] = max(scores_gen_gen.get(''.join(var_info), 0), gen_score)
            scores_con_gen[''.join(var_info)] = max(scores_con_gen.get(''.join(var_info), 0), gen_score)
            scores_gen_con[''.join(var_info)] = max(scores_gen_con.get(''.join(var_info), 0), con_score)
            scores_con_con[''.join(var_info)] = max(scores_con_con.get(''.join(var_info), 0), con_score)
        score_file.close()
        logging.info('Done loading spidex score file.')
        
    #ages = []
    #with open(ages_path) as ages_file:
        #for line in ages_file:
            #ages.append(float(line.strip()))
    
    genders = []
    with open(gender_path) as gender_file:
        for line in gender_file:
            genders.append(line.strip())

    c_c = []
    with open(cc_path) as cc_file:
        for line in cc_file:
            c_c.append((line.strip() == 'COPD_Case'))

    #inds = [i for i in range(len(ages)) if ages[i] < 46 or ages[i] >= 49]
    inds = range(len(c_c))
    # Get the rest of our stuff
    # List of variant info lines
    variant_lines = []
    # List of variant weight lines based on whatever scoring file we are using
    ggwl = []
    gcwl = []
    cgwl = []
    ccwl = []
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
        # Assume that func.refgene is 27th info field, otherwise look for it manually
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
                insert_variant(chrom, pos, ref, allele, func, genos, gene, variant_lines, matrix_lines, ggwl, gcwl, cgwl, ccwl, genes, scores_gen_gen, scores_gen_con, scores_con_gen, scores_con_con, which_alt=i+1)
        else:
                insert_variant(chrom, pos, ref, alt, func, genos, gene, variant_lines, matrix_lines, ggwl, gcwl, cgwl, ccwl, genes, scores_gen_gen, scores_gen_con, scores_con_gen, scores_con_con)

        vcf_line_count += 1
        if vcf_line_count % 100000 == 0:
            logging.info('Finished parsing %d variants' % vcf_line_count)
            logging.info('Currently %d genes in dictionary' % len(genes.keys()))
    

    assert len(c_c) == len(cols)
    logging.info('Writing to files')

    # Create out folder if necessary
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
    
    # Dump stuff to file
    gen_gen_variant_file = open(out_folder + '/variants_gen_gen.txt', 'w')
    gen_con_variant_file = open(out_folder + '/variants_gen_con.txt', 'w')
    con_gen_variant_file = open(out_folder + '/variants_con_gen.txt', 'w')
    con_con_variant_file = open(out_folder + '/variants_con_con.txt', 'w')
    gene_file = open(out_folder + '/genes.txt', 'w')
    patient_file = open(out_folder + '/patients.txt', 'w')
    exome_file = open(out_folder + '/exome.txt', 'w')

    # Write variants
    for gene in sorted(genes.keys()):
        for i, var_ind in enumerate(genes[gene]):
            info = '_'.join(variant_lines[var_ind][:4])
            annotation = variant_lines[var_ind][4]
            ggscore = ggwl[var_ind]
            gcscore = gcwl[var_ind]
            cgscore = cgwl[var_ind]
            ccscore = ccwl[var_ind]
            maf = float(sum([int(x) for x in matrix_lines[i] if x == '1' or x == '2'])) / (2.0*len(cols))
            gene_index = str(i+1)
            gen_gen_variant_file.write('\t'.join([gene, info, annotation, str(ggscore), str(maf), gene_index]) + '\n')
            gen_con_variant_file.write('\t'.join([gene, info, annotation, str(gcscore), str(maf), gene_index]) + '\n')
            con_gen_variant_file.write('\t'.join([gene, info, annotation, str(cgscore), str(maf), gene_index]) + '\n')
            con_con_variant_file.write('\t'.join([gene, info, annotation, str(ccscore), str(maf), gene_index]) + '\n')

    # Write genes
    for i, gene in enumerate(sorted(genes.keys())):
        gene_file.write('\t'.join([gene, str(i+1), str(len(genes[gene]))]) + '\n')

    # Write patients
    for i in inds:
        patient_id = cols[i]
        patient_index = i+1
        #age = int(ages[i])
        pheno = c_c[i]
        gen = genders[i]
        patient_file.write('\t'.join([patient_id, str(patient_index), str(int(pheno)), gen]) + '\n')

    # Write exomes
    for i, var in enumerate(variant_lines):
        info = '_'.join(var[:4])
        for ind in inds:
            if matrix_lines[i][ind] == '1' or matrix_lines[i][ind] == '2':
                exome_file.write('\t'.join([info, cols[ind], matrix_lines[i][ind]]) + '\n')
                
    gen_gen_variant_file.close()
    gen_con_variant_file.close()
    con_gen_variant_file.close()
    con_con_variant_file.close()
    gene_file.close()
    patient_file.close()
    exome_file.close()
    
    logging.info('Done writing!')

def insert_variant(chrom, pos, ref, alt, func, genos, gene, variant_lines, matrix_lines, ggwl, gcwl, cgwl, ccwl, genes, scores_gen_gen, scores_gen_con, scores_con_gen, scores_con_con, which_alt=1):
        if not ''.join([chrom, pos, ref, alt]) in scores_gen_gen: return
        ggwl.append(scores_gen_gen[''.join([chrom, pos, ref, alt])])
        gcwl.append(scores_gen_con[''.join([chrom, pos, ref, alt])])
        cgwl.append(scores_con_gen[''.join([chrom, pos, ref, alt])])
        ccwl.append(scores_con_con[''.join([chrom, pos, ref, alt])])
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
    #parser.add_argument('ages_path', metavar='AGES')
    parser.add_argument('gender_path', metavar='GENDERS')
    parser.add_argument('cc_path', metavar='CASE_OR_CONTROL')
    parser.add_argument('out_folder', metavar='OUT')
    parser.add_argument('--score_field', metavar='FIELD', help='Integer (0-indexed) indicating which field in the tab separated file contains the score we are'
            'interested in. Default is 4 for cadd. Use 5 for spidex z-scores', type=int, default=4)
    parser.add_argument('--pool_spidex', help='Spidex score file to optionally use for pooling with cadd scores')
    return parser.parse_args(args)

def main(args = sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
