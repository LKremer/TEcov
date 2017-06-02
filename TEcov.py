#!/usr/bin/env python3


'''
This pipeline can determine the repeat content of genes and their neighboring regions using bedtools.
It is similar to the following bedtools workflow:

# extract the genes and CDSs from our gene annotation GFF so we can use them to find flanking regions and remove CDSs
grep -P "\tgene\t" Cryptotermes_secundus.gff > Csec_2.6_genes.gff
grep -P "\tCDS\t" Cryptotermes_secundus.gff > Csec_2.6_cds.gff

# we need to know the scaffold lengths to get flanking regions
[CUSTOM_SCRIPT] Csec_genome.fa > Csec_scaffold_lengths.txt

# get 10kb flanking regions around genes
bedtools flank -i Csec_2.6_genes.gff -g Csec_scaffold_lengths.txt -b 10000 > Csec_2.6_gene_10kb_flanking.gff

# remove CDSs from flanking regions (like Janousek et al.)
bedtools subtract -a Csec_2.6_gene_10kb_flanking.gff -b Csec_2.6_cds.gff > Csec_2.6_gene_10kb_flanking_noCDS.gff

# find % coverage of a specific kind of repeat in flanking regions, e.g. for "Mariner" TEs:
# step 1: extract Mariner TEs from repeatmasker.out and convert to GFF format
[CUSTOM_SCRIPT] Csec_repeatmasker.out > Csec_Mariner_TEs.gff
# step 2: calculate Mariner TE coverage of e.g. flanking regions:
bedtools coverage -a Csec_Mariner_TEs.gff -b Csec_2.6_gene_10kb_flanking_noCDS.gff
'''


import pyfaidx
import argparse
import os
import subprocess
import shlex
import re
import csv
import math


def weighted_mean(data, min_len):
    '''
    e.g. data2 = [
      (20, 300),
      (100, 200),
      (225, 150),
    ]
    returns 173.18840579710144
    '''
    total_len = 0
    running_total = 0
    for cov_len, cov_percent in data:
        running_total += (cov_len * cov_percent) 
        total_len += cov_len
    if total_len < min_len:
        return 'NA'
    try:
        return (running_total / total_len)
    except ZeroDivisionError:
        return 0.0


def sanitize_header_for_SQL(header_rowdict):
    'remove special characters from the header to allow csv import into SQL'
    for fieldname, row_header in header_rowdict.items():
        clean_row_header = row_header.replace('%', '').replace('/', '_').replace('-', '_').replace('10kb_', '')
        header_rowdict[fieldname] = clean_row_header
    return header_rowdict


class Main(object):
    def __init__(self, species_gff, te_annotation, genome, flank_size,
                 jobname, min_basepairs, gene_gff_feat, cds_gff_feat):
        self.gene_gff_feat = gene_gff_feat.lower()
        self.cds_gff_feat = cds_gff_feat.lower()
        self.species_gff = species_gff
        self.te_annotation = te_annotation
        self.genome = genome
        self.flank_size_kb = str(flank_size / 1000)
        if self.flank_size_kb.endswith('.0'):
            self.flank_size_kb = self.flank_size_kb[:-2]
        self.flank_size_bp = flank_size
        self.jobname = jobname
        self.min_basepairs = min_basepairs
        self.te_gffs = {}
        os.mkdir(jobname)
        return


    def get_scaffold_length_file(self):
        '''
        Write a tab-separated file with the name and length of each scaffold.
        Required for 'bedtools flank'.
        '''
        self.scaffold_lengths = '{}/scaffold_lengths.tsv'.format(self.jobname)
        print('\nWriting the scaffold lengths of {} to '
              '{}'.format(self.genome, self.scaffold_lengths))
        with open(self.scaffold_lengths, 'w') as outf:
            scaffold_fa_entries = pyfaidx.Fasta(self.genome)  # raises PermissionError if no writing permissions
            for scaffold_name in scaffold_fa_entries.keys():
                scaffold_size_in_bp = len(scaffold_fa_entries[scaffold_name])
                outf.write('{}\t{}\n'.format(scaffold_name, scaffold_size_in_bp))
        return self.scaffold_lengths


    def split_gff_into_genes_and_cds(self):
        '''
        equivalent to grepping all "CDS" and all "gene" lines
        of the GFF into two separate GFF files.
        '''
        self.gene_gff = '{jobname}/genes.gff'.format(**self.__dict__)
        self.cds_gff = '{jobname}/cds.gff'.format(**self.__dict__)
        print('\nSplitting {species_gff} into {gene_gff} and {cds_gff}'.format(**self.__dict__))
        with open(self.species_gff, 'r') as in_gff, open(self.gene_gff, 'w') as g_gff, open(self.cds_gff, 'w') as c_gff:
            for line in in_gff:
                if '\t{}\t'.format(self.gene_gff_feat) in line.lower():
                    g_gff.write(line)
                elif '\t{}\t'.format(self.cds_gff_feat) in line.lower():
                    c_gff.write(line)
        #self.cds_gff = self.merge_isoform_cds(self.cds_gff)
        return self.gene_gff, self.cds_gff


    def determine_te_subfamilies(self):
        '''
        parses the repeatmasker.out to find out which TE-classes
        (e.g. LINE, SINE, DNA) and which families (e.g. gypsy, Pao,
        Ginger) there are
        '''
        print('\nDetermining the TE classes and families')
        te_classes = set()
        te_families = set()
        with open(self.te_annotation, 'r') as repeatmasker_out:
            in_header = True
            for row in repeatmasker_out:
                if in_header:
                    if not row.strip():
                        in_header = False  # reached an empty line, header is over!
                    continue
                te_class_family = row.split()[10]  # the TE class/family, e.g. LINE/Jockey
                if '/' in te_class_family:
                    te_class, te_family = te_class_family.split('/')
                    te_classes.add(te_class)
                te_families.add(te_class_family)

        te_families -= te_classes
        print('Found {} TE classes ({}) and {} TE families'.format(
            len(te_classes), ', '.join(te_classes), len(te_families)
        ))
        self.TE_classes = te_classes
        self.TE_families = te_families
        return
        

    def split_te_annotation(self):
        ''' turns repeatmasker .out file into a different GFF for each TE class '''
        files = {}  # dict containing output file handles for each TE class
        print('\nSplitting {} into GFF files with different TE classes:'.format(self.te_annotation))
        try:
            repeatmasker_out = open(self.te_annotation, 'r')

            te_categories = self.TE_classes | self.TE_families
            for TE_type in te_categories | {'other', 'all'}:
                file_name = '{}/{}-TEs.gff'.format(self.jobname, TE_type.replace('/', '-'))
                print('\t{:.<40}({})'.format(file_name, TE_type))
                file_object = open(file_name, 'w')
                files[TE_type] = file_object
                self.te_gffs[TE_type] = file_name

            # iterate over repeatmasker out and skip the header
            in_header = True
            for row in repeatmasker_out:
                if in_header:
                    if not row.strip():
                        in_header = False  # reached an empty line, header is over!
                    continue
                values = row.strip().split()
                percent = values[1]
                scaffold = values[4]
                start = values[5]
                stop = values[6]
                te_motif = values[9]  # which specific TE pattern was matched, e.g. "LTR_retrotransposon_1275"
                te_class_family = values[10]  # the TE class/family, e.g. LINE/Jockey
                # decide to which class a TE belongs
                if '/' in te_class_family:
                    # handling cases like 'LINE/Jockey'
                    te_class = te_class_family.split('/')[0]
                    te_family = te_class_family
                    assert te_class in self.TE_classes
                    assert te_family in self.TE_families
                else:
                    if te_class_family in self.TE_classes:
                        # handling cases like 'LINE'
                        te_class = te_class_family
                        assert te_class in self.TE_classes
                        te_family = 'uncategorized'
                    else:
                        # handling cases like 'mariner'
                        te_class = 'other'
                        te_family = te_class_family
                        assert te_family in self.TE_families
                comment = 'class={};motif={};broadclass={}\n'.format(
                    te_class_family, te_motif, te_class
                )
                gff_row = [scaffold, 'RepeatMasker', 'similarity',
                           start, stop, percent, '.', '.', comment]
                # write GFF line to the respective GFF file(s):
                files[te_class].write('\t'.join(gff_row))
                if te_family != 'uncategorized':
                    files[te_family].write('\t'.join(gff_row))
                files['all'].write('\t'.join(gff_row))  # one GFF that contains all TEs
        # make sure all files are closed since I cannot use "with" here :(
        finally:
            repeatmasker_out.close()
            for TE_class, file_object in files.items():
                file_object.close()
                assert os.stat(self.te_gffs[TE_class]).st_size > 0, 'The ' \
                    'gff for {} TEs is empty: {}'.format(
                        TE_class, self.te_gffs[TE_class])


    def get_flanking_regions(self):
        ''' bedtools flank -i csec_2.6_genes.gff -g csec_scaffold_lengths.txt -b 10000 > csec_2.6_gene_10kb_flanking.gff '''
        self.flank_with_cds = '{jobname}/{flank_size_kb}kb_flanking.gff'.format(**self.__dict__)
        cmd = 'bedtools flank -i {gene_gff} -g {scaffold_lengths} -b {flank_size_bp}'.format(**self.__dict__)
        print('\nGetting the {}kb flanking regions around genes with command\n\t{}'.format(self.flank_size_kb, cmd))
        with open(self.flank_with_cds, 'w') as flank_gff:
            proc = subprocess.Popen(shlex.split(cmd), stdout=flank_gff)
            proc.communicate()
        return self.flank_with_cds


    def remove_cds_from_flanking_regions(self):
        ''' bedtools subtract -a csec_2.6_gene_10kb_flanking.gff -b csec_2.6_cds.gff > csec_2.6_gene_10kb_flanking_noCDS.gff '''
        self.flank_no_cds = '{jobname}/{flank_size_kb}kb_flanking_no_cds.gff'.format(**self.__dict__)
        cmd = 'bedtools subtract -a {flank_with_cds} -b {cds_gff}'.format(**self.__dict__)
        print('\nRemoving CDS from flanking regions with command\n\t', cmd)
        with open(self.flank_no_cds, 'w') as flank_gff:
            proc = subprocess.Popen(shlex.split(cmd), stdout=flank_gff)
            proc.communicate()
        return self.flank_no_cds


    def get_promoter_regions(self):
        ''' bedtools flank -s -l 2000 -r 0 -i genes.gff -g scaffold_lengths.tsv > promotors.gff '''
        self.promotor_gff = '{jobname}/promotors.gff'.format(**self.__dict__)
        cmd = 'bedtools flank -s -l 2000 -r 0 -i {gene_gff} -g {scaffold_lengths}'.format(**self.__dict__)
        print('\nGetting the promotor regions (2k upstream) of genes with command\n\t{}'.format(cmd))
        with open(self.promotor_gff, 'w') as prom_gff:
            proc = subprocess.Popen(shlex.split(cmd), stdout=prom_gff)
            proc.communicate()
        return self.promotor_gff


    def get_te_coverages(self):
        ''' bedtools coverage -a csec_Mariner_TEs.gff -b csec_2.6_gene_10kb_flanking_noCDS.gff '''
        self.te_coverage_gffs = {}
        query_gffs = {
            '{}kb-flanks'.format(self.flank_size_kb) : self.flank_no_cds,
            'gene' : self.gene_gff,
            'promotor' : self.promotor_gff,
            #'cds' : self.cds_gff,  # THIS IS TOO HARD I GIVE UP
        }
        for gene_subset, gene_subset_gff in query_gffs.items():
            self.te_coverage_gffs[gene_subset] = {}
            for TE_class, TE_gff in self.te_gffs.items():
                coverage_fname = '{}/{}_{}-TE_coverage.gff'.format(
                    self.jobname, gene_subset, TE_class.replace('/', '-'))
                cmd = 'bedtools coverage -a {} -b {}'.format(TE_gff, gene_subset_gff)
                print('\nGetting {}-TE coverage of {} with command\n\t{} '
                      '> {}'.format(TE_class, gene_subset, cmd, coverage_fname))
                with open(coverage_fname, 'w') as cov_gff:
                    proc = subprocess.Popen(shlex.split(cmd), stdout=cov_gff)
                    proc.communicate()
                self.te_coverage_gffs[gene_subset][TE_class] = coverage_fname
        return self.te_coverage_gffs


    def make_summary_table(self):
        ''' parse all coverage gffs to get a % value for each TE '''
        print('\nGenerating TE% output summary table...')
        gene_id_pattern = re.compile(r'(?<=\w\=)[a-zA-Z]+\_[a-zA-Z]*\d+')
        gene2cov = {}

        csv_headers = set()
        for genic_region, dictionary in self.te_coverage_gffs.items():
            for te_class, coverage_file in dictionary.items():
                csv_header = '{}_{}-TE%'.format(genic_region, te_class)
                csv_headers.add(csv_header)
                print('Parsing bedtools coverage file {}'.format(coverage_file))
                with open(coverage_file, 'r') as cov_gff:
                    for line in cov_gff:
                        values = line.split('\t')
                        ids = values[8]
                        te_percent = float(values[12])
                        gene_id = re.search(gene_id_pattern, ids).group(0)
                        start, stop = int(values[3]), int(values[4])
                        feature_len = stop - start
                        if start >= stop:  # meaningless bedtools artifact?!
                            te_percent = 0.0
                            feature_len = 0

                        if gene_id not in gene2cov:
                            gene2cov[gene_id] = {}
                        if csv_header not in gene2cov[gene_id]:
                            gene2cov[gene_id][csv_header] = []
                        gene2cov[gene_id][csv_header].append(
                            (feature_len, te_percent)
                        )
        # now gene2cov has an entry for each gene that looks like:
        '''
        'Csec_G03224': {'10kb-flanks_DNA-TE%': [(9999, 0.0), (9999, 0.0)],
                        '10kb-flanks_LINE-TE%': [(9999, 0.4764), (9999, 0.3883)],
                        '10kb-flanks_LTR-TE%': [(9999, 0.0), (9999, 0.0)],
                        '10kb-flanks_all-TE%': [(9999, 0.572), (9999, 0.454)],
                        '10kb-flanks_other-TE%': [(9999, 0.2123), (9999, 0.2231)],
                        'gene_DNA-TE%': [(407, 0.0)],
                        'gene_LINE-TE%': [(407, 0.0)],
                        'gene_LTR-TE%': [(407, 0.0)],
                        'gene_all-TE%': [(407, 0.1813726)],
                        'gene_other-TE%': [(407, 0.1813726)]},
        '''

        out_tbl = '{}_summary.csv'.format(self.jobname)
        with open(out_tbl, 'w') as csv_obj:
            fieldnames = ['gene_id'] + sorted(csv_headers)
            writer = csv.DictWriter(csv_obj, fieldnames=fieldnames,
                                    restval='NA', delimiter='\t')
            header = dict(zip(writer.fieldnames, writer.fieldnames))
            writer.writerow(sanitize_header_for_SQL(header))
            #writer.writeheader()
            for gene_id, cov_dict in gene2cov.items():
                for key, percent_list in cov_dict.items():
                    if key.startswith('gene'):
                        min_bp = 0
                    else:  # shorter sequences will get 'NA' as TE percentage
                        min_bp = self.min_basepairs
                    cov_dict[key] = weighted_mean(percent_list, min_len=min_bp)
                cov_dict['gene_id'] = gene_id
                writer.writerow(cov_dict)
        print('\nWrote TE% summary table to', out_tbl)






def main(species_gff, te_annotation, genome, flank_size, jobname,
         min_basepairs, gene_gff_feat, cds_gff_feat):
    print('''   - Input Options -

   species_gff: {}
 te_annotation: {}
        genome: {}
    flank_size: {} bp
       jobname: {}
 min_basepairs: {} bp
 gene_gff_feat: {}
  cds_gff_feat: {}
    '''.format(species_gff, te_annotation, genome, flank_size, jobname,
               min_basepairs, gene_gff_feat, cds_gff_feat))

    m = Main(species_gff, te_annotation, genome, flank_size, jobname,
             min_basepairs, gene_gff_feat, cds_gff_feat)

    # we need to know the scaffold lengths to get flanking regions
    m.get_scaffold_length_file()

    # extract the genes and CDSs from the Csec gff so we can use them to find flanking regions and remove CDSs
    m.split_gff_into_genes_and_cds()

    # parse the TE annotation to determine which kinds of TEs there are
    m.determine_te_subfamilies()

    # split the TE annotation into one GFF for each TE type
    m.split_te_annotation()

    # get 10 kb flanking regions around genes
    m.get_flanking_regions()

    # remove CDSs from flanking regions (like Janousek et al)
    m.remove_cds_from_flanking_regions()

    # get 2kb upstream regions of genes (promotors)
    m.get_promoter_regions()

    # calculate the TE% of flanking regions and gene body for all genes
    m.get_te_coverages()

    # report TE coverage percentages for each gene in a summary table
    m.make_summary_table()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='A pipeline to determine the repeat content (e.g. content '
        'of transposable elements) of genes and their neighboring regions. '
        'Requires a gene annotation and a repeatmasker repeat annotation. '
        'Reports the repeat content (%bp covered) of genic regions (exons, '
        'introns, UTRs), promoter (upstream) regions and flanking (up- and '
        'downstream) regions.')
    parser.add_argument('species_gff', help='path to the species\' '
                        'genome annotation (.gff, only longest isoforms!)')
    parser.add_argument('te_annotation', help='path to the species\' '
                        'TE annotation (repeatmasker .out but NOT .gff)')
    parser.add_argument('genome', help='path to the species\' '
                        'genome (.fasta).')
    parser.add_argument('-f', '--flank_size', help='how long the '
                        'flanking regions should be (in bp, '
                        'default is 10.000 bp)', default=10000, type=int)
    parser.add_argument('-j', '--jobname', help='name a directory to '
                        'which all intermediate files will be dumped',
                        default='TEcov')
    parser.add_argument('-m', '--min_basepairs', help='percentages get '
                        'meaningless when the overlapped sequences are '
                        'very short. This value defines the minimun seq'
                        'uence length (in bp). Shorter sequences will '
                        'be "NA" (default: 200 bp)', default=200, type=int)
    parser.add_argument('--gene_gff_feat', default='gene', help='the name '
                        'of the GFF feature that denotes the location of '
                        'genes (default: "gene")')
    parser.add_argument('--cds_gff_feat', default='cds', help='the name '
                        'of the GFF feature that denotes the location of '
                        'CDSs (default: "cds")')
    args = parser.parse_args()
    main(**vars(args))


# TODO: 
# - add CDS TE%
#    needs a preparatory step where CDS are intersected with genes using bedtools
# - normalize by genome-wide TE content
