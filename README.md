TEcov, a pipeline to to determine the repeat content of genes and their neighboring regions.
=========

Requires a gene annotation in ".gff" format and a [RepeatMasker](http://www.repeatmasker.org/)
repeat annotation in ".out" format. Reports the repeat content (% of base pairs covered) of

- genic regions (exons, introns, untranslated regions)
- promoter regions (upstream regions)
- flanking regions (up- and downstream regions)

The results can be used to study the influence of repetetive content such as
transposable elements (TEs) on gene expression or gene family evolution.


Requirements
------------

TEcov requires [Python3.4+](https://www.python.org/downloads/) and the Python package pyfaidx:

`pip3 install pyfaidx`


Usage
------------

```
usage: TEcov.py [-h] [-f FLANK_SIZE] [-j JOBNAME] [-m MIN_BASEPAIRS]
                [--gene_gff_feat GENE_GFF_FEAT] [--cds_gff_feat CDS_GFF_FEAT]
                species_gff te_annotation genome

A pipeline to determine the repeat content (e.g. content of transposable
elements) of genes and their neighboring regions. Requires a gene annotation
and a repeatmasker repeat annotation. Reports the repeat content (%bp covered)
of genic regions (exons, introns, UTRs), promoter (upstream) regions and
flanking (up- and downstream) regions.

positional arguments:
  species_gff           path to the species' genome annotation (.gff, only
                        longest isoforms!)
  te_annotation         path to the species' TE annotation (repeatmasker .out
                        but NOT .gff)
  genome                path to the species' genome (.fasta).

optional arguments:
  -h, --help            show this help message and exit
  -f FLANK_SIZE, --flank_size FLANK_SIZE
                        how long the flanking regions should be (in bp,
                        default is 10.000 bp)
  -j JOBNAME, --jobname JOBNAME
                        name a directory to which all intermediate files will
                        be dumped
  -m MIN_BASEPAIRS, --min_basepairs MIN_BASEPAIRS
                        percentages get meaningless when the overlapped
                        sequences are very short. This value defines the
                        minimun sequence length (in bp). Shorter sequences
                        will be "NA" (default: 200 bp)
  --gene_gff_feat GENE_GFF_FEAT
                        the name of the GFF feature that denotes the location
                        of genes (default: "gene")
  --cds_gff_feat CDS_GFF_FEAT
                        the name of the GFF feature that denotes the location
                        of CDSs (default: "cds")
```
