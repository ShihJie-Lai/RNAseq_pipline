#!/usr/bin/perl -w
use Cwd;
use strict;
no strict 'refs';
use Encode qw/encode decode/;

require "/home/shihjielai/RNAseq_PIP/enrichment_analysis2.pl";
use Getopt::Long qw(GetOptions);

&enrich_ANA();

