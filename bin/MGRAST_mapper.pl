#!/usr/bin/env perl
use strict;
use warnings;

=pod

=head1 NAME

MGRAST_mapper.pl -- Download MG-RAST metagenomes and map reads to query

=head1 VERSION

This is version 0.0.1

=head1 USAGE

    MGRAST_mapper.pl [options]

=head1 REQUIRED ARGUMENTS

=over

=item -id <ids>

File containing a list of MG-RAST metagenome ids (1 ID per line).
Use '-' if from STDIN.

=for Euclid:
ids.type: input

=item -x <bt2_idx>

Index filename prefix (minus trailing .X.bt2).

=for Euclid:
bt2_idx.type: string

=back

=head1 OPTIONS

=over

=item -stage <stage>...

Processing stage for downloading files.
If multiple stages provided, will try each
in succession until entries are returned.

Default: stage.default

=for Euclid:
stage.type: num
stage.default: [150,100]


=item -b <bt2_params> | -bt2_params <bt2_params>

Parameters passed to bowtie2
(besides input & output params).

Default: bt2_params.default

=for Euclid:
bt2_params.type: string
bt2_params.default: '--local'


=item -o <output> | -output <output>

Output report file name ('-' = STDOUT).

Default: output.default

=for Euclid:
output.type: string
output.default: '-'


=item -prefix <prefix>

Output file prefix.

Default: prefix.default

=for Euclid:
prefix.type: string
prefix.default: ''


=item --debug [<log_level>]

Set the log level. Default is log_level.default but if you provide --debug,
then it is log_level.opt_default.

=for Euclid:
    log_level.type:        int
    log_level.default:     0
    log_level.opt_default: 1

=item --quiet

=item --version

=item --usage

=item --help

=item --man

Print the usual program information

=back

=head1 DESCRIPTION

=head2 General workflow of the pipeline:

For each metagenome ID:

=over

=item * download metagenome as temporary file using MG-RAST API

=item * collect general stats on metagenome

=item * map metagenome genome to query (index) with bowtie2

=item * created filtered bam file of hits (just reads that hit)

=back

Finally: write report on mappings


=head2 -stage flag

Multiple stages can be provided as
a backup in case the first stage doesn't
return any entries (some stages are optional).

Stages can be found in Appendix A of the
MG-RAST manual
(ftp://ftp.metagenomics.anl.gov/data/manual/mg-rast-tech-report-v3_r1.pdf).

=head1 AUTHOR

Nick Youngblut (ndy2@cornell.edu)

=head1 BUGS

There are undoubtedly serious bugs lurking somewhere in this code.
Bug reports and other feedback are most welcome.

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

#--- modules ---#
use Data::Dumper;
use Getopt::Euclid;
use JSON;
use LWP::UserAgent;
use URI::Escape;
use Data::Dumper;

use FindBin;
use lib "$FindBin::RealBin/../lib/";
use seqDB_tools qw/MGRAST_mapper_main/;
  

#--- main ---#
MGRAST_mapper_main(\%ARGV);
