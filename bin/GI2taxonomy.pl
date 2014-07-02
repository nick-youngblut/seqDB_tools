#!/usr/bin/env perl

=pod

=head1 NAME

GI2taxonomy.pl -- adding NCBI taxonomy to a table based on a column if GI numbers

=head1 SYNOPSIS

GI2taxonomy.pl [flags] < table.txt > table_wTaxIDs.txt

=head2 Required flags

=over

=item names

names.dmp file from NCBI (see blast2lca documentation).
[/var/seq_data/ncbi_db/taxonomy/names.dmp]

=item nodes

nodes.dmp file from NCBI (see blast2lca documentation)
[/var/seq_data/ncbi_db/taxonomy/nodes.dmp]

=item dict

(gi_taxid_nucl.bin|gi_taxid_prot.bin) file from NCBI (see blast2lca documentation)
[/var/seq_data/ncbi_db/taxonomy/gi_taxid_nucl.bin]

=back

=head2 Optional flags

=over

=item -column  <int>

Column number containing GIs (indexed by 1). [1]

=item -verbose  <bool>

Verbose output. [TRUE]

=item -h  <bool>

Print this help message & exit. [FALSE]

=back

=head2 For more information:

perldoc GI2taxonomy.pl

=head1 DESCRIPTION


It uses Bio::LITE::Taxonomy::NCBI 
and requires the following input files: 

=over

=item * a nodes.dmp file

=item * a names.dmp file

=item * a dict file (gi_taxid_nucl.bin|gi_taxid_prot.bin)

=back

No table header allowed!

Taxonomic levels appended to the table:
root,domain,phylum,class,order,family,genus,species,strain.
'NA' used if taxnomic level not found.

=head1 EXAMPLES

=head2 Basic usage:

GI2taxonomy.pl -names names.dmp -node nodes.dmp -dict gi_taxid_nucl.bin
< tbl.txt > tbl_wtaxonomy.txt

=head1 AUTHOR

Nick Youngblut <ndy2@cornell.edu>

=head1 AVAILABILITY

email me

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


#--- modules ---#
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use Bio::LITE::Taxonomy;
use Bio::LITE::Taxonomy::NCBI;

#--- args/flags ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my $verbose_b; 
my $names_loc = "/var/seq_data/ncbi_db/taxonomy/names.dmp";
my $nodes_loc = "/var/seq_data/ncbi_db/taxonomy/nodes.dmp";
my $dict_loc = "/var/seq_data/ncbi_db/taxonomy/gi_taxid_nucl.bin";

my $column = 1;
GetOptions(
	   "names=s" => \$names_loc,
	   "nodes=s" => \$nodes_loc,
	   "dict=s" => \$dict_loc,
	   "column=i" => \$column,
	   "verbose" => \$verbose_b,
	   "help|?" => \&pod2usage # Help
	  );

#--- I/O error ---#
$column--;

#--- MAIN ---#
# loading names & nodes files; making taxDB
print STDERR "Making taxonomy DB from names & nodes files...\n";
my $taxDB = Bio::LITE::Taxonomy::NCBI->new (
                                            db=>"NCBI",
                                            names=> $names_loc,
                                            nodes=>$nodes_loc,
                                            dict=>$dict_loc
                                           );                                        

# loading blast; getting taxonomy; filtering
print STDERR "Loading input table & finding taxonomy...\n";

my %n_filt = (  root => 0,
		domain => 0,
		env_sample => 0,
		n_levels => 0,
		classification => 0 );

print STDERR "Getting taxonomy from GI numbers...\n";
my @column_chk;
while(<>){
  chomp;
  next if /^\s*$/;
    
  # parsing 
  my @l = split /\t/;
  die "ERROR: cannot find column '$column' in line $.!\n"
    unless defined $l[$column];

  # sanity check
  if(defined $column_chk[0]){
    die "Input ERROR: line $. does not contain $column_chk[0] columns as the 1st line does!\n"
      unless $column_chk[0] == @l;
  }
  else{ $column_chk[0] = @l; }
  
  
  # getting GI
  (my $GI = $l[$column]) =~ s/^[^|]+\|//;	# [1] = gi number
  my $tax = $taxDB->get_taxonomy_from_gi($GI);	# @$[taxonomy]

  
  # appending to table
  if(! $tax){
    push @l, ('NA') x 9;
  }
  else{
    for my $i (0..8){
      if(defined $$tax[$i]){
	push @l, $$tax[$i];
      }
      else{
	push @l, 'NA';
      }
    }
  }

  # sanity check
  if(defined $column_chk[1]){
    die "Taxonomy addition ERROR: line $. does not contain $column_chk[1] columns as the 1st line does!\n"
      unless $column_chk[1] == @l;
  }
  else{ $column_chk[1] = @l; }
  
  # writing output
  print join("\t", @l), "\n";
}

