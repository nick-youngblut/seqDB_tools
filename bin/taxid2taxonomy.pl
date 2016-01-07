#!/usr/bin/env perl
use strict;
use warnings;

=pod

=head1 NAME

taxid2taxonomy.pl -- get taxonomies for a list of taxids.

=head1 SYNOPSIS

taxid2taxonomy.pl [options] < taxids.txt > taxonomies.txt 

=head2 Options

=over

=item <taxids.txt> 

A list of taxids (taxonomy IDs; 1 per line)

=item -attempts  <int>

Number of attempts at retrieving taxonomy info.
[default: 5]

=item -procs  <int>

Number of parallel processes. 
[default: 1]

=item -help

=back

=head1 DESCRIPTION

Convert taxids (taxonomy IDs) to taxonomies
(phylum, class, order, etc.).

BioPerl is used to query entrez for the taxonomies. 

The input can be a multi-column txt file, as long
as the first column is the taxid.

The output is a *.txt file: [taxid lineage]

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
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use Time::Local;
use List::Util qw/max/;
use List::MoreUtils qw/any/;
use Bio::DB::Taxonomy;
use Parallel::ForkManager;
use File::Path qw/rmtree/;


#--- I/O error ---#

pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my $procs = 1;
my $attempts = 5;
GetOptions(
	   "procs=i" => \$procs,
	   "attempts=i" => \$attempts,
	   "help|?" => \&pod2usage
	   );

#--- MAIN ---#

# loading taxids and calling entrez
my $pm = Parallel::ForkManager->new($procs);

$pm->run_on_finish(
                   sub {
                     my ($pid, $exit_code, $ident,
                         $exit_signal, $core_dump, $tax_r) = @_;
		     print join("\t", @$tax_r), "\n";                     
                     }
                   );

# header
my @levels = qw/superkingdom phylum class order family genus species/;
print join("\t", 'taxid', @levels), "\n";

while(<>){
  chomp;
  next if /^\s*$/;
  
  my @line = split /\t/;
  my $taxid = $line[0];
  
  $pm->start and next;
  my $tries = 0;
  while (1){
    $tries++;

    # get fastas & write
    my $tax_r = get_lineage($taxid, \@levels);
    
    # print if successful
    if ($tax_r){
      unshift(@$tax_r, $taxid);
      $pm->finish(0, $tax_r);  
      last;
    } 

    # too many tries to connect to DB
    if ($tries >= $attempts){
      warn "Exceeded $attempts tries to get taxonomy for $taxid\n";
      $pm->finish;      
      last;
    }
    else{
      warn "Failed query (retrying): $taxid\n";
      next;
    }
    $pm->finish;
  }
}
$pm->wait_all_children;
exit;

  
sub get_lineage{
# if lineage needed #
  my $id = shift;
  my $levels = shift || die "Provide array of taxonomy levels";
  my @levels = @$levels;
    
  my %lineage;
  my $db = new Bio::DB::Taxonomy(-source => 'entrez');
  my $node= $db->get_Taxonomy_Node(-taxonid => $id);

  return 0 unless $node;
  
  # get basal node #
  my $anc = $node;
  while(1){
    last unless $anc->ancestor();
    $anc = $anc->ancestor();
    my $rank = $anc->rank;
    if( grep(/$rank/i, @levels) ){
      $lineage{$rank} = $anc->scientific_name;
    }
  }
  
  my @desc = $db->get_all_Descendents($node);
  
  for my $child ( @desc ) {
    my $rank = $anc->rank;
    if( grep(/$rank/i, @levels) ){
      $lineage{$rank} = $anc->scientific_name;
    }
  }
  
  # checking for each level #
  map{ $lineage{$_} = "NA" unless $lineage{$_} } @levels;
  
  # return 
  return [@lineage{@levels}];
}



