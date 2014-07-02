#!/usr/bin/env perl

=pod

=head1 NAME

accession2genbankInfo.pl -- Get Genbank info associated with accession numbers

=head1 SYNOPSIS

accession2genbankInfo.pl [options] < accessions.txt > genbank_info.txt

=head2 options

=over

=item -column 

Column containing the accession number (index by 1). [1]

=item -help			

This help message

=back

=head2 For more information:

perldoc accession2genbankInfo.pl

=head1 DESCRIPTION

Provide a list of accession numbers
(1 number per line; can be 1st column of table)
and get metadata info associated with
the accession number.

Taxonomic classification data can be obtained from accession number.

=head2 Key for column names

=over

=item * 'sub_' = info directly related to the accession number

=item * 'man_' = info related to the manuscript associated with 
the accession number

=back

=head2 Reqirements:

Bioperl -> use Bio::DB::GenBank

=head1 EXAMPLES

=head2 Basic usage

accession2genbankInfo.pl < accessions.txt > genbank_info.txt

=head2 Using output from blastxml_parse.pl

cat blast_result.xml | blastxml_parse.pl | cut -f 2 | accession2genbankInfo.pl

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut



### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::DB::GenBank;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $header_b);
my $column = 1; 
my $batch = 100000;
GetOptions(
	   "column=i" => \$column,
	   "x" => \$header_b,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
$column--;

### MAIN
while(1){
  last if eof;
	
  my %lines;
  for(0..($batch-1)){  	# loading select number of lines at a time into memory
    # loading line
    my $line = <>;
    if($.==1 and ! $header_b){
      print;
      next;
    }
    last unless defined $line;
    
    chomp $line;
    next if $line =~ /^\s*$/;
    
    # loading %@
    my @line = split /\t/, $line;
    
    if ( defined $line[$column] && $line[$column] ne 'NA' ){
      push @{$lines{'acc'}{$line[$column]}}, \@line;
    }
    else{
      push @{$lines{'blank'}{'blank'}}, \@line;
    }
  }
  
  # getting accession info
  my $acc_r = get_accession_info( \%lines );
  
  # writing lines 
  foreach my $l ( @{$lines{'blank'}{'blank'}} ){
    print join("\t", @$l), "\n";
  }
  foreach ( keys %{$lines{'acc'}} ){
    foreach my $l ( @{$lines{'acc'}{$_}} ){
      print join("\t", @$l), "\n";
    }
  }
}


sub get_accession_info{
  my ($lines_r) = @_;
		
  my $gb = new Bio::DB::GenBank;
  my $seqio = $gb->get_Stream_by_acc( [keys %{$lines_r->{'acc'}}] );
  
  my @ref_q = qw/title authors location pubmed/;
  
  my %acc;
  while( my $seq =  $seqio->next_seq ) {	
    # seqfeature #
    my @line = (
		$seq->accession_number,
		join(";", reverse($seq->species->classification()) ),
		$seq->authority,
		$seq->desc
	       );
    
    # annotations #
    my @annotations = $seq->annotation->get_Annotations('reference');
    foreach my $val (@annotations){
      my $hash_r = $val->hash_tree;
      foreach my $q (@ref_q){
	if(exists $hash_r->{$q}){
	  push @line, $hash_r->{$q};
	}
	else{
	  push @line, "";
	}
      }
    }
    
    # checking for undefined variables #
    map{$_ = "" unless $_} @line;
        
    # adding info to lines of table
    (my $acc = $line[0]) =~ s/$/.1/ unless /\.[0-9]+$/;
    foreach (@{$lines_r->{'acc'}{$acc}}){
      push @$_, @line;
    }
    
  }
  
  #print Dumper $lines_r; exit;
}


