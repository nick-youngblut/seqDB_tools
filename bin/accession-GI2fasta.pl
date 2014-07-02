#!/usr/bin/env perl

=pod

=head1 NAME

accession-GI2fasta.pl -- get fasta from NCBI accession numbers

=head1 SYNOPSIS

accession-GI2fasta.pl [flags] < accessions.txt > sequences.fasta

=head2 Required flags

=over

NONE

=back

=head2 Optional flags

=over

=item -batch  <int>

Number of sequences to try and download at a time. [10]

=item -trial  <int>

Number of trials to attempt for downloading a batch. [5]

=item -desc  <bool>

Add sequence description to name in fasta? [TRUE]

=item -fork  <int>

Number of batches to process in parallel. [1]

=item -h	This help message

=back

=head2 For more information:

perldoc accession-GI2fasta.pl

=head1 DESCRIPTION

A wrapper for Bio::DB::GenBank.

Provide a list of accession numbers
(1st column if  tab-delimited table)

Accession numbers, versioned accession numbers,
or GI numbers can be used.

Downloading is not always successful, so
sequences are downloaded in batches. If
not all sequences in the batch are downloaded
successfully, the batch is re-tried. 

Any accesssion lacking a sequence will be skipped.

=head1 EXAMPLES

=head2 Accession:

echo "CAB02640" | accession-GI2fasta.pl > seqs.fasta

=head2 GI:

echo "405830" | accession-GI2fasta.pl > seqs.fasta

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/

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
use File::Spec;
use Bio::SeqIO;
use Bio::DB::GenBank;
use File::Temp qw/ tempfile /;
use Parallel::ForkManager;

#--- args/flags ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my $trial_limit = 5;
my $batch_size = 10;
my $threads = 0;
my $add_desc;
GetOptions(
	   "trial=i" => \$trial_limit,
	   "batch=i" => \$batch_size,
	   "fork=i" => \$threads,
	   "desc" => \$add_desc,   # add desc to id? [TRUE]
	   "help|?" => \&pod2usage # Help
	  );

#--- I/O error ---#


#--- MAIN ---#
my $acc_r = load_accessions($batch_size);
my $gb = new Bio::DB::GenBank;
my $pm = Parallel::ForkManager->new($threads);

my $batch_cnt = 0;
foreach my $batch_r (@$acc_r){
  # status #
  $batch_cnt++;
  print STDERR "Starting batch: $batch_cnt\n";

  $pm->start and next;

  my $trials = 0;
  while(1){
    # status #
    $trials++;
    print STDERR "Starting trial: $trials\n";

    # making temp file for writing
    my $tmp = File::Temp->new();
    my $filename = $tmp->filename;

    # streaming in sequences
    my $seqio = $gb->get_Stream_by_id( $batch_r );
    my $seq_out = Bio::SeqIO->new( -file => ">$filename", -format => 'fasta');
  
    # writing out sequences
    my $cnt = 0; # keeping track of number of successful writes 
    my $skip = 0;
    while( my $seq = $seqio->next_seq ){
      $seq->id( join("|", $seq->id, $seq->desc) ) unless $add_desc;

      if(defined $seq->seq()){
	$seq_out->write_seq($seq);
	$cnt++;
      }
      else{
	print STDERR "WARNING: '", $seq->id, "' has no sequence. Skipping!\n";
	$skip++;
	next;
      }
    }
    
    # checking to make sure all have been written
    if($cnt == scalar(@$batch_r) - $skip){ # all sequences appear to have been written out; writing from tempFile to STDOUT
      open IN, $filename or die $!;
      while(<IN>){ print; }
      close IN;
      last;
    }
    elsif( $trials >= $trial_limit){
      print STDERR "WARNING: exceded number of trials. Skipping\n";
      last;
    }
    else{      
      print STDERR "WARNING: only $cnt of ", scalar @$batch_r, " sequences were written. Retrying\n";
    }   
  }

  $pm->finish;
}
$pm->wait_all_children;


#--- Subroutines ---#
sub load_accessions{
  my $batch_size = shift;
  my @acc;
  my @tmp;
  while(<>){
    chomp;
    next if /^\s*$/;

    my @line = split /\t/;
    push @tmp, $line[0];

    if($. % $batch_size == 0 || eof){
      push @acc, [@tmp];      
      @tmp = ();
    }
    
  }
 # print Dumper @acc; exit;
  return \@acc;
}
