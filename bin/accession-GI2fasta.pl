#!/usr/bin/env perl

=pod

=head1 NAME

accession-GI2fasta.pl -- get fasta from NCBI accession numbers

=head1 SYNOPSIS

accession-GI2fasta.pl [options] < accessions.txt > sequences.fasta

=head2 Options

=over

=item -acc_col <int>

Column containing accessions. [1]

=item -name_col <int>

Column containing taxon names (if not using accession IDs).

=item -batch  <int>

Number of sequences to try and download at a time. [10]

=item -trial  <int>

Number of trials to attempt for downloading a batch. [5]

=item -desc  <bool>

Add sequence description to name in fasta? [TRUE]

=item -header <bool>

Header in accessions.txt file? [FALSE]

=item -fork  <int>

Number of batches to process in parallel. [1]

=item -output <str>

Output directory. [.]

=item -h	This help message

=back

=head2 For more information:

perldoc accession-GI2fasta.pl

=head1 DESCRIPTION

A wrapper for Bio::DB::GenBank.

Downloading is not always successful, so
sequences are downloaded in batches. If
not all sequences in the batch are downloaded
successfully, the batch is re-tried. 

Any accesssion lacking a sequence will be skipped.

=head2 Input

Provide a table of accession numbers (accessions.txt).
Accession numbers separated by commas will be 
written to the same fasta file.

A name column can be specified, which will 
be used to name all of the output fasta files.

Accession numbers, versioned accession numbers,
or GI numbers can be used.

=head2 Output

Each entry (row in the table) will be written
to a separate file. The file will be named
by the accession number(s) for the entry,
or if a name column is provided, it will be
named by the name column value.

=head3 Note

Existing & non-empty sequence files
of the same name will not be over-written.
You must delete old versions sequences.

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

my $acc_col = 1;
my $name_col ;
my $trial_limit = 5;
my $batch_size = 10;
my $threads = 0;
my $add_desc;
my $header;
my $outDir = '.';
GetOptions(
	   "acc_col=i" => \$acc_col,
	   "name_col=i" => \$name_col,
	   "trial=i" => \$trial_limit,
	   "batch=i" => \$batch_size,
	   "fork=i" => \$threads,
	   "header" => \$header,
	   "desc" => \$add_desc,   # add desc to id? [TRUE]
	   "output=s" => \$outDir,
	   "help|?" => \&pod2usage # Help
	  );

#--- I/O error ---#
$acc_col--;
$name_col--;
$outDir = File::Spec->rel2abs($outDir);
mkdir $outDir unless -d $outDir;


#--- MAIN ---#
my $acc_r = load_accessions($batch_size, $acc_col, $name_col, $header);
my $gb = new Bio::DB::GenBank;
my $pm = Parallel::ForkManager->new($threads);

#my $batch_cnt = 0;
foreach my $batch_cnt (sort keys %$acc_r){
  my $batch_r = $acc_r->{$batch_cnt};

  $pm->start and next;

  my $trials = 0;
  while(1){
    # status #
    $trials++;
    print STDERR "Starting batch->trial: $batch_cnt->$trials\n";

    # streaming in sequences
    my $cnt = 0;
    my $skip = 0;
    foreach my $taxon (keys %$batch_r){
      my $outFile = File::Spec->catfile($outDir, "$taxon.fna");
      if (-f $outFile and ! -z $outFile){
	# file exists & is non-empty
	$cnt++;
      }
      else{
	my $seqio = $gb->get_Stream_by_id($batch_r->{$taxon});
	my $seq_out = Bio::SeqIO->new( -file => ">$outFile", 
				       -format => 'fasta');      
	my ($cnt_x, $skip_x) = write_seqs($seqio, $seq_out, $add_desc);
	$cnt += $cnt_x;
	$skip += $skip_x;
      }
    }
      
    # checking to make sure all have been written
    if($cnt >= scalar(keys %$batch_r) - $skip){ # all sequences appear to have been written out
      last;
    }
    elsif($trials >= $trial_limit){
      print STDERR "WARNING: exceded number of trials. Skipping\n";
      last;
    }
    else{      
      print STDERR "WARNING: only $cnt of ", 
	scalar keys %$batch_r, 
	  " sequences were written. Retrying\n";
    }   
  }

  $pm->finish;
}
$pm->wait_all_children;


#--- Subroutines ---#
sub write_seqs{
  my $seqio = shift || die $!;
  my $seq_out = shift || die $!;
  my $add_desc = shift;

  # writing out sequences
  my $cnt = 0; # keeping track of number of successful writes 
  my $skip = 0;
  while (my $seq = $seqio->next_seq){
    $seq->id( join("|", $seq->id, $seq->desc)) unless $add_desc;

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
  return $cnt, $skip;
}


sub load_accessions{
  my $batch_size = shift;
  my $acc_col = shift;
  my $name_col = shift;
  my $header = shift;
  
  my %acc;
  my %tmp;
  my $batch = 0;
  while(<>){
    if (defined $header and $. == 1){
      next;
    }
    chomp;
    next if /^\s*$/;

    my $acc;
    my $name;
    my @line = split /\t/;

    # accession
    if (defined $line[$acc_col]){
      $acc = [split /\s*,\s*/, $line[$acc_col]]
    }
    else{ die "Error: accession column not found!"; }

    # name
    if (defined $name_col){
      if (defined $line[$name_col]){
	$name = $line[$name_col];
	$name =~ s/\s/_/g;
      } 
      else{ die "Error: name column not found!"; }
    }
    else{ $name = $acc; }
    

    #push @tmp, @tmp2;  #$line[0];
    
    $acc{$batch}{$name} = $acc;

    if($. % $batch_size == 0 || eof){
      $batch++;
    }
    
  }
  #print Dumper %acc; exit;
  return \%acc;
}
