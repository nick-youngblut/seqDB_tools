#!/usr/bin/env perl
use strict;
use warnings;

=pod

=head1 NAME

NCBIprokaryoteTableFilter.pl -- filter (& add taxonomy to) prokaryote.txt file for NCBI genome ftp site

=head1 VERSION

This is version 0.0.1

=head1 USAGE

    NCBIprokaryoteTableFilter.pl [options]

=head1 REQUIRED ARGUMENTS

=over

=item <prokaryotes.txt>

prokaryotes.txt file from NCBI genome ftp site. ('-' = STDIN)

=back

=head1 OPTIONS

=over

=item -w[[rite][_fastas]]

Write out fasta files for each organism? [FALSE]


=item -d[irectory] <directory>

Directory to write out all fastas (only needed with '-w').

Default: directory.default

=for Euclid:
directory.type: writable
directory.default: 'genomes'


=item -a[ttempts] <attempts>

Number of attempts at trying download a genome.

Default: attempts.default

=for Euclid:
attempts.type: int >= 1
attempts.default: 5


=item -t[hreads] <threads>

Number of threads to run jobs in parallel.

Default: threads.default

=for Euclid:
threads.type: int >= 0
threads.default: 1


=item --debug [<log_level>]

Set the log level. Default is log_level.default but if you provide --debug,
then it is log_level.opt_default.

=for Euclid:
    log_level.type:        int
    log_level.default:     0
    log_level.opt_default: 1

=item --version

=item --usage

=item --help

=item --man

Print the usual program information

=back

=head1 DESCRIPTION

This script helps get the prokaryotes.txt
table from NCBI genome ftp site in a usable
format. 

The filtered or unfiltered table and be used
to write out the fastas of all genomes in
the table that have accessions in the
'Chromosomes/INSDC' field.

=head2 Filtering:

* 'Status' field must contain 'complete' (caps-invariant)

* For entries with the same organism name: using the most recent version.

=head2 Taxonomy:

Using the 'TaxID' column to add taxonomy: from domain to strain.

=head1 EXAMPLES

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
use Pod::Usage;
use Data::Dumper;
use Getopt::Euclid;
use Time::Local;
use List::Util qw/max/;
use List::MoreUtils qw/any/;
use Bio::DB::Taxonomy;
use Parallel::ForkManager;
use File::Path qw/rmtree/;

#--- I/O error ---#

if (any {$_ eq '-h'} values %ARGV){
  print "Type '--help' for script info\n";
  exit;
}

#--- MAIN ---#
# just writing out genome fastas from accesion
if($ARGV{'-write_fastas'}){
  print STDERR "Writing out fastas of each genome",
    " (fastas from GenBank accession numbers)...\n";
  my $tbl_r = load_table_accession($ARGV{'<prokaryotes.txt>'});

  # output diretory
  my $outdir =  $ARGV{'-directory'};
  rmtree($outdir) if -d $outdir;
  mkdir $outdir or die $!;
  
  # forking & writing
  my $pm = Parallel::ForkManager->new($ARGV{-threads});
  my $tries = 0;
  foreach my $org (keys %$tbl_r){
    $pm->start and next;
    while (1){
      $tries++;

      # get fastas & write
      my $ret = write_fasta_from_accession($org, $tbl_r->{$org}, $outdir);	

      # success/fail
      if ($ret == 0){
	warn "Successful download: $org\n";
	last;
      }
      elsif ($tries >= $ARGV{'-attempts'}){
	warn "Exceeded tries to download: $org\n";
	last;
      }
      else{
	warn "Failed download (retrying): $org\n";
	next;
      }
    }
    $pm->finish;
  }
  $pm->wait_all_children;
  exit;
}
  
# filtering table
print STDERR "Loading table...\n";
my $tbl_r = load_table($ARGV{'<prokaryotes.txt>'});

print STDERR "Just keeping most recent entries for same organism...\n";
filter_by_most_recent($tbl_r);

print STDERR "Adding taxonomy to each entry...\n";
## edit header & write
edit_header($tbl_r->{header});

## forking
my $pm = Parallel::ForkManager->new($ARGV{-threads});
$pm->run_on_finish(
		   sub {
		     my ($pid, $exit_code, $ident, 
			 $exit_signal, $core_dump, $tbl_r) = @_;
		     write_table($tbl_r);
		     }
		   );

foreach my $name (keys %{$tbl_r->{body}}){
  $pm->start and next;
  my $tmp = add_taxonomy($tbl_r->{body}{$name}, 
			$tbl_r->{header}{'TaxID'});
  $pm->finish(0, $tmp);
}
$pm->wait_all_children;



#--- Subroutines ---#
sub write_fasta_from_accession{
  use Bio::SeqIO;
  use Bio::DB::GenBank;

  my ($org, $acc_r, $dir) = @_;
  print STDERR "Downloading $org\n";
  
  if($$acc_r[0] eq '-'){
    warn "'$org' does not have an accession. Skipping\n";
    return 0;
  }

  # fixing org name
  $org =~ s/'//g;
  $org =~ s/[. \/()=+\[\]:]+/_/g;
  $org =~ s/_$//;

  # making sure genome file starts empty
  open OUT, ">$dir/$org.fasta" or die $!;
  close OUT;
  
  # getting fasta and writing 
  my $gb = new Bio::DB::GenBank;
  my $ret_val = 0;
  foreach my $acc (@$acc_r){
    my $seqio = $gb->get_Stream_by_id( [$acc] );
    my $seq_out = Bio::SeqIO->new( -file => ">>$dir/$org.fasta", 
				   -format => 'fasta');
    
    while(my $seq = $seqio->next_seq){
      if(defined $seq->seq()){
	$seq_out->write_seq($seq);
      }
      else{
	warn $seq->id, "has no sequence. Skipping!\n";
	$ret_val = 1;
	next;
      }
    }
  }
  return $ret_val;
}

sub load_table_accession{
  my ($file) = @_;

  my $fh;
  $file eq '-' ? $fh = \*STDIN : 
    open $fh, $file or die $!;

  my %header;
  my %tbl;
  while(<$fh>){
    chomp;
    next if /^\s$/;

    my @l = split /\t/;

    # header
    unless( %header ){
      for my $i (0..$#l){
	$header{$l[$i]} = $i;
      }
      next;
    }
    
    # body
    my $chrom_col = exists $header{'Chromosomes/INSDC'} ? 
      $header{'Chromosomes/INSDC'} : 
	die "ERROR: 'Chromosomes/INSDC' field not found\n";
    my $org_col = exists $header{'#Organism/Name'} ? 
      $header{'#Organism/Name'} : 
	die "ERROR: '#Organism/Name' field not found\n";

    $tbl{$l[$org_col]} = [split /\s*,\s*/, $l[$chrom_col]];
  }

 #print Dumper %tbl; exit;
  return \%tbl;
}

sub edit_header{
# in: hash_ref
  my $header_r = shift;

  my @levels = qw/superkingdom phylum class order family genus species/;

  # adding taxonomic classifications to the header
  my $max_col = max( values %$header_r );
  for my $i (0..$#levels){
    $header_r->{$levels[$i]} = $max_col + $i + 1;
  }

  #writing header
  print join("\t", sort{$header_r->{$a} <=>
			  $header_r->{$b}} 
	     keys %$header_r), "\n";
}

sub write_table{
  my ($tbl_r) = @_;

  foreach my $Uid ( keys %$tbl_r ){
    print join("\t", @{$tbl_r->{$Uid}}), "\n";
  }
}

sub add_taxonomy{
  my ($tbl_r, $taxid_col) = @_;

  # getting taxonomy for each taxid
  foreach my $Uid ( keys %$tbl_r ){
    my $id = ${$tbl_r->{$Uid}}[$taxid_col];
    
    my $lineage_r = get_lineage($id);
    
    push @{$tbl_r->{$Uid}}, @$lineage_r;            
  }

return $tbl_r;
}

sub get_lineage{
# if lineage needed #
  my $id = shift;
  
  my @levels = qw/superkingdom phylum class order family genus species/;
  
  my %lineage;
  my $db = new Bio::DB::Taxonomy(-source => 'entrez');
  my $node= $db->get_Taxonomy_Node(-taxonid => $id);
  
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


sub filter_by_most_recent{
# just keeping the most recent assembly
  my ($tbl_r) = @_;

  foreach my $name (keys %{$tbl_r->{body}}){
    next if scalar keys %{$tbl_r->{body}{$name}} == 1; 

    # getting submission dates 
    my %dates;
    foreach my $Uid ( keys %{$tbl_r->{body}{$name}} ){
      my $rel_date_col = $tbl_r->{header}{'Release Date'};
          
      my $rel_date =  ${$tbl_r->{body}{$name}{$Uid}}[$rel_date_col];
      $rel_date =~ s|/||g;
      next unless $rel_date =~ /^\d+$/;
      $dates{$Uid} = $rel_date;
    }

    # deleting all but most recent
    foreach my $Uid (sort{$dates{$a} cmp $dates{$b}} keys %dates){
      delete $tbl_r->{body}{$name}{$Uid} 
	unless scalar keys %{$tbl_r->{body}{$name}} == 1;
    }
  }

}

sub load_table{
# loading prokaryotes.txt table
  my $file = shift;

  my $fh;
  $file eq '-' ? $fh = \*STDIN : 
    open $fh, $file or die $!;  
  
  my %tbl;
  while(<$fh>){
    chomp;
    next if /^\s*$/;
    
    my @l = split /\t/;

    # header
    unless( $tbl{header} ){
      for my $i (0..$#l){
	$tbl{header}{$l[$i]} = $i;
      }
      next;
    }

    # filtering
    ## must be complete genome
    my $status_col = $tbl{header}{Status};
    next if $l[$status_col] !~ /Complete/i;

    # body
    my $name_col = $tbl{header}{'#Organism/Name'};
    $tbl{body}{$l[$name_col]}{$.} = \@l;

  }
  
  close $fh or die $!;

  return \%tbl;
}

