#!/usr/bin/env perl
use strict;
use warnings;

=pod

=head1 NAME

entrez_get.pl -- wrapper for esearch & efetch

=head1 VERSION

This is version 0.0.1

=head1 USAGE

    entrez_get.pl [options]

=head1 REQUIRED ARGUMENTS

=item -db <db>

Entrez database.

=for Euclid:
db.type: string

=item -query <query>

Query statement.

=for Euclid:
query.type: string


=head1 OPTIONS

=over

=item -SRA[meta] <SRAmeta>

SRA metadata sqlite3 DB file from SRAdb R package.

Default: SRAmeta.default

=for Euclid:
SRAmeta.type: input
SRAmeta.default: "/var/seq_data/ncbi_db/SRA/SRAmetadb.sqlite"


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

The script uses the Entrez EUtilities scripts:
esearch and efetch to get metadata.

As of now, only the SRA database is supported.

ftp links for fastq files associated with the
samples are written along with the metadata.

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
use IPC::Cmd qw/can_run run/;
use Text::ParseWords;
use DBI;
use LWP::Simple;
use List::MoreUtils qw/all/;

#--- I/O error ---#
can_run('esearch') or die "ERROR: esearch no in \$PATH\n";
can_run('efetch') or die "ERROR: efetch no in \$PATH\n";


#--- MAIN ---#
# db connect
my $dbh = connect2db( $ARGV{-SRAmeta} );

# parsing entries
my $fh_pipe = call_esearch_efetch( \%ARGV );
my $entries_r = parse_efetch_text($fh_pipe, \%ARGV);


# determining how many are in SRA meta sqlite3 file
#count_entries($dbh, $entries_r, \%ARGV);
filter_entries($dbh, $entries_r, \%ARGV);

# getting read ftp
get_fastq_ftp($dbh, $entries_r, \%ARGV);
#get_sra_ftp($dbh, $sample_acc_r); 

# writing
## as csv
write_entries_csv($entries_r, \%ARGV);


# closing pipe
close $fh_pipe;
$dbh->disconnect or die "ERROR: could not disconnect from db\n";



#--- Subroutines ---#

=head2 write_entries_csv

writing entries as csv file

=cut

sub write_entries_csv{
  my $entries_r = shift or die "Provide entries\n";
  my $argv_r = shift;
  
  # making header
  ## getting header values
  my %header;
  foreach my $entry_id (keys %$entries_r){
    # getting keys and adding to header    
    foreach my $key ( keys %{$entries_r->{$entry_id}} ){
      foreach my $cat ( keys %{$entries_r->{$entry_id}{$key}} ){
	$header{$key}{$cat} = 1;
      }
    }
  }
  ## writing out header values
  my @header;
  foreach my $key (sort keys %header){
    foreach my $cat (sort keys %{$header{$key}}){
      push @header, $cat;
    }
  }
  print join("\t", 'entry_id', @header), "\n";


  # writing body
  foreach my $entry_id (keys %$entries_r){
    my @row;
    push @row, $entry_id;

    # loading values from each category 
    foreach my $key ( sort keys %header ){
      foreach my $cat ( sort keys %{$header{$key}} ){
	if(exists $entries_r->{$entry_id}{$key}{$cat}){
	  push @row, $entries_r->{$entry_id}{$key}{$cat};
	}
	else{
	  push @row, 'NA';
	}
      }
    }
    # printing line of output
    print join("\t", @row), "\n";
  }
}


=head2 get_fastq_ftp

Getting the ftp links associated with a biosample

=head3 OUT

edited entries hash

=cut

sub get_fastq_ftp{
  my $dbh = shift or die "Provide dbh\n";
  my $entries_r = shift or die "Provide entries\n";
  my $argv_r = shift;

  # foreach entry, using sample_accession (SRA) to get fastq file ftp table
  foreach my $entry_id (keys %$entries_r){

    # sample accession
    my $samp_acc;
    if( exists $entries_r->{$entry_id}{Identifiers}{SRA} ){
      $samp_acc = $entries_r->{$entry_id}{Identifiers}{SRA};
    }
    else{
      print STDERR "WARNING: no sample_accession for entry $entry_id\n";
      next;
    }

    # getting table from ENA
   # $samp_acc = "SRS450530";   # debug
    my $site = join("/", "http://www.ebi.ac.uk/ena/data/view/reports/sra/fastq_files",
		    $samp_acc);
    print STDERR "Getting: $site\n";
    my $text = get $site; 
    open  IN, '<', \$text or die $!;

    # parsing table
    my @header;
    while(<IN>){
      chomp;
      next if /^\s*$/;
      my @l = split /\t/;
      
      # header
      if($. == 1){
	@header = @l;
      }
      # body: adding values to entries_r
      else{
	for my $i (0..$#l){
	 # die "Logic error: entry_id $entry_id : line $.\n'$_'\n"
	 #   unless defined $header[$i];

	  next unless defined $header[$i];
	  $entries_r->{$entry_id}{fastq_files}{ $header[$i] } = $l[$i];	  
	}	
      }
    }

    # debug
#    print Dumper $entries_r->{$entry_id}; exit;

   # last if $argv_r->{'--debug'};
  }
}



=head2 filter_entries

filtering entries: only entries found in SRAmetaDB kept

=cut

sub filter_entries{
  my $dbh = shift or die "Provide dbh\n";
  my $entries_r = shift or die "Provide entries\n";


  # getting all SRA sample accessions
  my @SRA;
  foreach my $entry_id (keys %$entries_r){
    push @SRA, $entries_r->{$entry_id}{Identifiers}{SRA} if
      exists $entries_r->{$entry_id}{Identifiers}{SRA};
  }

  # preparing query
  ## quote
  map{ $_ = $dbh->quote($_) } @SRA; 

  ## sql
  my $cmd = "SELECT sample_accession from experiment ".
    "where sample_accession IN (".
      join(",", @SRA).
	")";

  # querying
  my $ret = $dbh->selectall_hashref($cmd, 'sample_accession') or die $dbh->err;

  
  # filtering entries
  ## status: pre-filter
  printf STDERR "Number of entries from query: %i\n",
    scalar keys %$entries_r;
  ## delete entries
  foreach my $entry_id (keys %$entries_r){
    unless(exists $entries_r->{$entry_id}{Identifiers}{SRA} and
	   exists $ret->{$entries_r->{$entry_id}{Identifiers}{SRA}} ){      

      delete $entries_r->{$entry_id};
    }
  }
  ## status: pre-filter
  printf STDERR "Number of entries with hits in SRAmetaDB: %i\n",
    scalar keys %$entries_r;
  
}


sub count_entries{
  my $dbh = shift or die "Provide dbh\n";
  my $entries_r = shift or die "Provide entries\n";


  # getting all SRA sample accessions
  my @SRA;
  foreach my $entry_id (keys %$entries_r){
    push @SRA, $entries_r->{$entry_id}{Identifiers}{SRA} if
      exists $entries_r->{$entry_id}{Identifiers}{SRA};
  }

  # preparing query
  ## quote
  map{ $_ = $dbh->quote($_) } @SRA; 

  ## sql
  my $cmd = "SELECT count(*) from experiment ".
    "where sample_accession IN (".
      join(",", @SRA).
	")";

  # querying
  my $ret = $dbh->selectall_arrayref($cmd) or die $dbh->err;

  # status
  printf STDERR "Number of entries from query: %i\n",
    scalar keys %$entries_r;
  printf STDERR "Number of entries with hits in SRAmetaDB: %i\n",
    $ret->[0][0];

#  return \@SRA;
}


=head2 parse_efetch_text

Parsing the output from efetch.

Assuming the output is in the text report format.

=cut

sub parse_efetch_text{
  my $fh_pipe = shift or die "Provide pipe\n";
  my $argv_r = shift;

  my %entries;
  my $sample_id = 0;
  while (my $line = <$fh_pipe>){
    chomp $line;

    # start of entry
    if($line =~ /^(\d+): (.+)/){
      $sample_id++;
      $entries{$sample_id}{Title}{Title} = $2;

      # status
      print STDERR "Processed $sample_id samples\n"
	if $sample_id % 100 == 0  
	  and not exists $argv_r->{'--quiet'};
    }
    # identifiers (nested)
    elsif($line =~ /^Identifiers: (.+)/){
      my %words = quotewords(": ", 0, split / *; */, $1);
      map{ $entries{$sample_id}{Identifiers}{$_} = $words{$_} } keys %words;
    }
    # other categories
    elsif($line =~ /^\S+: *.+/){
      my @tmp = quotewords(": ", 0, $line);      

      my %words;
      for(my $i=0; $i<=$#tmp; $i+=2){   # balanced number of values
        $words{$tmp[$i]} =$tmp[$i+1] if defined $tmp[$i+2];
      }

      map{ $entries{$sample_id}{$_} = $words{$_} } keys %words;
    }
    # attributes (nothing needed)
    elsif($line =~ /^Attributes/){
      next;
    }
    # attribute lines
    elsif($line =~ /^\s+\/(.+)/){
      my %words = quotewords("=", 0, $1);
      map{ $entries{$sample_id}{Attributes}{$_} = $words{$_} } keys %words;      
    }
    # skip blank lines and any others
    else{
      next;
    }

    # debug
    last if $argv_r->{'--debug'} and scalar keys %entries >= 300;
  }

  #print Dumper %entries; exit;

  return \%entries;
}


=head2 call_esearch_efetch

esearch | efetch  pipeline for querying

returning pipe file handle

=cut

sub call_esearch_efetch{
  my $argv_r = shift or die "Provide argv\n";
  
  my $cmd = join(" ", 
		 "esearch", 
		 "-db", $argv_r->{-db},
		 "-query", $argv_r->{-query},
		 "|",
		 "efetch |"
		 );


  open my $fh_pipe, $cmd or die $!;
  return $fh_pipe;
}


=head2 connect2db

connecting to sqlite db

=cut 

sub connect2db{
# connecting to CLdb
# $db_file = database file
  my $db_file = shift or die "Provide db file name\n";
  
  my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
  my $dbh = DBI->connect("dbi:SQLite:dbname=$db_file", '','', \%attr)
                or die " Can't connect to $db_file!\n";
  
  return $dbh;
}



=head2 get_sra_ftp

=cut

sub get_sra_ftp{
  my $dbh = shift or die "Provide dbh\n";
  my $sample_acc_r = shift or die "Provide entries\n";

  ## sql
  my $cmd = "SELECT sample_accession, experiment_accession, run_accession ".
    "from sra ".
      "where sample_accession IN (".
	join(",", @$sample_acc_r).
	  ")";
  
  # querying
  my $ret = $dbh->selectall_hashref($cmd, 'sample_accession') or die $dbh->err;

  # making ftp
  foreach my $samp_acc (keys %$ret){
    $ret->{$samp_acc}{ftp} = join("/", 
				'ftp://ftp-trace.ncbi.nlm.nih.gov',
				'sra',
				'sra-instant',
				'reads',
				'ByExp',
				'sra',
				substr( $ret->{$samp_acc}{experiment_accession}, 0, 3 ),
				substr( $ret->{$samp_acc}{experiment_accession}, 0, 6 ),
				$ret->{$samp_acc}{experiment_accession},			       
				$ret->{$samp_acc}{run_accession},
				$ret->{$samp_acc}{run_accession} . ".sra"
				);
				
  }
}
