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

=item -db <db> | -database <db>

SRAmetadb file (from SRAdb R package)

=head1 OPTIONS

=over

=item -sql <sql>

sql statement to get sample accessions. 
Some form of the sample accession must be returned by query.
Only the 1st column returned by the query will be used.

Default: sql.default

=for Euclid:
sql.type: string
sql.default: 'SELECT distinct(sample_accession) FROM sra_ft WHERE sra_ft MATCH "*soil*metagenom*"'


=item -nt <n_threads> | -num_threads <n_threads>

Number of threads to use (0 = no forking).

Default: n_threads.default

=for Euclid:
n_threads.type: int >= 0
n_threads.default: 0


=item --keep

Keep 'fastq_files' table if exists and append entries to the existing table.

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

Add/update a 'fastq_files' table
in the SRAmetaDb from R package SRAdb.

This table contains info on fastq
files associated with SRA biosamples.

The info is collected from the ENA at: 
http://www.ebi.ac.uk/ena/data/view/reports/sra/fastq_file/


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
use DBI;
use Parallel::ForkManager;
use LWP::Simple;
use Term::ProgressBar;

#--- I/O error ---#


#--- MAIN ---#
# db connect
my $dbh = connect2db( $ARGV{-db} );

# getting all sample accessions from db
my $samp_acc_r = get_db_samp_acc($dbh, \%ARGV);

# getting read ftp
my $tbl_r = get_fastq_ftp($samp_acc_r, \%ARGV);

# adding fastq_file entries to database
## getting list of tables
my @tables = $dbh->tables() ;
#map{ s/.+\."|"//g } @tables;

## creating new fastq_files table in db
### if table doesn't exist or user-defined overwrite
if(! (grep /"main"\."fastq_files"/, @tables)
   or ! $ARGV{'--keep'} ){

  print STDERR "Creating fastq_files table\n";
  create_fastq_files_table($dbh, \%ARGV);
}


## addng fastq_files entries to db
add_fastq_files_entries($dbh, $tbl_r, \%ARGV);


# closing 
$dbh->commit or die $dbh->err;
$dbh->disconnect or die $dbh->err;




#--- Subroutines ---#

=head2 add_fastq_files_entries

adding entires to fastq_files

=cut

sub add_fastq_files_entries{
  my $dbh = shift or die "Provide dbh\n";
  my $tbl_r = shift or die "Provide tbl\n";
  my $argv_r = shift;

  # getting table columns
  my $sth = $dbh->prepare("SELECT * from fastq_files");
  my @columns = @{$sth->{NAME}};
  my %columns;
  map{ $columns{$columns[$_]} = $_ } 0..$#columns;
  
  # preparing sql
  my $sql = join(" ", 
		 "INSERT INTO fastq_files(",
		 join(",", sort keys %columns),
		 ") VALUES(",
		 join(",", ('?') x scalar keys %columns),
		 ")"
		);
  

  $sth = $dbh->prepare($sql) or die $dbh->err;

  # adding entries
  my $add_cnt = 0;
  foreach my $samp_acc (keys %$tbl_r){
    foreach my $entry_id (keys %{$tbl_r->{$samp_acc}}){

      # adding sample accession and entry id to entry
      $tbl_r->{$samp_acc}{$entry_id}{sample_accession} = $samp_acc;
      $tbl_r->{$samp_acc}{$entry_id}{entry_id} = $entry_id;

      # binding params to sth
      my $col_cnt = 0;
      foreach my $column (sort keys %columns){
	$col_cnt++;
	my $val = exists $tbl_r->{$samp_acc}{$entry_id}{$column} ?
	  $tbl_r->{$samp_acc}{$entry_id}{$column} : undef;
       
	$sth->bind_param($col_cnt, $val) or die $dbh->err;
      }

      # calling with bound params
      $sth->execute() or die $dbh->err;
      $add_cnt++;
    }
  }


  # status
  printf STDERR "Number of entries added to fastq_files table: $add_cnt\n";
    
}


=head2 create_fastq_files_table

creating fastq_files table in sqlite3 database

=cut

sub create_fastq_files_table{
  my $dbh = shift or die "Provide dbh\n";
  my $argv_r = shift;

  # status
  print STDERR "Creating a new 'fastq_files' table\n"
    unless $argv_r->{'--quiet'};

  # dropping table if exists
  my $sql = "DROP TABLE IF EXISTS fastq_files";
  $dbh->do($sql) or die $dbh->err;
  

  # creating new table
  $sql = <<HERE;
CREATE TABLE
fastq_files
(
sample_accession	TEXT,
entry_id                NUM,
analysis		TEXT,
file_size		NUM,
library_source		TEXT,
ru			TEXT,
instrument_model	TEXT,
sample			TEXT,
library_selection	TEXT,
file_name		TEXT,
organism		TEXT,
library_name		TEXT,
instrument_platform	TEXT,
run_base_count		NUM,
run_read_count		NUM,
experiment		TEXT,
library_strategy	TEXT,
study			TEXT,
library_layout		TEXT,
md5			TEXT,
ftp			TEXT
)
HERE

  $dbh->do($sql) or die $dbh->err;


  # creating index
  $sql = <<HERE;
CREATE INDEX fastq_files_sample_acc_idx 
ON fastq_files (sample_accession)
HERE
  
  $dbh->do($sql) or die $dbh->err;


  # commit
  $dbh->commit() or die $dbh->err;
}

=head2 get_fastq_ftp

Getting the ftp links associated with a biosample

=cut

sub get_fastq_ftp{
  my $samp_acc_r = shift or die "Provide samp_acc_r\n";
  my $argv_r = shift;

  # status
  print STDERR "Getting info from ENA\n" unless $argv_r->{'--quiet'};

  # forking initilize
  my $pm = Parallel::ForkManager->new($argv_r->{-num_threads});

  my %tbl;
  $pm->run_on_finish(
		     sub{
		       my ($pid, $exit_code, 
			   $ident, $exit_signal, 
			   $core_dump, $ret_r) = @_;
		      
		       foreach my $samp_acc (keys %$ret_r){			 
			 $tbl{$samp_acc} = $ret_r->{$samp_acc};
		       }
		     }
		    );

  # status
  my $prog = Term::ProgressBar->new(scalar @$samp_acc_r);
  $prog->minor(0);

  my ($cnt, $next_update) = (0,0);  # status
  foreach my $samp_acc (@$samp_acc_r){
    # status
    $cnt++;
    $next_update = $prog->update($cnt) if $cnt > $next_update;
#    print STDERR "Number of sample accessions processed: $cnt\n"
#      if $cnt % 100 == 0 and ! $argv_r->{'--quiet'};    
   # last if $cnt >= 2000;   # debug

    # forking start
    $pm->start and next;

    my $site = join("/", "http://www.ebi.ac.uk/ena/data/view/reports/sra/fastq_files",
		    $samp_acc->[0]);
    my $text = get $site; 
    open IN, '<', \$text or die $!;

    my @header;
    my %ret;
    while(<IN>){
      chomp;
      next if /^\s*$/;
      my @l = split /\t/;
      
      # header
      if($. == 1){
	@header = @l;
	map{ tr/A-Z /a-z_/ } @header;
      }
      # body: adding values to entries_r
      else{
	for my $i (0..$#l){	  
	  next unless defined $header[$i];
	  $ret{$samp_acc->[0]}{$.-1}{$header[$i]} = $l[$i];
	}	
      }
    }
    close IN or die $!;
    
    # forking end
    $pm->finish(0, \%ret );
  }
  $pm->wait_all_children;


  # return
#  print Dumper %tbl; exit;
  return \%tbl;
}


=head2 get_db_samp_acc

Getting sample accessions from SRAmetadb

=cut

sub get_db_samp_acc{
  my $dbh = shift or die "Provide dbh\n";
  my $argv_r = shift;


  my $sql = $argv_r->{-sql};

  # debug
#  $sql .= " limit 200";

  # query
  print STDERR "sql: '$sql'\n\n" unless $argv_r->{'--quiet'};
  my $ret = $dbh->selectall_arrayref($sql) or die $dbh->err;

  # status
  printf STDERR "Number of sample accessions to query: %i\n\n",
    scalar @$ret;


  return $ret;
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


