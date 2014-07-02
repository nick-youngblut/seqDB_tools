#!/usr/bin/env perl
use strict;
use warnings;

=pod

=head1 NAME

SRAmetadb_add_efetch_res.pl -- add table containing results from esearch & efetch query with sample accessions

=head1 VERSION

This is version 0.0.1

=head1 USAGE

    SRAmetadb_add_efetch_res.pl [options]

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
sql.default: 'SELECT distinct(sample_accession) FROM sra_ft WHERE sra_ft MATCH "*soil*metagenom*" and sample_accession != "Multiplex"'


=item -nt <n_threads> | -num_threads <n_threads>

Number of threads to use (0 = no forking).

Default: n_threads.default

=for Euclid:
n_threads.type: int >= 0
n_threads.default: 0


=item -ma <max_attributes> | -max_attributes <max_attributes>

Max number of attributes allowed in an entry (skipped otherwise). 
This limits the number of fields created in the efetch table.

Default: max_attributes.default

=for Euclid:
max_attributes.type: int >= 0
max_attributes.default: 15


=item --lat_long_edit

Try to get latitude and longitude values into decimal format
and into their corresponding fields.

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

Add/update a 'efetch' table
in the SRAmetaDb from R package SRAdb.

This table contains all info returned from
queries with esearch and efetch with
SRA biosample accessions.

The accessions are pulled from the SRAmetadb
database using the '-sql' flag.

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
use Term::ProgressBar;
use Text::ParseWords;
use Regexp::Common;
use Text::Unidecode;

#--- I/O error ---#


#--- MAIN ---#
# db connect
my $dbh = connect2db( $ARGV{-db} );

# getting all sample accessions from db
my $samp_acc_r = get_db_samp_acc($dbh, \%ARGV);

# calling esearch/efetch
my $tbl_r = esearch_efetch($samp_acc_r, \%ARGV);


# adding efetch entries to database
## getting list of tables
my @tables = $dbh->tables() ;
map{ s/.+\."|"//g } @tables;

## creating new efetch table in db
### if table doesn't exist or user-defined overwrite
if(! (grep /"main"\."efetch"/, @tables)
   or ! $ARGV{'--keep'} ){

  ## getting list of all attributes
  my $fields_r = get_all_needed_fields($tbl_r, \%ARGV);

  ## creating efetch table and deleting old one
  create_efetch_table($dbh, $fields_r, \%ARGV);
}


## addng efetch entries to db
add_efetch_entries($dbh, $tbl_r, \%ARGV);


# closing 
$dbh->commit or die $dbh->err;
$dbh->disconnect or die $dbh->err;



#--- Subroutines ---#

=head2 add_efetch_entries

adding entires to efetch table

=cut

sub add_efetch_entries{
  my $dbh = shift or die "Provide dbh\n";
  my $tbl_r = shift or die "Provide tbl\n";
  my $argv_r = shift;

  # status
  print STDERR "Adding entries to efetch table\n"
    unless $argv_r->{'--quiet'};
  

  # getting table columns
  my $sth = $dbh->prepare("SELECT * from efetch");
  my @columns = @{$sth->{NAME}};
  my %columns;
  map{ $columns{$columns[$_]} = $_ } 0..$#columns;
  

  # preparing sql
  my $sql = join(" ", 
		 "INSERT INTO efetch(",
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
  printf STDERR "Number of entries added to efetch table: $add_cnt\n";
}


=head2 create_efetch_table

creating fastq_files table in sqlite3 database

=cut

sub create_efetch_table{
  my $dbh = shift or die "Provide dbh\n";
  my $fields_r = shift or die "Provide fields_r hashref\n";
  my $argv_r = shift;

  # status
  print STDERR "Creating a new 'efetch' table\n"
    unless $argv_r->{'--quiet'};

  # dropping table if exists
  my $sql = "DROP TABLE IF EXISTS efetch";
  $dbh->do($sql) or die $dbh->err;
  

  # creating new table
  my @columns;
  foreach my $k (sort keys %$fields_r){
    push @columns, join("\t", $k, $fields_r->{$k});
  }

  $sql = join("\n", 
	      "CREATE TABLE efetch (",
	      join(",\n", @columns), 
	      ")");


  $dbh->do($sql) or die "ERROR: could not create table:\n$sql\n";


  # creating index
  $sql = <<HERE;
CREATE INDEX efetch_sample_acc_idx 
ON efetch (sample_accession)
HERE
  
  $dbh->do($sql) or die $dbh->err;


  # commit
  $dbh->commit() or die $dbh->err;
}



=head2 get_all_needed_fields

Getting all needed fields for creating efetch table.

Combining redundant fields.

=cut

sub get_all_needed_fields{
  my $tbl_r = shift or die "Provide tbl_r hashref\n";
  my $argv_r = shift;


  # getting all fields
  my %fields;
  foreach my $samp_acc (keys %$tbl_r){
    $fields{sample_accession} = 'TEXT';

    foreach my $entry_id (keys %{$tbl_r->{$samp_acc}}){
      $fields{entry_id} = 'NUM';

      
      # skipping if number of attributes exceeds max
      next if scalar keys %{$tbl_r->{$samp_acc}{$entry_id}}
	> $argv_r->{-max_attributes};

      # getting entry
      foreach my $cat (keys %{$tbl_r->{$samp_acc}{$entry_id}} ){
	next if $cat =~ /^\s*$/;
	unless(exists $fields{$cat}){
	  # is field text or numeric?
	  $fields{$cat} = $cat =~ /%RE{num}{real}/ ? 'NUM' : 'TEXT';
	}		
      }
    }
  }

  # filtering fields
  foreach my $field (keys %fields){
    delete $fields{$field} if
      $field =~ /^[atgc]+$/ or
	$field =~ /.+_tag_primer_sequences$/;	
  }


#  print Dumper %fields; exit;
#  print Dumper %$tbl_r; exit;
  return \%fields;
}


=head2 esearch_efetch

esearch | efetch  pipeline for querying

parallelizing querying.

=cut

sub esearch_efetch{
  my $samp_acc_r = shift or die "Provide sample_acc_r\n";
  my $argv_r = shift or die "Provide argv\n";

  # forking initialize
  my $pm = Parallel::ForkManager->new($argv_r->{-num_threads});

  ## run on finish
  my %tbl;
  $pm->run_on_finish(
                     sub{
                       my ($pid, $exit_code,
                           $ident, $exit_signal,
                           $core_dump, $ret_r) = @_;

		       if(defined $ret_r->[0] and defined $ret_r->[1]){			 
			 $tbl{ $ret_r->[0] } = $ret_r->[1];
		       }
                     }
                    );


  # status
  my $prog = Term::ProgressBar->new(scalar @$samp_acc_r);
  $prog->minor(0);
  my ($cnt, $next_update) = (0,0);


  ## calling entrez utilz
  foreach my $samp_acc (@$samp_acc_r){
    $cnt++;
    $next_update = $prog->update($cnt) if $cnt > $next_update;

    # sanity check
    next unless defined $samp_acc->[0];

    # forking
    $pm->start and next;

    # calling esearch/efetch
    my $entries_r = call_esearch_efetch($samp_acc->[0], $argv_r);

    # end fork
    $pm->finish(0, [$samp_acc->[0], $entries_r]);
  }
  $pm->wait_all_children;


  # return hash
  return \%tbl;
}


=head2 call_esearch_efetch

esearch | efetch  pipeline for querying

returning hashref of parsed output

=cut

sub call_esearch_efetch{
  my $samp_acc = shift or die "Provide samp_acc string\n";
  my $argv_r = shift or die "Provide argv\n";

  # query
  my $cmd = join(" ",
                 "esearch",
                 "-db biosample",
                 "-query $samp_acc",
                 "|",
                 "efetch |"
                 );

  open my $fh_pipe, $cmd or die $!;

  # parsing output
  my $entries_r = parse_efetch_text($fh_pipe, $argv_r);

  # adding sample accession to each entry
  map{ $entries_r->{$_}{sample_accession} = $samp_acc } keys %$entries_r;

 
  # close/return
  close $fh_pipe or die $!; 
  return $entries_r;
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
      $entries{$sample_id}{title} = $2;
    }
    # identifiers (nested)
    elsif($line =~ /^Identifiers: (.+)/){
      my %words = quotewords(": ", 0, split / *; */, $1);
      
      foreach my $k (keys %words){
	(my $cat = $k) =~ s/\W/_/g;
	$cat =~ tr/A-Z/a-z/;
	$entries{$sample_id}{$cat} = $words{$k};
      }
    }
    # other categories
    elsif($line =~ /^\S+: *.+/){
      my @tmp = quotewords(": ", 0, $line);

      my %words;
      for(my $i=0; $i<=$#tmp; $i+=2){   # balanced number of values
        $words{$tmp[$i]} =$tmp[$i+1] if defined $tmp[$i+2];
      }

      foreach my $k (keys %words){
	(my $cat = $k) =~ s/\W/_/g;
	$cat =~ tr/A-Z/a-z/;
	$entries{$sample_id}{$cat} = $words{$k};
      }
    }
    # attributes (nothing needed)
    elsif($line =~ /^Attributes/){
      next;
    }
    # attribute lines
    elsif($line =~ /^\s+\/(.+)/){
      my %words = quotewords("=", 0, $1);

      foreach my $k (keys %words){
	(my $cat = $k) =~ s/\W/_/g;
	$cat =~ tr/A-Z/a-z/;
	
	# combining certain attributes (reducing repetitive fields)
	($cat, $words{$k}) = edit_lat_long($cat, $words{$k})
	  if $argv_r->{'--lat_long_edit'};

	# adding to entries
	if(ref $cat eq "ARRAY"){
	  for my $i (0..$#$cat){
	    die "Internal error\n" unless defined $words{$k}->[$i];
	    $entries{$sample_id}{$cat->[$i]} = $words{$k}->[$i];
	  }
	}
	else{
	  $entries{$sample_id}{$cat} = $words{$k};
	}
      }
    }
    # skip blank lines and any others
    else{
      next;
    }

    # debug
    last if $argv_r->{'--debug'} and scalar keys %entries >= 300;
  }

#  print Dumper %entries; exit;
  return \%entries;
}


=head2 edit_lat_long

Try to properly format latitude and longitude

=cut

sub edit_lat_long{
  my $field = shift;
  my $value = shift;

  return $field, $value unless 
    defined $field and defined $value;

  # general formating of lat-long values
  if( $field =~ /latitude/ or $field =~ /longitude/
      or $field =~ /lat.+long*/ ){
    $value =~ s/ +[NSEW]//gi;  # removing any directional characters    
    $value =~ s/ +DD//;
    
    # unicode decode
    $value = unidecode($value);

    # more editing
    $value =~ s/[)']+$//g if defined $value;

    # if nothing return 
    return $value unless defined $value;        
  }
  
  # latitude and longitude
  if( $field =~ /latitude.+longitude/ or
      $field =~ /lat_+long*/ ){
    $value =~ s/^\s+//;
    $value =~ s/\s+$//;

    # split to lat & long
    my @tmp = split /\s*[ _,]\s*/, $value;

    # conversion if needed
    @tmp = map{ DMS_to_decimal($_) } @tmp;
    map{ s/.*?(-*[\d.]+).*/$1/ } @tmp;

    # return
    if(scalar @tmp == 2){   # if correctly split lat-long
      return ['latitude', 'longitude'], \@tmp; 
    }
    else{
      return $field, $value;
    }
  }
  elsif( $field =~ /latitude/ ){
    # conversion if needed
    $value = DMS_to_decimal($value);
    $value =~ s/.*?(-*[\d.]+).*/$1/;

    return 'latitude', $value;
  }
  elsif( $field =~ /longitude/ ){
    # conversion if needed
    $value = DMS_to_decimal($value);
    $value =~ s/.*?(-*[\d.]+).*/$1/;

    return 'longitude', $value;
  }
}


=head2 DMS_to_decimal

converting lat-long from degree-min-sec to decimal

=cut

sub DMS_to_decimal{
  my $value = shift or die "Provide lat-long, lat, or long string\n";

  # return if already in decimal format
  return $value if $value =~ /[\d.]+/;  

#  my $test = "19°25<U+0092>|85°04<U+0092>";
#  my $test2 = "50°56′02″"; 
#  $test = unidecode($test);
#  $test2 = unidecode($test2); 
#  print Dumper $test, $test2; exit;
#  $value = $test2;

  # initialize
  my $degree;
  my ($min,$sec) = (0,0);

  # degree
  if( $value =~ /(\d+)A*deg/ ){
    $degree = $1;
  }
  else{   # can't find degree, returning value
    return $value;
  }

  # min
  if($value =~ /(\d+)a2/){
    $min = $1;
  }
  elsif($value =~ /deg(\d+)/){   # assuming numeric after degree is min
    $min = $1;
  }

  # sec
  if($value =~ /(\d+)a3/){
    $sec = $1;
  }
  elsif($value =~ /a2(\d+)/){   # assumning numeric after min is sec
    $sec = $1;   
  }

  # calc decimal
  my $dec = $degree + ($min / 60) + ($sec / 3600);

  #print Dumper $dec; exit;
  return $dec;
}



=head2 get_db_samp_acc

Getting sample accessions from SRAmetadb

=cut

sub get_db_samp_acc{
  my $dbh = shift or die "Provide dbh\n";
  my $argv_r = shift;


  my $sql = $argv_r->{-sql};

  # debug
#  $sql .= " limit 100";

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


