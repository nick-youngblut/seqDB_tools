#!/usr/bin/env perl
use strict;
use warnings;

=pod

=head1 NAME

SRAmetadb_edit_lat_long.pl -- editing latitude and longitude values in SRAmetadb (or other sqlite3 DB)

=head1 VERSION

This is version 0.0.1

=head1 USAGE

    SRAmetadb_edit_lat_long.pl [options]

=head1 REQUIRED ARGUMENTS

=item -db <db> | -database <db>

SRAmetadb file (from SRAdb R package)

=head1 OPTIONS

=over

=item -sql <sql>

sql statement to get lat-long values.

Default: sql.default

=for Euclid:
sql.type: string
sql.default: 'SELECT sample_accession, latitude_and_longitude, latitude, longitude from efetch'

=item -table <table_name>

The name of the database table to update values.

Default: table_name.default

=for Euclid:
table_name.type: string
table_name.default: 'efetch'


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

Add/update a latitude and longitude columns
in the SRAmetaDb from R package SRAdb.
The script should also work for any other
database with a table containing 3 columns
for latitude, longitude, and latitude_and_longitude.

The '-sql' and '-table' flags define which table(s)
is is used to select the needed columns 
(latitude, longitude, and latitude_and_longitude)
and which table to update the values of those columns.

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
use Text::ParseWords;
use Regexp::Common;
use Text::Unidecode;
use Encode;
use List::MoreUtils qw/none lastidx/;

#--- I/O error ---#


#--- MAIN ---#
# db connect
my $dbh = connect2db( $ARGV{-db} );

# getting all sample accessions from db
my $entries_r = get_db_entries($dbh, \%ARGV);


# edit latitude and longitude
## trying to get all entries in to lat and long fields
$ARGV{entry_edit_cnt} = 0;
foreach my $samp_acc (keys %$entries_r){
  edit_lat_long(
		$entries_r->{$samp_acc},
		\%ARGV
	       );
}

# status
printf STDERR "Number of lat-long values edited: %i\n\n",
  $ARGV{entry_edit_cnt};


# updating entries
update_db_entries($dbh, $entries_r, \%ARGV );


# closing 
$dbh->commit or die $dbh->err;
$dbh->disconnect or die $dbh->err;



#--- Subroutines ---#

sub debug_entries{
  my $entries = shift or die "Provide hashref of entries\n";
  foreach my $k (keys %$entries_r){
    my @row;
    if(defined $entries_r->{$k}{latitude} and defined $entries_r->{$k}{longitude}){
      map{ $_ = '' unless defined $_ } values %{$entries_r->{$k}};
      map{ push @row, $entries_r->{$k}{$_} } sort keys %{$entries_r->{$k}};
      print join("|", @row), "\n";
    }
  }

}

=head2 update_db_entries

Updating db entries. 

=cut

sub update_db_entries{
  my $dbh = shift or die "Provide dbh\n";
  my $entries_r = shift or die "Provide hashref of entries\n";
  my $argv_r = shift;

  my $sql = join(" ", 
		 "UPDATE ",
		 $argv_r->{-table},
		 "SET latitude = ?, longitude = ?",
		 "WHERE sample_accession = ?");

  my $sth = $dbh->prepare($sql) or die $dbh->err;


  # loading entriesq
  foreach my $samp_acc ( keys %$entries_r ){
    # bind params
    $sth->bind_param(1, $entries_r->{$samp_acc}{latitude})
      or die $dbh->err;
    $sth->bind_param(2, $entries_r->{$samp_acc}{longitude})
      or die $dbh->err;
    $sth->bind_param(3, $entries_r->{$samp_acc}{sample_accession})
      or die $dbh->err;

    # execute
    $sth->execute or die $dbh->err;
  }

  print STDERR "Entry updating complete.\n";
}



=head2 edit_lat_long

Top-level subroutine for editing lat-long.
Trying to parse human semi-readable lat-long.
Lat-long can include unicode.

=cut

sub edit_lat_long{
  my $entry_r = shift or die "Provide hashref of entry\n";
  my $argv_r = shift;

  # debug
  if($argv_r->{'--debug'}){
    $entry_r = { 
		sample_accession => 'TEST',
#		latitude_and_longitude => '33: 56: 34.0 S; 60: 33: 59.9 W',
#		latitude_and_longitude => "29Á 2\' 16.332 95Á 16\' 0.948",
		latitude_and_longitude => "82 N 62 W",
		latitude => undef,
		longitude => undef
	       };
  }

  
  # defined 
  ## if lat,long defined: trying to get from lat_and_long
  if(defined $entry_r->{latitude_and_longitude} and
     ! defined $entry_r->{latitude} and
     ! defined $entry_r->{longitude} ){
     
    # value edit
    ## colons
    $entry_r->{latitude_and_longitude} =~ s/([NS]):/$1;/;

    ## (N:W) style
    if( $entry_r->{latitude_and_longitude} =~ 
	/\(([NS]):([EW])\)\s*(-*\d+\.*\d*)[ _:;\/](-*\d+\.*\d*)/){
      $entry_r->{latitude} = join(" ", $3, $1);
      $entry_r->{longitude} = join(" ", $4, $2);
    }
    # simple numeric style
    # elsif($entry_r->{latitude_and_longitude} =~ 
    # 	  /^\D*?(-*\d+\.*\d*)\D*?[ \t:;,\/]\D*?(-*\d+\.*\d*)\D*?([NSns])*$/){

    #   my @regex_vals = ($1, $2);

    #   # getting direction
    #   my ($lat_dir, $long_dir) = ('','');
    #   if( $entry_r->{latitude_and_longitude} =~ /\d\D*([NSns])/ ){
    # 	$lat_dir = $1;
    #   }
    #   if( $entry_r->{latitude_and_longitude} =~ /\d\D*([EWew])/ ){
    # 	$long_dir = $1;
    #   }



    #   map{ $_ = '' unless defined $_ } @regex_vals;
    #   $entry_r->{latitude} = join(" ", @regex_vals[0..1]);
    #   $entry_r->{longitude} = join(" ", @regex_vals[2..3]);
    #   }

    ## lat-long DMS style
    elsif( $entry_r->{latitude_and_longitude} =~
	   /^\D*?(-*[\d.]+)\D*[ ,:;_\/]\D*?(-*[\d.]+)\D*$/ 
	   or $entry_r->{latitude_and_longitude} =~ 
	   /^\D*?(-*[\d.]+)\D*?(-*[\d.]+)\D*[ ,:;_\/]\D*?(-*[\d.]+)\D*?(-*[\d.]+)\D*$/ 
	   or $entry_r->{latitude_and_longitude} =~ 
	   /^\D*?(-*[\d.]+)\D*?(-*[\d.]+)\D*?(-*[\d.]+)\D*[ ,:;_\/]\D*?(-*[\d.]+)\D*?(-*[\d.]+)\D*?(-*[\d.]+)\D*$/ 
	 ){
	     
      # saving regex
      my @regex_vals = ($1, $2, $3, $4, $5, $6);
      
      # splitting regex so each side has even number of defined values
      my $idx = lastidx { defined $_ } @regex_vals;
      $idx += 1 if ($idx+1) % 2 != 0;  # round up if odd
      my $mididx = int($idx / 2); 
      @regex_vals = ( [@regex_vals[0..$mididx]],
		      [@regex_vals[$mididx+1..$idx]] );
      ## adding values not present
      foreach my $x (@regex_vals){
	for my $i (0..2){
	  $x->[$i] = '' unless defined $x->[$i];
	}
      }

      # getting direction
      my ($lat_dir, $long_dir) = ('','');
      if( $entry_r->{latitude_and_longitude} =~ /\d\D*([NSns])/ ){
	$lat_dir = $1;
      }
      if( $entry_r->{latitude_and_longitude} =~ /\d\D*([EWew])/ ){
	$long_dir = $1;
      }

      #map{ $_ = '' unless defined $_ } @regex_vals;
      $entry_r->{latitude} = join(" ", @{$regex_vals[0]}, $lat_dir);
      $entry_r->{longitude} = join(" ", @{$regex_vals[1]}, $long_dir);

    }
    ## undef if entry is not numeric
    elsif( $entry_r->{latitude_and_longitude} =~ /^\D*$/ ){
      $entry_r->{latitude} = undef;
      $entry_r->{longitude} = undef;
    }
    ## undef if entry just has 1 numeric 
    elsif( $entry_r->{latitude_and_longitude} =~ /^\D*$RE{num}{real}\D*$/) {
      $entry_r->{latitude} = undef;
      $entry_r->{longitude} = undef;
    }
    ## split failed; don't know why
    else{      
      printf STDERR "WARNING: don't know how to parse lat-long value '%s' in sample '%s'. ",
	 $entry_r->{latitude_and_longitude}, $entry_r->{sample_accession};
      print STDERR "Making NULL\n";
      $entry_r->{latitude} = undef;
      $entry_r->{longitude} = undef;
    }
  }

  # debug split
  if($argv_r->{'--debug'}){
    printf STDERR "split-lat: %s\n", $entry_r->{latitude} if defined $entry_r->{latitude};
    printf STDERR "split-long: %s\n", $entry_r->{longitude} if defined $entry_r->{longitude};
  }

  ## if lat or long defined
  my ($latitude, $longitude);
  if(defined $entry_r->{latitude} or
     defined $entry_r->{longitude} ){

    # getting lat-long from just lat or long
    ($latitude,$longitude) = find_lat_or_long( $entry_r->{latitude}, 
					       'longitude')
      unless defined $entry_r->{longitude};

    ($latitude,$longitude) = find_lat_or_long( $entry_r->{longitude}, 
					       'latitude')
      unless defined $entry_r->{latitude};

    # interpreting lat or long
    $latitude = interpret_lat_long( 
				   $entry_r->{latitude},
				   dir => 'latitude'
				  );
    $longitude = interpret_lat_long(  
				    $entry_r->{longitude}, 
				    dir => 'longitude'
				   );    
  }
  ## none defined
  elsif( none{ defined $entry_r->{$_} } 
	 qw/latitude_and_longitude latitude longitude/){
    
  }
  ## else ?
  else{
    printf STDERR "WARNING: could not determine decimal lat-long for %s\n",
      $entry_r->{sample_accession};
    map{ $entry_r->{$_} = undef } qw/latitude_and_longitude latitude longitude/;
  }


  # reassigning lat and long
  ## counting edits
  if( ! defined $latitude or ! defined $longitude ){
    
  }
  elsif( (! defined $entry_r->{latitude} and defined $latitude) or
      (! defined $entry_r->{longitude} and defined $longitude) ){
    $argv_r->{entry_edit_cnt}++;
  }
  elsif( $entry_r->{latitude} ne $latitude or
	 $entry_r->{longitude} ne $longitude ){
    $argv_r->{entry_edit_cnt}++;
  }
  ## making change
  $entry_r->{latitude} = $latitude;
  $entry_r->{longitude} = $longitude;


  # final check of lat-long
  foreach my $cat (qw/latitude longitude/){
    next unless defined $entry_r->{$cat};   # skipping undef values
 
    # is decimal?
    unless($entry_r->{$cat} =~ /^-*\d+\.*\d*$/){
      printf STDERR "ERROR: $cat value is not in decimal form! Value = '%s'\n",
	$entry_r->{$cat};
      exit(1);
    }
    # is realistic value?
    if( $entry_r->{$cat} > 500 ){
      printf STDERR "WARNING: $cat value does not appear to be correct (Value = '%s'). Making NULL\n",
	$entry_r->{$cat};
      $entry_r->{$cat} = undef;
    }
  }

  # debug
  if( $argv_r->{'--debug'} ){
    print Dumper "final entry: ", $entry_r; exit;
  }
}


=head2 interpret_lat_long

Interpreting lat and long.
Many contain unicode and many need to
be converted from DMS format to decimal.

=cut 

sub interpret_lat_long{
  my $value = shift or die "Provide latitude or longitude value\n";
  my $kwargs = @_;

  my $decimal;
  # just real number? assumed decimal
  if($value =~ /^\D*?(-*\d+\.*\d*)\D*$/){
    $decimal = $1;
    $decimal *= -1 if $value =~ /\d+\D*[SWsw]/;
  }
  ## decimal with direction?
#  elsif($value =~ /^\D*?(-*\d+\.*\d*)\D*?([NSEWnsew])\D*$/){
#    my ($decimal, $dir) = ($1, $2);
#    $decimal *= -1 if $dir =~ /[SWsw]/;
#  }
  # is DMS format?
  elsif($value =~ /\D*?(-*\d+\.*\d*)\D+?(-*\d+\.*\d*)(\D+?-*\d+\.*\d*)*/){
    my ($degrees, $minutes, $seconds) = ($1, $2, $3);
    
    # getting direction if available
    my $dir;
    if($value =~ /\d\D*([NSEWnsew])/ ){
      $dir = $1;
    }
    $decimal = DMS_to_decimal($degrees, $minutes, $seconds, $dir);
  }
  # don't know what to do
  else{    
    warn "Don't know how to parse lat or long value: '$value'\n";
  }

#  print Dumper "ret: $decimal"; exit;
  return $decimal;
}

=head2 DMS_to_decimal

converting lat-long from degree-min-sec to decimal

=cut

sub DMS_to_decimal{
  my $degrees = shift;
  defined $degrees or die "Provide degrees\n";
  my $minutes = shift;
  $minutes = 0 unless defined $minutes;
  my $seconds = shift;
  $seconds = 0 unless defined $seconds;
  my $direction = shift;
  $direction = 'N' unless defined $direction;


  # checking that degrees, minutes, seconds are numeric
  foreach my $x ($degrees, $minutes, $seconds){
    die "ERROR: value '$_' does not contain a real number\n"
      unless $x =~ /-*[\d.]+/;
    $x =~ s/\D*?(-*\d+\.*\d*)\D*/$1/;
  }


  # calc decimal
  my $decimal = $degrees + ($minutes / 60) + ($seconds / 3600);

  $decimal = -$decimal if $direction =~ /^\s*[SWsw]/;

#  print Dumper $decimal; exit;
  return $decimal;
}


=head2 find_lat_or_long

Either lat or long missing.
Checking for missing value defined value.

=cut

sub find_lat_or_long{
  my $present = shift or die "Provide defined lat or long\n";
  my $missing = shift or die "Provide missing: latitude or longitude\n";
  
  my ($latitude, $longitude);
  if( $missing eq 'latitude' ){
    if( $present =~ /^\s*\([NS]:[EW]\) +(-*[\d.]+):+(-*[\d.]+)/ ){
      $latitude = $1;
      $longitude = $2;
    }
    else{
      die "Dont know how to parse: '$present'\n";
    }
  }
  elsif( $missing eq 'longitude' ){
      die "Dont know how to parse: '$present'\n";
  }
  else{
    die "ERROR: cannot determine missing variable: '$missing'\n";
  }

  return $latitude, $longitude;
}


=head2 get_db_entries

Getting entries from SRAmetadb

=cut

sub get_db_entries{
  my $dbh = shift or die "Provide dbh\n";
  my $argv_r = shift;


  # getting column names
  my $sth = $dbh->prepare( $argv_r->{-sql} ) or die $dbh->err;

  # checking column names
  my $col_names = $sth->{NAME_lc_hash};
  die "ERROR: 'sample_accession' must be returned by sql\n"
    unless grep /sample_accession/, keys %$col_names; 
  die "ERROR: 'latitude_and_longitude', 'latitude', or 'longitude' must be returned by sql\n"
    unless grep /(latitude_and_longitude|latitude|longitude)/, keys %$col_names; 


  # query
  $sth->execute() or die $dbh->err;
  
  # get columns
  my $ret = $sth->fetchall_hashref('sample_accession') or die $dbh->err;


  # status
  printf STDERR "Number of entries returned from query: %i\n\n",
    scalar keys %$ret;


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


