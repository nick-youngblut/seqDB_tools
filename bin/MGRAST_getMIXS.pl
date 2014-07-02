#!/usr/bin/env perl
use strict;
use warnings;

=pod

=head1 NAME

MGRAST_getMIXS.pl -- get MIXS metadata in csv format from MG-RAST API RESTFUL submission.

=head1 VERSION

This is version 0.0.1

=head1 USAGE

    MGRAST_getMIXS.pl [options]

=head1 REQUIRED ARGUMENTS

=head1 OPTIONS

=over

=item -metadata <metadata>...

>=1 query string for any metadata field

=item -limit <limit>

Limit the number of returned entries.

Default: limit.default

=for Euclid:
limit.type: int >= 1
limit.default: 10000

=item -direction <direction>

Sort direction (asc = ascending order; desc = descending order).

Default: direction.default

=for Euclid:
direction.type: string, direction eq 'asc' || direction eq 'desc'
direction.default: 'asc'


=item -status <status>

Permissions status of metagenome.
(both = all data; public = just public data; private = just private).

Default: status.default

=for Euclid:
status.type: string, status eq 'both' || status eq 'public' || status eq 'private'
status.default: 'public'


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

This is a simple script to call
the MG-RAST API and download metagenome
MIXS metadata. The output is in tab-delimited format.

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
use Data::Dumper;
use Getopt::Euclid;
use JSON;
use LWP::UserAgent;
use URI::Escape;
use Data::Dumper;




#--- define params ---#


#--- MAIN ---#
# get a user agent
my $ua = LWP::UserAgent->new;
$ua->agent("MyClient/0.1 ");


# constructing url
my $url = "http://api.metagenomics.anl.gov/1/metagenome?status=public&verbosity=mixs";
## metadata
if( exists $ARGV{-metadata} ){
  $url = join("&metadata=", $url, @{$ARGV{-metadata}} );
}
## limit
$url .= "&limit=" . $ARGV{-limit};
## direction
$url .= "&direction=" . $ARGV{-direction};
## status
$url .= "&status=" . $ARGV{-status};



# response
print STDERR "url: '$url'\n" unless $ARGV{'--quiet'};
my $response = $ua->get($url);
## unless success
my $mixs_r;
if( ! $response->is_success ){
  die "ERROR: unsuccessful url: '$url'\n";
}
else{
  my $jsonO = new JSON;
  $mixs_r = $jsonO->decode( $response->content );
}


# status
if(exists $mixs_r->{total_count}){
  printf STDERR "Number of entries returned by query: %i\n",
    $mixs_r->{total_count};
}

# mgrast_mixs to csv
mgrast_mixs_to_csv($mixs_r);



#--- Subroutines ---#

sub mgrast_mixs_to_csv{
  my $mixs_r = shift or die "Provide mixs ref\n";
  my $argv_r = shift;

  # input check
  die "ERROR: no 'data' key in mixs ref\n"
    unless exists $mixs_r->{data};

  # collect all unique fields
  my %fields;
  foreach my $entry_r ( @{$mixs_r->{data}} ){
    map{ $fields{$_} = 1 unless exists $fields{$_} } keys %$entry_r;
  }
  my @fields = sort keys %fields;

  # writing header
  print join("\t", @fields),"\n";

  # writing table
  foreach my $entry_r ( @{$mixs_r->{data}} ){
    my @row;
    foreach my $field (@fields){
      if(exists $entry_r->{$field}){
	$entry_r->{$field} = 'NA' if $entry_r->{$field} =~ /^\s*$/;
	$entry_r->{$field} =~ s/\t/ /g;  # making sure no tabs in field value
	push @row, $entry_r->{$field};
      }
      else{
	push @row, 'NA';
      }
    }
    print join("\t", @row),"\n";
  }
}


