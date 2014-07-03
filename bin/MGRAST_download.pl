#!/usr/bin/env perl
use strict;
use warnings;

=pod

=head1 NAME

MGRAST_download.pl -- General file downloading with MG-RAST API

=head1 VERSION

This is version 0.0.1

=head1 USAGE

    MGRAST_download.pl [options]

=head1 REQUIRED ARGUMENTS

=over

=item -id <ids>...

>=1 metagenome ID

=for Euclid:
ids.type: str

=back

=head1 OPTIONS

=over

=item -stage <stage>...

Processing stage for downloading files.
If multiple stages provided, will try each
in succession until entries are returned.

Default: stage.default

=for Euclid:
stage.type: num
stage.default: [150,100]


=item -prefix <prefix>

Output file prefix.

Default: prefix.default

=for Euclid:
prefix.type: string
prefix.default: ''


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

This is a simple script that uses
the MG-RAST API to download 
metagenome nucleotide files.

=head2 -stage flag

Multiple stages can be provided as
a backup in case the first stage doesn't
return any entries (some stages are optional).

Stages can be found in Appendix A of the
MG-RAST manual
(ftp://ftp.metagenomics.anl.gov/data/manual/mg-rast-tech-report-v3_r1.pdf).

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


#--- Input check ---#
map{ die "Don't recognize id format: " . $_ 
       unless $_ =~ /\d\d\d\d\d\d\d.\d/ } @{$ARGV{-id}};

#--- MAIN ---#
# get a user agent
my $ua = LWP::UserAgent->new;
$ua->agent("MyClient/0.1 ");

# each metagenome ID
foreach my $id (@{$ARGV{-id}}){
  # getting file entries from first stage, checking next stage if nothing
  my $file_entries;
  foreach my $stage (@{$ARGV{-stage}}){
    $file_entries = get_files_in_stage($id, $stage);
    if( defined $file_entries ){
      last;  # keeping these file entries
    }
    else{
      print STDERR "WARNING: no entries found for id->$id, stage->$stage.\n";      
    }
  }
    
  # downloading fastx files
  get_fastx_files($file_entries, \%ARGV);
}


#--- Subroutines ---#
sub get_fastx_files{
  my $file_entries = shift or die "Provide hashref of entries\n";
  my $argv_r = shift or die "Provide argv\n";

  foreach my $entry (@$file_entries){
    # passed fastx files
    if(exists $entry->{file_type} and 
       $entry->{file_type} =~ /^f(na|asta|astq|q)$/ and
       exists $entry->{file_name} and
       $entry->{file_name} !~ /removed\.f(na|asta|astq|q)(\.gz)*/ and
       exists $entry->{id} and 
       exists $entry->{file_id} ){
      
      # unpacking
      my $url = exists $entry->{url} ? 
	$entry->{url} : die "ERROR: no url for ".$entry->{id}.".\n";
      
      # making filename
      my $filename = join("_",
			  $entry->{id}, 
			  $entry->{file_id}, 
			  $entry->{file_name}
			 );
      $filename = join("_", $ARGV{-prefix}, $filename) if 
	defined $ARGV{-prefix} and $ARGV{-prefix} ne '';

      # response
      print STDERR "Download url: '$url'\n" unless $argv_r->{'--quiet'};
      my $response = $ua->get($url, ":content_file" => $filename);

      ## unless success
      if($response->is_success){
	print STDERR "Successful download to: '$filename'\n";
      }
      else{
	warn "WARNING: unsuccessful url: '$url'\n";
      }
    }
  }
}


sub get_files_in_stage{
  my $id = shift or die "Provide id\n";
  my $stage = shift or die "Provide stage\n";
  my $argv_r = shift;

  my $url = 'http://api.metagenomics.anl.gov/1/download/';
  $url .= $id;
  $url .= "?stage=$stage";

  # response
  print STDERR "Stage url: '$url'\n" unless $argv_r->{'--quiet'};
  my $response = $ua->get($url);
  ## unless success
  my $file_entries;
  if( ! $response->is_success ){
    warn "WARNING: unsuccessful url: '$url'\n";    
  }
  else{
    my $jsonO = new JSON;
    $file_entries = $jsonO->decode( $response->content );
  }

#  print Dumper $file_entries; exit;
  return $file_entries;
}


