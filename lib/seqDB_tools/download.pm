package seqDB_tools::download;

use 5.006;
use strict;
use warnings;

=head1 NAME

download -- subroutines for database info downloading

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

=head1 EXPORT


=head1 SUBROUTINES/METHODS

=cut

use base 'Exporter';
our @EXPORT_OK = '';
use Carp qw/ croak confess carp /;
use Data::Dumper;
use File::Temp;
use JSON;


=head2 MGRAST_api_download

Main for downloading via mgrast api

=head3 IN

id -- metagenome id <string>
argv_r -- argv

=cut

push @EXPORT_OK, 'MGRAST_api_download';

sub MGRAST_api_download{
  my $id = shift or confess "Provide metagenome id\n";
  my $argv_r = shift or confess "Provide getopt::euclide argv\n";
  my %opts = @_; 

  # get a user agent
  my $ua = LWP::UserAgent->new;
  $ua->agent("MyClient/0.1 ");

  # getting file entries from first stage, checking next stage if nothing
  my $file_entries;
  exists $argv_r->{-stage} or confess "Provide arg: -stage\n";
  foreach my $stage (@{$argv_r->{-stage}}){
    $file_entries = get_MGRAST_stage_files($id, 
					   $stage,
					   $ua);
    if( defined $file_entries ){
      last;  # keeping these file entries
    }
    else{
      print STDERR "WARNING: no entries found for id->$id, stage->$stage.\n";
    }
  }


  # downloading fastx files
  my $files_r = download_nuc_files($file_entries, 
				   $ua,
				   $argv_r, 
				   %opts);

  # returning hashref of downloaded file info
  return $files_r;
}


=head2 download_nuc_files

Getting fast[aq] files from mg-rast.
Writing to a file.

=cut

sub download_nuc_files{
  my $file_entries = shift or confess "Provide hashref of entries\n";
  my $ua = shift or confess "Provide a user agent\n";
  my $argv_r = shift or confess "Provide argv\n";
  my %opts = @_;


  my %file_info;
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
      $filename = join("_", $argv_r->{-prefix}, $filename) if
	defined $argv_r->{-prefix} and $argv_r->{-prefix} ne '';

      ## persistent or temporary file?
      ### temp
      my $tmpdiro;
      if( exists $opts{-temp} ){ 	
	$tmpdiro = File::Temp->newdir();
	my $tmpdirname = $tmpdiro->dirname;
	$filename = File::Spec->catfile($tmpdirname, $filename);
      }
      ### persist
      elsif( exists $argv_r->{-prefix} and $argv_r->{-prefix} ne '' ){
	$filename = join("_", $argv_r->{-prefix}, $filename);
      }

      # response
      print STDERR "Download url: '$url'\n" unless $argv_r->{'--quiet'};
      my $response = $ua->get($url, ":content_file" => $filename);

      ## unless success
      if($response->is_success){
	print STDERR "Successful download to: '$filename'\n";
	
	# saving file info
	$file_info{$filename}{tmpdir} = $tmpdiro;	
	map{ $file_info{ $filename }{ $_ } = $entry->{$_} }
	  keys %$entry;
      }
      else{
        warn "WARNING: unsuccessful url: '$url'\n";
      }
    }
  }
  
  # return
  return \%file_info;
}


=head2 get_MGRAST_stage_files

getting all files in MGRAST stage

=cut

sub get_MGRAST_stage_files{
  my $id = shift or confess "Provide id\n";
  my $stage = shift or confess "Provide stage\n";
  my $ua = shift or confess "Provide a user agent\n";
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


=head1 AUTHOR

Nick Youngblut, C<< <ndy2 at cornell.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-download at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=download>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc download


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=download>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/download>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/download>

=item * Search CPAN

L<http://search.cpan.org/dist/download/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Nick Youngblut.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of download
