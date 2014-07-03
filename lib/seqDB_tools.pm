package seqDB_tools;

use 5.006;
use strict;
use warnings;

=head1 NAME

seqDB_tools - tools for getting and analyzing data from nucleotide databases

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
use LWP::UserAgent;
use URI::Escape;
use IPC::Cmd qw/run can_run/;

use seqDB_tools::download qw/MGRAST_api_download/;


=head2 MGRAST_mapper_main

main subroutine for downloading from mgrast

=head3 IN

$argv_r -- Getopt::Euclid %ARGV

=cut

push @EXPORT_OK, 'MGRAST_mapper_main';

sub MGRAST_mapper_main{
  my $argv_r = shift or confess "Provide getopt::euclid argv\n";

  # getting ids
  my $ids_r = parse_ids($argv_r->{-id});

  # check for bowtie2 
  can_run('bowtie2') or confess "ERROR: 'bowtie2' is not in your \$PATH\n";

  # each metagenome ID
  my %report;
  foreach my $id (@$ids_r){ 
    # downloading file(s)
    my $tmp_files_r = MGRAST_api_download($id, $argv_r, -temp => 1);

    # mapping
    call_bowtie2($tmp_files_r, \%report, $argv_r);
  }

  # writing out report
  write_report( \%report, $argv_r );

}

=head2 parse_ids

Parsing file containing ids. 1 id per line

=cut

sub parse_ids{
  my $id_opt = shift or confess "Provide -id arg\n";
  
  my $infh;
  $id_opt eq '-' ? $infh = \*STDIN : 
    open $infh, $id_opt or confess $!;

  my @ids;
  while(<$infh>){
    next if /^\s*$/;
    s/^\s*(\S+)\s*$/$1/;

    ## id format
    die "Don't recognize id format: " . $_
      unless $_ =~ /\d\d\d\d\d\d\d.\d/;

    push @ids, $_;
  }
  close $infh or die $!;

  return \@ids;
}

=head2 write_report

Writing report on mappings of each metagenome.

=cut

sub write_report{
  my $report_r = shift or confess "Provide hashref of report\n";
  my $argv_r = shift;

  # output file handle
  my $outfh;
  if( defined $argv_r and 
      exists $argv_r->{-output} and 
      $argv_r->{-output} ne '-'){
    open $outfh, ">", $argv_r->{-output} or confess $!;
  }
  else{
    $outfh = \*STDOUT;
  }

  # getting all possible columns
  my %header;
  foreach my $filename (keys %$report_r){
    foreach my $cat1 (keys %{$report_r->{$filename}}){
      foreach my $cat2 (keys %{$report_r->{$filename}{$cat1}}){
	$header{$cat1}{$cat2} = 1 unless exists $header{$cat2};
      }
    }
  }
  my @header;
  foreach my $cat1 (sort {$b cmp $a} keys %header){
    push @header, (sort keys %{$header{$cat1}});
  }
  print join("\t", 'bam_file', @header), "\n";

  # writing
  foreach my $filename (keys %$report_r){
    my @row;
    foreach my $cat1 (keys %{$report_r->{$filename}}){      
      foreach my $col ( sort keys $header{$cat1} ){
	if(exists $report_r->{$filename}{$cat1}{$col}){
	  push @row, $report_r->{$filename}{$cat1}{$col};
	}
	else{
	  #print Dumper "col: $col"; 
	  #print Dumper $report_r->{$filename}{$cat1}; exit;
	  push @row, 'NA';
	}
      }
    }
    # writing row
    print $outfh join("\t", $filename, @row), "\n";
  }
}


=head2 call_bowtie2

Calling bowtie2 on a set of files

=cut

sub call_bowtie2{
  my $files_r = shift or confess "Provide hashref of files\n";
  my $report_r = shift or confess "Provide report hashref\n";
  my $argv_r = shift;

  # input check
  my $bt2_params = exists $argv_r->{-bt2_params} ? 
    $argv_r->{-bt2_params} : '';
  my $idx = exists $argv_r->{-x} ? $argv_r->{-x} : '';


  #print Dumper $files_r; exit;  

  foreach my $filename (keys %$files_r){
    my $entry = $files_r->{$filename};

    # output file
    (my $basename = $entry->{file_name}) =~ s/\.[^.]+$//;
    my $outfile = join("_",
		   $entry->{id},
		   $entry->{file_id},
		   $basename . ".bam"		    
		  );
    ## adding prefix if needed
    $outfile = join("_", $argv_r->{-prefix}, $filename) if
      defined $argv_r->{-prefix} and $argv_r->{-prefix} ne ''; 

    # adding '-f' param if fasta (fna) file
    exists $entry->{file_type} or confess "Cannot find file_type for entry\n";
    my $fna_flag = $entry->{file_type} =~ /^(fna|fasta)$/i ? '-f' : '';

    # cmd
    my $cmd = join(" ", 
		   "bowtie2", 
		   "-U", $filename,
		   "-x", $idx, 
		   $bt2_params,
		   $fna_flag,
		   "|", 
		   "samtools view -b -S -F 0x0004 -",
		   ">", $outfile
		   );
		
    ## call
    print STDERR $cmd,"\n" unless $argv_r->{'--quiet'};
    my( $success, $error_message, $full_buf, 
	$stdout_buf, $stderr_buf ) =
	  run( command => $cmd, 
	       verbose => $argv_r->{'--debug'} );    

    ## parsing stderr
    $report_r->{$outfile}{bowtie2} = parse_bowtie2_stderr($stderr_buf);

    ## adding entry info
    map{ $report_r->{$outfile}{metagenome}{$_} = $entry->{$_} 
	   unless $_ eq 'tmpdir' } keys %$entry;
  }

  #print Dumper %$report_r; exit;
}

=head2 parse_bowtie2_stderr

Parsing bowtie2 stderr to get stats on mapping

=cut

sub parse_bowtie2_stderr{
  my $stderr_buf = shift or confess "Provide stderr buffer\n";
 
  my %res;
  foreach my $l (@$stderr_buf){
    my @ll = split /(\n|\r)/, $l;

    foreach my $ll (@ll){
      chomp $ll;
      next if $ll =~ /^\s*$/;

      if( $ll =~ /^\s*(\d+) +\(([\d.]+)%\) +([^;]+)/ ){
	next unless defined $3;
	my ($count, $percent, $cat) = ($1, $2, $3);
	$cat =~ s/ /_/g;

	$res{"$cat\__count"} = $count;
	$res{"$cat\__percent"} = $percent;
      }
    }
  }

  return \%res;
}



=head1 AUTHOR

Nick Youngblut, C<< <ndy2 at cornell.edu> >>

=head1 BUGS



=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc seqDB_tools


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=seqDB_tools>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/seqDB_tools>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/seqDB_tools>

=item * Search CPAN

L<http://search.cpan.org/dist/seqDB_tools/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Nick Youngblut.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of seqDB_tools
