#!/usr/bin/env perl
use strict;
use warnings;

=pod

=head1 NAME

ENA_sample_xmlParse.pl -- European Nucleotide Archive ERS/DRS/SRS accession attribute parser

=head1 VERSION

This is version 0.0.1

=head1 USAGE

    ENA_sample_xmlParse.pl [options]

=head1 REQUIRED ARGUMENTS

=over

=item <ERS-DRS-SRS_accession_list>

A list of ERS/DRS/SRS accessions (1 per line; '-' if from STDIN)

=back

=head1 OPTIONS

=over

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

ENA attributes of samples can only be obtained
from the ENA Browser by parsing xml of the 
ERS/DRS/SRS accession (ie., 'secondary_sample_accession')
for the sample.

Output written to STDOU

=head1 EXAMPLES

=head2 1 secondary sample accession

    echo 'ERS326120' | perl ENA_sample_xmlParse.pl - | less

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
use XML::Parser;

#--- I/O error ---#
my $fh;
$ARGV{'<ERS-DRS-SRS_accession_list>'} eq "-" ? $fh = \*STDIN :
  open $fh, $ARGV{'<ERS-DRS-SRS_accession_list>'} or die $!;

#--- MAIN ---#
# parsing accession list
my $acc_r = parse_accession_list($fh);

# foreach accession: parsing xml
my %res;
foreach my $acc (@$acc_r){
  my $parser = new XML::Parser( Style => 'Tree');
  my $url = "http://www.ebi.ac.uk/ena/data/view/$acc&display=xml";
  print STDERR "curl '$url'\n";
  my $curl = `curl '$url'`;
  my $parsed = $parser->parse($curl);

  # value table
  my $tbl_r =  ${$parsed->[1]}[4];
  
  # center = [0]
  # SUBMITTER_ID = [4][14?]
  # TITLE = [8][2]
  # SAMPLE_NAME = [12]
  # DESCRIPTION = [16]
  # SAMPLE_ATTRIBUTES = [24] (name = 23)
  
  $res{$acc} = parse_attributes($tbl_r);
}

# getting all tags used
my %all_tags;
foreach my $acc (keys %res){
  map{ $all_tags{$_} = 1 } keys %{$res{$acc}};
}

# writing header
my @tags = sort keys %all_tags;
print join("\t", 'secondary_sample_accession', @tags), "\n";

# writing output
foreach my $acc (keys %res){
  # '' if undefined values
  map{ $_ = '' unless defined $_ } @{$res{$acc}}{ keys %all_tags }; 
  # output
  print join("\t", $acc, @{$res{$acc}}{ @tags }), "\n";
}



#--- Subroutines ---#
sub parse_attributes{
  my $tbl_r = shift;

  my %attribs;
  # sample attribute parsing
  for my $i (0..$#$tbl_r){
    #print Dumper $tbl_r->[$i];
    if( $tbl_r->[$i] eq "SAMPLE_ATTRIBUTES" ){
      my $sample_att_r = $tbl_r->[$i+1];
      #print Dumper @$sample_att_r; exit;
      for my $ii (0..$#$sample_att_r){
	if( $sample_att_r->[$ii] eq "SAMPLE_ATTRIBUTE" ){
	  my ($tag, $value);
	  
	  my $list_r = $sample_att_r->[$ii+1];
	  foreach my $iii (0..$#$list_r){
	    if( $list_r->[$iii] eq "TAG" ){
	      $tag = $list_r->[$iii+1][2];
          }
	    elsif( $list_r->[$iii] eq "VALUE" ){
	      $value = $list_r->[$iii+1][2];
	    }
	  }

	  $attribs{$tag} = $value;
	}
      }
    }
  }
  
  return \%attribs;
}

sub parse_accession_list{
  # parsing the accession number list provided via a filehandle
  my ($fh) = shift;

  my @acc;
  while(<$fh>){
    chomp;
    next if /^\s*$/;
    
    push @acc, $_;
  }

  return \@acc;
}
