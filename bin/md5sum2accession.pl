#!/usr/bin/env perl

=pod

=head1 NAME

md5sum2accession.pl -- get accession numbers for md5sum values

=head1 SYNOPSIS

md5sum2accession.pl [flags] < md5sum.txt > md5sum-acc.txt

=head2 Required flags

=over

NONE

=back

=head2 Optional flags

=over

=item -column  <int>

Column containing the md5sum values (indexed by 1). [2]

=item -help  <bool>

Print this help message & exit. [FALSE]

=back

=head2 For more information:

perldoc md5sum2accession.pl

=head1 DESCRIPTION

Get accession numbers for md5sum values.

=head1 EXAMPLES

=head2 Basic usage:

md5sum2accession.pl < md5sum.txt > md5sum-acc.txt

=head1 AUTHOR

Nick Youngblut <ndy2@cornell.edu>

=head1 AVAILABILITY

email me

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


#--- modules ---#
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use JSON;
use LWP::UserAgent;
use URI::Escape;


#--- args/flags ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose_b);
my $column = 2;
my $batch = 100000;
GetOptions(
	"column=i" => \$column,
	"verbose" => \$verbose_b,
	"batch=i" => \$batch,
	"help|?" => \&pod2usage # Help
	);

#--- I/O error ---#
$column--;

#--- setting defaults ---#

#--- MAIN ---#
# get a user agent
my $ua = LWP::UserAgent->new;
$ua->agent("$0 ");

# foreach line; get md5sum & then get accession 
## header
my $l = <>;
print $l;

while(1){		# loading 100000 lines into memory at once
	last if eof;

	my %lines;
	for(0..($batch-1)){		
		my $l = <>;
		last unless defined $l;
		
		chomp $l;
		my @l = split /\t/, $l;
		die "ERROR: cannot find column: ", $column +1, " in line $.\n"
			unless defined $l[$column];

		push @{$lines{$l[$column]}}, \@l;
	}
	
	last unless %lines;
	
	# getting accessions & writing
	foreach my $md5sum (keys %lines){
		my $biom = get_accession($md5sum); 	
		if( exists $biom->{'data'}[0]{'accession'}){
			foreach (@{$lines{$md5sum}}){
				push @$_, $biom->{'data'}[0]{'accession'};
			}
		}
		else{
			print STDERR "WARNING: no accession number found for line $.\n";
			foreach (@{$lines{$md5sum}}){
				push @$_, "NA";
			}
		}
		foreach (@{$lines{$md5sum}}){
			print join("\t", @$_), "\n";
		}
	}	
}



#--- Subroutines ---#
sub get_accession{
	my ($md5sum) = @_;
	
	my $base_url = "http://api.metagenomics.anl.gov/1/m5nr/md5/$md5sum";
	my $esc = join("?", "", "source=GenBank");
				
	my $url = $base_url.uri_escape($esc);

	print STDERR "Retrieving $url\n";
	my $response = $ua->get($url);
	my $content = $response->content;
	die "Error with HTTP request:  ". $response->status_line."\n".$content unless $response->is_success ;
	
	my $json = new JSON;
	my $biom = $json->decode( $content );
	return $biom;
}

