#!/usr/bin/env perl

=pod

=head1 NAME

RAST_tblFromTarBall.pl -- pull out 'FIG-ID/Features/peg/tbl' file from RAST tarball

=head1 SYNOPSIS

RAST_tblFromTarBall.pl [flags] *tar

=head2 Required flags

NONE

=over

=back

=head2 Optional flags

=over

=item -f 	Number of files to process in parallel. [1]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc RAST_tblFromTarBall.pl

=head1 DESCRIPTION

Pull out the 'FIG-ID/Features/peg/tbl' file, which
has the RAST FIG-PEG IDs.

=head1 EXAMPLES

=head2 Basic usage (4 files processed in parallel):

RAST_tblFromTarBall.pl -f 4 *tar

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

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
use File::Spec;
use Archive::Tar;
use File::Path qw/rmtree/;
use Parallel::ForkManager;

#--- args/flags ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose_b);
my $fork = 0;
GetOptions(
	"fork=i" => \$fork,
	"verbose" => \$verbose_b,
	"help|?" => \&pod2usage # Help
	);

#--- I/O error ---#
my $pm = new Parallel::ForkManager($fork);

#--- MAIN ---#
foreach my $infile (@ARGV){
	my $pid = $pm->start and next; 
	
	die "ERROR: cannot find '$infile'\n"
		unless -e $infile;
	
	my $tar = Archive::Tar->new($infile);
	
	my @tarfiles = $tar->list_files();
	
	(my $tmp = $infile) =~ s/\.(tar|tgz)$//i;
	
	my $tbl_file = "$tmp/Features/peg/tbl";
	die "ERROR: cannot find '$tbl_file'"
		unless grep /$tbl_file$/, @tarfiles;
	
	# extracting file to folder #
	$tar->extract_file($tbl_file);
	print STDERR "'$tbl_file' extracted to '$tmp\_tbl'\n";

	rename($tbl_file, "$tmp\_tbl") or die $!;
	rmtree($tmp) or die $!;
	
	$pm->finish
	}
$pm->wait_all_children;


#--- Subroutines ---#



