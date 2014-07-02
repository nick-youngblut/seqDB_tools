#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $username, $password, $sub_in);
GetOptions(
	   "username=s" => \$username,
	   "password=s" => \$password,
	   "submission=s" => \$sub_in,
	   "xample" => \&write_example,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );
	   

### I/O error & defaults
die " ERROR: provide a username and password!\n" if ! $username || ! $password;
die " ERROR: provide a job submission table!\n" if ! $sub_in;
die " ERROR: $sub_in doesn't exist!\n" if ! -e $sub_in;

#my $jobfh = make_job_fh($sub_in);

### MAIN
my $sub_tbl_r = load_submission_table($sub_in);

foreach my $job (sort {$a <=> $b} keys %$sub_tbl_r){
	my @line;
	foreach my $class (keys %{$$sub_tbl_r{$job}}){
		push(@line, join(" ", "--$class", $$sub_tbl_r{$job}{$class}));
		}
	my $cmd = join(" ", "svr_submit_RAST_job", "--user $username", "--passwd $password", @line);
	print STDERR "$cmd\n";
	system($cmd);
	
	}

#close $jobfh;


### Subroutines
sub make_job_fh{
	my $sub_in = shift;
	
	(my $job_name = $sub_in) =~ s/\.[^\.]+$|$/_jobid.txt/;
	open my $jobfh, ">$job_name" or die $!;
	
	print STDERR " Job IDs will be written to '$job_name\n";
	
	return $jobfh;
	}

sub load_submission_table{
# loading a submission table for batch job submission to RAST #
	my ($sub_in) = @_;
	
	open IN, $sub_in or die $!;

	my %sub_tbl;
	my %header;
	while(<IN>){
		chomp;
		s/#.+//g;
		next if /^\s*$/;
		my @line = split /\t/;

		if($. == 1){			# header
			for my $i (0..$#line){
				$header{$line[$i]} = $i;
				}
			}
		else{					# loading body 
			# fasta | genbank #
			if(exists $header{"fasta"} && exists $header{"genbank"}){
				die " ERROR: provide file names under 'fasta' OR 'genbank' column header\n";
				}
			if(exists $header{"fasta"}){
				$sub_tbl{$.}{"fasta"} = $line[$header{"fasta"}];
				}
			elsif(exists $header{"genbank"}){
				$sub_tbl{$.}{"genbank"} = $line[$header{"genbank"}];			
				}
			else{
				die " ERROR: provide file names under 'fasta' or 'genbank' column header\n";
				}
				
			# domain #
			if(exists $header{"domain"}){
				$sub_tbl{$.}{"domain"} = $line[$header{"domain"}];
				}
			else{
				die " ERROR: provide domain [Archaea | Bacteria ] under 'domain' column header\n";
				}
			
			# bioname #
			if(exists $header{"genus"} && exists $header{"species"} && exists $header{"strain"}){
				$sub_tbl{$.}{"bioname"} = join("", '"', 
					join(" ", $line[$header{"genus"}], $line[$header{"species"}], 
						$line[$header{"strain"}]), '"'); 
				}
			else{
				print STDERR " WARNING: no bioname information provided. Using default\n";
				}
				
			# genetic code #
			if(exists $header{"genetic_code"}){
				$sub_tbl{$.}{"genetic_code"} = $line[$header{"genetic_code"}];
				}
			else{ 
				die " ERROR: provide a genetic code number [11 | 4] under 'genetic_code' column name\n";
				}
			
			# gene caller #
			if(exists $header{"gene_caller"}){
				$sub_tbl{$.}{"gene_caller"} = $line[$header{"gene_caller"}];
				}
			else{
				die " ERROR: provide a gene caller under 'gene_caller' column name\n";
				}
			}
		}

	close IN;
	
		#print Dumper %sub_tbl; exit;
	return \%sub_tbl;
	}

sub write_example{
 my $tbl = <<HERE;
fasta	domain	genus	species	strain	genetic_code	gene_caller
file1.fna	Archaea	Methanosarcina	mazei	Go1	11	rast
file2.fna	Archaea	Methanosarcina	mazei	TMA	11	rast
file3.fna	Archaea	Methanosarcina	mazei	LYC	11	rast

HERE
	print $tbl;
    exit(1);
	}


__END__

=pod

=head1 NAME

RAST_batch_submit_job.pl -- submit multiple jobs to RAST

=head1 SYNOPSIS

=head2 Batch submission

RAST_batch_submit_job.pl -u -p -s

=head2 Write example submission table

RAST_batch_submit_job.pl -x > example.txt

=head2 options

=over

=item -u 	Username for RAST server

=item -p 	Password for RAST server

=item -s 	Submission table (*txt)

=item -x 	Write an example submission table

=item -h	This help message

=back

=head2 For more information:

perldoc RAST_batch_submit_job.pl

=head1 DESCRIPTION

Batch upload of multiple assemblies for RAST annotation.

=head2 Submission table format guidelines

=over

=item *

The table requires a header.

=item *

Table header labels must match svr_submit_RAST_job options (case-sensitive),
except for --bioname, which has been split into 'genus', 'species', & 'strain' columns.

=item *

If an option is not required, remove the column from the table.

=item *

A boolean option just requires a column header.

=item *

The 'fasta' or 'genbank' columns require the assembly file name with the full path from
the current working directory.

=back

=head1 EXAMPLES

=head2 Example table

RAST_batch_submit_job.pl -x > submission.txt

=head2 Batch submission

RAST_batch_submit_job.pl -u USERNAME -p PASSWD -s submission.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

