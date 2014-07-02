#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use DBI;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_in, $header_bool, $overlap);
my ($coding_in, $mapping_in, $filter_in);
my $db_name = "MGRAST_db";
my $m5nr_repeat = 3;
my $sleep_sec = 120;					# 2 min between repeat 
GetOptions(
	   "db=s" => \$database_in,
	   "coding=s" => \$coding_in,
	   "mapping=s" => \$mapping_in,
	   "filter=s" => \$filter_in,
	   "name=s" => \$db_name,			# new sqlite3 db name
	   "overlap" => \$overlap, 			# partial overlap or complete overlap of blast hits? [complete]
	   "x" => \$header_bool,			# header? [TRUE]
	   "repeat=i" => \$m5nr_repeat, 	# if nothing returned, repeat X times
	   "sleep=i" => \$sleep_sec, 		# sleep time between repeats
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide either an sqlite database or 3 MG-RAST tables: coding, mapping, & md5 filter\n"
	unless $database_in || ($coding_in && $mapping_in && $filter_in);

### MAIN
# making database #
unless($database_in){
	print STDERR "### making an sqlite3 DB ###\n" unless $verbose;
	make_tables($db_name);
	add_CDS_table($coding_in, $db_name);
	add_mapping_table($mapping_in, $db_name);
	add_no_mapped($db_name);
	add_filter_table($filter_in, $db_name);
	$database_in = $db_name;
	}

# loading database & blast #
my $dbh = connect2db($database_in);
my $blast_res_r = load_blast();

# querying sqlite3 db #
my $md5_r = query_db($blast_res_r, $dbh, $overlap);

# querying m5nr #
my %m5nr_res; my $sleep_bool = 0;
my $m5nr_res_r = call_m5nr($md5_r, $m5nr_repeat, \%m5nr_res, $sleep_sec, $sleep_bool);

# writing out appended blast table #
write_m5nr_res_table($m5nr_res_r, $blast_res_r, $header_bool);

$dbh->disconnect();

### subroutines
sub add_filter_table{
	my ($filter_in, $db_name) = @_;
	
	
	open IN, $filter_in or die $!;
	my $pipe = "| sqlite3 $db_name '.imp \"/dev/stdin\" \"Filter\"' ";
	print STDERR $pipe, "\n" unless $verbose;
	open PIPE, $pipe or die $!;

	my @last;	
	while(<IN>){	
		chomp; 
		s/\t[\d.]+%.*$//; 

		print PIPE join("|", split /\t/), "\n";
		}
	close IN;
	close PIPE;
	}

sub add_no_mapped{
# adding all 'singletons' to cluster table #
	my ($db_name) = @_;
	
	my $dbh = connect2db($db_name);
	
	# getting CDS_ids w/ no cluster #	
	my $q = "select CDS.cds_id from CDS left outer join cluster on CDS.CDS_id = cluster.cds_id where cluster.cluster_id is null";
	my $res =  $dbh->selectall_arrayref($q);
	
	$dbh->disconnect();
	
	
	# adding to cluster table #
	my $pipe = "| sqlite3 $db_name '.imp \"/dev/stdin\" \"Cluster\"' ";
	print STDERR $pipe, "\n" unless $verbose;
	open PIPE, $pipe or die $!;
	
	foreach my $r (@$res){
		foreach my $rr (@$r){
			print PIPE join("|", $rr, $rr), "\n";
			}
		}	
	close PIPE;
	}

sub add_mapping_table{
# adding cluster table #
	my ($mapping_in, $db_name) = @_;
	
	
	open IN, $mapping_in or die $!;
	my $pipe = "| sqlite3 $db_name '.imp \"/dev/stdin\" \"Cluster\"' ";
	print STDERR $pipe, "\n" unless $verbose;
	open PIPE, $pipe or die $!;

	my @last;	
	while(<IN>){	
		chomp; 
		s/\t[\d.]+%.*$//; 
		
		my @line=split /[\t,]/; 
		for my $i (1..$#line){ 
			print PIPE join("|", @line[($i, 0)]), "\n";
			}
		}
	close IN;
	close PIPE;
	}

sub add_CDS_table{
# adding CDS table from sequence data (genes called) #
	my ($coding_in, $db_name) = @_;
	
	open IN, $coding_in or die $!;
	my $pipe = "| sqlite3 $db_name '.imp \"/dev/stdin\" \"CDS\"' ";
	print STDERR $pipe, "\n" unless $verbose;
	open PIPE, $pipe or die $!;

	my @last;	
	while(<IN>){
		chomp;
		next if /^\s*$/;
		
		if(/^>/){
			s/>//; 
			@last=split /_/;
			} 
			
		else{
			print PIPE join("|", join("_", @last), 
				  join("_", @last[0..3]), @last[4..$#last], $_), "\n";
			}
		}
	close IN;
	close PIPE;
	
	}

sub make_tables{
	my $db_name = shift;

# making tables in sqlite db #
my $sql = <<HERE;
BEGIN TRANSACTION;

DROP TABLE IF EXISTS CDS;

CREATE TABLE CDS (
CDS_id	TEXT	PRIMARY KEY,
scaffold_id	TEXT	NOT NULL,
start	INTEGER	NOT NULL,
end	INTEGER	NOT NULL,
strand	TEXT	NOT NULL,
seq_AA	TEXT
);


DROP TABLE IF EXISTS Cluster;

CREATE TABLE Cluster (
CDS_id	TEXT	NOT NULL,
Cluster_id	TEXT	NOT NULL
);


DROP TABLE IF EXISTS Filter;

CREATE TABLE Filter (
query_id	TEXT,
hit_id	TEXT	NOT NULL,
percentage_identity	REAL	NOT NULL,
alignment_length	INTEGER	NOT NUll,
number_mismatches	INTEGER	NOT NUll,
number_gaps_openings	INTEGER	NOT NUll,
q_start	INTEGER	NOT NUll,
q_end	INTEGER	NOT NUll,
s_start	INTEGER	NOT NUll,
s_end	INTEGER	NOT NUll,
evalue	TEXT	NOT NULL,
bit_score	INTEGER	NOT NUll
);

COMMIT;

HERE
	
	unlink $db_name or die $! if -e $db_name;
	open PIPE, " | sqlite3 $db_name" or die $!;
	print PIPE $sql;
	close PIPE;
	
	}

sub write_m5nr_res_table{
	my ($m5nr_res_r, $blast_res_r, $header_bool) = @_;

	print STDERR "### appending M5nr info to blast table ###\n" unless $verbose;

	# getting all uniq categories #
	my %cat;
	foreach my $subject (keys %$m5nr_res_r){			
		foreach my $query (keys %{$m5nr_res_r->{$subject}}){
			foreach my $row (keys %{$m5nr_res_r->{$subject}{$query}}){
				foreach my $md5 (keys %{$m5nr_res_r->{$subject}{$query}{$row}}){
					foreach my $source (keys %{$m5nr_res_r->{$subject}{$query}{$row}{$md5}}){
						foreach my $cat (keys %{$m5nr_res_r->{$subject}{$query}{$row}{$md5}{$source}}){
							next if $cat eq "NA";
							$cat{$cat} = 1;
							}
						}
					}
				}
			}
		}
	my @cat = sort keys %cat;
	
	# writing header #
	print join("\t", qw/qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore m5d source/, 
			@cat), "\n" unless $header_bool;
	
	# writing table #
	foreach my $subject (keys %$m5nr_res_r){			
		foreach my $query (keys %{$m5nr_res_r->{$subject}}){
			foreach my $row (keys %{$m5nr_res_r->{$subject}{$query}}){
				foreach my $md5 (keys %{$m5nr_res_r->{$subject}{$query}{$row}}){
					foreach my $source (keys %{$m5nr_res_r->{$subject}{$query}{$row}{$md5}}){
						my @vals;
						foreach my $cat (@cat){
							if(exists $m5nr_res_r->{$subject}{$query}{$row}{$md5}{$source}{$cat}) {
								push @vals, $m5nr_res_r->{$subject}{$query}{$row}{$md5}{$source}{$cat};
								}
							else{
								push @vals, "NA";
								}
							}
						
						#print Dumper @vals;
					
						# writing blast line w/ appended info #
						print join("\t", 
							@{$blast_res_r->{$subject}{$query}{$row}},
							$md5, $source, @vals), "\n";
							
						}
					}
				}
			}
		}

		#print Dumper %$m5nr_res_r; exit;
	}

sub parse_m5nr_cgi{
# http://api.metagenomics.anl.gov/m5nr/m5nr_rest.cgi/ output #
	my ($res) = @_;
	$res =~ s/^\[|]$|\n|\r//g;
	my @parts = split /{|}/, $res;
	
	my %m5nr;
	foreach my $pp (@parts){
		next unless $pp;
		next if $pp eq ",";
		
		my @pp = split /","/, $pp;			# "tag:value"
		my %tmp;
		my $source;
		my $md5; 
		foreach my $ppp (@pp){
			$ppp =~ s/"//g;
			my @ppp = split /:/, $ppp;
			if($ppp[0] eq "source"){
				$source = $ppp[1];
				}
			elsif($ppp[0] eq "md5"){
				$md5 = $ppp[1];
				}
			else{
				$tmp{$ppp[0]} = $ppp[1];
				}
			}
			
		$m5nr{$md5}{$source} = \%tmp;
		}
		#print Dumper %m5nr; exit;
	return \%m5nr;
	}
	
sub parse_m5nr_api{
# http://api.metagenomics.anl.gov/m5nr/md5/ output #
	my ($res) = @_;
	
	my @init = split /\[|\]/, $res;	
	my @parts = split /{|}/, $init[1];
	
	my %m5nr;
	foreach my $pp (@parts){
		next unless $pp;
		next if $pp eq ",";
		
		my @pp = split /"*,"/, $pp;			# parsing to: "tag:value"
		
		my %tmp;
		my $source;
		my $md5; 
		foreach my $ppp (@pp){
			$ppp =~ s/"//g;
			my @ppp = split /:/, $ppp;
			
			if($ppp[0] eq "source"){
				$source = $ppp[1];
				}
			elsif($ppp[0] eq "md5"){
				$md5 = $ppp[1];
				}
			else{
				$tmp{$ppp[0]} = $ppp[1];
				}
			}
			
		$m5nr{$md5}{$source} = \%tmp;
		}
		
		#print Dumper %m5nr; exit;
	return \%m5nr;
	}

sub call_m5nr{
# querying m5nr database by using RESTFUL api #
	my ($md5_r, $m5nr_repeat, $m5nr_res_r, $sleep_sec, $sleep_bool) = @_;
	
	print STDERR "### Calling M5nr with RESTFUL ###\n" unless $verbose;
	print STDERR "#-- Repeat countdown: $m5nr_repeat --#\n" if $sleep_bool && ! $verbose;
	
	foreach my $subject (keys %$md5_r){
		foreach my $query (keys %{$md5_r->{$subject}}){
			foreach my $row (keys %{$md5_r->{$subject}{$query}}){
				if( ${$md5_r->{$subject}{$query}{$row}}[0] ne "NA" ){			# if blast hit overlaps w/ m5d

					# skipping if m5nr annotations already found in past repeat #
					next if exists $m5nr_res_r->{$subject}{$query}{$row} && 
						! exists $m5nr_res_r->{$subject}{$query}{$row}{"NA"};
					
					# using RESTFUL app #
					my $cmd = join("", "lynx -dump \"http://api.metagenomics.anl.gov/m5nr/md5/",
						join(";", @{$md5_r->{$subject}{$query}{$row}} ), "\"");
					print STDERR "$cmd\n" unless $verbose;
					my $res = `$cmd`;
				
					# checking & loading m5nr annotations #
					if($res && $res =~ /\{.*\[.*\{.*\}.*\].*\}/){							# if RESTFUL app returns something
						$m5nr_res_r->{$subject}{$query}{$row} = parse_m5nr_api($res);	
						}
					elsif($m5nr_repeat == 0){			# all attempts failed
						print STDERR " WARNING: all attempts for query \"$cmd\" have failed! Using 'NA'\n";
						$m5nr_res_r->{$subject}{$query}{$row}{"NA"}{"NA"}{"NA"} = "NA";
						}
					else{							# no hit 
						$m5nr_res_r->{$subject}{$query}{$row}{"NA"}{"NA"}{"NA"} = "NA";
						$sleep_bool = 1;
						}
					}
				else{		# if no blast hit
					$m5nr_res_r->{$subject}{$query}{$row}{"NA"}{"NA"}{"NA"} = "NA";
					}
				}
			}
		}
	
	# repeating this run; recursive subroutine #
	$m5nr_repeat--;
	if($m5nr_repeat > 0 && $sleep_bool){
		print STDERR "...Sleeping for $sleep_sec seconds before next repeat of m5nr api queries\n";
		sleep $sleep_sec;
		call_m5nr($md5_r, $m5nr_repeat, $m5nr_res_r, $sleep_sec, $sleep_bool);
		}
		
		#print Dumper %$m5nr_res_r, "here"; exit;
	return $m5nr_res_r;		# cds => m5d => source => cat => value
	}

sub query_db{
# querying sqlite3 db of MG-RAST download files #
	my ($blast_res_r, $dbh, $overlap) = @_;

	print STDERR "### Querying the sqlite3 DB ###\n" unless $verbose;
	
	my $sth;
	if($overlap){
		$sth = $dbh->prepare(q{SELECT hit_id FROM CDS LEFT OUTER JOIN cluster ON CDS.CDS_id = Cluster.CDS_id LEFT OUTER JOIN filter ON cluster.cluster_id = filter.query_id WHERE CDS.scaffold_id = ? and (CDS.start BETWEEN ? and ? OR CDS.end BETWEEN ? OR ?)});
		}
	else{
		$sth = $dbh->prepare(q{SELECT hit_id FROM CDS LEFT OUTER JOIN cluster ON CDS.CDS_id = Cluster.CDS_id LEFT OUTER JOIN filter ON cluster.cluster_id = filter.query_id WHERE CDS.scaffold_id = ? and CDS.start <= ? and CDS.end >= ?});
		}
	
	my %md5;
		#my $cnt = 0;
	foreach my $subject (keys %$blast_res_r){		# scaffold=1, start=8, stop=9
		foreach my $query (keys %{$blast_res_r->{$subject}}){
			foreach my $row (keys %{$blast_res_r->{$subject}{$query}}){
				# status on query #
				print STDERR "db_query= query \"$query, subject ${$blast_res_r->{$subject}{$query}{$row}}[1], blast_subject_start ${$blast_res_r->{$subject}{$query}{$row}}[8], blast_subject_end ${$blast_res_r->{$subject}{$query}{$row}}[9]\"\n" 
					unless $verbose;
				$sth->bind_param(1, ${$blast_res_r->{$subject}{$query}{$row}}[1]);
				
				my $blast_start = ${$blast_res_r->{$subject}{$query}{$row}}[8];
				my $blast_end = ${$blast_res_r->{$subject}{$query}{$row}}[9];
				
				# flipping blast start & end if other direction #
				if($blast_start > $blast_end){
					my $tmp = $blast_start;
					$blast_start = $blast_end;
					$blast_end = $tmp;
					}
				
				# binding parameters #
				$sth->bind_param(2, $blast_start);
				if($overlap){			# BETWEEN: 5 '?'
					$sth->bind_param(3, $blast_start);
					$sth->bind_param(4, $blast_end);
					$sth->bind_param(5, $blast_end);
					}
				else{					# just 3 values needed
					$sth->bind_param(3, $blast_end);
					}
				$sth->execute();	
				while( my $r = $sth->fetchrow_arrayref()){						
					my @new_r = @$r;
					if($new_r[0]){
						$md5{$subject}{$query}{$row} = \@new_r;			# subject=>query=>blast = m5d
						}
					else{
						$md5{$subject}{$query}{$row} = [("NA")];
						}
					}			
				}
			}

			#$cnt++;
			#last if $cnt >= 2;		#<-- debuggin
		}
	$sth->finish();
		
		#print Dumper %md5; exit;
	return \%md5;		# CDS_id => m5d_id
	}

sub connect2db{
	my ($database_in) = @_;
	my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
	my $dbh = DBI->connect("dbi:SQLite:dbname=$database_in", '','', \%attr) 
		or die " Can't connect to $database_in: $DBI::errstr";
	return $dbh;
	}

sub load_blast{
# loading blast results #
# subject should have m5d query names #
	my %blast_res;
	while(<>){
		chomp;
		next if /^\s*$/;
		
		my @line = split /\t/;
		$blast_res{$line[1]}{$line[0]}{$.} = \@line;
		}
		#print Dumper %blast_res; exit;
	return \%blast_res;
	}



__END__

=pod

=head1 NAME

MGRAST_blast2md5.pl -- add md5 info to a tab-delimited blast table

=head1 SYNOPSIS

=head2 Creating, then querying an sqlite3 DB from MG-RAST files

MGRAST_blast2md5.pl -c *350.genecalling.coding.faa -m *550.cluster.aa90.mapping -f *.650.superblat.sims.filter < blastn.txt > blastn_md5.txt

=head2 Querying an existing sqlite3 DB

MGRAST_blast2md5.pl -db MGRAST_db < blastn.txt > blastn_md5.txt

=head2 Required flags

=over

=item -coding

*350.genecalling.coding.faa file from MG-RAST

=item -mapping

*550.cluster.aa90.mapping file from MG-RAST

=item -filter

*.650.superblat.sims.filter file from MG-RAST

=item -db

sqlite3 database created by MGRAST_blast2md5.pl

=back

=head2 Optional flags

=over

=item -name

Database name (if creating a sqlite3 database with '-c' '-m' '-f'). [MGRAST_db]

=item -overlap

Complete overlap required for blast hits to md5 gene calls? [TRUE] 

=item -repeat

Number of times to repeat M5nr query with API. [3]

=item -sleep

Number of seconds to sleep between repeat M5nr queries (if repeats needed). [120]

=item -x

Add a header to the blast table? [TRUE]

=item -v	Verbose output. [TRUE]

=item -h	This help message

=back

=head2 For more information:

perldoc MGRAST_blast2md5.pl

=head1 DESCRIPTION

Appends md5 database info (MG-RAST) onto a blast table.
Helps for blasting against contigs/scaffolds annotated by
MG-RAST. The blast hits are duplicated so that only one M5nr
'source' are appended onto the end of the row.

Blast hits must fall into a CDS called by MG-RAST in order to 
return M5nr annotation information.

"NA" means no md5 was found.

The input blast table is assumed to be in '-outfmt 6' format.

=head1 EXAMPLES

=head2 Creating, then querying an sqlite3 DB from MG-RAST files

MGRAST_blast2md5.pl -c 4514930.3.350.genecalling.coding.faa -m 4514930.3.550.cluster.aa90.mapping -f 4514930.3.650.superblat.sims.filter < blastn.txt > blastn_md5.txt

=head2 Querying an existing sqlite3 DB

MGRAST_blast2md5.pl -db MGRAST_db < blastn.txt > blastn_md5.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/annotation/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

