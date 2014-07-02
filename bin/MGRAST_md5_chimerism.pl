#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use List::Util qw/max min sum/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, @mapping_in, @lca_in, @scaffolds_in);
my $tax_level_cut = 5;
GetOptions(
	   "mapping=s{,}" => \@mapping_in,
	   "lca=s{,}" => \@lca_in,
	   "scaffolds=s{,}" => \@scaffolds_in,
	   "taxonomy=i" => \$tax_level_cut,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: the number of scaffold, mapping, and lca files must match!\n"
	unless scalar @mapping_in == scalar @lca_in && 
		scalar @mapping_in == scalar @scaffolds_in;
		
print STDERR "### taxonomic level cutoff: $tax_level_cut\n";

### MAIN
my %lens;
for my $i (0..$#scaffolds_in){
	die " ERROR: $scaffolds_in[$i] not found!\n" unless -e $scaffolds_in[$i];
	get_scaffold_lengths($scaffolds_in[$i], $i, \%lens);
	}

my %mapping;
for my $i (0..$#mapping_in){
	die " ERROR: $mapping_in[$i] not found!\n" unless -e $mapping_in[$i];
	load_mapping(\%mapping, $mapping_in[$i], $i);
	}

for my $i (0..$#lca_in){
	die " ERROR: $lca_in[$i] not found!\n" unless -e $lca_in[$i];
	my $lca_r = load_lca($lca_in[$i], $i, \%mapping);
	taxonomic_distance($lca_r, $lca_in[$i], $i, $tax_level_cut, \%lens);
	}
	
### Subroutines
sub get_scaffold_lengths{
# getting the lengths of each scaffold #
	my ($infile, $f_order, $lens_r) = @_;
	
	# status #
	print STDERR "...Loading scaffold file: $infile\n";
	
	# loading #
	open IN, $infile or die $!;
	my $seq_name;
	my $seq_len = 0;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		if(/^>/){
			$lens_r->{$f_order}{$seq_name} = $seq_len unless ! $seq_name;
			($seq_name = $_) =~ s/^>|_\[cov.+$//g;
			$seq_len = 0;
			}
		else{
			$seq_len += length($_);
			}
		$lens_r->{$f_order}{$seq_name} = $seq_len if eof IN;
		}
	close IN;
	
		#print Dumper %$lens_r; exit;
	}

sub load_mapping{
# loading AA-clustering-90% mg-rast file #
	my ($mapping_r, $infile, $f_order) = @_;
	
	# status #
	print STDERR "...loading mapping file: $infile\n";
	
	# loading #
	open IN, $infile or die $!;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @line = split /\t/;

		for my $i (1..($#line - 1)){
			$mapping_r->{$f_order}{$line[0]}{$line[$i]}++;
			die " ERROR: $infile -> $line[0] -> $line[$i] is not unique!\n" 
				if $mapping_r->{$f_order}{$line[0]}{$line[$i]} > 1;
			}
		}
	close IN;
	
		#print Dumper %$mapping_r; exit;
	return $mapping_r; 
	}

sub load_lca{
# loading mg-rast lca file #
	my ($infile, $f_order, $mapping_r) = @_;

	# status #
	print STDERR "...loading lca file: $infile\n";

	# loading #
	open IN, $infile or die $!;

	my %lca;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		
		# parsing #
		my @line = split /\t/;
		next if $line[5] =~ /^unclassified/i;
		#next if $line[5] =~ /^other/i;
		my @tax = split /;/, $line[5];						# taxonomy
		
		# loading hash #
		if($line[1] =~ /aa90/){ 			# expanding aa90 #
			if(exists $mapping_r->{$f_order}{$line[1]}){
				#print Dumper $mapping_r->{$f_order}{$line[1]}; exit;
				foreach my $scaf (keys %{$mapping_r->{$f_order}{$line[1]}}){
					my @scaf = split /_\[cov=\d+\]_/, $scaf;			# scaffold location
					#$lca{$scaf[0]}{$scaf[1]} = \@tax;
					push(@{$lca{$scaf[0]}}, \@tax);
					}
				}
			else{
				die " ERROR: $infile -> $line[1] not found in mapping file!\n";
				}
			}
		else{						# if not an aa90 cluster #
			my @scaf = split /_\[cov=\d+\]_/, $line[1];			# scaffold location
			#$lca{$scaf[0]}{$scaf[1]} = \@tax;			# scaffold_name=>AA-pos=>@taxonomy
			push(@{$lca{$scaf[0]}}, \@tax);				# scaffold_name => [AA#, taxonomy]
			}
		 
		}
	close IN;

		#print Dumper %lca; exit;
	return \%lca;
	}

sub taxonomic_distance{
# determining taxonomic distance among annotations from the same scaffold #
	my ($lca_r, $lca_in, $f_order, $tax_level_cut, $lens_r) = @_;
	
	# status #
	print STDERR "...calculating taxonomic distances\n";
	
	# taxonomic distance among annotations in a scaffold #
	foreach my $scaf (keys %$lca_r){
		my @tax_dists;
		my @min_depths;
		for my $i (0..$#{$lca_r->{$scaf}}){			# foreach annotation in scaffold
			for my $ii (0..$#{$lca_r->{$scaf}}){
				next if $i <= $ii;		# lower triangle (no-diag)
					#print Dumper $i, $ii;
				my ($dist, $min_depth) = calc_tax_distance(${$lca_r->{$scaf}}[$i], 
										${$lca_r->{$scaf}}[$ii],
										$tax_level_cut);
				push(@tax_dists, $dist);
				push(@min_depths, $min_depth);
				}
			}
		
		# checking for scaffold length existence #
		die " ERROR: cannot find '$scaf' in scaffold lengths\n" 
			unless exists $lens_r->{$f_order}{$scaf};
	
		# writting out table #
		my $scaf_tax_dist = sum_multi_tax_dist(\@tax_dists);
		my $mean_min_depth;
		if(@min_depths){
			$mean_min_depth = sum(@min_depths)/ scalar @min_depths; 
			}
		else{
			$mean_min_depth = "NA";
			}
		
		print join("\t", $lca_in, $scaf, $lens_r->{$f_order}{$scaf}, 
			scalar @{$lca_r->{$scaf}}, $scaf_tax_dist, $mean_min_depth), "\n";
		}
	}

sub sum_multi_tax_dist{
# getting sum of taxonomic distances weighted by occurance #
	my ($tax_dists_r) = @_;
	
	my %occur;
	map{ $occur{$_}++ } @$tax_dists_r;		# occurance of taxonomic distance
	
	my $sums = 0;
	map{ $sums += ($_ * $occur{$_} / scalar @$tax_dists_r) } keys %occur;		# dist_value * occurance / N-dist_values
	
	return $sums;
	}

sub calc_tax_distance{
# getting taxonomic distance between to taxonomy arrays #
	my ($tax1_r, $tax2_r, $tax_level_cut) = @_;
	
	
	# getting min taxonomic depth of comparison #
	my $min_depth = min(max_tax_depth($tax1_r), max_tax_depth($tax2_r));
	
	# reduing min_depth if greater than taxonomic cutoff (e.g. min_depth=species, but cutoff=genus) #
	$min_depth = $tax_level_cut if $tax_level_cut && $tax_level_cut < $min_depth;
	
	# comparing taxonomies #
	my $dist = 0;
	for (my $i=$min_depth; $i>=0; $i--){		# countdown from max taxonomic depth
			#print Dumper $$tax1_r[$i], $$tax2_r[$i];
		return($dist, $min_depth) if $$tax1_r[$i] eq $$tax2_r[$i]; 		# if taxonomies are the same 
		$dist++;				# number of taxonomic levels that differ between annotations
		}
	return $dist, $min_depth;							# if taxonomies array completely different
	

	# subroutines #
	sub max_tax_depth{
		# max taxonomic depth resolved (no "-") #
		my $arr_r = shift;
		for my $i (0..$#$arr_r){
			return ($i - 1) if $$arr_r[$i] eq "-";		# returning length of array up to "-"
			}
		return $#$arr_r; 				# returning length of array -1
		}
	}

__END__

=pod

=head1 NAME

MGRAST_md5_chimerism.pl -- checking for metagenome chimerism based on mg-rast annotations

=head1 SYNOPSIS

MGRAST_md5_chimerism.pl -scaf -map -lca > output.txt

=head2 options

=over

=item -scaf 

"*.preprocess.passed.fna" file(s) from MG-RAST

=item -map

"*cluster.aa90.mapping" file(s) from MG-RAST

=item -lca

"*superblat.expand.lca" files(s) from MG-RAST

=item -tax

Taxonomic level cutoff for comparing annotations. 0-indexed. (Genus is default). [5]

=item -h	This help message

=back

=head2 For more information:

perldoc MGRAST_md5_chimerism.pl

=head1 DESCRIPTION

Check for chimerism in scaffolds by determining the taxonomic distances 
among annotations on the same scaffold.

=head2 Taxonomic distance calculation:

=head3 Pairwise comparisons of annotations in a scaffold:

The number of bifurcations to the LCA with max possible bifurcations
being the max taxonomic depth of the least resolved annotation. 

=head3 Summary of pairwise comparisons

Taxonomic distance = distance_value * fraction of occurance

=head2 If multiple of each file type give (e.g. multiple genomes):

The order of files must match for '-scaf', '-map', & '-lca'

=head2 Output columns:

=over

=item * LCA file name

=item * Scaffold ID

=item * Scaffold length

=item * Number of annotations on scaffold

=item * Taxonomic distance

=item * Mean of max taxonomic resolutions for comparisons between annotations

=back

=head1 EXAMPLES

=head2 General Usage

MGRAST_md5_chimerism.pl -l 4514931.3.650.superblat.expand.lca -map 4514931.3.550.cluster.aa90.mapping -scaf 4514931.3.100.preprocess.passed.fna > md5_chimera.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/annotation/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

