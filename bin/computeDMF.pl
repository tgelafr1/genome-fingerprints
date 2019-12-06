use strict;
use List::Util qw(min);
use List::Util qw(max);
use List::Util qw(any);
use List::Util qw(sum);
use List::MoreUtils qw(uniq);
my $version = '181217';
####
#
# This software computes a genome fingerprint for a single personal genome.
# The method is described in:
#    Glusman G, Mauldin DE, Hood LE, Robinson M. Ultrafast Comparison of Personal
#    Genomes via Precomputed Genome Fingerprints. Front Genet. 2017 Sep 26;8:136. doi:
#    10.3389/fgene.2017.00136. eCollection 2017. PubMed PMID: 29018478; PubMed Central
#    PMCID: PMC5623000.
# 
# Copyright 2017 by Gustavo Glusman, Institute for Systems Biology, Seattle, WA, USA.
# It is provided by the Institute for Systems Biology as open source software,
# free for non-commercial use.
# See the accompanying LICENSE file for information about the governing license.
#
####
#
# Accepted input formats include VCF, BCF, RCF (ISB's range call format), and Complete Genomics' var and masterVar formats.
# The first parameter is an 'id' for the job; output files will use this as base. You can include in it a path to where you want the output files to be located.
# The second parameter is the input file (the genome as VCF, RCF, var or masterVar).
# The third (optional) parameter is the number of consecutive SNPs to use (integer)
# The fourth (optional) parameter is the distance metric to use between SNPs.
# The fifth (optional) parameter is the format of the input file: 'vcf', 'bcf', 'rcf', 'var' or 'masterVar'. Defaults to 'vcf'.
# The sixth (optional) parameter is the fingerprint size. Multiple sizes can be specified, comma-delimited. Defaults to including several sizes.
# The seventh (optional) parameter is the distance between consecutive SNVs that are considered 'too close'. Default is 20.
# The eight (optional) parameter is a bed file specifying regions of interest to be included in the analysis. For example, one could specify the definition of exome segments to compute an exome-compatible fingerprint from whole-genome data. This is available only for VCF and RCF input.
#
####
#
# Examples of usage:
#   computeDMF.pl myGenome vcfs/myGenome.vcf.gz
#   computeDMF.pl fingerprints/anotherGenome vcfs/aGenome.vcf.gz vcf 5,20,120 20 exomeRegions.bed
#
####

# Set default values for missing inputs (2 is the default SNP # to use, mean
# the default distance metric, vls is Ls to try, 20 is the default C,
# and vcf is the default format)
my($id, $file, $k, $dist, $format, $L, $C, $bedmask) = @ARGV;
my @vls = (5, 7, 11, 13, 17, 19, 20, 40, 50, 80, 100, 120, 200);
@vls = split /,/, $L if $L;	
$C ||= 20;
$format ||= 'vcf';
$dist ||= 'mean';
$k ||= 2;
my @keys = makekeys($k);

# Sanitize inputs.
die if $file =~ /[\s\;]/;
die if $bedmask =~ /[\s\;]/;

# Preparation.
my $rawFileExt = 'out';
my $closeFileExt = 'out.close';
my $normFileExt = 'outn';

my($prevChrom,$prevprevChrom, $prevStart, $prevprevStart, $prevKey,$prevprevKey, %count, %close, %binary, $snvPairs);
my $cat = 'cat';
if ($file =~ /\.gz$/) {
	$cat = 'gunzip -c';
} elsif ($file =~ /\.bz2$/) {
	$cat = 'bzcat';
} elsif ($file =~ /\.bcf$/) {
	$cat = 'bcftools view';
}
my %filter = (
	'vcf' => "grep -v ^\# | grep -v \"0[\/\|]0\"",
	'bcf' => "grep -v ^\# | grep -v \"0[\/\|]0\"",
	'rcf' => 'grep -v ^\#',
	'masterVar' => 'grep -v no-call',
	'var' => 'grep -v ref | grep -v no-call');
my %fields = (
	'vcf' => '1,2,4,5', 'bcf' => '1,2,4,5', 'rcf' => '1,2,5,6',
	'masterVar' => '3,4,8,9,10', 'var' => '4,5,8,9');

# Process input file.
if ($bedmask && ($format eq 'vcf' || $format eq 'rcf')) {
	open INF, "bedtools intersect -a $file -b $bedmask | $filter{$format} | cut -f$fields{$format} |";
} else {
	open INF, "$cat $file | $filter{$format} | cut -f$fields{$format} |";
}

my @prevChroms;
my @prevStarts;
my @prevKeys;

while (<INF>) {
	chomp;
	my($chrom, $start, $ref, $var, $othervar) = split /\t/;

	# Focus the analysis on autosomes only, excluding sex chromosomes, mitochondrial chromosome, alternative haplotypes, etc.
	next unless $chrom =~ /^(chr)?\d+$/;
	
	# Interpret othervar statement from masterVar.
	$var = $othervar if $othervar && $var eq $ref;
	
	# Pay attention only to SNVs.
	next unless $var =~ /^[ACGT]$/i && $ref =~ /^[ACGT]$/i && uc $var ne uc $ref;
	
	# Compute the key for the current SNV.
	my $key = uc "$ref$var";
	my @D;
	my $d;
	if (sameChroms($k, $chrom, @prevChroms)) {

		@D = getDs(@prevStarts, $start);

		next if any {$_ < 0} @D;

		if ($dist eq "mean") {
		 	$d = int(sum(@D)/@D);
		} elsif($dist eq "max") {
			$d = max(@D);
		} elsif($dist eq "min") {
			$d = min(@D);
		} else {
			print("Bad distance function");
			die;
		}	

		# Compute the key for the pair.
		my $pairKey = join "", @prevKeys, $key;
		
		# Add to table, segregating by $C, by pair key and by reduced distance.
		if ($d < $C) {
		 	$close{$pairKey}[$d]++;
		} else {
		 	$binary{$pairKey}[$d % 2]++;
		 	$count{$_}{$pairKey}[$d % $_]++ foreach @vls;
		}
		$snvPairs++;
	}
	
	# Store info on current SNV for next round, remove info on old SNV
	@prevChroms = updateList($k, @prevChroms, $chrom);
	@prevStarts = updateList($k, @prevStarts, $start);
	@prevKeys = updateList($k, @prevKeys, $key);
}
close INF;

my $bin = join("", map {$binary{$_}[1]>$binary{$_}[0] ? 1 : 0} @keys);

my @headers = (
	['#software-version', $version],
	['#source', $file],
	['#k', $k],
	['#distanceMetric', $dist],
	['#format', $format],
	['#SNVpairs', $snvPairs],
	['#vectorLengths', join("\t", @vls)],
	['#tooCloseCutoff', $C],
	['#created', `date`],
	);
my $header = join("\n", map {join("\t", @{$_})} @headers);

# Output main fingerprint table.
open OUTF, ">$id.$rawFileExt";
print OUTF $header;
print OUTF "#binary\t$bin\n";
foreach my $vl (@vls) {
	foreach my $key (@keys) {
		print OUTF join("\t", $vl, $key, map {$count{$vl}{$key}[$_] || 0} (0..$vl-1)), "\n";
	}
}
close OUTF;

# Output secondary (short distance) fingerprint table.
open OUTF, ">$id.$closeFileExt";
print OUTF $header;
foreach my $key (@keys) {
	print OUTF join("\t", $key, map {$close{$key}[$_] || 0} (0..$C-1)), "\n";
}
close OUTF;

# Normalize the fingerprints.
foreach my $vl (@vls) {
	# Normalize fingerprint per reduced distance.
	foreach my $col (0..$vl-1) {
		my @v = ();
		push @v, $count{$vl}{$_}[$col] foreach @keys;
		my($avg, $std) = avgstd(\@v);
		$std ||= 1;
		$count{$vl}{$_}[$col] = ($count{$vl}{$_}[$col]-$avg)/$std foreach @keys;
	}
	# Normalize fingerprint per SNV pair key.
	foreach my $sig (@keys) {
		my($avg, $std) = avgstd($count{$vl}{$sig});
		$std ||= 1;
		$_ = ($_-$avg)/$std foreach @{$count{$vl}{$sig}};
	}
}

# Output normalized fingerprint table.
open OUTF, ">$id.$normFileExt";
print OUTF $header;
print OUTF "#binary\t$bin\n";
foreach my $vl (@vls) {
	foreach my $key (@keys) {
		print OUTF join("\t", $vl, $key, map {sprintf("%.3f", $_)} @{$count{$vl}{$key}}), "\n";
	}
}
close OUTF;

# Compress output files.
`gzip -f $id.$rawFileExt; gzip -f $id.$closeFileExt; gzip -f $id.$normFileExt`;


###
sub avgstd {
	my($values) = @_;
	my($sum, $devsqsum);

	my $n = scalar @$values;
	return unless $n>1;
	$sum += $_ foreach @$values;
	my $avg = $sum / $n;
	$devsqsum += ($_-$avg)**2 foreach @$values;
	my $std = sqrt($devsqsum/($n-1));
	return $avg, $std;
}

# A function that makes the key values based on how many consecutive pairs of
# SNPs are being used.
#
# Args: k	The number of consecutive SNPs to use.
#
# Returns: A list of character values corresponding to the keys.
sub makekeys {
	my @vars = ('AC', 'AT', 'AG',
				'CA', 'CT', 'CG',
				'TA', 'TC', 'TG',
				'GA', 'GC', 'GT');
	my $k = $_[0];
	my @oldkeys = @vars;
	my @keys;

	for (my $i = 1; $i < $k; $i++) {
		@keys = ();
		for my $a (@vars) {
			for my $b (@oldkeys) {
				push @keys, $a.$b;
			}
		}
		@oldkeys = @keys;
	}
	return @keys;
}

# A function to get pairwise distances between consecutive elements of a list.
# 
# Args: prevStarts		A list of all the start values of the k-1 previous SNPs.
# 		start			The start value of the current SNP.
# 
# Returns: A list of k-1 values corresponding to the distances between
#		   consecutive elements of the input items.
sub getDs {
	my $start = shift;
	my @prevStarts = @_;
	my @D;
	my $currStart = $start;

	for my $s (@prevStarts) {
		push @D, $s - $currStart - 1;
		$currStart = $s;
	}

	return @D;
}

# Check if all SNPs to be analyzed are on the same chromosome.
#
# Args: k	The number consecutive SNPs.
#		_	The array of SNP chromosome values.
#
# Return: A boolean indicating whether or not all the chromosomes are the same.
sub sameChroms {
	my $k = shift;
	my $array_length = @_;
	if ($k > $array_length) {
		return 0;
	}

	my @uniq_array = uniq @_;
	my $uniq_length = @uniq_array;

	$array_length != 1 && $uniq_length == 1 ? return 1 : return 0;	
}

# A function to pop the first element off a list if it's too long (has more than
# k elements)
#
# Args: k		The number of elements to keep in the list.
#		_		The old list.
#
# Returns: The new list.
sub updateList{
	my $k = shift;
	my $l = @_;

	if ($l >= $k) {
		shift;
	}
	return @_;
}