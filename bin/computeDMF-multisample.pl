#!/usr/bin/env perl
$|=1;
use strict;
use List::Util qw(min);
use List::Util qw(max);
use List::Util qw(any);
use List::Util qw(sum);
use List::MoreUtils qw(uniq);
my $version = '181217';

my($dir, $k, $dist, $vls, $out) = @ARGV;

# First parameter: the directory where the vcfs/bcfs are located. [cwd]
# The second (optional) parameter is the number of consecutive SNPs to use (integer)
# The third (optional) parameter is the distance metric to use between SNPs.
# Fourth parameter (optional): the comma-delimited list of fingerprint lengths to compute. [120]
# Fifth parameter (optional): the directory where output should be saved. [DFM_multisample_outdir]

unless (length($dir) && -e $dir) {
	print "Usage: computeDMF-multisample.pl inputDirectory [consecutiveSNPs] [distanceMetric] [fingerprintLengths] [outputDirectory]\n";
	print "       consecutiveSNPs defaults to 2; distanceMetric defaults to 'mean'; fingerprintLengths defaults to 120; outputDirectory defaults to DMF_multisample_outdir\n";
	print "Examples: computeDMF-multisample.pl dirWithVCFs\n";
	print "			 computeDMF-multisample.pl dirWithVCFs 3 min\n";
	print "          computeDMF-multisample.pl dirWithVCFs 2 mean 20,100\n";
	print "          computeDMF-multisample.pl dirWithBCFs 2 mean 0 desiredOutputDir\n";
	exit;
}

my @vls = (120);
@vls = split /,/, $vls if $vls;
my $tooCloseCutoff = 20;
my $computeBinary = 0;
my $computeClose = 0;
my $numberOfPartitions = 100;

# Set defaults
$out ||= "DMF_multisample_outdir";
$dist ||= 'mean';
$k ||= 2;
mkdir $out, 0700;

# Make possible base pair keys, dependent on how many consecutive SNPs are being
# used.
my @keys = makekeys($k);

FILE: foreach my $file (slicedirlist($dir, "[bv]cf")) {
	# Not modified by us, this code just reads in the files.
	my $cat;
	if ($file =~ /\.gz$/) {
		$cat = 'gunzip -c';
	} elsif ($file =~ /\.bcf$/) {
		$cat = 'bcftools view';
	} elsif ($file =~ /\.vcf$/) {
		$cat = 'cat';
	} elsif ($file =~ /\.bz2$/) {
		$cat = 'bzcat';
	} else {
		print "skipping file: $file\n";
		next;
	}
	
	open F, "$cat $dir/$file |";
	my @samples;
	while (<F>) {
		next if /^##/;
		if (/^#CHROM/) {
			chomp;
			(undef, undef, undef, undef, undef, undef, undef, undef, undef, @samples) = split /\t/;
			#print "Found ", scalar @samples, " samples\n";
			last;
		}
	}
	
	my($currentChromosome, $outdir, $lines, @prevStarts, @prevKeys, @close, @count, @binary, @snvPairs);

	# While the file is open...
	while (<F>) {
		chomp;

		# Read in line information, make sure it is well formed, create 
		# appropriate output directory.
		my($chrom, $start, $rsid, $ref, $alt, $qual, $filter, $infostring, $format, @obs) = split /\t/;
		$chrom = "chr$chrom" unless $chrom =~ /^chr/;
		if (-e "$out/$chrom") {
			if ($chrom ne $currentChromosome && defined($currentChromosome)) {
				close F;
				next FILE;
			}
		} elsif ($currentChromosome) {
			die "Unexpected chromosome $chrom when working on $currentChromosome\n";
		} else {
			$outdir = "$out/$chrom";
			mkdir $outdir, 0700;
			$currentChromosome = $chrom;
			print "(", scalar @samples, " samples) $chrom";
		}
		
		next if $rsid =~ /CNV/;
		next unless $alt =~ /^[ACGT]$/io && $ref =~ /^[ACGT]$/io;
		my $key = uc "$ref$alt";
		
		# For each sample in the file...
		# This is where the bulk of our implementation appears!
		foreach my $i (0..$#obs) {

			my $gt = $obs[$i];
			next if $gt =~ /0.0/ || $gt =~ /\./;

			my @D;
			my $d;
			# Only compute the metric if you've already read in k keys.
			if (enoughKeys($k, $i, @prevKeys)) {

				# Get the distances between each pair of SNPs
				@D = getDs($prevStarts[$i], $start);

				# If any of the distances are negative, stop, somethings wrong.
				next if any{$_ < 0} @D;

				# Calculate the final distance metric based on user input.
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

				# Create the key for the set of SNPs by concatenating all of the
				# individual keys
				my $pairKey = join "", @{$prevKeys[$i]}, $key;

				# If the distance is too close, add to "close" matrix
				if ($d<$tooCloseCutoff) {
					$close[$i]{$pairKey}[$d]++ if $computeClose;
				
				# Else add to fingerprint
				} else {
					$binary[$i]{$pairKey}[$d % 2]++ if $computeBinary;
					foreach my $vl (@vls) {
						$count[$i]{$vl}{$pairKey}[$d % $vl]++;
					}
				}
			}

			# Update lists containing information about previous SNPs within 
			# a range of k.
			$prevStarts[$i] = updateList($k, $start, $prevStarts[$i]);
			$prevKeys[$i] = updateList($k, $key, $prevKeys[$i]);

			$snvPairs[$i]++;
		}
		
		# Debugging output
		$lines++;
		#last if $lines>10000;
		print "." unless $lines % 10000;
		print " " unless $lines % 1000000;
	}
	close F;
	
	# This code was not modified by us! It writes the saved fingerprint to the
	# output file.
	foreach my $i (0..$#samples) {
		my $sample = $samples[$i];
		my $partition = partition($sample);
		mkdir "$outdir/$partition", 0700;
		
		if ($computeBinary) {
			open OUTF, ">$outdir/$partition/$sample.binary";
			print OUTF "#source\t$file\n";
			print OUTF "#tooCloseCutoff\t$tooCloseCutoff\n";
			print OUTF "#SNVpairs\t$snvPairs[$i]\n";
			foreach my $key (@keys) {
				print OUTF join("\t", $key, map {$binary[$i]{$key}[$_] || 0} (0,1)), "\n";
			}
			close OUTF;
		}
		
		open OUTF, ">$outdir/$partition/$sample.out";
		print OUTF "#source\t$file\n";
		print OUTF "#tooCloseCutoff\t$tooCloseCutoff\n";
		print OUTF "#SNVpairs\t$snvPairs[$i]\n";
		foreach my $vl (@vls) {
			foreach my $key (@keys) {
				print OUTF join("\t", $vl, $key, map {$count[$i]{$vl}{$key}[$_] || 0} (0..$vl-1)), "\n";
			}
		}
		close OUTF;
		
		if ($computeClose) {
			open OUTF, ">$outdir/$partition/$sample.close";
			print OUTF "#source\t$file\n";
			print OUTF "#tooCloseCutoff\t$tooCloseCutoff\n";
			print OUTF "#SNVpairs\t$snvPairs[$i]\n";
			foreach my $key (@keys) {
				print OUTF join("\t", $key, map {$close[$i]{$key}[$_] || 0} (0..$tooCloseCutoff-1)), "\n";
			}
			close OUTF;
		}
		
		`gzip -f $outdir/$partition/$sample.*`;
		print "o" unless ($i+1) % 1000;
		print " " unless ($i+1) % 10000;
	}
	print "\n";
}

# Not written by us.
sub slicedirlist {
	my($dir, $pat) = @_;
	opendir (DIR, $dir);
	my @files = grep /$pat/, readdir DIR;
	closedir DIR;
	return @files;
}

# Not written by us.
sub partition {
	my($name) = @_;
	
	my $v = $name;
	$v =~ s/\D//g;
	unless (length($v)) {
		$v += ord() foreach split //, $name;
	}
	$v = $v % $numberOfPartitions;
	$v = "0$v" if length($v)<2;
	return $v;
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
# Args: prevStartRef	A reference to a list of all the start values of the
#						k-1 previous SNPs
# 		start			The start value of the current SNP.
# 
# Returns: A list of k-1 values corresponding to the distances between
#		   consecutive elements of the input items.
sub getDs {
	my $prevStartRef = shift;
	my $start = shift;
	my @prevStarts = @{$prevStartRef};
	my @D;
	my $currStart = $start;

	for my $s (reverse(@prevStarts)) {
		push @D, $currStart - $s - 1;
		$currStart = $s;
	}

	return @D;
}

# A function to pop the first element off a list and push a new element onto the
# end of a list. These lists keep track of the k SNP elements we are 
# computing across. The function won't remove elements if the list is less than
# k elements long.
#
# Args: k		The number of elements to keep in the list.
#		new 	The new element to add to the list.
#		ref		A reference to the old list.
#
# Returns: A reference to the new list.
sub updateList {
	my $k = shift;
	my $new = shift;
	my $ref = shift;

	if (not defined($ref)) {
		my $out = [$new];
		return $out;
	}

	my @old = @{$ref};
	my $l = @old;

	if ($l >= $k-1) {
		shift @old;
	}

	push @old, $new;
	return \@old;
}

# A function to check if enough keys have been collected to begin calculating
# distance metrics.
#
# Args: k	The number of consecutive SNPs being evaluated.
#		i	The index of the sample currently being evaluated.
#		_	A list of references to previous keys collected. Essentially an
#		array of arrays which stores all of previous keys extracted for each
#		sample.
#
# Returns: a boolean indicating whether or not the key list is long enough to
#          proceed.
sub enoughKeys {
	my $k = shift;
	my $i = shift;
	my $ref = @_[$i];
	if (not defined($ref)) {
		return 0;
	}
	my @prev = @{$ref};
	my $l = @prev;
	return $l + 1 >= $k;
}