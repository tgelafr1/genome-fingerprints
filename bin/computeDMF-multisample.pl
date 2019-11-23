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
	print "Usage: computeDMF-multisample.pl inputDirectory [fingerprintLengths] [outputDirectory]\n";
	print "       fingerprintLengths defaults to 120; outputDirectory defaults to DMF_multisample_outdir\n";
	print "Examples: computeDMF-multisample.pl dirWithVCFs\n";
	print "          computeDMF-multisample.pl dirWithVCFs 20,100\n";
	print "          computeDMF-multisample.pl dirWithBCFs 0 desiredOutputDir\n";
	exit;
}

my @vls = (120);
@vls = split /,/, $vls if $vls;
my $tooCloseCutoff = 20;
my $computeBinary = 0;
my $computeClose = 0;
my $numberOfPartitions = 100;
$out ||= "DMF_multisample_outdir";
$dist ||= 'mean';
$k ||= 2;
mkdir $out, 0700;
my @keys = makekeys($k);

FILE: foreach my $file (slicedirlist($dir, "[bv]cf")) {
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
	while (<F>) {
		chomp;

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
		
		foreach my $i (0..$#obs) {

			my $gt = $obs[$i];
			next if $gt =~ /0.0/ || $gt =~ /\./;

			my @D;
			my $d;
			if (enoughKeys($k, $i, @prevKeys)) {
				@D = getDs($prevStarts[$i], $start);

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

				my $pairKey = join "", @{$prevKeys[$i]}, $key;
				if ($d<$tooCloseCutoff) {
					$close[$i]{$pairKey}[$d]++ if $computeClose;
				} else {
					$binary[$i]{$pairKey}[$d % 2]++ if $computeBinary;
					foreach my $vl (@vls) {
						$count[$i]{$vl}{$pairKey}[$d % $vl]++;
					}
				}
			}
			# @prevChroms = updateList($k, @prevChroms, $chrom);
			$prevStarts[$i] = updateList($k, $start, $prevStarts[$i]);
			$prevKeys[$i] = updateList($k, $key, $prevKeys[$i]);
			$snvPairs[$i]++;
		}
		
		$lines++;
		#last if $lines>10000;
		print "." unless $lines % 10000;
		print " " unless $lines % 1000000;
	}
	close F;
	
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

sub slicedirlist {
	my($dir, $pat) = @_;
	opendir (DIR, $dir);
	my @files = grep /$pat/, readdir DIR;
	closedir DIR;
	return @files;
}

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

sub updateList {
	my @poop = @_;
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