#!/usr/bin/env perl 

use strict;
use warnings;
use utf8;

my ($input,$data,$versusT) = @ARGV;
defined($versusT) or die "USAGE: $0 input data 0 | 1\n";

my ($tt,$tb,$te,$mt,$mb,$me);
open(FILE,$input) or die "$0: Cannot open $input : $!\n";
while (<FILE>) {
	next if (/^#/);
	chomp;
	if (/TemperatureTotal=(.*$)/) {
		$tt = $1;
		next;
	}

	if (/TemperatureBegin=(.*$)/) {
		$tb = $1;
		next;
	}

	if (/TemperatureEnd=(.*$)/) {
		$te = $1;
		next;
	}

	if (/MuTotal=(.*$)/) {
		$mt = $1;
		next;
	}

	if (/MuBegin=(.*$)/) {
		$mb = $1;
		next;
	}

	if (/MuEnd=(.*$)/) {
		$me = $1;
		next;
	}
}

close(FILE);

(defined($te) and defined($tb) and defined($tt)
and defined($mt)) or
die "$0: One or more labels not found in $input\n";

my $ts = ($te-$tb)/$tt;
my $ms = ($me-$mb)/$mt;
my @data;
readData(\@data,$data,$tt,$mt);

if ($versusT) {
	printVersusT($tb,$ts,$tt,\@data);
} else {
	printVersusMu($mb,$ms,$mt,$tt,\@data);
}

sub printVersusT
{
	my ($tb,$ts,$tt,$data) = @_;
	for (my $i = 0; $i < $tt; ++$i) {
		my $t = sprintf("%.3f", $tb + $ts*$i);
		print "$t ";
		for (my $j = 0; $j < $mt; ++$j) {
			print "$data->[$i+$j*$tt] ";
		}

		print "\n";
	}
}

sub printVersusMu
{
	my ($mb,$ms,$mt,$tt,$data) = @_;
	for (my $i = 0; $i < $mt; ++$i) {
		my $mu = sprintf("%.3f", $mb + $ms*$i);
		print "$mu ";
		for (my $j = 0; $j < $tt; ++$j) {
			print "$data->[$j+$i*$tt] ";
		}

		print "\n";
	}
}

sub readData
{
	my ($a,$file,$rows,$cols) = @_;
	my $row = 0;
	open(FILE,$file) or die "$0: Cannot open $file : $!\n";
	while (<FILE>) {
		last if (/^$rows $cols$/);
	}

	while (<FILE>) {
		next if (/^#/);
		my @temp = split;
		my $n = scalar(@temp);
		next unless ($n == $cols);
		for (my $i = 0; $i < $cols; ++$i) {
			$a->[$row + $i*$rows] = $temp[$i];
		}

		$row++;
	}

	close(FILE);

	die "$0: Too many rows in $file\n" if ($row>$rows);
	die "$0: Too few rows in $file $row $rows\n" if ($row<$rows);
}

