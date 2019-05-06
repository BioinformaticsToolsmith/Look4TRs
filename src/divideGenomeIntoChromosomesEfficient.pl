#! /usr/bin/perl -w

#
# Author: Hani Zakaria Girgis, PhD
#
# This program read one FASTA file that includes multiple sequences
# Write each sequence to a separate file.
#

use strict;

# DM6
# my $file      = '/Users/zakarota/Data/Genomes/Dm6/dm6.fa';
# my $directory = '/Users/zakarota/Data/Genomes/Dm6/Fa';

# my $file      = '/Users/zakarota/Data/Genomes/GlycineMax/Gmax_109.fa';
# my $directory = '/Users/zakarota/Data/Genomes/GlycineMax/Fa';

#my $file = '/Users/zakarota/Data/Genomes/TGACv1/Triticum_aestivum.TGACv1.dna.toplevel.fa';
#my $directory = '/Users/zakarota/Data/Genomes/TGACv1/Fa';


die "Usage: $0 file output_directory" unless $#ARGV==1; 

my ($file, $directory) = @ARGV;

drive();

sub drive {
	open( IN, '<', $file ) or die "Cannot open $file: $!\n";
	my $seq  = "";
	my $info = "";
	while (<IN>) {
		if ( $_ =~ m/^>/ ) {
			if ( $seq ne "" ) {
				writeFastaFile( $info, $seq );
				$seq = "";
			}
			$info = $_;
		}
		else {
			$seq = $seq . $_;
		}
	}

	# Last sequence
	writeFastaFile( $info, $seq );
	close(IN) || die $!;
}

sub writeFastaFile {
	my ( $info, $seq ) = @_;

	$info =~ m/\>(.+)\s?/;
	my @t = split( /\s/, $1 );
	my $seqFile = "$directory/$t[0].fa";

	print "Writting to $seqFile\n";

	open( OUT, '>', $seqFile ) or die "Cannot open $seqFile: $!\n";

	# Commented out on 3/19/2019
	# Now, the header is the same as the file name
	#print OUT $info;
	print OUT ">$t[0]\n";
	print OUT $seq;

	close(OUT);
}
