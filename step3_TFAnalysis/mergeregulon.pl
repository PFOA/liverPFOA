use strict;
use warnings;

my $reg=shift;

open IN,$reg;
while(<IN>){
	chomp;
	my $gene=`basename  $_ "_regulon.txt"`;
	chomp($gene);
	open F,$_;
	while(<F>){
		print "$gene\t$_";
	}
}
