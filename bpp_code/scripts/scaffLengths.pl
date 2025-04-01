# libraries
use strict; # help
use warnings;

## input parameters
my $fasta = $ARGV[0];
#my $prefix=$ARGV[1];


## merge scaffolds in same chromosome
open('fa', '<', $fasta) or die "can't open $fasta $!\n";
while (my $line = <fa>) {
    chomp($line);
    ## determine scaffold name
    if ($line =~ /^\>/) {
        my @fields = split(" ", $line);
        my $chrname = $fields[0];
        print "$chrname\t";
    }
    else {
        my $len = length $line;
        print "$len\n";
    }
}
