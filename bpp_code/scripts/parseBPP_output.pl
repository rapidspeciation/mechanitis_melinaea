## load libraries ******
use strict; # help you to keep track of variable names and scope

#/ arguments
my $bpp = $ARGV[0]; # bpp .out file to read

#/ read bpp oputput file and get info
if ($bpp =~ /.gz$/) { open(file, "gunzip -c $bpp |") || die "can't open pipe to $bpp"; }
else { open(file, $bpp) or die "can't open $bpp $!\n"; }

my $startRead = 0;
my @trees;

while (my $line = <file>) {
    chomp($line);
    # check if reach the part of the file with the trees and once you do change $startRead state [to active]
    if ( $line =~ /^(\Q(A)\E)/ ) { $startRead = 1; }
    else {
        if ( $startRead == 1) {
            # reach section (B)
            if ( $line =~ /^(\Q(B)\E)/ )  { $startRead = 0; }
            # still reading trees 
            else {
                # remove leading spaces, unnecessary spaces in tree and final ;
                $line=~ s/^\s+//;
                $line=~ s/,\s/,/g;
                $line=~ s/;//g;
                $line =~ s/ +/ /;
                chomp $line;
                # print "$line\n";
                # empty line
                if ( $line =~ /^\s/) { next; }
                # non-empty line
                else {
                    my @vals = split('\s', $line);
                    my $tre = $vals[4];
                    my $cnt = $vals[0];
                    my $ln = "$tre\t$cnt";
                    push(@trees, $ln);
                } 
            }
        }
    }
}

# last element
my $l_el = scalar @trees - 2;
my @extracted_elements = @trees[0..$l_el];
# print scalar @trees, "\n";
foreach (@extracted_elements) {
    print "$bpp\t";
    print "$_\n";
}
# print join("\n", @extracted_elements),"\n";