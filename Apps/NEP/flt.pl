#/usr/bin/perl -w

use strict;

while(<>) {
    my ($i, $j, $k, $l, $m, $lon, $lat) = split;
    $lon += 360.;
    print "$i $j $k $l $m    $lon   $lat\n";
}
