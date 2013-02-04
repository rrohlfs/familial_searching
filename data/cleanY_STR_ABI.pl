#!/usr/local/bin/perl
#this script should read in the Y_STR_ABI file and print out files for each pop group with alleles as plain ints

#set constants
$nloci = 16;
$datafile = "Y_STR_ABI.txt";
@popnames = ("African_American", "Asian", "Asian_Indian", "Caucasian", "Chinese", "Filipino", "Hispanic", "Japanese", "Malay", "Native_American", "Sub-sharan_African", "Thai", "Vietnamese");

#read in ref haps, saving them appropriately
open(RHS, $datafile);
$headerline = <RHS>;

while ($rawline = <RHS>) {
    chomp($rawline);
    ($curpop, @curhap) = split(/\s/, $rawline);
    
    #track what alleles have been observed
    for $i (0..$#curhap) {
	$newallele = 1;
	for $j (0..$#{$alleles[$i]}) {
	    if($alleles[$i][$j] == $curhap[$i]) {
		$newallele = 0;
	    }
	}
	if($newallele == 1) {
	    push(@{$alleles[$i]}, $curhap[$i]);
	}
    }

    #save hap in right place
    push(@{$refHaps{$curpop}}, ${@curhap});
}

close(RHS);

#code alleles as plain ints

for $curloci (0..$nloci) {
    $curalleles = sort($[@alleles[$curloci]);
    for $i (0..$#curalleles) {
	$curallelehash{$curalleles[$i]} = $i;
    }

    for $curpop (@popnames) {
	for $indivi (0..$#[$refHaps{$curpop}]) {
	    codeRefHaps{$curpop}[$curloci][$indivi] = $curallelehash{$refHaps{$curpop}[$curloci][$indivi]};
    }
    
}




#write allele freqs to files
#for $j (0..$#locinames) { 
#    $locusname = $locinames[$j];
#    $locusfile = $cumdirname.'/'.$locusname.'.dat';
#    open(OUT, ">$locusfile");
    
#    print(OUT "$locusname\t$subgroupnames\n");
#    for $k (0..($nalleles[$j]-1)) {
#	print(OUT "-1");
#	for $i (0..(2*$#popnames+1)) {
#	    print(OUT "\t$AFs[$i][$j][$k]");
#	}
#	print(OUT "\n");
#    }
#    close(OUT);
#}
