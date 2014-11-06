#!/usr/bin/perl -w
$pathExe=shift;
$alignmentfileName=shift;
$treefileName=shift;
$alignmentLength=shift;
$Name=shift;
###################################
#./1-exe.pl . /Users/lindadib/UNIL/source/COEV-MODEL/COINH/RUN/alignment.fasta /Users/lindadib/UNIL/source/COEV-MODEL/COINH/RUN/tree.txt
system "mkdir RESULTS$Name/";
print "EXECUTE Coev for every PAIR OF POSITIONS (Dib et al., 2014)\n";
####################################################################################### 
$index1=1;
while($index1< $alignmentLength){
	$index2=$index1+1;
        while($index2<= $alignmentLength){
	
	$outfile="RESULTS$Name/out-".$index1."_".$index2.".log";
	system "$pathExe/coev  -method ml -tree $treefileName -align $alignmentfileName -out $outfile -cols $index1 $index2";
	$index2=$index2+1;
	}	
	$index1=$index1+1;
}
exit(0);