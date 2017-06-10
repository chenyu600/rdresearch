#!/usr/bin/perl -w
use strict;
use File::Basename;
my $help="Usage: perl $0 <OT|out.bam>\n";

die "Usage:\nperl $0 <OT|path2out.bam> \nFor Example:\nperl $0 pe.bam\n" unless @ARGV>0;

#my ($sam,$out)=@ARGV[0,1];
#my ($out)=@ARGV;
#my ($out)=@ARGV;
my $out=shift;
my $out_file = basename $out;
my($filename, $dirs, $suffix) = fileparse($out);
my $err_file = "$dirs/abnormal.sam";
my $samtools="/path/samtools";
my $fai="/path/ucsc.hg19.fasta.fai";

open OT,"| $samtools view -bS - > $out" or die $!;
#open EOT,"| $samtools view -bS - > $err_file" or die $!;
open EOT,"> $err_file" or die $!;

#open OT,">$out/pe.filter.sam" or die $!;
#if($sam=~/\.gz$/){open SM,"gunzip -cd $sam|" or die $!;}else{open SM,"<$sam" or die $!;}
#A80HNWABXX:5:6:16084:133982#ATCACGAT    147     chrY    6633    0       90M     =       6565    -158    AAGCCATTTCCTTCCTTCCTTCATTCCTTCCTTCCTCCCTCCCTCCCTCCTTCCCTCCCTCCCTTTTTTTTTTCAGGGTCTTGCTCTGTC       FF?DFBFF=EEEBED:ADFDFF?FBFFDEFDDBDD;@DD;?EE?EFG?GGC=@EE3ADD5EFE>EFGGFGGGGGGGGGEGFGGGGGFGGG      XT:A:R  NM:i:0  SM:i:0  AM:i:0  X0:i:2  X1:i:0  XM:i:0  XO:i:0   XG:i:0  MD:Z:90 XA:Z:chrY,-6633,90M,0;
#while(my $line=<>){
while(<>){
    chomp;
    if (/^\@/){
        print OT "$_\n";
        print EOT "$_\n";
        next;
    }
# FCC3U20ACXX:5:1101:14786:2241#GGCAACAG  83      chr21   45105704        60      100M    =       45105550        -254    TTTTCGTGAATATGGACCCCACTGGGGGCGTTGATGCCAAACGCTGAGAAACAAGATTGAAAGAACTTTACTGAAAACCTCTGATAAGAAATAATAGTCC    ccccacbddccdca^\ccccdccbagghihhiihhiiiihgiiiiiiiihfiiiiiiiiiiihhhihfhghiiihgfhiiiiiiihhgggegeeeeebbb    NM:i:0  AS:i:100        XS:i:0  RG:Z:140318_I712_FCC3U20ACXX_L5
    my @t = (split /\t/);
    my ($qname, $flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual, @tags) = @t;
    my $field_num = scalar @t;
    my ($cigar_len, $seq_len, $qua_len) = (0,0,0);
    print EOT "$_\tXE:Z:field_num_lt_11\n" if $field_num < 11;
    next if $field_num < 11;

    print EOT "$_\tXE:Z:flag_not_num\n" if $flag !~ /^\d+$/;
    next if $flag !~ /^\d+$/;

    print EOT "$_\tXE:Z:pos_not_num\n"  if $pos !~ /^\d+$/;
    next if $pos !~ /^\d+$/;

    print EOT "$_\tXE:Z:bases_not_stored\n" if $seq eq '*';
    next if $seq eq '*';
#    print EOT "$_\tXE:Z:pos_eq0\n" if $pos == 0;
#    next if $pos == 0;

    print EOT "$_\tXE:Z:mapq_of_umaped_reads_ne0\n" if  ($flag & 4)==4 && $mapq !=0;
    next if ($flag & 4)==4 && $mapq !=0;

#    map{$cigar_len += $_}(split /\D/,$cigar);
    $cigar_len = length($cigar);
    $seq_len = length($seq);
    $qua_len = length($qual);

    print EOT "$_\tXE:Z:seq_len_lt10\n" if $seq_len <= 10;
    next if $seq_len <= 10;

    print EOT "$_\tXE:Z:cigar_len_eq0\n" if $cigar_len == 0;
    next if $cigar_len == 0;

    print EOT "$_\tXE:Z:seq_len_ne_qua_len\n" if $seq_len != $qua_len;
    next if $seq_len != $qua_len;

 #   print EOT "seq_len ne cigar_len\t$_\n" if $seq_len != $cigar_len;
 #   next if $seq_len != $cigar_len;;

    print OT "$_\n";
}
close OT;
close EOT;
