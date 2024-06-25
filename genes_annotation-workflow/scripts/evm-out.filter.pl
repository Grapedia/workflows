#! /usr/bin/env perl
my $file = shift;

my $count_abinitio=0;
my $mark=0;
my %supports;
my %details;
my %order;
my %key;
open(FILE,"<$file") or die;
while(<FILE>){
    chomp $_;
    if($_ =~ /^# EVM/){
       $_ =~ /(\d+)-(\d+)/;
       $key = $1."-".$2;
       $mark++;
       $order{$key} = $mark;
       $details{$key} .= $_."\n";
    }elsif($_ !~ /^$/){
       $details{$key} .= $_."\n";
       #-----------
       # ab initio
       #-----------
       if($_ =~ /AUGUSTUS/){
          if($supports{$key} !~ /AUGUSTUS/){
             $supports{$key} .= "AUGUSTUS;";
          }
       }
       if($_ =~ /geneid_v1.4/){
          if($supports{$key} !~ /geneid_v1.4/){
             $supports{$key} .= "geneid_v1.4;";
          }
       }
       if($_ =~ /maker/){
          if($supports{$key} !~ /maker/){
             $supports{$key} .= "maker;";
          }
       }
       if($_ =~ /Liftoff/){
          if($supports{$key} !~ /Liftoff/){
             $supports{$key} .= "Liftoff;";
          }
       }
       if($_ =~ /GlimmerHMM/){
          if($supports{$key} !~ /GlimmerHMM/){
             $supports{$key} .= "GlimmerHMM;";
          }
       }
       #-----------
       # evidence data
       #-----------
       if($_ =~ /exonerate/){
          if($supports{$key} !~ /exonerate;/){
             $supports{$key} .= "exonerate;";
          }
       }
       if($_ =~ /PsiCLASS_RNAseq_stranded/){
          if($supports{$key} !~ /PsiCLASS_RNAseq_stranded;/){
             $supports{$key} .= "PsiCLASS_RNAseq_stranded;";
          }
       }
       if($_ =~ /PsiCLASS_RNAseq_unstranded/){
          if($supports{$key} !~ /PsiCLASS_RNAseq_unstranded;/){
             $supports{$key} .= "PsiCLASS_RNAseq_unstranded;";
          }
       }
    }
}
close(FILE);

open(EVIDENCEDATA,">tmp/evm.evidencedata_only.out") or die;
open(DETAILS,">tmp/evm.out.details") or die;
open(TWO_ABINITIO,">tmp/evm.at_least_2_ABINITIO.out") or die;
open(ONE_ABINITIO,">tmp/evm.1_ABINITIO.out") or die;

foreach $ccc (sort {$order{$a} <=> $order{$b}} keys %order){
     $count_abinitio=0;
     if($supports{$ccc} =~ "AUGUSTUS;"){
       $count_abinitio++;
     }
     if ($supports{$ccc} =~ "geneid_v1.4;"){
       $count_abinitio++;
     }
     if ($supports{$ccc} =~ "maker;"){
       $count_abinitio++;
     }
     if ($supports{$ccc} =~ "Liftoff;"){
       $count_abinitio++;
     }
     if ($supports{$ccc} =~ "GlimmerHMM;"){
       $count_abinitio++;
     }
     # print STDOUT "--------------------";
     # print STDOUT "\n";
     # print STDOUT $count_abinitio;
     # print STDOUT "\n";
     # print STDOUT $details{$ccc};
     # print STDOUT "\n";
     # print STDOUT "--------------------";
     if($count_abinitio >= 2){
          print TWO_ABINITIO "$details{$ccc}\n";
     } elsif ($count_abinitio == 1 ) {
          print ONE_ABINITIO "$details{$ccc}\n";
     } else {
          print EVIDENCEDATA "$details{$ccc}\n";
     }
     print DETAILS "$ccc\t$supports{$ccc}\n";
}
close(EVIDENCEDATA);
close(DETAILS);
close(TWO_ABINITIO);
close(ONE_ABINITIO);
