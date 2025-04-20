#!/usr/bin/env -S perl -w

my $debug = 1;

my $randomSeed   = 0;
my $badnessMinOverlap = 4; # as suggested in the paper
my $maxIteration = 60000;
my $maxTolIter   = 20000;
my $tolC         = 100;
my %outIters     = ();
my $badnessEval  = "";

my $usage = "Usage: saddle_pair.pl [options] <genePairedPrimerFile> <outPrefix>\n";
$usage   .= "           -rand       : random seed (default: $randomSeed)\n";
$usage   .= "           -minOverlap : minimum reversecomplement to compute badness (default: $badnessMinOverlap)\n";
$usage   .= "           -maxIt      : maximum iteration (default: $maxIteration)\n";
$usage   .= "           -maxTolIt   : maximum tolerance iteration (default: $maxTolIter)\n";
$usage   .= "           -tolC       : factor C for tolerance threshold (default: $tolC)\n";
$usage   .= "           -outIter    : iterations to output answers (default: none)\n";
$usage   .= "           -badness    : no iteration and compute badness of iteration outputs (default: no)\n";

# optional parameter
my @arg_idx=(0..@ARGV-1);
for my $i (0..@ARGV-1) {
	if ($ARGV[$i] eq '-rand') {
		$randomSeed = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}elsif ($ARGV[$i] eq '-minOverlap') {
		$badnessMinOverlap = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}elsif ($ARGV[$i] eq '-maxIt') {
		$maxIteration = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}elsif ($ARGV[$i] eq '-maxTolIt') {
		$maxTolIter = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}elsif ($ARGV[$i] eq '-tolC') {
		$tolC = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}elsif ($ARGV[$i] eq '-outIter') {
		$outIters{$ARGV[$i+1]} = 1;
		delete @arg_idx[$i,$i+1];
	}elsif ($ARGV[$i] eq '-badness') {
		$badnessEval = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}
}
my @new_arg;
for (@arg_idx) { push(@new_arg,$ARGV[$_]) if (defined($_)); }
@ARGV=@new_arg;

if(length($badnessEval)>0){
    open(FILE,"<$badnessEval");
    my @strings = ();
    while(<FILE>){
        my @t=split;
        shift @t;
        push @strings, (shift @t);
        push @strings, (shift @t);
    }
    close FILE;

    my %nullCache = ();
    print "".badness_total(\@strings, $badnessMinOverlap, \%nullCache)."\n";
    
    exit 0;
}

# regular parameters
my $primerFile = shift or die $usage;
my $outPrefix  = shift or die $usage;

# init random seed
srand($randomSeed);

# read file
my %genePrimerPairs = (); # $genePrimerPair{gene} = [[p11, p12], [p21, p22], ...]

open(FILE,"<$primerFile");
while(<FILE>){
    chomp;
    my @t=split;
    
    # very simple format checking, assuming 3 tokens
    if(3 != @t){
        print STDERR "not 3 tokens in one line: $_\n";
        exit 1;
    }
    
    my $g = shift @t;
    push @{$genePrimerPairs{$g}}, [$t[0], $t[1]];
}
close FILE;

# cache for pair badness answers
my %pairBadnessCache = (); # $pairBadness{p1, p2} = badness

# randomized initial set and its badness
my %answerSet = (); # $answerSet{gene} = index of answer
my $answerSetRef = \%answerSet;

for my $g (sort keys %genePrimerPairs){
    $answerSet{$g} = int(rand(@{$genePrimerPairs{$g}}));
}

# iteration 0
my $iteration = 0;
my $current_badness = badness_total(collect_primers_pair(\%genePrimerPairs,$answerSetRef), $badnessMinOverlap, \%pairBadnessCache);

# iteration output, if specified
if(exists $outIters{$iteration}){
    open(iterOUT, ">$outPrefix.$iteration");
    for my $g (sort keys %answerSet){
        my @primers = @{${$genePrimerPairs{$g}}[$answerSetRef->{$g}]};
        print iterOUT "$g\t$primers[0]\t$primers[1]\n";
    }
    close(iterOUT);
}

open(scoreOUT,">$outPrefix.score");
print scoreOUT "$iteration\t$current_badness\n";

#### SADDLE iteration

my $stop = 0;
while(not $stop){
    # alter one primer pair, randomly
    my $newAnswerSetRef = answerSet_alterOne(\%genePrimerPairs,$answerSetRef);

    # badness
    my $new_badness = badness_total(collect_primers_pair(\%genePrimerPairs,$newAnswerSetRef), $badnessMinOverlap, \%pairBadnessCache);

    # decide accept or not
    if( $new_badness < $current_badness){
        # accept
        $current_badness = $new_badness;
        $answerSetRef = $newAnswerSetRef;
        
    }else{
        if($iteration <= $maxTolIter){
            my $toleranceTh = exp( ($current_badness - $new_badness)/($tolC) );
            
            if(rand() < $toleranceTh){
                # accept
                $current_badness = $new_badness;
                $answerSetRef = $newAnswerSetRef;
            }
        }
    }

    $iteration++;
    print scoreOUT "$iteration\t$current_badness\n";
    
    # iteration output, if specified
    if(exists $outIters{$iteration}){
        open(iterOUT, ">$outPrefix.$iteration");
        for my $g (sort keys %answerSet){
            my @primers = @{${$genePrimerPairs{$g}}[$answerSetRef->{$g}]};
            print iterOUT "$g\t$primers[0]\t$primers[1]\n";
        }
        close(iterOUT);
    }

    # stop check
    if($iteration >= $maxIteration){
        $stop = 1;
    }
}

#### SADDLE iteration, END

# subroutine: answer set manipulation, alter one primer pair
sub answerSet_alterOne {
    my $genePrimerPairsRef = shift;
    my $answerSetRef = shift;
    
    my %newAnswerSet = %{$answerSetRef}; # should be a clone

    # randomly pick a gene
    my @genes = sort keys %{$genePrimerPairsRef};
    my $geneAltered = rand(@genes);

    # next generation, pick another primer
    my $oldPrimerIdx = $newAnswerSet{$genes[$geneAltered]};
    my $newPrimerIdx = 0;
    if(1==@{$genePrimerPairsRef->{$genes[$geneAltered]}}){
        # no need to random pick
    }else{
        while($newPrimerIdx == $newAnswerSet{$genes[$geneAltered]}){
            $newPrimerIdx = rand(@{$genePrimerPairsRef->{$genes[$geneAltered]}});
        }
    }
    $newAnswerSet{$genes[$geneAltered]} = $newPrimerIdx;
    
    return \%newAnswerSet;
}

# subroutine: collect all primers
sub collect_primers_pair {
    my $genePrimerPairsRef = shift;
    my $answerSetRef = shift;

    my @strings = ();

    for my $g (sort keys %{$genePrimerPairsRef}){
        my @primers = @{${$genePrimerPairsRef->{$g}}[$answerSetRef->{$g}]};
        push @strings, $primers[0], $primers[1];
    }
    
    return \@strings;
}

# subroutines for badness
sub badness_total {
    my $arrRef = shift;
    my $minOverlap = shift;
    my $cacheRef = shift;

    my $ans = 0;
    for(my $i=0; $i<@{$arrRef}; $i++){
        for(my $j=0; $j<=$i; $j++){
            $ans += badness_pair($$arrRef[$i], $$arrRef[$j], $minOverlap, $cacheRef);
        }
    }
    
    return $ans;
}

sub badness_pair {
    my $p1 = uc(shift);
    my $p2 = uc(shift);
    my $minOverlap = shift;
    my $cacheRef = shift;
    
    if(exists $cacheRef->{$p1, $p2}){
        return $cacheRef->{$p1, $p2};
    }else{
        my $ans = badness_pair_compute($p1, $p2, $minOverlap);
        
        $cacheRef->{$p1, $p2} = $ans;
        $cacheRef->{$p2, $p1} = $ans;
        return $ans;       
    }
}

sub badness_pair_compute {
    my $p1 = shift;
    my $p2 = shift;
    my $minOverlap = shift;

    # p2 reverse complement
    $p2 = reverse($p2);
    $p2 =~ tr/ACGT/TGCA/;
    
    # establish min-mer hash for overlap search
    my @testArr;
    my %wordHash = (); # $wordHash{min-mer}{i} = pos's
    @testArr = ($p1, $p2);
    for(my $i=0; $i<@testArr; $i++){
        for(my $j=0; $j+$minOverlap<=length($testArr[$i]); $j++){
            push @{$wordHash{substr($testArr[$i], $j, $minOverlap)}{$i}}, $j;
        }
    }
    # establish a dot matrix based on min-mer matches
    my %dotHash = (); # $dotHash{x, y} = 1 if an min-mer match
    for my $minmer (keys %wordHash){
        my $cnt = keys %{$wordHash{$minmer}};
        if($cnt==2){
            for(my $i=0; $i<@{$wordHash{$minmer}{0}}; $i++){
                for(my $j=0; $j<@{$wordHash{$minmer}{1}}; $j++){
                    $dotHash{ $wordHash{$minmer}{0}[$i], $wordHash{$minmer}{1}[$j] } = 1;
                }
            }
        }
    }
    # generate all diagonal
    my %diagonals = (); # $diagonals{x, y, len} = 1;
    
    my $nDots = keys %dotHash;
    while($nDots){
        my @dots = sort {
                my @arrA=split(/$;/,$a);
                my @arrB=split(/$;/,$b);
                $arrA[0]<=>$arrB[0] || $arrA[1]<=>$arrB[1]
            } keys %dotHash;
        my @xy0 = split(/$;/, $dots[0]);

        my @xy = @xy0;
        delete $dotHash{$xy[0], $xy[1]};
        my $len = 4;

        while(exists $dotHash{$xy[0]+1, $xy[1]+1}){
            $xy[0]++;
            $xy[1]++;
            delete $dotHash{$xy[0], $xy[1]};
            $len++;
        }
        $diagonals{$xy0[0], $xy0[1], $len} = 1;

        $nDots = keys %dotHash;
    }
    # for each diagonal
    my $ans = 0;
    for my $diag (keys %diagonals){
        my @diagInfo = split(/$;/, $diag);
        my $subRC = substr($p1, $diagInfo[0], $diagInfo[2]);
        
        # count GC's
        my $gc = () = $subRC =~ /[GC]/g;
        # d1 and d2
        my $d1 = length($p1) - ($diagInfo[0]+$diagInfo[2]);
        my $d2 = length($p2) - ($diagInfo[1]+$diagInfo[2]);
        
        $ans += ((2**$diagInfo[2])*(2**$gc)) / (($d1+1)*($d2+1));
    }
    
    return $ans;
}
