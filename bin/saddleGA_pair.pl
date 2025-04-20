#!/usr/bin/env -S perl -w

use AI::Genetic::Pro;

$| = 1;

my $randomSeed   = 0;
my $badnessMinOverlap = 4; # as suggested in the paper
my $maxIteration = 600;
my $numReport    = 5;

# GA parameters
my $ga_population = 200;
my $ga_crossover  = 0.9;
my $ga_mutation   = 0.05;

my $usage = "Usage: saddleGA_pair.pl [options] <genePairedPrimerFile> <outPrefix>\n";
$usage   .= "           -rand       : random seed (default: $randomSeed)\n";
$usage   .= "           -minOverlap : minimum reversecomplement to compute badness (default: $badnessMinOverlap)\n";
$usage   .= "           -maxIt      : maximum iteration (default: $maxIteration)\n";
$usage   .= "           -report     : number of top answers to report (default: $numReport)\n";
$usage   .= "           -GApopulation : GA parameter population (default: $ga_population)\n";
$usage   .= "           -GAcrossover  : GA parameter crossover (default: $ga_crossover)\n";
$usage   .= "           -GAmutation   : GA parameter mutation (default: $ga_mutation)\n";

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
	}elsif ($ARGV[$i] eq '-report') {
		$numReport = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}elsif ($ARGV[$i] eq '-GApopulation') {
		$ga_population = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}elsif ($ARGV[$i] eq '-GAcrossover') {
		$ga_crossover = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}elsif ($ARGV[$i] eq '-GAmutation') {
		$ga_mutation = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}
}
my @new_arg;
for (@arg_idx) { push(@new_arg,$ARGV[$_]) if (defined($_)); }
@ARGV=@new_arg;

# regular parameters
my $primerFile = shift or die $usage;
my $outPrefix  = shift or die $usage;

# init random seed
srand($randomSeed);

# read file
my %genePrimerPairs = (); # $genePrimerPair{gene} = [[p11, p12], [p21, p22], ...]
my $genePrimerPairsRef = \%genePrimerPairs;

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


#### SADDLE GA iteration

my $ga = AI::Genetic::Pro->new(
    -fitness         => \&fitness,
    -type            => 'rangevector',
    -population      => $ga_population,
    -corssover       => $ga_crossover,
    -mutation        => $ga_mutation,
    -preserve        => 10,
    -variable_length => 0,
    -parents         => 2,
    -selection       => [ 'RouletteBasic' ],
    -strategy        => [ 'Points', 2 ],
    -native          => 1,
    -history         => 1,
# it seems that the current MCE implementation is not feasible for rangevector, so 0
    -mce             => 0,
    -workers         => 3,
);

# GA init
my @initArray = ();
for my $g (sort keys %{$genePrimerPairsRef}){
    push @initArray, [0, @{$genePrimerPairsRef->{$g}}-1];
}

$ga->init(\@initArray);

my $iteration = 0;
my ($max, $mean, $min) = $ga->getAvgFitness();
$max = 1/$max;
$min = 1/$min;

open(scoreOUT,">$outPrefix.score");
scoreOUT->autoflush(1);
print scoreOUT "$iteration $max $min\n";

while($iteration<$maxIteration){
    $ga->evolve(1);
    $iteration++;
    ($max, $mean, $min) = $ga->getAvgFitness();
    $max = 1/$max;
    $min = 1/$min;
    
    print scoreOUT "$iteration $max $min\n";
}
close scoreOUT;

# GA output
my @bests = $ga->getFittest($numReport, 1); # 1 for uniq answers

for(my $i=0; $i<@bests; $i++){
    my @indexes = $ga->as_array($bests[$i]);
    my @genes   = sort keys %{$genePrimerPairsRef};

    open(iterOUT,">$outPrefix.best$i");
    for(my $i=0; $i<@genes; $i++){
        my @primers = @{${$genePrimerPairsRef->{$genes[$i]}}[$indexes[$i]]};
        print iterOUT "$genes[$i]\t$primers[0]\t$primers[1]\n";
    }
    close iterOUT;
}
    
#### SADDLE GA iteration, END

sub fitness {
    my $ga = shift;
    my $chromosome = shift;
    
    # return inverse of basness
    my @indexes = $ga->as_array($chromosome);
    my @genes   = sort keys %{$genePrimerPairsRef};

    my @strings = ();
    for(my $i=0; $i<@genes; $i++){
        my @primers = @{${$genePrimerPairsRef->{$genes[$i]}}[$indexes[$i]]};
        push @strings, $primers[0], $primers[1];
    }
    
    my $badness = badness_total(\@strings, $badnessMinOverlap, \%pairBadnessCache);
    
    return 1/$badness;
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
