#!/usr/bin/env -S perl -w

use MCE::Map;
use List::Util qw(shuffle);

use File::Spec;
use File::Path qw(make_path remove_tree);
use File::Basename;

use IPC::Open2;

my $baseDir = dirname(__FILE__);
if($baseDir!~/^\//){
    $baseDir = File::Spec->rel2abs( $baseDir ) ;
}

my $debug = 0;

$| = 1;

# group division
my $nDiv = 2;
my $ncpu = 2;
my $numReport    = 5;

# GA1 parameters
my $randomSeed   = 0;

my $ga1_iteration = 100;
my $ga1_population = 200;
my $ga1_preserve   = 1;
my $ga1_crossover  = 0.9;
my $ga1_mutation   = 0.05;

# GA2: the underlying SADDLE
my $badnessMinOverlap = 4; # as suggested in the paper

# GA2 parameters
my $ga2_iteration = 150;
my $ga2_population = 1000;
my $ga2_crossover  = 0.60;
my $ga2_mutation   = 0.20;

my $usage = "Usage: saddleGA_pair_div.pl [options] <genePairedPrimerFile> <outPrefix>\n";
$usage   .= "    Group division:\n";
$usage   .= "        -nDiv       : divide targets into this number of groups (default: $nDiv)\n";
$usage   .= "        -cpu           : GA parameter mutation (default: $ga1_mutation)\n";
$usage   .= "        -report     : number of top answers to report (default: $numReport)\n";
$usage   .= "        -maxIt      : maximum iteration (default: $ga1_iteration)\n";
$usage   .= "        -GApopulation1 : GA parameter population (default: $ga1_population)\n";
$usage   .= "        -GApreserve1   : GA parameter preserve (default: $ga1_preserve)\n";
$usage   .= "        -GAcrossover1  : GA parameter crossover (default: $ga1_crossover)\n";
$usage   .= "        -GAmutation1   : GA parameter mutation (default: $ga1_mutation)\n";
$usage   .= "    SADDLE GA options:\n";
$usage   .= "        -minOverlap : minimum reversecomplement to compute badness (default: $badnessMinOverlap)\n";
$usage   .= "        -testIt     : test iteration (default: $ga2_iteration)\n";
$usage   .= "        -GApopulation2 : GA parameter population (default: $ga2_population)\n";
$usage   .= "        -GAcrossover2  : GA parameter crossover (default: $ga2_crossover)\n";
$usage   .= "        -GAmutation2   : GA parameter mutation (default: $ga2_mutation)\n";
$usage   .= "    Others:\n";
$usage   .= "        -rand       : random seed (default: $randomSeed)\n";
$usage   .= "        -cpu        : number of CPUs for computation (default: $ncpu)\n";

# optional parameter
my @arg_idx=(0..@ARGV-1);
for my $i (0..@ARGV-1) {
	if ($ARGV[$i] eq '-rand') {
		$randomSeed = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}elsif ($ARGV[$i] eq '-cpu') {
		$ncpu = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
# Group division parameters
	}elsif ($ARGV[$i] eq '-nDiv') {
		$nDiv = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}elsif ($ARGV[$i] eq '-report') {
		$numReport = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}elsif ($ARGV[$i] eq '-maxIt') {
		$ga1_iteration = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}elsif ($ARGV[$i] eq '-GApopulation1') {
		$ga1_population = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}elsif ($ARGV[$i] eq '-GApreserve1') {
		$ga1_preserve = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}elsif ($ARGV[$i] eq '-GAcrossover1') {
		$ga1_crossover = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}elsif ($ARGV[$i] eq '-GAmutation1') {
		$ga1_mutation = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
# SADDLE GA parameters
	}elsif ($ARGV[$i] eq '-minOverlap') {
		$badnessMinOverlap = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}elsif ($ARGV[$i] eq '-testIt') {
		$ga2_iteration = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}elsif ($ARGV[$i] eq '-GApopulation2') {
		$ga2_population = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}elsif ($ARGV[$i] eq '-GAcrossover2') {
		$ga2_crossover = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}elsif ($ARGV[$i] eq '-GAmutation2') {
		$ga2_mutation = $ARGV[$i+1];
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

#### Group division GA iteration

# MCE::Map init
MCE::Map->init({
    max_workers => $ncpu,
    chunk_size  => 1
});

# fitness function reference
my $fitnessFuncRef = \&fitness;

# initialization, let populationRef fitnessRef fitnessAnsRef represent one generation
# populationRef: array ref to individuals, where a individual is an ref to an array of permutation
# fitnessRef   : array ref to an array of computed fitness scores of this generation
# fitnessAnsRef{individual} = [ [target, left, right], ... ]
#        (index in population)   ^^^^^^^^^^^^^^^^^^^best answer of target (in some chunk)
my @initArray = sort keys %genePrimerPairs;
my ($populationRef, $fitnessScoreRef, $fitnessAnsRef) = initialize_population($ga1_population, 
                                                                              \@initArray, 
                                                                              $fitnessFuncRef);

my $iteration = 0;
my ($max, $min) = getBestWorstScore($fitnessScoreRef);
$max = 1/$max;
$min = 1/$min;

open(scoreOUT,">$outPrefix.score");
scoreOUT->autoflush(1);

# GA parameters
print scoreOUT "GA related parameters:\n";
print scoreOUT "  -GApopulation : $ga1_population\n";
print scoreOUT "  -GAcrossover  : $ga1_crossover\n";
print scoreOUT "  -GAmutation   : $ga1_mutation\n";

# iterations
print scoreOUT "\n";
print scoreOUT "Iteration best and worst\n";
print scoreOUT "$iteration $max $min\n";

while($iteration<$ga1_iteration){
    ($populationRef, $fitnessScoreRef, $fitnessAnsRef) = evolve($populationRef,
                                                                $fitnessScoreRef,
                                                                $fitnessAnsRef,
                                                                1,
                                                                $ga1_preserve,
                                                                $ga1_crossover,
                                                                $ga1_mutation);
    $iteration++;

    ($max, $min) = getBestWorstScore($fitnessScoreRef);
    $max = 1/$max;
    $min = 1/$min;
    
    print scoreOUT "$iteration $max $min\n";
}
close scoreOUT;

# GA output
# my @bests = $ga->getFittest($numReport, 1); # 1 for uniq answers
# 
# for(my $i=0; $i<@bests; $i++){
#     my @indexes = $ga->as_array($bests[$i]);
#     my @genes   = sort keys %{$genePrimerPairsRef};
#
#     open(iterOUT,">$outPrefix.best$i");
#     for(my $i=0; $i<@genes; $i++){
#         my @primers = @{${$genePrimerPairsRef->{$genes[$i]}}[$indexes[$i]]};
#         print iterOUT "$genes[$i]\t$primers[0]\t$primers[1]\n";
#     }
#     close iterOUT;
# }

#### Group division GA iteration, END

sub fitness {
    my $inArrRef = shift;

    my $chromosome = $inArrRef->[0]; # this should be an array ref
    my $answerRef  = $inArrRef->[1]; # this should be an array ref of
                                     # refs of target-left-right triples

    print "answer size [1]: ".@$answerRef."\n" if $debug;

    # build target answer hash
    my %targetAnswerHash = ();
    for(my $i=0;$i<@$answerRef;$i++){
        my @triple = @{$$answerRef[$i]};
        push @{$targetAnswerHash{$triple[0]}}, $triple[1], $triple[2];
    }
    
    my @genes = @{$chromosome};
    my @geneChunks = split_array(\@genes, $nDiv);
    
    my $badness = 0;
    my @ansArr = ();
    for(my $i=0; $i<@geneChunks; $i++){
        my $ansRef = call_saddleGA($geneChunks[$i],
                                   \%targetAnswerHash,
                                   $badnessMinOverlap,
                                   $ga2_iteration,
                                   $ga2_population,
                                   $ga2_crossover,
                                   $ga2_mutation);
        # first item, best badness
        $ansArr[0] += (shift @$ansRef);
        # rest triples, target left/right primers
        while(@$ansRef){
            my $target = shift @$ansRef;
            my $leftP  = shift @$ansRef;
            my $rightP = shift @$ansRef;
            push @ansArr, $target, $leftP, $rightP;
        }
    }
    print "answer size [2]: ".(@ansArr-1)."\n" if $debug;
    # turn total badness into fitness
    $ansArr[0] = 1/$ansArr[0];
    
    return \@ansArr;
}

#### sub routines

sub call_saddleGA {
    my ($geneChunkRef,$targetAnsRef,$minOverlap,$maxIt,$GApopulation,$GAcrossover,$GAmutation) = @_;

    my $randomSeed=int(1000000*rand());
    
    my $cmd  = "$baseDir/saddleGA_pair.pl";
       $cmd .= " -rand $randomSeed";
       $cmd .= " -minOverlap $minOverlap";
       $cmd .= " -maxIt $maxIt";
       $cmd .= " -GApopulation $GApopulation";
       $cmd .= " -GAcrossover $GAcrossover";
       $cmd .= " -GAmutation $GAmutation";
       $cmd .= " -report 1";
       $cmd .= " -as_child 1";
       $cmd .= " dummy dummy";
    
    my $pid = open2(my $child_out, my $child_in, $cmd);
    
    # subset genePrimerPairs hash
    my %subsetHash = map { $_ => $genePrimerPairs{$_} } @{$geneChunkRef};
    # output subset hash to child
    for my $target (keys %subsetHash){
        for my $arrRef (@{$subsetHash{$target}}){
            print $child_in "$target @{$arrRef}\n";
        }
    }
    # subset targetAnsRef
    if((scalar keys %{$targetAnsRef}) > 0){
        my %subsetAnswer = map { $_ => $targetAnsRef->{$_} } @{$geneChunkRef};
        # minor check, SHOULD not happen
        if((scalar keys %subsetHash) != (scalar keys %subsetAnswer)){
            die "chunk answer not enough\n";
        }
        for my $target (keys %subsetAnswer){
            print $child_in "!INJECT!$target\t$subsetAnswer{$target}[0]\t$subsetAnswer{$target}[1]\n";
        }
    }
    close $child_in;

    my @childSays;    
    while (<$child_out>) {         # Read child's STDOUT
        chomp $_;
        s/^\s+|\s+$//g;

        push @childSays, $_ if length($_)>0;
    }
    close $child_out;
    
    my @ans;
    my @t=split(/\s+/, $childSays[-(@{$geneChunkRef})-1]);
    push @ans, $t[1];
    for(my $i=@childSays-(@{$geneChunkRef}); $i<@childSays; $i++){
        my @t=split(/\s+/, $childSays[$i]);
        push @ans, @t;
    }
    
    return \@ans;
}

sub split_array {
    my ($array_ref, $num_chunks) = @_;
    my $size = @{$array_ref};
    my $chunk_size = int($size / $num_chunks);
    my $remainder = $size % $num_chunks;
    
    my @chunks;
    my $start = 0;
    
    for my $i (0 .. $num_chunks - 1) {
        my $end = $start + $chunk_size - 1;
        $end++ if $i < $remainder;  # Distribute extra items evenly
        
        push @chunks, [@$array_ref[$start .. $end]];
        
        $start = $end + 1;
    }
    
    return @chunks;
}

#### GA implementation for permutation encoding

sub getBestWorstScore{
    my $fitnessScoreRef = shift;
    
    my %fitnessHash = ();
    for(my $i=0;$i<@$fitnessScoreRef;$i++){
        $fitnessHash{$i} = $$fitnessScoreRef[$i];
    }
    my @sortedIdx = sort {$fitnessHash{$b}<=>$fitnessHash{$a}} keys %fitnessHash;
    
    return ($fitnessHash{$sortedIdx[0]}, $fitnessHash{$sortedIdx[-1]});
}

sub get_best{
    my $populationRef   = shift;
    my $fitnessScoreRef = shift;
    my $bestN           = shift;

    my %fitnessHash = ();
    for(my $i=0;$i<@$fitnessScoreRef;$i++){
        $fitnessHash{$i} = $$fitnessScoreRef[$i];
    }
    my @sortedIdx = sort {$fitnessHash{$b}<=>$fitnessHash{$a}} keys %fitnessHash;

    my @ansList;
    for(my $i=0;$i<$bestN;$i++){
        push @ansList, $$populationRef[$sortedIdx[$i]];
    }

    return \@ansList;
}

sub evolve {
    if($debug){
        print "@_\n";
    }
    
    my $populationRef   = shift;
    my $fitnessScoreRef = shift;
    my $fitnessAnsRef   = shift;
    my $GENERATIONS     = shift;
    my $preserve        = shift;
    my $crossoverRate   = shift;
    my $mutationRate    = shift;

    # checking
    if($preserve >= @$populationRef){
        die "preserve number >= population size\n";
    }
    
    # every generation
    for my $gen (1..$GENERATIONS) {
    
        my @new_population;
        my %new_fitnessAnsHash = ();

        # fitness hash and indexes from best to worst    
        my %fitnessHash = ();
        for(my $i=0;$i<@$fitnessScoreRef;$i++){
            $fitnessHash{$i} = $$fitnessScoreRef[$i];
        }
        my @sortedIdx = sort {$fitnessHash{$b}<=>$fitnessHash{$a}} keys %fitnessHash;
        
        if($debug){
            print "parent BEST: ".(1/$fitnessHash{$sortedIdx[0]})."\n";
        }

        # preserve best
        for(my $i=0; $i<$preserve || $i<int((1-$crossoverRate)*@$populationRef); $i++){
            push @new_population, $$populationRef[$sortedIdx[$i]];
            if($i>=$preserve && rand()<$mutationRate){
                # non-preserved ones might be mutated
                swap_mutation($new_population[-1]);
            }
            # answer part
            $new_fitnessAnsHash{$i} = $fitnessAnsRef->{$sortedIdx[$i]};
        }

        # establish a cumulative array of fitness
        my @cFitness = (0);
        for($i=0;$i<@sortedIdx;$i++){
            push @cFitness, ($fitnessHash{$sortedIdx[$i]} + $cFitness[-1]);
        }
        my $totalFitness = $cFitness[-1];

        # fill rest slots by Roulette selection + crossover rate + mutation rate
        while(@new_population < @$populationRef){
            # pick parent 1
            my $randSel = rand($totalFitness);
            my $parent1Idx = binary_search(\@cFitness, $randSel);
            # pick parent 2
            my $parent2Idx = $parent1Idx;
            while($parent2Idx!=$parent1Idx){
                $randSel = rand($totalFitness);
                $parent2Idx = binary_search(\@cFitness, $randSel);
            }
            # crossover
            my ($new_chr, $new_ans) = order_crossover($$populationRef[$sortedIdx[$parent1Idx]],
                                                      $$populationRef[$sortedIdx[$parent2Idx]],
                                                      $fitnessAnsRef->{$sortedIdx[$parent1Idx]},
                                                      $fitnessAnsRef->{$sortedIdx[$parent2Idx]});
            push @new_population, $new_chr;
            $new_fitnessAnsHash{@new_population-1} = $new_ans;
            # mutation
            if(rand() < $mutationRate){
                swap_mutation($new_population[-1]);
            }

            # answer part, give nothing
        }
        
        $populationRef = \@new_population;
        $fitnessAnsRef = \%new_fitnessAnsHash;
        ($fitnessScoreRef, $fitnessAnsRef) = population_fitness($populationRef, 
                                                                $fitnessFuncRef,
                                                                $fitnessAnsRef);
    }
    
    return ($populationRef, $fitnessScoreRef, $fitnessAnsRef);
}

sub binary_search {
    my ($arr, $target) = @_;
    my ($low, $high) = (0, $#$arr);

    while ($low <= $high) {
        my $mid = int(($low + $high) / 2);

        if ($arr->[$mid]<=$target && $target<$arr->[$mid+1]) {
            return $mid;
        } elsif ($arr->[$mid+1] <= $target) {
            $low = $mid + 1;  # Search right half
        } else {
            $high = $mid;  # Search left half
        }
    }
    return -1;  # Not found
}

sub population_fitness {
    my $populationRef = shift;
    my $fitnessFuncRef = shift;
    my $fitnessAnsRef = shift;

    my @tmpKeys = keys %{$fitnessAnsRef};
    print "answer ref [4]: $fitnessAnsRef ".(keys %{$fitnessAnsRef}).": @tmpKeys\n" if $debug;

    # build input array of arrays based on populationRef and fitnessAnsRef
    my @inputArray = ();
    for(my $i=0;$i<@$populationRef;$i++){
        # first element: population
        $inputArray[$i] = [ $populationRef->[$i] ];
        # second element: previous best answer, if any
        if(exists $fitnessAnsRef->{$i}){
            push @{$inputArray[$i]}, $fitnessAnsRef->{$i};
        }else{
            my @empty = ();
            push @{$inputArray[$i]}, \@empty;
        }
        print "answer size [6]: $i ".(@{$inputArray[$i]->[1]})."\n" if $debug;
    }

    # MCE::Map
    my @fitnessArrays = mce_map {
        $fitnessFuncRef->($_);
    } \@inputArray;

    # build fitnessdScores and fitnessAns as definition
    my @fitnessScores;
    my %fitnessAnsHash = ();
    for(my $i=0;$i<@fitnessArrays;$i++){
        # first element in an array
        push @fitnessScores, (shift @{$fitnessArrays[$i]});
        # rest elements, ith answer
        my @answerArray = ();
        for(my $j=0;$j<@{$fitnessArrays[$i]};$j+=3){
            my @triple = @{$fitnessArrays[$i]}[$j..($j+2)];
            push @answerArray, \@triple;
        }
        $fitnessAnsHash{$i} = \@answerArray;
        print "answer size [3]: ".(@{$fitnessArrays[$i]})." ".(@answerArray)." ".(@{$fitnessAnsHash{$i}})."\n" if $debug;
    }
    
    return (\@fitnessScores, \%fitnessAnsHash);
}

sub initialize_population {
    my ($size, $elementsRef, $fitnessFuncRef) = @_;
    
    # generate randomized population
    my @population;
    for (1 .. $size) {
        push @population, [shuffle(@{$elementsRef})];  # Random permutation
    }
    
    # compute fitness for all
    my %emptyHash = ();
    my ($fitnessScoreRef, $fitnessAnsRef) = population_fitness(\@population, $fitnessFuncRef , \%emptyHash);
    
    return (\@population, $fitnessScoreRef, $fitnessAnsRef);
}

# OX crossover
sub order_crossover {
    my ($parent1, $parent2, $ans1, $ans2) = @_;
    my @child = @$parent1;
    my @new_ans = ();

    # build hash for the two answers
    my %ansHash1 = ();
    my %ansHash2 = ();
    
    for(my $i=0;$i<@$ans1;$i++){
        my @triple = @{$ans1->[$i]};
        $ansHash1{$triple[0]} = \@triple;
    }
    for(my $i=0;$i<@$ans2;$i++){
        my @triple = @{$ans2->[$i]};
        $ansHash2{$triple[0]} = \@triple;
    }

    # Select a random subset from parent1
    my $start = int(rand(@$parent1));
    my $end   = $start + int(rand(@$parent1 - $start));
    
    my @selected = @child[$start..$end];

    my %selectedHash = ();
    for my $g (@selected){
        $selectedHash{$g} = 1;
    }
    for my $g (@$parent1){
        if(exists $selectedHash{$g}){
            push @new_ans, $ansHash1{$g};
        }else{
            push @new_ans, $ansHash2{$g};
        }
    }
    
    # Fill remaining positions using parent2's order
    my %used = map { $_ => 1 } @selected;
    my @remaining = grep { !$used{$_} } @$parent2;

    @child = (@remaining[0..$start-1], @selected, @remaining[$start..$#remaining]);

    return (\@child, \@new_ans);
}

# in-place swap, no need to change its fitness answer 
sub swap_mutation {
    my $chromosome = shift;
    my $i = int(rand(@$chromosome));
    my $j = int(rand(@$chromosome));

    @$chromosome[$i, $j] = @$chromosome[$j, $i];  # Swap two elements
}
