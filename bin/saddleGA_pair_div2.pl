#!/usr/bin/env -S perl -w

use MCE::Map;
use List::Util qw(shuffle sum);

use File::Spec;
use File::Path qw(make_path remove_tree);
use File::Basename;

use JSON;

use IPC::Open2;

my $baseDir = dirname(__FILE__);
if($baseDir!~/^\//){
    $baseDir = File::Spec->rel2abs( $baseDir ) ;
}

my $debug = 1;
my $debug_sig = 0;

$| = 1;

# group division
my $nDiv = 2; # fixed
my $ncpu = 2;
my $numReport    = 5;

# GA1 parameters
my $randomSeed   = 0;

my $ga1_iteration = 100;
my $ga1_population = 200;
my $ga1_preserve   = 1;
my $ga1_crossover  = 0.9;
my $ga1_mutation   = 0.05;

my $startStateFile = "";

# GA2: the underlying SADDLE
my $badnessMinOverlap = 4; # as suggested in the paper

# GA2 parameters
my $ga2_iteration = 150;
my $ga2_population = 1000;
my $ga2_crossover  = 0.60;
my $ga2_mutation   = 0.20;

my $usage = "Usage: saddleGA_pair_div2.pl [options] <genePairedPrimerFile> <outPrefix>\n";
$usage   .= "    Group division:\n";
$usage   .= "        -report     : number of top answers to report (default: $numReport)\n";
$usage   .= "        -maxIt      : maximum iteration (default: $ga1_iteration)\n";
$usage   .= "        -GApopulation1 : GA parameter population (default: $ga1_population)\n";
$usage   .= "        -GApreserve1   : GA parameter preserve (default: $ga1_preserve)\n";
$usage   .= "        -GAcrossover1  : GA parameter crossover (default: $ga1_crossover)\n";
$usage   .= "        -GAmutation1   : GA parameter mutation (default: $ga1_mutation)\n";
$usage   .= "        -startState    : load saved GA state (default: \"$startStateFile\")\n";
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
	}elsif ($ARGV[$i] eq '-startState') {
		$startStateFile = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
# Group division parameters
	}elsif ($ARGV[$i] eq '-report') {
		$numReport = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}elsif ($ARGV[$i] eq '-maxIt') {
		$ga1_iteration = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}elsif ($ARGV[$i] eq '-GApreserve1') {
		$ga1_preserve = $ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}elsif ($ARGV[$i] eq '-GApopulation1') {
		$ga1_population = $ARGV[$i+1];
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

#### target combination GA iteration

# GA parameters
open(scoreOUT,">$outPrefix.score");
scoreOUT->autoflush(1);

print scoreOUT "GA related parameters:\n";
print scoreOUT "Level 1:\n";
print scoreOUT "    -maxIt         : $ga1_iteration\n";
print scoreOUT "    -GApopulation1 : $ga1_population\n";
print scoreOUT "    -GApreserve1   : $ga1_preserve\n";
print scoreOUT "    -GAcrossover1  : $ga1_crossover\n";
print scoreOUT "    -GAmutation1   : $ga1_mutation\n";
print scoreOUT "Level 2:\n";
print scoreOUT "    -testIt        : $ga2_iteration\n";
print scoreOUT "    -GApopulation2 : $ga2_population\n";
print scoreOUT "    -GAcrossover2  : $ga2_crossover\n";
print scoreOUT "    -GAmutation2   : $ga2_mutation\n";
print scoreOUT "\n";

# fitness function reference
my $fitnessFuncRef = \&fitness;

# initialization, let populationRef fitnessRef fitnessAnsRef represent one generation
# populationRef: array ref to individuals, where a individual is an ref to an array of combination
# fitnessRef   : array ref to an array of computed fitness scores of this generation
# fitnessAnsRef{individual} = [ [target, left, right], ... ]
#        (index in population)   ^^^^^^^^^^^^^^^^^^^best answer of target (in some chunk)
my @targetArray = sort keys %genePrimerPairs; # for quick reference
my @initArray = sort keys %genePrimerPairs;
# make first half of the initArray "selected"
for(my $i=0;$i<@initArray;$i++){
    if($i<int(@initArray/2)){
        $initArray[$i]=0;
    }else{
        $initArray[$i]=1;
    }
}

my ($populationRef, $fitnessScoreRef, $fitnessAnsRef);

if(length($startStateFile)){
    if($debug){
        print STDERR "state read start\n";
    }
    # read from a state JSON file
    open(stateIN,"<$startStateFile");
    my $json_text = do { local $/; <stateIN> };
    close stateIN;
    
    my $data = decode_json($json_text);
    my @population = ();
    my @fitnessScores = ();
    my %fitnessAnswer = ();
    # population
    for my $i (sort {$a<=>$b} keys %{$data->{"population"}}){
        my @arr = @{$data->{"population"}->{$i}};
        push @population, \@arr;
    }
    # fitness
    for my $i (sort {$a<=>$b} keys %{$data->{"fitness"}}){
        $fitnessScores[$i] = $data->{"fitness"}->{$i};
    }
    # answer
    for my $i (sort {$a<=>$b} keys %{$data->{"answer"}}){
        my $arrRef = $data->{"answer"}->{$i};
        $fitnessAnswer{$i} = $arrRef;
    }

    if($debug){
        print STDERR "state read end\n";
    }
    
    ($populationRef, $fitnessScoreRef, $fitnessAnsRef) = (\@population, \@fitnessScores, \%fitnessAnswer);
}else{
    if($debug){
        print STDERR "population initialization start\n";
    }
    # regular initialization
    ($populationRef, $fitnessScoreRef, $fitnessAnsRef) = initialize_population($ga1_population, 
                                                                               \@initArray, 
                                                                               $fitnessFuncRef);
    if($debug){
        print STDERR "population initialization end\n";
    }
}

if($debug_sig){
    print STDERR "population out:\n";
    for(my $i=0; $i<@$populationRef; $i++){
        print STDERR "".sprintf("%03d",$i)." ".chromosome_signature($$populationRef[$i])."\n";
    }
}

my $iteration = 0;
my ($max, $min) = getBestWorstScore($fitnessScoreRef);
$max = 1/$max;
$min = 1/$min;


# iterations
print scoreOUT "Iteration best and worst\n";
print scoreOUT "$iteration $max $min\n";

while($iteration<$ga1_iteration){
    if($debug){
        print STDERR "iteartion ".($iteration+1)."\n";
    }
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

# GA state output, should include
# (i) current population
# (ii) best primer selection inside each individual
# so that this can be loaded and continued
my %jsonHash = ();
for(my $i=0; $i<@$populationRef; $i++){
    $jsonHash{"population"}{$i} = $$populationRef[$i];
    $jsonHash{"fitness"}{$i}    = $$fitnessScoreRef[$i];
    $jsonHash{"answer"}{$i}     = $$fitnessAnsRef{$i};
}
open(stateOUT, ">$outPrefix.state");
print stateOUT to_json(\%jsonHash, { pretty => 1 })."\n";
close stateOUT;

# GA report output, should include
# (i) best current population
# (ii) best primer selection inside each individual
report_bestN($populationRef, $fitnessScoreRef, $fitnessAnsRef, $outPrefix, $numReport);

#### target combination GA iteration, END

sub fitness {
    my $inArrRef = shift;

    my $chromosome = $inArrRef->[0]; # this should be an array ref
    my $answerRef  = $inArrRef->[1]; # this should be an array ref of
                                     # refs of target-left-right triples

    # build target answer hash
    my %targetAnswerHash = ();
    for(my $i=0;$i<@$answerRef;$i++){
        my @triple = @{$$answerRef[$i]};
        push @{$targetAnswerHash{$triple[0]}}, $triple[1], $triple[2];
    }
    
    my @genes = @{$chromosome}; # 0's and 1's in the array
    my @array0 = ();
    my @array1 = ();
    for(my $i=0;$i<@genes;$i++){
        if($genes[$i]){ # array1
            push @array1, $targetArray[$i];
        }else{          # array0
            push @array0, $targetArray[$i];
        }
    }
    my @geneChunks = (\@array0, \@array1);
    
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

#### GA implementation for combination encoding

sub getBestWorstScore{
    my $fitnessScoreRef = shift;
    
    my %fitnessHash = ();
    for(my $i=0;$i<@$fitnessScoreRef;$i++){
        $fitnessHash{$i} = $$fitnessScoreRef[$i];
    }
    my @sortedIdx = sort {$fitnessHash{$b}<=>$fitnessHash{$a}} keys %fitnessHash;
    
    return ($fitnessHash{$sortedIdx[0]}, $fitnessHash{$sortedIdx[-1]});
}

sub report_bestN{
    my $populationRef   = shift;
    my $fitnessScoreRef = shift;
    my $fitnessAnsRef   = shift;
    my $outPrefix       = shift;
    my $bestN           = shift;

    my %fitnessHash = ();
    for(my $i=0;$i<@$fitnessScoreRef;$i++){
        $fitnessHash{$i} = $$fitnessScoreRef[$i];
    }
    my @sortedIdx = sort {$fitnessHash{$b}<=>$fitnessHash{$a}} keys %fitnessHash;

    # for each of the best individual
    for(my $i=0;$i<$bestN;$i++){
        open(iterOUT,">$outPrefix.best$i");
        
        # build target-answer map
        my %targetAnswer = ();
        for my $tripleRef (@{$fitnessAnsRef->{$sortedIdx[$i]}}){
            $targetAnswer{$$tripleRef[0]} = $tripleRef;
        }
        # build group-target map
        my %groupTarget = ();
        for(my $j=0; $j<@{$$populationRef[$sortedIdx[$i]]}; $j++){
            $groupTarget{${$$populationRef[$sortedIdx[$i]]}[$j]}{$targetArray[$j]}=1;
        }
        # output
        for my $group (sort {$a<=>$b} keys %groupTarget){
            for my $target (sort keys %{$groupTarget{$group}}){
                print iterOUT "$group\t@{$targetAnswer{$target}}\n";
            }
        }
        
        close iterOUT;
    }
}

sub evolve {
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
    
        if($debug_sig){
            print STDERR "population in:\n";
            for(my $i=0; $i<@$populationRef; $i++){
                print STDERR "".sprintf("%03d",$i)." ".chromosome_signature($$populationRef[$i])."\n";
            }
        }
        my @new_population = ();
        my %new_fitnessAnsHash = ();
        
        my %populationPool = ();

        # fitness hash and indexes from best to worst    
        my %fitnessHash = ();
        for(my $i=0;$i<@$fitnessScoreRef;$i++){
            $fitnessHash{$i} = $$fitnessScoreRef[$i];
        }
        my @sortedIdx = sort {$fitnessHash{$b}<=>$fitnessHash{$a}} keys %fitnessHash;
        
        # best
        if($debug){
            print STDERR "iteartion ".($iteration+1).": best preserve\n";
        }
        while(@new_population<$preserve || @new_population<int((1-$crossoverRate)*@$populationRef)){
            my @new_chr = @{$$populationRef[$sortedIdx[@new_population]]};
            
            # no mutation if preserved
            if(@new_population<$preserve){
                # populationPool update
                $populationPool{chromosome_signature(\@new_chr)}=1;
                if($debug_sig){
                    print STDERR "sig pool preserv add: $sortedIdx[@new_population] : ".chromosome_signature(\@new_chr)."\n";
                }
                # add answer part
                $new_fitnessAnsHash{@new_population} = $fitnessAnsRef->{$sortedIdx[@new_population]};
                # add to new population
                push @new_population, \@new_chr;
            }else{
                # non-preserved, see if need mutation
                if(rand()<$mutationRate){
                    swap_mutation(\@new_chr);
                }
                # signature check to avoid duplicates
                my $signature = chromosome_signature(\@new_chr);
                if(not exists $populationPool{$signature}){
                    # populationPool update
                    $populationPool{$signature}=1;
                    if($debug_sig){
                        print STDERR "sig pool best add: $sortedIdx[@new_population] : $signature\n";
                    }
                    # add answer part
                    $new_fitnessAnsHash{@new_population} = $fitnessAnsRef->{$sortedIdx[@new_population]};
                    # add to new population
                    push @new_population, \@new_chr;
                }
            }
        }

        # establish a cumulative array of fitness
        my @cFitness = (0);
        for(my $i=0;$i<@sortedIdx;$i++){
            push @cFitness, ($fitnessHash{$sortedIdx[$i]} + $cFitness[-1]);
        }
        my $totalFitness = $cFitness[-1];

        # fill rest slots by Roulette selection + crossover rate + mutation rate
        if($debug){
            print STDERR "iteartion ".($iteration+1).": crossover\n";
        }
        while(@new_population < @$populationRef){
            # pick parent 1
            my $randSel = rand($totalFitness);
            my $parent1Idx = binary_search(\@cFitness, $randSel);
            # pick parent 2
            my $parent2Idx = $parent1Idx;
            while($parent2Idx==$parent1Idx){
                $randSel = rand($totalFitness);
                $parent2Idx = binary_search(\@cFitness, $randSel);
            }
            
            my $sigFail = 0; # this is for debug
            
            # crossover
            my ($new_chr, $new_ans) = shrinking_crossover($$populationRef[$sortedIdx[$parent1Idx]],
                                                          $$populationRef[$sortedIdx[$parent2Idx]],
                                                          $fitnessAnsRef->{$sortedIdx[$parent1Idx]},
                                                          $fitnessAnsRef->{$sortedIdx[$parent2Idx]});
            # mutation
            if(rand() < $mutationRate){
                swap_mutation($new_chr);
            }
            # signature check to avoid duplicates
            my $signature = chromosome_signature($new_chr);
            if($debug_sig){
                print STDERR "iteartion ".($iteration+1).": signature check\n";
            }
            if(not exists $populationPool{$signature}){
                if($debug_sig){
                    print STDERR "iteartion ".($iteration+1).": signature pass\n";
                }
                # populationPool update
                $populationPool{$signature}=1;
                if($debug_sig){
                    print STDERR "sig pool add: $signature\n";
                }
                # add answer part
                $new_fitnessAnsHash{@new_population} = $new_ans;
                # add to new population
                push @new_population, $new_chr;
            }elsif($debug_sig){
                $sigFail++;
                
                if($sigFail>=10){
                    print STDERR "cFitness: @cFitness\n";
                    print STDERR "last random sel val: $randSel\n";
                    print STDERR "last random sel idx: $parent1Idx $parent2Idx\n";
                    print STDERR "prnt1: ".chromosome_signature($$populationRef[$sortedIdx[$parent1Idx]])."\n";
                    print STDERR "prnt2: ".chromosome_signature($$populationRef[$sortedIdx[$parent2Idx]])."\n";
                    print STDERR "child: $signature\n";
                    print STDERR "current pool\n";
                    for(my $i=0; $i<@new_population; $i++){
                        print STDERR "".sprintf("%03d",$i)." ".chromosome_signature($new_population[$i]);
                        print STDERR "<==" if chromosome_signature($new_population[$i]) eq $signature;
                        print STDERR "\n";
                    }
                    
                    exit 1;
                }
            }
        }
        
        if($debug){
            print STDERR "iteartion ".($iteration+1).": fitness\n";
        }
        $populationRef = \@new_population;
        $fitnessAnsRef = \%new_fitnessAnsHash;
        ($fitnessScoreRef, $fitnessAnsRef) = population_fitness($populationRef, 
                                                                $fitnessFuncRef,
                                                                $fitnessAnsRef);
        if($debug){
            print STDERR "iteartion ".($iteration+1).": done\n";
        }
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
    }

    # MCE::Map init
    MCE::Map->init({
        max_workers => $ncpu,
        chunk_size  => 1
    });

    # MCE::Map
    my @fitnessArrays = mce_map {
        $fitnessFuncRef->($_);
    } \@inputArray;
    
    # MCE::Map manual shutdown
    MCE::Map->finish;

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
    }
    
    return (\@fitnessScores, \%fitnessAnsHash);
}

sub initialize_population {
    my ($size, $elementsRef, $fitnessFuncRef) = @_;
    
    my %populationPool = ();
    
    # generate randomized population
    my @population;
    while(@population < $size) {
        my @new_chr = shuffle(@{$elementsRef}); # Random permutation
        my $signature = chromosome_signature(\@new_chr);
        # signature check to avoid duplicates
        if(not exists $populationPool{$signature}){
            push @population, \@new_chr;
            $populationPool{$signature} = 1;
            if($debug_sig){
                print STDERR "sig pool add: $signature\n";
            }
        }
    }
    
    # compute fitness for all
    my %emptyHash = ();
    my ($fitnessScoreRef, $fitnessAnsRef) = population_fitness(\@population, $fitnessFuncRef , \%emptyHash);
    
    return (\@population, $fitnessScoreRef, $fitnessAnsRef);
}

sub chromosome_signature {
    my $chromosome = shift;
    return "@$chromosome";
}

# SX crossover
# REF: Chen, JS., Hou, JL.  (2006).  A Combination Genetic Algorithm with
# Applications on Portfolio Optimization.  In: Ali, M., Dapoigny, R.  (eds)
# Advances in Applied Artificial Intelligence.  IEA/AIE 2006.  Lecture Notes
# in Computer Science(), vol 4031.  Springer, Berlin, Heidelberg. 
# https://doi.org/10.1007/11779568_23
#
# procedure SX(C1, C2, s, t)
#     until (C1 and C2 have the same number of 1’s between s & t)
#         t = t-1
#     enduntil
#     for (i = s to t)
#         Temp = C1 [i]
#         C1 [i] = C2 [i]
#         C2 [i] = Temp
#     endfor
# endprocedure
sub shrinking_crossover {
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
    my $start = 0;
    my $end   = 0;
    my $sum1  = 0;
    my $sum2  = 0;
    # ensure exchange, avoid (i) no 1's (ii) all 1's
    while( $sum1==0 || $sum2==0 || 
           $sum1==($end-$start+1) || $sum2==($end-$start+1)){
        $start = int(rand(@$parent1));
        $end   = $start + int(rand(@$parent1 - $start));
        $sum1  = sum(@{$parent1}[$start..$end]);
        $sum2  = sum(@{$parent2}[$start..$end]);
        # the shrinking loop, this ensures sum1==sum2
        while($sum1 != $sum2){
            $end--;
            # force redo if empty range
            if($end<$start){
                $sum1=0;
                last;
            }
            # use differential
            $sum1 -= $$parent1[$end+1];
            $sum2 -= $$parent2[$end+1];
        }
    }
    
    # NOTE: in this practice, one child was generated based on crossover of two parents
    # will take parent1 as the template, with parent2 in the start-end region
    for(my $i=$start; $i<=$end; $i++){
        $child[$i] = $$parent2[$i];
    }

    for(my $i=0;$i<@targetArray;$i++){
        if($start<=$i && $i<=$end){
            push @new_ans, $ansHash2{$targetArray[$i]};
        }else{
            push @new_ans, $ansHash1{$targetArray[$i]};
        }
    }
    
    return (\@child, \@new_ans);
}

# in-place swap, no need to change its fitness answer 
sub swap_mutation {
    my $chromosome = shift;
    my $index0 = int(rand(@$chromosome));
    my $index1 = int(rand(@$chromosome));

    while($$chromosome[$index0]!=0){
        $index0 = int(rand(@$chromosome));
    }
    while($$chromosome[$index1]!=1){
        $index1 = int(rand(@$chromosome));
    }

    if($debug_sig){
        print STDERR "swap before: ".chromosome_signature($chromosome)."\n";
    }
    @$chromosome[$index0, $index1] = @$chromosome[$index1, $index0];  # Swap two elements
    if($debug_sig){
        print STDERR "swap after : ".chromosome_signature($chromosome)."\n";
    }
}
