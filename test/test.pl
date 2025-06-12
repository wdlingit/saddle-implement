#!/usr/bin/env -S perl -w

use MCE::Map;
use List::Util qw(shuffle);

# Define fitness function
sub fitness {
    my $chromosome = shift; # should be an array ref
    
    my $ans = 0;
    for(my $i=0; $i<@{$chromosome}; $i++){
        $ans += ((@{$chromosome}-$i)*$$chromosome[$i]);
    }

    return $ans;
}

# fitness function reference
my $fitnessFuncRef = \&fitness;

# template input array
my @items = 1..100;

MCE::Map->init({
    max_workers => 35,
    chunk_size  => 1
});

# initialization, let populationRef and fitnessRef represent one generation
my ($populationRef, $fitnessScoreRef) = initialize_population(1000, \@items, $fitnessFuncRef);

($populationRef, $fitnessScoreRef) = evolve($populationRef,
                                            $fitnessScoreRef,
                                            2000,
                                            10,
                                            0.8,
                                            0.1);

$bestsRef = get_best($populationRef,$fitnessScoreRef,1);
print "@{$bestsRef->[0]}\n";

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
            push @new_population, order_crossover($$populationRef[$sortedIdx[$parent1Idx]],
                                                  $$populationRef[$sortedIdx[$parent2Idx]]);
            # mutation
            if(rand() < $mutationRate){
                swap_mutation($new_population[-1]);
            }
        }
        
        $populationRef = \@new_population;
        $fitnessScoreRef = population_fitness($populationRef, $fitnessFuncRef);
    }
    
    return ($populationRef, $fitnessScoreRef);
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
    
    my @fitnessScores = mce_map {
        $fitnessFuncRef->($_);
    } $populationRef;
    
    return \@fitnessScores;
}

sub initialize_population {
    my ($size, $elementsRef, $fitnessFuncRef) = @_;
    
    # generate randomized population
    my @population;
    for (1 .. $size) {
        push @population, [shuffle(@{$elementsRef})];  # Random permutation
    }
    
    # compute fitness for all
    my $fitnessScoreRef = population_fitness(\@population, $fitnessFuncRef);
    
    return (\@population, $fitnessScoreRef);
}

# OX crossover
sub order_crossover {
    my ($parent1, $parent2) = @_;
    my @child = @$parent1;  

    # Select a random subset from parent1
    my $start = int(rand(@$parent1));
    my $end   = $start + int(rand(@$parent1 - $start));
    
    my @selected = @child[$start..$end];
    
    # Fill remaining positions using parent2's order
    my %used = map { $_ => 1 } @selected;
    my @remaining = grep { !$used{$_} } @$parent2;

    @child = (@remaining[0..$start-1], @selected, @remaining[$start..$#remaining]);

    return \@child;
}

# in-place swap
sub swap_mutation {
    my $chromosome = shift;
    my $i = int(rand(@$chromosome));
    my $j = int(rand(@$chromosome));

    @$chromosome[$i, $j] = @$chromosome[$j, $i];  # Swap two elements
}
