# saddle-implement

Implement of the iteration part of the SADDLE framework (PMID: 35410464) for multiplex PCR primer design. Two versions in this repository.

NOTE: In this implementation, _badness_ is computed directly following Eq1 and Eq2 from the SADDLE paper. The _hash_ described in the paper was not included. Insteadly, this implementation caches _badness_ of all ever-computed primer pairs. When computing _badness_ of one given primer pair, min-mer lookup table and dot plot approaches were applied. This should be faster than simple substring comparisons.

## saddle_pair.pl

Perl implementation of the iteration part as described in the paper. The only exception is that the C(g) function in Eq4 was replaced by a constant, which can be specified by `-tolC`. You may use `-maxTolIt` to specify the iteration to stop the tolerance mechanism.

```
saddle-implement/bin$ ./saddle_pair.pl
Usage: saddle_pair.pl [options] <genePairedPrimerFile> <outPrefix>
           -rand       : random seed (default: 0)
           -minOverlap : minimum reversecomplement to compute badness (default: 4)
           -maxIt      : maximum iteration (default: 60000)
           -maxTolIt   : maximum tolerance iteration (default: 20000)
           -tolC       : factor C for tolerance threshold (default: 100)
           -outIter    : iterations to output answers (default: none)
           -badness    : no iteration and compute badness of iteration outputs (default: no)
           -startWith  : start iterations with primer pairs assigned in the given file (default: "")
```

The input file `genePairedPrimerFile` should be a three column tab-delimited text file. In this input file, each line represents one primer pair by three tokens: (i) target ID, (ii) left primer, and (iii) right primer. The selection would be for primer pairs but not for individual primers. `outPrefix` specifies the output prefix of a number of output files:
1. `<outPrefix>.score`: badness of each iteration
2. `<outPrefix>.n`: primer pair combination of iteration `n`, specified by `-outIter`. This option can be applied multiple times. You may assign one such file to option `-badness` for confirming its badness.

Under WSL ubuntu of my desktop PC, this script (with default parameters) took less than 30 minutes for 96 targets.

## saddleGA_pair.pl

Genetic Algorithm version of the SADDLE framework. In this implementation, GA _fitness_ is defined as inverse of _badness_. The tolerance was not explicitly implemented becasue _mutation_ of GA can be considered as a kind of tolerance. This script requires the `AI::Genetic::Pro` module. You may use `cpan` to install it. In case you are using EasyBuild, we have an EB file. Please contact me directly.

```
saddle-implement/bin$ ./saddleGA_pair.pl
Usage: saddleGA_pair.pl [options] <genePairedPrimerFile> <outPrefix>
           -rand       : random seed (default: 0)
           -minOverlap : minimum reversecomplement to compute badness (default: 4)
           -maxIt      : maximum iteration (default: 600)
           -report     : number of top answers to report (default: 5)
           -GApopulation : GA parameter population (default: 200)
           -GAcrossover  : GA parameter crossover (default: 0.9)
           -GAmutation   : GA parameter mutation (default: 0.05)
```

Note that it seems that the option `-rand` is affecting the initial generation only. This might because `AI::Genetic::Pro` is generating random numbers by some specific module.

The input file `genePairedPrimerFile` is as described for script `saddle_pair.pl`. Other files with name prefix `outPrefix` are:
1. `<outPrefix>.score`: GA parameters at the beginning, and the best and the worst badness for each generation/iteration.
2. `<outPrefix>.bestN`: primer pair combination of best N individuals. This number can be specified by `-report`. These files can be applied to options `-badness` and `-startWith` of `saddle_pair.pl`.
3. `<outPrefix>.GA`: GA state file. This file should be loadable by the `AI::Genetic::Pro` module.

Under WSL ubuntu of my desktop PC, this script (with default parameters) took about one hour for 96 targets.

## saddleGA_pair_divN.pl

Two-layered Genetic Algorithm implementation that divides targets into groups of specified sizes (the first layer GA) and reduce of all within-group _badness_ (the second layer GA, i.e. `saddleGA_pair.pl`). Would be suitable if there are too many targets for single multiplex PCR reaction or if small size mulplex PCR reactions were desired.

```
saddle-implement/bin$ ./saddleGA_pair_divN.pl
Usage: saddleGA_pair_divN.pl [options] <genePairedPrimerFile> <outPrefix>
    Group division:
        -group      : if a number, divide targets into this number of groups
                      if a file of group name-target numbers, divide targets
                      into specified groups with specified numbers (default: 2)
        -report     : number of top answers to report (default: 5)
        -maxIt      : maximum iteration (default: 100)
        -GApopulation1 : GA parameter population (default: 200)
        -GApreserve1   : GA parameter preserve (default: 1)
        -GAcrossover1  : GA parameter crossover (default: 0.9)
        -GAmutation1   : GA parameter mutation (default: 0.05)
        -startState    : load saved GA state (default: "")
    SADDLE GA options:
        -minOverlap : minimum reversecomplement to compute badness (default: 4)
        -testIt     : test iteration (default: 150)
        -GApopulation2 : GA parameter population (default: 1000)
        -GAcrossover2  : GA parameter crossover (default: 0.6)
        -GAmutation2   : GA parameter mutation (default: 0.2)
    Others:
        -rand       : random seed (default: 0)
        -cpu        : number of CPUs for computation (default: 2)
```

The input file `genePairedPrimerFile` is as described for script `saddle_pair.pl`. Other files with name prefix `outPrefix` are:
1. `<outPrefix>.score`: GA parameters at the beginning, and the best and the worst badness for each generation/iteration.
2. `<outPrefix>.bestN`: Text files showing group composition and primer pair combination of best N individuals. This number can be specified by `-report`. These files can be applied to options `-badness` of `saddle_pair.pl`. See below for an example.
3. `<outPrefix>.state`: GA state file. This file can be loaded by using the `-startState` option.

### Example of using the state file

The `.state` file is for storing GA state so that `saddleGA_pair_divN.pl` can be executed multiple times with different GA parameters and passing the _state of GA_. In so doing, the best solutions can be kept evolving along all the executions. Here is an example. File `groups.txt` assigned four groups A/B/C/D for dividing 96 targets.

```
$ cat groups.txt
A 24
B 24
C 24
D 24
```

The first execution of `saddleGA_pair_divN.pl` (stage0) iterated the first layer GA with 20 generations (`-maxIt 20`). After the 20 generations, we can see the best and worst _badness_ scores in the population in file `stage0.score`. The population was saved in file `stage0.state`
```
$ saddleGA_pair_divN.pl -group groups.txt -cpu 35 -maxIt 20 -GApopulation1 200 -GApreserve1 10 -GAcrossover1 0.6 -GAmutation1 0.1 -testIt 10 -GApopulation2 100 genePrimers_220a.txt stage0
(output deleted)

$ tail -n 2 stage0.score
19 21054.1035539059 28231.4290461235
20 21010.8588862839 28114.5916906677

$ wc -l stage0.state
116209 stage0.state
```

The GA state was passed to the next execution (stage1) by assigning the `-startState` option and the same `-group` setting. In this execution, we assigned the second layer GA to iterate 20 generations. In file `stage1.score`, it was observed that the best and worst _badness_ scores at iteration 0 is exactly the same with those at the 20th iteration of stage0. This showed that the _state_ of the GA was actully passed from stage0 to stage1. In so doing, we can apply different GA parameters to the two layers in a series executions of `saddleGA_pair_divN.pl` to achieve better performance then monotone GA setting. 
```
$ saddleGA_pair_divN.pl -startState stage0.state -group groups.txt -cpu 35 -maxIt 20 -GApopulation1 200 -GApreserve1 10 -GAcrossover1 0.6 -GAmutation1 0.1 -testIt 20 -GApopulation2 100 genePrimers_220a.txt stage1
(output deleted)

$ head -17 stage1.score
GA related parameters:
Level 1:
    -maxIt         : 20
    -GApopulation1 : 200
    -GApreserve1   : 10
    -GAcrossover1  : 0.6
    -GAmutation1   : 0.1
Level 2:
    -testIt        : 20
    -GApopulation2 : 100
    -GAcrossover2  : 0.6
    -GAmutation2   : 0.2

Iteration best and worst
0 21010.8588862839 28114.5916906677
1 20478.9256562984 26327.8371593806
2 20154.4374091816 26380.2592372486
```

### Example of using the bestN file

The following perl one-liner helps to calculate within group badness scores and total badness scores in a `.bestN` file. Remember to adjust the path to `saddle_pair.pl`.
```
$ cat stage4.best0 | perl -MIPC::Open2 -ne 'chomp; @t=split; push @{$hash{shift @t}}, join("\t",@t); if(eof){ for $k (sort keys %hash){ open2(COUT,CIN,"/path/to/saddle_pair.pl -badness /dev/stdin"); print CIN join("\n",@{$hash{$k}})."\n"; close CIN; while($outline=<COUT>){ chomp $outline; if(length($outline)>0){ $sum+=$outline; print "$k\t$outline\n"} } } print "SUM\t$sum\n"; }'
A       3468.85676282468
B       3872.91467408912
C       6626.79507158898
D       3014.22188743838
SUM     16982.7883959412
```
