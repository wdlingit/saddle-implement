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

Two-layered Genetic Algorithm implementation that divides targets into groups of specified sizes and reduce sum of all within-group _badness_. Would be suitable if there are too many targets for single multiplex PCR reaction or if small size mulplex PCR reactions were desired.

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

### Example of using the bestN file

```
$ cat stage4.best0 | perl -MIPC::Open2 -ne 'chomp; @t=split; push @{$hash{shift @t}}, join("\t",@t); if(eof){ for $k (sort keys %hash){ open2(COUT,CIN,"/path/to/saddle_pair.pl -badness /dev/stdin"); print CIN join("\n",@{$hash{$k}})."\n"; close CIN; while($outline=<COUT>){ chomp $outline; if(length($outline)>0){ $sum+=$outline; print "$k\t$outline\n"} } } print "SUM\t$sum\n"; }'
A       3468.85676282468
B       3872.91467408912
C       6626.79507158898
D       3014.22188743838
SUM     16982.7883959412
```
