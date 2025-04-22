# saddle-implement

Implement of the iteration part of the SADDLE framework (PMID: 35410464) for multiplex PCR primer design. Two versions in this repository.

## `saddle_pair.pl`

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

Under WSL ubuntu of my desktop PC, this script took less than 30 minutes for 96 targets.

## `saddleGA_pair.pl`

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
2. `<outPrefix>.bestN`: primer pair combination of best N individuals. This number can be specified by `-report`. These files can be applied to options `-badness` and `startWith` of `saddle_pair.pl`.
3. `<outPrefix>.GA`: GA state file. This file should be loadable by the `AI::Genetic::Pro` module.

Under WSL ubuntu of my desktop PC, this script (with default parameters) took about one hour for 96 targets.
