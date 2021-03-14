# LinearPartition: Linear-Time Approximation of RNA Folding Partition Function and Base Pairing Probabilities

This repository contains the C++ source code for the LinearPartition project, the first linear-time partition function and base pair probabilities calculation algorithm/software for RNA secondary structures.

[LinearPartition: linear-time approximation of RNA folding partition function and base-pairing probabilities](https://academic.oup.com/bioinformatics/article/36/Supplement_1/i258/5870487). Bioinformatics, Volume 36, Issue Supplement_1, July 2020, Pages i258–i267. ISMB 2020

He Zhang, Liang Zhang, David Mathews, Liang Huang*

\* corresponding author

Web server: http://linearfold.org/partition


## Dependencies
gcc 4.8.5 or above; 
python2.7

## To Compile
```
make
```

## To Run
LinearPartition can be run with:
```
echo SEQUENCE | ./linearpartition [OPTIONS]

OR

cat SEQ_OR_FASTA_FILE | ./linearpartition [OPTIONS]
```
Both FASTA format and pure-sequence format are supported for input.

OPTIONS:
```
-b BEAM_SIZE
```
The beam size (default 100). Use 0 for infinite beam.
```
-V
```
Switches LinearPartition-C (by default) to LinearPartition-V.
```
--verbose
```
Prints out beamsize, Log Partition Coefficient or free energy of ensemble (-V mode) and runtime information. (default False)
```
--sharpturn
```
Enable sharpturn. (default False)
```
-o FILE_NAME
```
Outputs base pairing probability matrix to a file with user specified name. (default False)
```
-r FILE_NAME
```
Output base pairing probability matrix to a file with user specified name (overwrite if the file exists). (default False)
```
--prefix PREFIX_NAME
```
Outputs base pairing probability matrices to files with user specified prefix. (default False)
```
-p
```
Partition function calculation only. (default False)
```
-c THRESHOLD
```
Only output base pair probability larger than user specified threshold between 0 and 1. (DEFAULT=0.0)
```

--dumpforest
```
dump forest (all nodes with inside [and outside] log partition functions but no hyperedges) for downstream tasks such as sampling and accessibility (DEFAULT=None)

```
--MEA
```
get MEA structure, (DEFAULT=FALSE)

```
--gamma
```
set MEA gamma, (DEFAULT=3.0)

```
--MEA_output
```
output MEA structure to a file with user specified name (rewrite if the file exists) (DEFAULT=FALSE)


## Example: Run Predict
```
cat testseq | ./linearpartition -V --prefix testseq_output
UGAGUUCUCGAUCUCUAAAAUCG
Free Energy of Ensemble: -1.96 kcal/mol
Outputing base pairing probability matrix to testseq_output_1...
Done!
AAAACGGUCCUUAUCAGGACCAAACA
Free Energy of Ensemble: -9.41 kcal/mol
Outputing base pairing probability matrix to testseq_output_2...
Done!
AUUCUUGCUUCAACAGUGUUUGAACGGAAU
Free Energy of Ensemble: -7.72 kcal/mol
Outputing base pairing probability matrix to testseq_output_3...
Done!
UCGGCCACAAACACACAAUCUACUGUUGGUCGA
Free Energy of Ensemble: -9.09 kcal/mol
Outputing base pairing probability matrix to testseq_output_4...
Done!
GUUUUUAUCUUACACACGCUUGUGUAAGAUAGUUA
Free Energy of Ensemble: -13.58 kcal/mol
Outputing base pairing probability matrix to testseq_output_5...
Done!

echo GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGCCCUGGGUUCAAAUCCCAGCGAGUCCACCA | ./linearpartition -o output
GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGCCCUGGGUUCAAAUCCCAGCGAGUCCACCA
Log Partition Coefficient: 15.88268
Partition Function Calculation Time: 0.004343 seconds.
Base Pairing Probabilities Calculation Time: 0.003293 seconds.
Outputing base pairing probability matrix to output...
Done!
```

## Example: Run Partition Function Calculation Only
```
echo GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGCCCUGGGUUCAAAUCCCAGCGAGUCCACCA | ./linearpartition -V -p --verbose
beam size: 100
GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGCCCUGGGUUCAAAUCCCAGCGAGUCCACCA
Free Energy of Ensemble: -32.14 kcal/mol
Partition Function Calculation Time: 0.005509 seconds.
```

## Example: Run Prediction and Print MEA structure
```
echo GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGCCCUGGGUUCAAAUCCCAGCGAGUCCACCA | ./linearpartition -V --MEA
GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGCCCUGGGUUCAAAUCCCAGCGAGUCCACCA
Free Energy of Ensemble: -32.14 kcal/mol
(((((((..((((.......))))((((((((...)))))))).(((((.......))))))))))))....
```



## Example: Run Prediction and output ThreshKnot structure in bpseq format
```
echo GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGCCCUGGGUUCAAAUCCCAGCGAGUCCACCA | ./linearpartition -V --ThreshKnot --ThreshKnot_threshold 0.3
GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGCCCUGGGUUCAAAUCCCAGCGAGUCCACCA
Free Energy of Ensemble: -32.14 kcal/mol
1 G 68
2 G 67
3 G 66
4 C 65
5 U 64
6 C 63
7 G 62
8 U 0
9 A 0
10 G 24
11 A 23
12 U 22
13 C 21
14 A 0
15 G 0
16 C 0
17 G 0
18 G 0
19 U 0
20 A 0
21 G 13
22 A 12
23 U 11
24 C 10
25 G 43
26 C 42
27 U 41
28 U 40
29 C 39
30 C 38
31 U 37
32 U 36
33 C 0
34 G 0
35 C 0
36 A 32
37 A 31
38 G 30
39 G 29
40 A 28
41 A 27
42 G 26
43 C 25
44 C 0
45 C 61
46 U 60
47 G 59
48 G 58
49 G 57
50 U 0
51 U 0
52 C 0
53 A 0
54 A 0
55 A 0
56 U 0
57 C 49
58 C 48
59 C 47
60 A 46
61 G 45
62 C 7
63 G 6
64 A 5
65 G 4
66 U 3
67 C 2
68 C 1
69 A 0
70 C 0
71 C 0
72 A 0
```



References
-------------

Liang Zhang, He Zhang, David H Mathews, and Liang Huang. Threshknot: Thresholded probknot for improved RNA secondary structure prediction. arXiv preprint arXiv:1912.12796.