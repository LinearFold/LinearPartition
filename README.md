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
