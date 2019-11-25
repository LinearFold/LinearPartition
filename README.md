# LinearPartition: Linear-Time Approximation of RNA Folding Partition Function and Base Pairing Probabilities

This repository contains the C++ source code for the LinearPartition project, the first linear-time partition function and base pair probabilities calculation algorithm/software for RNA secondary structures.

LinearPartition: Linear-Time Approximation of RNA Folding Partition Function and Base Pairing Probabilities

He Zhang, Liang Zhang, David Mathews, Liang Huang*

\* corresponding author


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
Prints out beamsize, Log Partition Coefficient or free energy of ensumble (-V mode) and runtime information. (default False)
```
--sharpturn
```
Enable sharpturn. (default False)
```
-o
```
Outputs base pairing probability matrix to a file with user specified name. (default False)
```
--prefix
```
Outputs base pairing probability matrices to files with user specified prefix name. (default False)
```
-p
```
Partition function calculation only. (default False)

## Example Run Predict
```
cat testseq | ./linearpartition -V --prefix testseq_output
UGAGUUCUCGAUCUCUAAAAUCG
Free Energy of Ensumble: -1.96 kcal/mol
Outputing base pairing probability matrix to testseq_output_1...
Done!
AAAACGGUCCUUAUCAGGACCAAACA
Free Energy of Ensumble: -9.41 kcal/mol
Outputing base pairing probability matrix to testseq_output_2...
Done!
AUUCUUGCUUCAACAGUGUUUGAACGGAAU
Free Energy of Ensumble: -7.72 kcal/mol
Outputing base pairing probability matrix to testseq_output_3...
Done!
UCGGCCACAAACACACAAUCUACUGUUGGUCGA
Free Energy of Ensumble: -9.09 kcal/mol
Outputing base pairing probability matrix to testseq_output_4...
Done!
GUUUUUAUCUUACACACGCUUGUGUAAGAUAGUUA
Free Energy of Ensumble: -13.58 kcal/mol
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

## Example Run partition function calculation only
```
echo GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGCCCUGGGUUCAAAUCCCAGCGAGUCCACCA | ./linearpartition -V -p --verbose
beam size: 100
GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGCCCUGGGUUCAAAUCCCAGCGAGUCCACCA
Free Energy of Ensumble: -32.14 kcal/mol
Partition Function Calculation Time: 0.005509 seconds.
```
