# Multi-circular dynamic time warping alignment

## Introduction

In this repository you will find some functions that extend the capabilities of dynamic time warping alignment from the `dtw` package. 

One of its limitations is that in pairwise local alignments (open end, open begin or both), the gaps can only be placed at the query sequence and not the reference sequence. This can be an issue if one tries to align data that contains, for example, a gap at the beginning of the query sequence and a gap at the end of the reference sequence (see next example):

```
Type:       1              2           3           4
Supported:  No             Yes         Yes         Yes
Ref:        TEACUP-----    CUPBOARD    CUPBOARD    CUPBOARD
               |||         |||            |||||       ||||
Query:      ---CUPBOARD    CUP-----    ---BOARD    ---BOAR-
```

To be able to generate the type 1 alignments, we use a progressive elongation of the length of the query sequence to test all possible second and third type of alignments (see next example):

```
Length:     1           2           3           4           5           6        
Edit:       1           2           0           4           5           6       
Ref:    CUPBOARD    CUPBOARD    CUPBOARD    CUPBOARD    CUPBOARD    CUPBOARD               
        *           **          |||         ****        *****       ******         
Query:  P-------    UP------    CUP-----    ACUP----    EACUP---    TEACUP--
```

It is clear from the previous example that the 3rd alignment is the correct one. The only thing that remains is to add the TEA part. This same process is applied with a continuos time-series data. To decide which is the best alignment the normalized distance between sequences is used. The function `pairwiseUnivariateAlignment` tests all possible type1 alignments (sliding the query from the right and left sides), a local alignment with gaps at both sides of the query sequence, and a global alignment (no gaps in any sequence). After all possible alignments are tested, the best one is given as an output.



## Dependencies

* `dtw`  (required)
* `ggplot2` (recommended)
* `gridExtra` (recommended)

## Contents

* Functions
    * `pairwiseUnivariateAlignment`:
    * `multiUnivariateAlignmentProfile`:
    * `circularizeUnivariateSequence`:
    * `rotateSequence`:
    
* Demo: contains examples of all of the above functions. It uses as data a sequence generated from the `sin` function.

## Considerations

## In progress