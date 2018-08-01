# Multi-circular dynamic time warping alignment

## Introduction

In this repository you will find some functions that extend the capabilities of the `dtw` package. 

One of the limitations of `dtw` is that in pairwise local alignments (open end, open begin or both), gaps can only be placed at the query sequence and not the reference sequence. This can be an issue if one tries to align sequences in which the beginning of one sequence aligns with the end of the other sequence (see example below). The following examples, although depicted as discrete sequences, should give an idea of the rationale behind.

```
Type:       1              2           3           4
Supported:  No             Yes         Yes         Yes
Ref:        TEACUP-----    CUPBOARD    CUPBOARD    CUPBOARD
               |||         |||            |||||       ||||
Query:      ---CUPBOARD    CUP-----    ---BOARD    ---BOAR-
```

To be able to generate the type 1 alignments, we use a progressive elongation of the length of the query sequence to test all possible type 2 (see example below):

```
Length:     1           2           3           4           5           6        
Edit:       1           2           0           4           5           6       
Ref:    CUPBOARD    CUPBOARD    CUPBOARD    CUPBOARD    CUPBOARD    CUPBOARD               
        *           **          |||         ****        *****       ******         
Query:  P-------    UP------    CUP-----    ACUP----    EACUP---    TEACUP--
```

It is clear from the previous example that the 3rd alignment is the correct one. The only thing that remains is to add the TEA part at the beginning of the CUP query sequence. This same process is applied with continuous time-series data. To decide which is the best alignment the normalized distance between sequences is used. The function `pairwiseUnivariateAlignment` tests all possible type 1 alignments (sliding the query from both the right and left sides), a local alignment with gaps at both sides of the query sequence, and a global alignment (no gaps in any sequence). After all possible alignments are tested, the best one is given as an output.

When can this type of alignment be useful? This type of alignment can be used if you are not able to measure the whole sequence of a time-series. You have several measurements that together show the complete sequence, using the overlaps between those measurements alignments can be made to get the whole picture of the sequence (see example below).

```
Complete sequence   THEQUICKBROWNFOXJUMPSOVERTHELAZYDOG
Fragment1           THEQUICKBRO
                             ||
Fragment2                    ROWNFOXJUM
                                   ||||
Fragment3                          XJUMPSOVERTHELA    
                                              ||||
Fragment4                                     HELAZYDOG
```

## Contents

* Functions
    * `pairwiseUnivariateAlignment`:
    * `multiUnivariateAlignmentProfile`:
    * `circularizeUnivariateSequence`:
    * `rotateSequence`:
    
* Demo: contains examples of all of the above functions. It uses as data a sequence generated from the `sin` function.

## Considerations

## In progress

## Dependencies

* `dtw` (required)
* `ggplot2` (recommended)
* `gridExtra` (recommended)
