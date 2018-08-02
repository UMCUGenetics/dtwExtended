# dtwExtended

## Introduction

Dynamic time warping can be used to align time-series of continuous data. In R, the `dtw` package provides with a set of very useful tools to apply it. In this repository you will find some functions that extend the capabilities of this package. 

### Gaps both in reference and query sequence

One of the limitations of `dtw` is that in pairwise local alignments (open end, open begin or both), gaps can only be placed at the query sequence and not the reference sequence. This can be an issue if one tries to align sequences in which the beginning of one sequence aligns with the end of the other sequence (see example below). The following examples, although depicted as discrete sequences, should give an idea of the rationale behind.

```
Type:       1              2           3           4
Supported:  No             Yes         Yes         Yes
Ref:        TEACUP-----    CUPBOARD    CUPBOARD    CUPBOARD
               |||         |||            |||||       ||||
Query:      ---CUPBOARD    CUP-----    ---BOARD    ---BOAR-
```

### Sliding the query sequence approach

To be able to generate the type 1 alignments, we use a progressive elongation of the length of the query sequence to test all possible type 2 (see example below):

```
Length:     1           2           3           4           5           6        
Edit:       1           2           0           4           5           6       
Ref:    CUPBOARD    CUPBOARD    CUPBOARD    CUPBOARD    CUPBOARD    CUPBOARD               
        *           **          |||         ****        *****       ******         
Query:  P-------    UP------    CUP-----    ACUP----    EACUP---    TEACUP--
```

It is clear from the previous example that the 3rd alignment is the correct one. The only thing that remains is to add the TEA part at the beginning of the CUP query sequence. This same process is applied with continuous time-series data. To decide which is the best alignment the normalized distance between sequences is used. The function `pairwiseUnivariateAlignment` tests all possible type 1 alignments (sliding the query from both the right and left sides), a local alignment with gaps at both sides of the query sequence, and a global alignment (no gaps in any sequence). After all possible alignments are tested, the best one is given as an output.

### What can it be used for?

This type of alignment can be used if you are not able to measure the whole sequence of a time-series. You have several measurements that together show the complete sequence, using the overlaps between those measurements alignments can be made to get the whole picture of the sequence (see example below).

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

Another utility is to circularize a sequence. For example, time-series that involve circadian type of data my want to align their data without extending to more than one cycle (see example below). This can be done using the `circularizeUnivariateSequence` function.

```
First alignments
Cyclic sequence   THEQUICKBROWNFOXJUMPSOVERTHEQUICKBROWFOX
Fragment1         THEQUICKBRO
                           ||
Fragment2                  ROWNFOXJUM
                                 ||||
Fragment3                        XJUMPSOVERTHEQUICK    
                   ||||||||||
Fragment4          HEQUICKBROWNFOX                  (aligned to fragment1)
Outcome           THEQUICKBROWNFOXJUMPSOVERTHEQUICK

Circularization
                                THEQUICKBROWNFOXJUMPSOVERTHEQUICK
                                ||||||||
2nd half               JUMPSOVERTHEQUICK
Outcome                         THEQUICKBROWNFOXJUMPSOVER (aligned and duplicated parts removed)
```

## Contents

* Functions
    * `pairwiseUnivariateAlignment`: generates a local or global alignment based on the minimum normalized distance
    * `multiUnivariateAlignmentProfile`: generates a profile alignment based on multiple pairwise alignments
    * `circularizeUnivariateSequence`: reduces to 1 cycle a cyclic sequence
    * `rotateSequence`: sets the starting point of a circularized sequence
    
* Demo: contains examples of all of the above functions. It uses as data a sequence generated from the `sin` function.

## Considerations

* Profile value calculation: in `pairwiseUnivariateAlignment`, the `uni` column contains the `mean` of the two sequences. One should consider if this approach fits with ones requirements. This is specially true when using `multiUnivariateAlignmentProfile`, since if a sequence is the result of 4 previous sequences and it is aligned with an 'individual' sequence, both values will have the same weight. It could be argued that an 80%-20% weight should be applied in this case.

## In progress

* Multivariate alignment
* Alignment of independent data to the alignment

## Dependencies

* `dtw` (required)
* `ggplot2` (recommended)
* `gridExtra` (recommended)

