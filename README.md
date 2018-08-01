# Multi-circular dynamic time warping alignment

## Introduction

In this repository you will find some functions that extend the capabilities of dynamic time warping alignment from the `dtw` package. One of its limitations is that in pairwise local alignments (open end, open begin or both), the gaps can only be placed at the query sequence and not the reference sequence. 
<<<<<<< HEAD

This can be an issue if one tries to align data that contains, for example, a gap at the beginning of the query sequence and a gap at the end of the reference sequence:

```
Ref:   TEACUP-----
          |||
Query: ---CUPBOARD
```

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

