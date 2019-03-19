R package: InteRact
================
Guillaume Voisinne
2019 - 03 - 18

[![Travis-CI Build Status](https://travis-ci.org/VoisinneG/InteRact.svg?branch=master)](https://travis-ci.org/VoisinneG/InteRact)

InteRact
========

Analysis of affinity purification/tandem MS data with R

-   flexible : no strict input format
-   publication-ready figures

More info on the [website](https://voisinneg/github.io/InteRact)

Install
-------

The package is available from github:

    devtools::install_github("VoisinneG/InteRact")

A shiny based GUI can be accessed [here](https://voisinneg.shinyapps.io/interact/)

Use
---

Load example protein group dataset and run `InteRact`:

``` r
data("proteinGroups_Cbl")
names(proteinGroups_Cbl)[1:10]
```

    ##  [1] "Protein.IDs"                   "Majority.protein.IDs"         
    ##  [3] "Peptide.counts..all."          "Peptide.counts..razor.unique."
    ##  [5] "Peptide.counts..unique."       "Protein.names"                
    ##  [7] "Gene.names"                    "Fasta.headers"                
    ##  [9] "Number.of.proteins"            "Peptides"

``` r
res <- InteRact(proteinGroups_Cbl, bait_gene_name = "Cbl")
```

    ## Contaminant proteins discarded
    ## Proteins with no gene name available discarded
    ## Number of theoretically observable peptides unavailable : used MW instead
    ## Merge protein groups associated to the same gene name (sum of intensities) 
    ## Rescale median intensity across conditions
    ## Replace missing values and perform interactome analysis for 1 replicates
    ## Nrep=1
    ## Averaging 1 interactomes

Identify specific interactors

``` r
res <- identify_interactors(res, 
                            p_val_thresh = 0.001, 
                            fold_change_thresh = 3, 
                            n_success_min = 2, 
                            consecutive_success = TRUE)
print(res$interactor)
```

    ## [1] "Sh3kbp1"  "Ubash3a"  "Crkl"     "Pik3r1"   "StrepTag"

Show summary table

``` r
sum_tbl <- summary_table(res)
head(sum_tbl[, c("bait", "names", "max_fold_change" ,"max_stoichio")])
```

    ##   bait    names max_fold_change max_stoichio
    ## 1  Cbl      Cbl       4338.4207   1.00000000
    ## 2  Cbl StrepTag       1130.9978   1.85184340
    ## 3  Cbl     Crkl        881.4774   0.63190341
    ## 4  Cbl  Ubash3a        667.1752   0.08654619
    ## 5  Cbl  Sh3kbp1        468.3102   0.06527946
    ## 6  Cbl   Pik3r1        836.4604   0.05197326
