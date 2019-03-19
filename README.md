R package: InteRact
================
Guillaume Voisinne
2019 - 03 - 19

[![Travis-CI Build Status](https://travis-ci.org/VoisinneG/InteRact.svg?branch=master)](https://travis-ci.org/VoisinneG/InteRact)

InteRact
========

Analysis of affinity purification/tandem MS data with R

-   flexible : no strict input format
-   publication-ready figures

More info on the [website](https://voisinneg.github.io/InteRact)

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
                            n_success_min = 1, 
                            consecutive_success = TRUE)
print(res$interactor)
```

    ##  [1] "Cbl"      "Mccc1"    "Sh3kbp1"  "Ubash3a"  "Crkl"     "Pik3r1"  
    ##  [7] "Ywhah"    "Ywhag"    "Pik3ca"   "Pik3r2"   "Inpp5d"   "StrepTag"
    ## [13] "Pik3cd"   "Grap"     "Sdha"     "Tbce"     "Pccb"     "Rbm25"   
    ## [19] "Myl6"     "Unc13d"

Show summary table

``` r
sum_tbl <- summary_table(res)
head(sum_tbl[, c("bait", "names", "max_fold_change" ,"max_stoichio", "is_interactor")])
```

    ##   bait    names max_fold_change max_stoichio is_interactor
    ## 1  Cbl      Cbl      4309.24604    1.0000000             1
    ## 2  Cbl StrepTag      1095.88799    1.8517938             1
    ## 3  Cbl     Pccb        22.03694    1.2265462             1
    ## 4  Cbl     Crkl      1342.62339    0.6320675             1
    ## 5  Cbl    Ywhah       124.87974    0.2295337             1
    ## 6  Cbl    Ywhag       150.97868    0.2230696             1
