# EMC-92 Classifier
This is an implementation of the EMC-92 gene classifier based on publication by Kuiper et al.
http://www.nature.com/leu/journal/v26/n11/full/leu2012127a.html
and accompanying source code in sourcescript.txt

# Function

`emc.92.ensg`: Function to invoke the EMC-92 model on Ensembl gene expression data.

# Input parameters

`test.eset`: An ExpressionSet formatted _test_ dataset whose features are Ensembl genes.

`test.anno`: A data frame holding clinical annotation data corresponding to the _test_ dataset (default = `NULL`).  Should have columns ID (matching the sample names in test.eset), time, and event.  

`already.log2.transformed`: A boolean flag denoting
whether or not the dataset has already been log2
transformed (default = `FALSE`)

`train.eset`: An ExpressionSet formatted _training_ dataset whose features are Ensembl genes (default = `NULL`).

`train.anno`: A data frame holding clinical annotation data corresponding to the _training_ dataset (default = `NULL`).  Should have columns ID (matching the sample names in train.eset), time, and event.  

`threshold`:  A threshold to use to discriminate between high- and low-risk cases (default = `NULL`).

If threshold and train.eset are NULL, the default published threshold will be used.
If threshold is NULL and train.est and train.anno are non-null, the threshold
will be optimized based on the training data.

# To run
For more details, please see test code in `test.model.R`

Note that `test.model.R` pulls the expression and clinical annotations from Synapse.
If you do not want to do so, set `use.synapse` to `FALSE`.  And make sure that `data.file`
and `anno.file` are set to the locations of the MMRF expression _count_ data and clinical data,
respectively.

This script runs the EMC-92 classifier using the published threshold on all of the data, 
a fitted threshold on the test data, and the published threshold on the test data.  The Kaplan-Meier
curves and Cox hazard ratios for these three cases are output to the file `emc92-kms.pdf`.


```
source("test.model.R")
```

