# Test of EMC-92, UAMS-17, and UAMS-70 classifier
This is an implementation of the gene classifiers based on publication by Kuiper et al.
http://www.nature.com/leu/journal/v26/n11/full/leu2012127a.html
and accompanying source code in sourcescript.txt

The EMC-92 interface is well documented.  The UAMS-17 and UAMS-70 are not, but follow the
EMC-92 interface closely.

# To run
For more details, please see test code in `test.models.R`

Note that `test.models.R` pulls the expression and clinical annotations from Synapse.
If you do not want to do so, set `use.synapse` to `FALSE`.  And make sure that `data.file`
and `anno.file` are set to the locations of the MMRF expression _count_ data and clinical data,
respectively.

This script runs the EMC-92, UAMS-17, and UAMS-70 classifiers using the published threshold on all of the data, 
a fitted threshold on the test data, and the published threshold on the test data.  The Kaplan-Meier
curves and Cox hazard ratios for these three cases are output to the files `{emc92,uams17,uams70}-kms.pdf`.


```
source("test.models.R")
```

