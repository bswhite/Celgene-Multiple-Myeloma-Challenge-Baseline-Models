# UAMS-5 Classifier
This is an implementation of the UAMS-5 gene classifier based on publication by Heuck et al.
http://www.nature.com/leu/journal/v28/n12/full/leu2014232a.html

# Input parameters

`uams.5`: Function for calling the 5 gene model

`eset`: An ExpressionSet formatted dataset (needs to be Affy U133A or U133 Plus 2.0)

`already.log2.transformed`: A boolean flag denoting
whether or not the dataset has already been log2
transformed (default = `FALSE`)


# To Run
For more details, please see test code in `test.model.R`


```
source('uams-5.R')
output <- uams.5(eset)
```
