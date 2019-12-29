# UAMS-17 Classifier
This is an implementation of the UAMS-17 gene classifier based on publication and
patent by Shaughnessy et al.
http://www.bloodjournal.org/content/109/6/2276
https://www.google.com/patents/US20080187930

# Input parameters

`uams.17`: Function for calling the 17 gene model

`eset`: An ExpressionSet formatted dataset (needs to be Affy U133A or U133 Plus 2.0)

`already.log2.transformed`: A boolean flag denoting
whether or not the dataset has already been log2
transformed (default = `FALSE`)


# To Run
For more details, please see test code in `test.model.R`


```
source('uams-17.R')
output <- uams.17(eset)
```
