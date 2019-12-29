source("classifier-base.R")
source("uams-70.R")

coefficients <- uams.70.get.coefficients()
mapping <- uams.70.probe.to.entrez.mapping()
classifier.df <- translate.probe.coefficients.to.gene.coefficients(coefficients, mapping)


ensg.mapping <- uams.70.probe.to.ensg.mapping()
ensg.classifier.df <- translate.probe.coefficients.to.gene.coefficients(coefficients, ensg.mapping)

coefficients.uams17 <- uams.17.get.coefficients()
mapping.uams17 <- uams.17.probe.to.entrez.mapping()
classifier.uams17.df <- translate.probe.coefficients.to.gene.coefficients(coefficients.uams17, mapping.uams17)


