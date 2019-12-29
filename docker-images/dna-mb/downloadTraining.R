library(synapseClient)
synapseLogin()
directoryPath = "~/GitHub/Celgene-Multiple-Myeloma-Challenge/docker-images/dna-mb/trainingData"
unlink(directoryPath,recursive=T)
chal_data_table = synQuery("SELECT name,id FROM file WHERE parentId==\'syn7222257\'")
sapply(chal_data_table$file.id, function(x) { synGet(x, downloadLocation=directoryPath)})


directoryPath = "~/GitHub/Celgene-Multiple-Myeloma-Challenge/docker-images/dna-mb/trainingData"
unlink(directoryPath,recursive=T)
chal_data_table = synQuery("SELECT name,id FROM file WHERE parentId==\'syn7250060\'")
sapply(chal_data_table$file.id, function(x) { synGet(x, downloadLocation=directoryPath)})

#directoryPath = "~/GitHub/Celgene-Multiple-Myeloma-Challenge/docker-images/docker-dna-ml/testData"
#unlink(directoryPath,recursive=T)
#chal_data_table = synQuery("SELECT name,id FROM file WHERE parentId==\'syn7222262\'")
#sapply(chal_data_table$file.id, function(x) { synGet(x, downloadLocation=directoryPath)})

#directoryPath = "~/GitHub/Celgene-Multiple-Myeloma-Challenge/docker-images/docker-dna-ml/training-vcf"
#unlink(directoryPath,recursive=T)
#chal_data_table = synQuery("SELECT name,id FROM file WHERE parentId==\'syn9891644\'")
#sapply(chal_data_table$file.id, function(x) { synGet(x, downloadLocation=directoryPath)})
