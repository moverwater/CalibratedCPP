rates = read.csv("bd_rates.csv")

rho = 0.1

taxa = paste0("leaf_",1:100,collapse = ",")

clade1 = list(Age = 1.5,Taxa = paste0("leaf_",as.character(45:55),collapse=","))
clade2 = list(Age = 1.2,Taxa = paste0("leaf_",as.character(45:51),collapse=","))
clade3 = list(Age = 2,Taxa = paste0("leaf_",as.character(87:90),collapse=","))

cladeCalibrations = list(clade1,clade2,clade3)

stemAge = 3.0

for (i in 1:nrow(rates)){
  tree = simCladeCalibratedAges(birthRate = rates$birthRates[i], 
                                deathRate = rates$deathRates[i], 
                                samplingProbability = rho, 
                                taxa = taxa,
                                cladeCalibrations = cladeCalibrations,
                                stemAge = stemAge)
  writeLines(text = tree, paste0("./trees/tree_",i,".txt"))
}
