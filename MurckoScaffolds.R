library(data.table)
library(ggplot2)

library(foreach)
library(doParallel)
library(itertools)
library(purrr)

library(reticulate)

cl <- parallel::makeCluster(6)
registerDoParallel(cl)


ChEMBLstructures <- fread("zcat ./ChEMBL25dumps/ChEMBLid2SMILES.tsv.gz",
                          col.names = c("SMILES","ChEMBLid","molegno"))

tic <- Sys.time()

scaffoldsList <- foreach(SMILES = ichunk(ChEMBLstructures[,unique(SMILES)], chunkSize = 10000),
        .packages = c("data.table","reticulate", "purrr")) %dopar% {
  
  use_python("/usr/bin/python3")
  
  rdkit <- import("rdkit.Chem.Scaffolds.MurckoScaffold")
  
  safelyMapped <- map(SMILES ,
      
      #These are cannonical smiles
      safely(~ data.table(structure = .x,
                          scaffold = rdkit$MurckoScaffoldSmiles(.x)), otherwise = NA) )
  
  return(rbindlist( map(safelyMapped[map_lgl(safelyMapped, ~ is.data.table(.x$result))],  ~.x$result) ))
}

ChEMBLstructures[rbindlist(scaffoldsList) ,
                       scaffold := scaffold
                       , on =.(SMILES = structure)]

fwrite(ChEMBLstructures,
       './ChEMBL25dumps/ChEMBLid2SMILESplusMurckoScaffolds.tsv.gz',
       sep = "\t", col.names=TRUE, quote=FALSE)

toc <- Sys.time()

print(toc-tic)

