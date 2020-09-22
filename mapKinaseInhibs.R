library(purrr)
library(data.table)

# Downloaded from https://www.ppu.mrc.ac.uk/list-clinically-approved-kinase-inhibitors
# on 6th Feb 2020
# 'Small molecular inhibitors of protein kinases approved for clinical use by the FDA'
kinaseDT <- fread("~/kinaseInhib.raw", sep = "|" , header = FALSE)[V1 %like% "(ib|il|mus)\\b"]

kinaseDT[, drug := str_extract(V1, "\\b\\w+(ib|il|mus)\\b")]


kinaseMapDT <- kinaseDT$drug %>%
  map(~{print(.x); 
        resDT <- operations$compound_annotations(annotations = .x);
        fwrite(resDT, file = "~/kinaseInhibMatches.tsv",
               sep = "\t", quote = FALSE, append = TRUE);
        return(resDT)})


kinaseMapDT[!sapply(kinaseMapDT, is.data.table)] # all map!


kinaseInhibDT <- rbindlist(kinaseMapDT)[,.(annotation,nexcompound_id)][!duplicated(toupper(annotation))]
kinaseInhibDT[,.N,by = annotation] # 1-2-1 mapping

ChemblCmpdsDT <- operations$available_compounds(count = TRUE, datasources = "%ChEMBL%")

ChemblCmpdsActive10uMDT <- operations$available_compounds(count = TRUE,
                                                          nm_upper_bound = 10000,
                                                          datasources = "%ChEMBL%")

ChemblCmpdsDT[ChemblCmpdsActive10uMDT,
              positiveEntries := i.n_datapoints, 
              on = "nexcompound_id"]

ChemblCmpdsDT[is.na(positiveEntries), positiveEntries := 0]
ChemblCmpdsDT[, negativeEntries := (n_datapoints- positiveEntries) ]

setkey(ChemblCmpdsDT, nexcompound_id)


ChemblCmpdsDT[,fracPositive := positiveEntries/(1+n_datapoints)]

ChemblCmpdsDT[, qplot(fracPositive*20000)]

ChemblCmpdsDT[, sum(positiveEntries)]/ChemblCmpdsDT[, sum(negativeEntries)]

ChemblCmpdsDT[, sum(n_datapoints)]


operations$available_compounds


kinaseInhibCounts <- ChemblCmpdsDT[kinaseInhibDT$nexcompound_id, , nomatch = 0]

kinaseInhibCounts[,qplot(n_datapoints)] +
  scale_x_log10() +
  theme_bw()

library(reldist)

kinaseInhibCounts[,gini(n_datapoints)]
