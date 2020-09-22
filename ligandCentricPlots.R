library(data.table)

# These are all single protein targets in ChEMBL
bioactivityDT <- fread("zcat ./ChEMBL25dumps/ChEMBLstandardisedBioactivity.tsv.gz",
                       col.names = c("molregno","ChEMBLid","operator","value","unit","measurement","TID","name","species"))

ChemblID2Murko <- fread("zcat ./ChEMBL25dumps/ChEMBLid2SMILESplusMurckoScaffolds.tsv.gz")

ChEMBLid2ClinicalStatus <- fread("zcat ./ChEMBL25dumps/ChEMBLclinicalStatus.tsv.gz",
                                 col.names = c("molregno","ChEMBLid","prefName","therapetuicFlag","maxPhase","molType","firstApproved"))

ChEMBLid2ClinicalStatus[!is.na(firstApproved), decadeApproved := str_replace_all(firstApproved, "\\d$","0's")]

bioactivityDT[ChemblID2Murko, scaffold := i.scaffold, on = "ChEMBLid"]
bioactivityDT[ChEMBLid2ClinicalStatus, maxPhase := i.maxPhase, on = "ChEMBLid"]
bioactivityDT[ChEMBLid2ClinicalStatus, decadeApproved := i.decadeApproved, on = "ChEMBLid"]

nCompoundAnnotationsPerTarget <- bioactivityDT[, .(nAnnot = .N,
                                                   nDistinct = unique(.SD)[,.N]),
                                               by = TID,
                                               .SDcols = c("ChEMBLid")]

#This actually takes about five minutes
system.time(nDistinctTargetAnnotationsPerCompound <- bioactivityDT[, .(nAnnot = .N,
                                                           nDistinct = unique(.SD)[,.N]),
                                                       by = .(ChEMBLid,maxPhase,decadeApproved),
                                                       .SDcols = c("TID")])

nDistinctTargetAnnotationsPerCompound %>%
  ggplot(aes(y = nAnnot, x = phase, fill = preclinical)) +
  
  geom_hline(yintercept = 10^c(1,2,3,4), linetype = "twodash", colour = "#4d4d4d") +
  
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_jitter(data = plotDT[maxPhase > 0], colour = "#4d4d4d", height = 0) +
  
  
  geom_boxplot(data = nDistinctTargetAnnotationsPerCompound[maxPhase == 0], outlier.colour = "#4d4d4d") +
  geom_boxplot(data = nDistinctTargetAnnotationsPerCompound[maxPhase > 0],outlier.shape = NA, alpha = 0.8) +
  
  scale_y_log10(breaks = 10^as.vector(sapply(0:5, function(base) base + c(0,0.39794,0.69897))),
                labels = c('1','3','5',
                           '10','25','50',
                           '100','250','500',
                           '1000','2500','5000',
                           '10,000','25,000','50,000',
                           '100,000','250,000','500,000')) +
  scale_x_discrete(labels=my_xlab) +
  scale_fill_manual(values = c("#9999FF","#DCDCDC")) +
  theme_bw() +
  guides(fill=FALSE) +
  #facet_wrap(~type, nrow = 2) +
  labs(y = "Datapoint Count",x = NULL)






#This takes about 3 minutes after setDTthreads(6)
system.time(nDistinctTargetAnnotationsPerScaffold <- bioactivityDT[, .(nAnnot = .N,
                                                           nDistinct = unique(.SD)[,.N],
                                                           maxPhase = .SD[,max(maxPhase)]),
                                                       by = scaffold,
                                                       .SDcols = c("TID")])

nDistinctTargetAnnotationsPerScaffold[maxPhase == 0, .N, by = scaffold]
nDistinctTargetAnnotationsPerScaffold[maxPhase == 0,.N,by = nDistinct][order(N)]

nannotFreq <- nDistinctTargetAnnotationsPerScaffold[maxPhase == 0,.N,by = nAnnot]

setorderv(nannotFreq, "nAnnot", order = 1)

nannotFreq[, prop := N/sum(N)]
nannotFreq[, cumProp := cumsum(prop)]
# Roughly a third are N=1, a third 2 <= N <= 6 and a third the rest


clinicalLevels <-  c("Phase I","Phase II","Phase III","Phase IV", "Preclinical or earlier\n(1 annotation)", "Preclinical or earlier\n(2-5 annotations)","Preclinical or earlier\n(6+ annotations)")
nDistinctTargetAnnotationsPerScaffold[maxPhase != 0, phase := factor(clinicalLevels, levels = clinicalLevels[c(5,6,7,1,2,3,4)])[maxPhase]]

nDistinctTargetAnnotationsPerScaffold[maxPhase ==0, phase := "Preclinical or earlier\n(6+ annotations)"]
nDistinctTargetAnnotationsPerScaffold[maxPhase ==0 & nAnnot == 1, phase := "Preclinical or earlier\n(1 annotation)"]
nDistinctTargetAnnotationsPerScaffold[maxPhase ==0 & nAnnot >=2 & nAnnot <= 5, phase := "Preclinical or earlier\n(2-5 annotations)"]

nDistinctTargetAnnotationsPerScaffold[, preclinical := FALSE]
nDistinctTargetAnnotationsPerScaffold[maxPhase ==0, preclinical := TRUE]

library(xtable)
nDistinctTargetAnnotationsPerScaffold[,.(`Annotation Count` = sum(nAnnot), `Scaffold Count` = .N),by= .(Status = phase)][order(`Annotation Count`)] %>%
  xtable() %>%
  print(type = "html", file = "SupplementaryTable1.html")
  

# prepare a special xlab with the number of obs for each group
my_xlab <- paste(levels(nDistinctTargetAnnotationsPerScaffold$phase),"\n[N=",table(nDistinctTargetAnnotationsPerScaffold$phase),"]",sep="")


fread("numberOfPerScaffoldAnnotations.tsv", col.names = c("Collection","Annotation Sum"))


plotDT <- melt(nDistinctTargetAnnotationsPerScaffold,
     id.vars = c("preclinical","phase","maxPhase"),
     measure.vars = c("nAnnot","nDistinct")) 

plotDT[variable == "nDistinct", type := "Distinct Single-Protein Targets Annotated"]
plotDT[variable == "nAnnot", type := "All (Possibly Duplicate) Annotations Recored"]




pMurkoAnnot <- plotDT[variable == "nAnnot"] %>%
  ggplot(aes(y = value, x = phase, fill = preclinical)) +
  
  geom_hline(yintercept = 10^c(1,2,3,4,5), linetype = "twodash", colour = "#4d4d4d") +
  
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_jitter(data = plotDT[maxPhase > 0][variable == "nAnnot"], colour = "#4d4d4d", height = 0) +
  
  
  geom_boxplot(data = plotDT[maxPhase == 0][variable == "nAnnot"], outlier.colour = "#4d4d4d") +
  geom_boxplot(data = plotDT[maxPhase > 0][variable == "nAnnot"],outlier.shape = NA, alpha = 0.8) +
  
  scale_y_log10(breaks = 10^as.vector(sapply(0:5, function(base) base + c(0,0.39794,0.69897))),
                labels = c('1','3','5',
                           '10','25','50',
                           '100','250','500',
                           '1000','2500','5000',
                           '10,000','25,000','50,000',
                           '100,000','250,000','500,000')) +
  scale_x_discrete(labels=my_xlab) +
  scale_fill_manual(values = c("#9999FF","#DCDCDC")) +
  theme_bw(base_size = 15) +
  guides(fill=FALSE) +
  #facet_wrap(~type, nrow = 2) +
  labs(y = "Number Of Annotations Per Murcko Scaffold\n(Possibly Duplicate)", x = NULL,
       title = "ChEMBL25 annotations recorded for distinct chemotypes (Murcko scaffolds), stratified by clinical status",
       subtitle = "There are more datapoints for scaffolds relating to drugs that have completed Phase IV (far right)\nthan the bottom two-thirds least annotated scaffolds (leftmost two boxes)")



nDistinctTargetAnnotationsPerCompound[maxPhase == 4 & !is.na(decadeApproved)] %>%
  ggplot(aes(x = decadeApproved, y = nAnnot)) +
  geom_boxplot() +
  scale_y_log10() +
  geom_jitter(height = 0)




decadesBackwards <- nDistinctTargetAnnotationsPerCompound[!is.na(decadeApproved),sort(unique(decadeApproved), decreasing = FALSE)]

nDistinctTargetAnnotationsPerCompound[,decade := factor(decadeApproved, levels = decadesBackwards)]


nPhaseIVcompounds <- nDistinctTargetAnnotationsPerCompound[!is.na(decadeApproved) & maxPhase==4, .N]


pNmarketed <- nDistinctTargetAnnotationsPerCompound[!is.na(decade), .N, by = decade] %>%
  ggplot(aes(x = decade, y = N)) +
  geom_bar(stat="identity", fill = "#DCDCDC", colour = "black", alpha = 0.8) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0,400,by=50)) +
  guides(fill=FALSE) +
  labs(x="Decade of Approval", y = "Number of Compounds Approved",
       title = glue("Decade of approval in ChEMBL25 for {nPhaseIVcompounds} marketed drugs trialed in Phase IV"),
       subtitle = "Most maketed drugs documented are more recent than the 1980's")


nAnnotByDecade <-nDistinctTargetAnnotationsPerCompound[!is.na(decade), .(count = sum(nAnnot)), by = decade] 

setorder(nAnnotByDecade, decade)

nAnnotByDecade[,percentage := 100*cumsum(count)/sum(count)]

pCumulativeProportion <- nAnnotByDecade %>%
  ggplot(aes(x = decade, y = percentage)) +
  geom_bar(stat="identity", fill = "#9999FF", colour = "black" , alpha = 0.8) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0,100,by = 10)) +
  guides(fill=FALSE) +
  labs(x="Decade of Approval", y = "Cumulative Percentage Of Compound Annotations\n(Possibly Duplicated)",
       title = glue("Proportion of data in ChEMBL25 attibuted to {nPhaseIVcompounds} marketed drugs trialed in Phase IV"),
       subtitle = "Half of maketed-drug data in ChEMBL is more recent than 1990, with little from the 'Golden Age' of the 1950's and 1960's\nIt would seem reasonable to suggest that the trend for these molecules is representative of the larger database")


pCumulativeProportionMod <- pCumulativeProportion +  labs(title = NULL, subtitle = NULL)

pNmarketedMod <- pNmarketed +  labs(x = NULL, subtitle = NULL, title = NULL)

library(patchwork)

pNmarketedMod/pCumulativeProportionMod
# Plot cumulative
# https://gist.github.com/zachcp/f2429fc17cf6c59d0967
