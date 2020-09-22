

library(data.table)
library(ggplot2)

# p(M+) = p(M+|D+)p(D+) + p(M+|D-)p(D-)

# p(D+|M+)=p(M+|D+)p(D+)/p(M+)

`p(D+)` = c("0.1%","1%","5%","10%")

bayesianModelTable <- cross_df(
  list(
    positiveDataPC = factor(`p(D+)`, levels = `p(D+)`),
    TPR = seq(0,1,0.01),
    FPR = c(0.001,0.005,0.01,0.025,0.05,0.1,0.2)
  )
)

setDT(bayesianModelTable)


positiveDataGivenPositiveModel <- function(TPR, FPR, positiveDataPercentage){
  
  positiveDataFrac <- as.numeric(sub("%", "", positiveDataPercentage))/100
  
  TPR*positiveDataFrac/(TPR*positiveDataFrac + FPR*(1-positiveDataFrac))
}

convert2percentage <- function(fraction){
  
  glue("{100*fraction}%")
}

bayesianModelTable[, `p(D+|M+)` := positiveDataGivenPositiveModel(TPR,FPR,positiveDataPC), by = .I]

bayesianModelTable[, facetField := glue("Activity on {as.numeric(sub('%', '', positiveDataPC))/100*1000} out of 1,000 proteins: {positiveDataPC}", .envir = .SD)]

bayesianModelTable[,facetField := factor(facetField,
                                         levels = c("Activity on 1 out of 1,000 proteins: 0.1%",
                                                    "Activity on 10 out of 1,000 proteins: 1%",
                                                    "Activity on 50 out of 1,000 proteins: 5%",
                                                    "Activity on 100 out of 1,000 proteins: 10%"))]


library(ggrepel)

pBayesAnnot <- bayesianModelTable %>%
  ggplot(aes(x = TPR, y = `p(D+|M+)`)) +
  geom_hline(yintercept = 0.75, linetype = "dotted", colour = "#4d4d4d") +
  geom_line(aes(group = as.factor(FPR), colour = FPR), size = 3,alpha = 0.8) +
  geom_label_repel(data = bayesianModelTable[TPR == 0.75], aes(label = glue("{100*FPR}%"))) +
  facet_wrap(~ facetField, nrow = 2) +
  theme_bw(base_size = 15) +
  scale_x_continuous(breaks = seq(0,1,0.1), labels = glue("{100*seq(0,1,0.1)}%")) +
  scale_colour_gradient(low = "#DCDCDC", high = "#9999FF", 
                        breaks = c(seq(0,0.08,0.025),seq(0.1,0.2,0.05)),
                        labels = c("0%","2.5%","5%","7.5%","10%","15%","20%"),
                        name = "P(M+|D-)\nFalse Positive Rate") +
  labs(y = "P(D+|M+)\nProbability compound is in reality active given model predicts active",
       x = "P(M+|D+)\nTrue Positive Rate (TPR)",
       subtitle = "Improvements to False Positive Rate can have a dramatic effect on the reliability of model predictions for rare events",
       title = "Model reliability explored using Bayes rule: a function of TPR, FPR and the background activity rate") +
  guides(fill=FALSE)


pBayesAnnot + labs(subtitle = NULL)


