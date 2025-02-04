# Import libraries
library(mia)
library(TreeSummarizedExperiment)
library(MultiAssayExperiment)

# Load data
data("FMT_data",  package = "anansi")
data("dictionary", package = "anansi")

# Convert to (Tree)SummarizedExperiment objects
metab_se <- SummarizedExperiment(assays = SimpleList(conc = as.matrix(FMT_metab)))
KO_tse <- TreeSummarizedExperiment(assays = SimpleList(counts = as.matrix(FMT_KOs)))

# Select functions that are represented in the dictionary
keep <- row.names(KO_tse) %in% sort(unique(unlist(anansi_dic)))
KO_tse <- KO_tse[keep, ]

# Remove features with less than 10% prevalence
KO_tse <- subsetByPrevalent(KO_tse,
                            assay.type = "counts",
                            prevalence = 0.1)

# Perform a centered log-ratio transformation on the functional counts assay
KO_tse <- transformAssay(KO_tse,
                         assay.type = "counts",
                         method = "clr",
                         pseudocount = TRUE)

# Prepare colData
coldata <- FMT_metadata
rownames(coldata) <- coldata$Sample_ID
coldata <- coldata[match(colnames(KO_tse), rownames(coldata)), ]

# Combine experiments into MultiAssayExperiment object
mae <- MultiAssayExperiment(
  experiments = ExperimentList(metabolites = metab_se, functions = KO_tse),
  colData = coldata
)

# Perform anansi analysis
out <- getAnansi(mae,
                 experiment1 = "metabolites", experiment2 = "functions",
                 assay.type1 = "conc", assay.type2 = "clr",
                 formula = ~ Legend, translate = TRUE,
                 X_translation = KO_translation,
                 Y_translation = cpd_translation)

# Select significant interactions
out <- out[out$model_full_q.values < 0.1, ]

colour_by = "type"
out["colour"] <- out[["type"]]

signif.threshold <- 0.1
out["signif"] <- out["model_disjointed_Legend_p.values"] < signif.threshold
out["signif"] <- factor(out[["signif"]], levels = c(TRUE, FALSE))


library(ggplot2)

# if colour_by
# if show.signif

p <- ggplot(data = out, aes(x = r.values, y = feature_X, fill = colour, alpha = signif)) + 
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red")+
  geom_point(shape = 21, size = 3) + 
  ggforce::facet_col(~ feature_Y, space = "free", scales = "free_y") +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.2), expand = c(0, 0)) +
  scale_y_discrete(limits = rev, position = "right") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 1/3),
      paste("Disjointed association\np <", signif.threshold)) +
  theme_bw() +
  labs(x = "Pearson's \u03c1", fill = colour_by) +
  theme(axis.title.y = element_blank())

