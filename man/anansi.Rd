% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anansi.R
\name{anansi}
\alias{anansi}
\title{Calculate an association network}
\usage{
anansi(
  web,
  method = "pearson",
  groups = NULL,
  metadata = NULL,
  formula = ~1,
  adjust.method = "BH",
  modeltype = "lm",
  resampling = F,
  locality = F,
  reff = NULL,
  verbose = T,
  diff_cor = T,
  ignore_dictionary = F
)
}
\arguments{
\item{web}{An \code{anansiWeb} object, containing two tables with omics data and a dictionary that links them. See \code{weaveWebFromTables()} for how to weave a web.}

\item{method}{Correlation method. \code{method = "pearson"} is the default value. The alternatives to be passed to \code{cor()} are "spearman" and "kendall".}

\item{groups}{A vector of the column names of categorical values in the metadata object to control which groups should be assessed for simple correlations. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.}

\item{metadata}{A vector or data.frame of categorical or continuous value necessary for differential correlations. Typically a state or treatment score. If no argument provided, anansi will let you know and still to regular correlations according to your dictionary.}

\item{formula}{A formula object. Used to assess differential associations.}

\item{adjust.method}{Method to adjust p-values for multiple comparisons. \code{adjust.method = "BH"} is the default value. See \code{p.adjust()} in the base R \code{stats} package.}

\item{modeltype}{A string, either "lm" or "lmer" depending on the type of model that should be ran, or "propr" in the case of within-composition associations..}

\item{resampling}{A boolean. For p-value adjustment. Toggles the resampling of p-values to help reduce the influence of dependence of p-values. Will take more time on large datasets.}

\item{locality}{A boolean. For p-value adjustment. Toggles whether to prefer sampling from p-values from a comparison that shares an x or y feature. In a nutshell, considers the local p-value landscape when more important when correcting for FDR.}

\item{reff}{A categorical vector typically depicting a shared ID between samples. Only for mixed effect models.}

\item{verbose}{A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.}

\item{diff_cor}{A boolean. Toggles whether to compute differential correlations. Default is \code{TRUE}.}

\item{ignore_dictionary}{A boolean. Default is FALSE. If set to TRUE, regular all vs all associations will be tested regardless of the dictionary.}
}
\value{
A list of lists containing correlation coefficients, p-values and q-values for all operations.
}
\description{
This is the main workspider function in the anansi package. It manages the individual functionalities of anansi, including correlation analysis, correlation by group and differential correlation.
}
\examples{
\dontrun{
#Load example data:

data(dictionary)
data(FMT_data)

#Clean and prepare the example data.
#In the example dataset, the metabolites are already cleaned.

KOs   <- floor(FMT_KOs)
KOs   <- apply(KOs,c(1,2),function(x) as.numeric(as.character(x)))
KOs   <- KOs[apply(KOs == 0, 1, sum) <= (ncol(KOs) * 0.90), ]

KOs   <- KOs[row.names(KOs) \%in\% sort(unique(unlist(anansi_dic))),]

#CLR-transform.

KOs.exp = clr_c(KOs)

#Make sure that columns are features and rows are samples.

t1 = t(FMT_metab)
t2 = t(KOs.exp)

#Run anansi pipeline.

web        = weaveWebFromTables(tableY     = t1,
                                tableX     = t2,
                                dictionary = anansi_dic)

anansi_out = anansi(web     = web,
                    method  = "pearson",
                    groups  = FMT_metadata$Legend,
                    adjust.method = "BH",
                    verbose = TRUE)

results    = spinToWide(anansi_output = anansi_out, translate = T,
                        Y_translation = anansi::cpd_translation,
                        X_translation = anansi::KO_translation)

#To recreate the long plot:
library(ggplot2)

anansiLong <- spinToLong(anansi_output = anansi_out, translate = T,
                         Y_translation = anansi::cpd_translation,
                         X_translation = anansi::KO_translation)

#Now it's ready to be plugged into ggplot2, though let's clean up a bit more.

#Only consider interactions where the entire model fits well enough.
anansiLong <- anansiLong[anansiLong$model_full_q.values < 0.1,]



ggplot(data = anansiLong,
       aes(x      = r.values,
           y      = feature_X,
           fill   = type,
           alpha  = model_disjointed_p.values < 0.05)) +

#Make a vertical dashed red line at x = 0
geom_vline(xintercept = 0, linetype = "dashed", colour = "red")+

#Points show  raw correlation coefficients
geom_point(shape = 21, size = 3) +

#facet per compound
ggforce::facet_col(~feature_Y, space = "free", scales = "free_y") +

#fix the scales, labels, theme and other layout
scale_y_discrete(limits = rev, position = "right") +
scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 1/3)) +
scale_fill_manual(values = c("Young yFMT" = "#2166ac",
                             "Aged oFMT"  = "#b2182b",
                             "Aged yFMT"  = "#ef8a62",
                             "All"        = "gray"))+
theme_bw() +
ylab("") +
xlab("Pearson's rho")

#See also ?spinToPlots
}

}
