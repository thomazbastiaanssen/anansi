formula = ~ group * magnitude


# which.ind   <- apply(web@dictionary, 2, which)
#
# split.index <- split.data.frame(which.ind, f = which.ind[,2])
library(anansi)
library(tidyverse)
library(deleuze)

#load anansi dictionary and complementary human-readable names for KEGG compounds and orthologues
data(dictionary)

#load example data + metadata from FMT Aging study
data(FMT_data)

KOs.exp <- dclr(FMT_KOs)

t1 <- t(FMT_metab); t2<- t(KOs.exp)

web <- weaveWebFromTables(tableY = t1, tableX = t2, dictionary = anansi_dic)

web@dictionary <- matrix(TRUE, nrow = 4, ncol = 2000)
web@dictionary <- replicate(2000, sample(c(T, T, T, T, F), size = 4, replace = FALSE))
web@dictionary <- matrix(TRUE, nrow = 4, ncol = 2000)

web@tableY     <- matrix(rnorm(80*4),    ncol = 4)
web@tableX     <- matrix(rnorm(80*2000), ncol = 2000)

metadata = data.frame(group  = sample(letters[1:3],
                                      size = nrow(web@tableX), replace = TRUE),
                      magnitude = runif(n = nrow(web@tableX)))


formula <- ~ group * magnitude


a <- anansiDiffCor2(web = web, metadata = metadata,
               formula = formula, reff = NULL,
               modeltype = NULL, verbose = T)


# a$modelfit$full@p.values %>%                 hist(breaks = 20)
a$disjointed$group@p.values %>%              hist(breaks = 20)
a$disjointed$magnitude@p.values %>%          hist(breaks = 20)
a$disjointed$`group:magnitude`@p.values %>%  hist(breaks = 20)
a$emergent$group@p.values %>%                hist(breaks = 20)
a$emergent$magnitude@p.values %>%            hist(breaks = 20)
a$emergent$`group:magnitude`@p.values %>%    hist(breaks = 20)


hist(sort(a$emergent$group@p.values) + rev(sort(a$disjointed$group@p.values)))
hist(sort(a$emergent$magnitude@p.values) + rev(sort(a$disjointed$magnitude@p.values)))
hist(sort(a$emergent$`group:magnitude`@p.values) + rev(sort(a$disjointed$`group:magnitude`@p.values)))



a$modelfit$full@estimates %>%                hist(breaks = 20)
a$disjointed$group@estimates %>%             hist(breaks = 20)
a$disjointed$magnitude@estimates %>%         hist(breaks = 20)
a$disjointed$`group:magnitude`@estimates %>% hist(breaks = 20)
a$emergent$group@estimates %>%               hist(breaks = 20)
a$emergent$magnitude@estimates %>%           hist(breaks = 20)
a$emergent$`group:magnitude`@estimates %>%   hist(breaks = 20)

plot(a$modelfit$full@p.values,
a$modelfit$full@estimates)

plot(a$disjointed$group@p.values,
a$disjointed$group@estimates)

plot(a$disjointed$magnitude@p.values,
a$disjointed$magnitude@estimates)

plot(a$disjointed$`group:magnitude`@p.values,
a$disjointed$`group:magnitude`@estimates)

plot(a$emergent$group@p.values,
a$emergent$group@estimates)

plot(a$emergent$magnitude@p.values,
a$emergent$magnitude@estimates)

plot(a$emergent$`group:magnitude`@p.values,
a$emergent$`group:magnitude`@estimates)




a = anansiDiffCor2(web = web, metadata = metadata,
                   formula = formula, reff = NULL,
                   modeltype = NULL, verbose = T)

b = anansiDiffCor(web = web, metadata = metadata,
                    formula = formula, reff = NULL,
                    modeltype = "lm", verbose = T)
do.call(rbind, list(
  rbind(a$emergent$magnitude@estimates %>% t() %>%  data.frame() %>% add_column(model = "a", id = 1:ncol(a$modelfit$full@estimates)),
        b$emergent$magnitude@estimates %>% t() %>%  data.frame() %>% add_column(model = "b", id = 1:ncol(a$modelfit$full@estimates))
  ) %>% add_column(param = "emergent"),
  rbind(a$disjointed$magnitude@estimates %>% t() %>%  data.frame() %>% add_column(model = "a", id = 1:ncol(a$modelfit$full@estimates)),
        b$disjointed$magnitude@estimates %>% t() %>%  data.frame() %>% add_column(model = "b", id = 1:ncol(a$modelfit$full@estimates))
  ) %>% add_column(param = "disjointed"),
  rbind(a$modelfit$full@estimates %>% t() %>%  data.frame() %>% add_column(model = "a", id = 1:ncol(a$modelfit$full@estimates)),
        b$modelfit$full@estimates %>% t() %>%  data.frame() %>% add_column(model = "b", id = 1:ncol(a$modelfit$full@estimates))
  ) %>% add_column(param = "full")

)) %>%
  pivot_longer(!c(model, id, param)) %>%

  #filter(param != "emergent") %>%

  ggplot() +
  aes(x = value, fill = model) +
  geom_histogram(bins = 40) +
  #scale_x_continuous(limits = c(0, 1))+
    facet_grid(param+model~name, scales = "free_x") +
    theme_bw()

do.call(rbind, list(
  rbind(a$emergent$magnitude@p.values %>% t() %>%  data.frame() %>% add_column(model = "a", id = 1:ncol(a$modelfit$full@estimates)),
        b$emergent$magnitude@p.values %>% t() %>%  data.frame() %>% add_column(model = "b", id = 1:ncol(a$modelfit$full@estimates))
  ) %>% add_column(param = "emergent"),
  rbind(a$disjointed$magnitude@p.values %>% t() %>%  data.frame() %>% add_column(model = "a", id = 1:ncol(a$modelfit$full@estimates)),
        b$disjointed$magnitude@p.values %>% t() %>%  data.frame() %>% add_column(model = "b", id = 1:ncol(a$modelfit$full@estimates))
  ) %>% add_column(param = "disjointed"),
  rbind(a$modelfit$full@p.values %>% t() %>%  data.frame() %>% add_column(model = "a", id = 1:ncol(a$modelfit$full@estimates)),
        b$modelfit$full@p.values %>% t() %>%  data.frame() %>% add_column(model = "b", id = 1:ncol(a$modelfit$full@estimates))
  ) %>% add_column(param = "full")

)) %>%
  pivot_longer(!c(model, id, param)) %>%

  #filter(param != "emergent") %>%

  ggplot() +
  aes(x = value, fill = model) +
  geom_histogram(bins = 40) +
  scale_x_continuous(limits = c(0, 1))+
  facet_grid(param+model~name, scales = "free") +
  theme_bw() + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

anansiDiffCor(web = web, metadata = metadata,
              formula = formula, reff = NULL,
              modeltype = "lm", verbose = T)



sum((lm(web@tableY[,2] ~ 1, data = lm.metadata) %>% .$residuals)^2)

anova(lm(web@tableY[,2] ~ 1, data = lm.metadata),
      lm(web@tableY[,2] ~ K22373 + group + magnitude + K22373:group + K22373:magnitude, data = cbind(web@tableX, lm.metadata)))


#compute f-values from p-vals & r^2

microbenchmark::microbenchmark(naive = replicate(
summary(lm(web@tableY[,sample(size = 1,
                        x = 1:ncol(web@tableY))] ~ web@tableX[,sample(size = 1,
                                                                      x = 1:ncol(web@tableX))])),
  n = nrow(which.ind)),
  new = lapply(split.index, FUN = function(y) {

    #Compute P-value and R-squared through F-statistic
    p_Rsq_from_resids(y = y)}
  ), times = 100)

waldo::compare(f_lm() %>% lapply(., unname), f_ls() %>% lapply(., unname))
#Can we compute goodness of fit p-values without ever running summary, and hopefully with very few anovas?


   #summary_to_PR
  #if(nrow(y) == 1)
  #{print(summary(fit)); return(summary(fit)$fstatistic)}
  #else {return(lapply(summary(fit), function(x) x$fstatistic))}
  #fit


#Strategy to subset models

anansi_ls_new <- function() anansiDiffCor2(web = web, metadata = metadata,
                                formula = formula, reff = NULL,
                                modeltype = "lm", verbose = T)

anansi_lm_cur <- function() anansiDiffCor(web = web, metadata = metadata,
                               formula = formula, reff = NULL,
                               modeltype = "lm", verbose = T)



microbenchmark::microbenchmark(

  anansi_ls_new(),
  anansi_lm_cur(),

  times = 10

)




microbenchmark::microbenchmark(

    anansi1 = anansi(web = web,          #Generated above
           method   = "pearson",    #Define the type of correlation used
           formula  = formula,     #Compare associations between treatments
           metadata = metadata, #With data referred to in the formula as column
           adjust.method = "BH",    #Apply the Benjamini-Hochberg procedure for FDR
           verbose  = T             #To let you know what's happening
    ),
    anansi2 = anansi2(web      = web,          #Generated above
                     method   = "pearson",    #Define the type of correlation used
                     formula  = formula,     #Compare associations between treatments
                     metadata = metadata, #With data referred to in the formula as column
                     adjust.method = "BH",    #Apply the Benjamini-Hochberg procedure for FDR
                     verbose  = T             #To let you know what's happening
    ),

  times = 10

)




