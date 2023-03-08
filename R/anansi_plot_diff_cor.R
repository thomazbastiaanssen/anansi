#Used for manuscript figure
# library(tidyverse)
# library(patchwork)
#
#
# set.seed(1)
# a = data.frame(Y  = c(1:10, (1:10)*3/2, 10:1) + runif(30, -1, 1),
#                X  = rep(1:10,  3) + runif(30, -1, 1),
#                Phenotype = rep(LETTERS[1:3], each = 10),
#                plot = "Disjointed Association")
#
# b = data.frame(Y  = c(1:10, (1:10)*3/2, (1:10)*5/4) + c(runif(20, -1, 1), runif(10, -7, 7)),
#                X  = rep(1:10,  3)  + c(runif(20, -1, 1), runif(10, -7, 7)),
#                Phenotype = rep(LETTERS[1:3], each = 10),
#                plot = "Emergent Association")
#
# assoc = rbind(a, b) %>%
#   ggplot() +
#   aes(x = X, y = Y, fill = Phenotype, colour = Phenotype) +
#   geom_point(size = 3, shape = 21, colour = "black") +
#   geom_smooth(method = "lm", se = F, linetype = "dashed", show.legend = F) +
#   facet_wrap(~plot) +
#   theme_bw() +
#   ggtitle("Between data sets")
#
# c = data.frame(Y  = rep(c(3, 3, 8), each = 10) + runif(30, -1, 1),
#                Phenotype = rep(LETTERS[1:3], each = 10),
#                plot = "Disjointed Proportionality")
#
# d = data.frame(Y  = rep(3, 30) + c(runif(20, -1, 1), runif(10, -7, 7)),
#                Phenotype = rep(LETTERS[1:3], each = 10),
#                plot = "Emergent Proportionality")
#
#
#
# prop = rbind(c, d) %>%
#   ggplot() +
#   aes(y = Y, x = Phenotype, fill = Phenotype, colour = Phenotype) +
#   geom_point(size = 3, shape = 21, colour = "black") +
#   stat_summary(fun="mean", geom="segment", aes(xend=..x.. - 0.25, yend=..y..), size = 1, show.legend = F) +
#   stat_summary(fun="mean", geom="segment", aes(xend=..x.. + 0.25, yend=..y..), size = 1, show.legend = F) +
#   facet_wrap(~plot) +
#   theme_bw() +
#   ylab("log(Y/X)") +
#   xlab("Phenotype") +
#   ggtitle("Within data set")
#
#
#
# assoc / prop +
#   plot_layout(guides = 'collect') +
#   plot_annotation(tag_levels = "A") &
#   theme(legend.background = element_rect(colour = "black"))
#
