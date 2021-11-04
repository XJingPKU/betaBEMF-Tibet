# beta diversity vs. spatial turnover in ecosystem functions
# XJ
# 2021.10.07

rm(list = ls())

# load libraries
library(tidyverse)
library(TeachingDemos)
library(cowplot)

char2seed("betadiversity")

# load data
df.tibet <- read.csv("./data/data_processing/tibet_distance_matrix.csv")


###########################################################
# Part 1: Mantel analysis
###########################################################

# 1.1 sensitivity analysis for beta-diversity indicies
names(df.tibet)
id.vars <- c("plant.sor", "plant.repl", "plant.rich",
             "bac.sor", "bac.repl", "bac.rich", "bac.BC",
             "fun.sor", "fun.repl", "fun.rich", "fun.BC")
id.mfs <- c("EMF", "EMF.plant", "EMF.soil")
res.beta <- NULL
for (i in id.mfs) {
  for (j in id.vars){
    temp <- df.tibet[c(j, i)]
    names(temp) <- c("y", "x")
    res <- ecodist::mantel(temp$y ~ temp$x, nperm = 9999, mrank = TRUE)
    res <- data.frame(t(res))
    res$idEF <- i
    res$id <- j
    res.beta <- rbind(res.beta, res)
  }
}
res.beta.clean <- res.beta %>% 
  mutate(idEF = factor(idEF,
                       levels = c("EMF", "EMF.plant", "EMF.soil"),
                       labels = c("Overall", "Plant mediated", "Soil mediated")),
         id = gsub("\\.", "_", id),
         Organism = gsub("_.*$", "", id),
         beta_diversity = gsub(".*_", "", id),
         Organism = factor(Organism,
                           levels = c("plant", "bac", "fun"),
                           labels = c("Plant", "Bacteria", "Fungi")),
         beta_diversity = factor(beta_diversity,
                                 levels = c("sor", "repl", "rich", "BC"),
                                 labels = c("Sorensen", "Replacement",
                                            "Richness difference", "Bray-Curtis"))
  ) %>% 
  select(idEF, Organism, beta_diversity, mantelr, pval1, pval2, pval3, llim.2.5., ulim.97.5.)
write.csv(res.beta.clean, "./outs/TableS3.csv")

# 1.2 partial Mantel analysis
# inspect the bivariate associations
p1 <- df.tibet %>% 
  select(plant.sor, bac.BC, fun.BC, geo, clim, pH) %>% 
  pivot_longer(geo:pH, names_to = "abiotic.variable", values_to = "abiotic.value") %>% 
  pivot_longer(plant.sor:fun.BC, names_to = "biotic.variable", values_to = "biotic.value") %>% 
  mutate(biotic.variable = factor(biotic.variable,
                                  levels = c("plant.sor", "bac.BC", "fun.BC"),
                                  labels = c("Plant", "Bacteria", "Fungi")),
         abiotic.variable = factor(abiotic.variable,
                                   levels = c("geo", "clim", "pH"),
                                   labels = c("Geographic distance", 
                                              "Climatic distance",
                                              "pH distance"))) %>% 
  ggplot(aes(abiotic.value, biotic.value)) +
  geom_point(size = 0.2, aes(color = biotic.variable)) +
  geom_smooth(method = "lm", se = FALSE, size = 0.3, color = "black") +
  facet_grid(biotic.variable ~ abiotic.variable, 
             scales = "free_x", switch = "both") +
  labs(y = expression(paste(beta, "-diversity")),
       x = NULL) +
  scale_y_continuous(breaks = c(0, 0.5, 1.0), limits = c(0, 1.1)) +
  scale_x_continuous(breaks = scales::pretty_breaks(3), limits = c(0, NA)) +
  scale_color_manual(values = c("#8FC93E", "#EF7318", "#105E7F")) +
  theme_bw(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = 'none')

# 1.2.1 Geo + Clim + pH --> beta-diversity
partialMantel1 <- function(data) {
  temp <- data
  names(temp) <- c("y", "geo", "clim", "pH")
  geo1 <- ecodist::mantel(y ~ geo, data = temp, nperm = 9999, mrank = TRUE)
  geo2 <- ecodist::mantel(y ~ geo + clim, data = temp, nperm = 9999, mrank = TRUE)
  geo3 <- ecodist::mantel(y ~ geo + pH, data = temp, nperm = 9999, mrank = TRUE)
  geo4 <- ecodist::mantel(y ~ geo + clim + pH, data = temp, nperm = 9999, mrank = TRUE)
  clim1 <- ecodist::mantel(y ~ clim, data = temp, nperm = 9999, mrank = TRUE)
  clim2 <- ecodist::mantel(y ~ clim + geo, data = temp, nperm = 9999, mrank = TRUE)
  clim3 <- ecodist::mantel(y ~ clim + pH, data = temp, nperm = 9999, mrank = TRUE)
  clim4 <- ecodist::mantel(y ~ clim + geo + pH, data = temp, nperm = 9999, mrank = TRUE)
  pH1 <- ecodist::mantel(y ~ pH, data = temp, nperm = 9999, mrank = TRUE)
  pH2 <- ecodist::mantel(y ~ pH + geo, data = temp, nperm = 9999, mrank = TRUE)
  pH3 <- ecodist::mantel(y ~ pH + clim, data = temp, nperm = 9999, mrank = TRUE)
  pH4 <- ecodist::mantel(y ~ pH + geo + clim, data = temp, nperm = 9999, mrank = TRUE)
  res <- rbind(geo1, geo2, geo3, geo4, clim1, clim2, clim3, clim4,
               pH1, pH2, pH3, pH4)
  res <- data.frame(res)
  res$id <- c("geo1", "geo2", "geo3", "geo4", "clim1", "clim2", "clim3", "clim4",
              "pH1", "pH2", "pH3", "pH4")
  return(res)
}

plant.part <- partialMantel1(df.tibet[c("plant.sor", "geo", "clim", "pH")])
bac.part <- partialMantel1(df.tibet[c("bac.BC", "geo", "clim", "pH")])
fun.part <- partialMantel1(df.tibet[c("fun.BC", "geo", "clim", "pH")])


p2 <- rbind(plant.part, bac.part, fun.part) %>% 
  data.frame() %>% 
  mutate(variable = rep(c("Plant", "Bacteria", "Fungi"), each = 12)) %>% 
  mutate(variable = factor(variable,
                           levels = c("Plant", "Bacteria", "Fungi")),
         id = factor(id,
                     levels = c("geo1", "geo2", "geo3", "geo4", "clim1", "clim2", "clim3", "clim4",
                                "pH1", "pH2", "pH3", "pH4"),
                     labels = c("Geo", "Geo | Clim", "Geo | pH", "Geo | Clim+pH",
                                "Clim", "Clim | Geo", "Clim | pH", "Clim | Geo+pH",
                                "pH", "pH | Geo", "pH | Clim", "pH | Geo+Clim"))) %>% 
  filter(!id %in% c("Geo", "Clim", "pH")) %>%
  ggplot(aes(forcats::fct_rev(id), mantelr)) +
  geom_errorbar(aes(ymin = llim.2.5., ymax = ulim.97.5.), 
                color = "black", alpha = 0.66,
                width = 0.15) +
  geom_hline(yintercept = 0, color = "gray", lty = 2) +
  geom_vline(xintercept = 3.5, color = "gray") +
  geom_vline(xintercept = 6.5, color = "gray") +
  geom_point(size = 3, aes(color = variable)) +
  facet_grid(~ variable, switch = "x") +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6)) +
  scale_color_manual(values = c("#8FC93E", "#EF7318", "#105E7F")) +
  labs(x = NULL, y = "Mantel coefficient (rho)") +
  coord_flip() +
  theme_bw(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

plot_grid(p1, p2, ncol = 1, labels = c("a)", "b)"), rel_heights = c(1/2, 1/2))
ggsave("./outs/tibet_beta-diversity.pdf", width = 6.5, height = 9.5)  # Figure 1

# Simple Mantel tests
rbind(plant.part, bac.part, fun.part) %>% 
  data.frame() %>% 
  mutate(variable = rep(c("Plant", "Bacteria", "Fungi"), each = 12)) %>% 
  mutate(variable = factor(variable,
                           levels = c("Plant", "Bacteria", "Fungi")),
         id = factor(id,
                     levels = c("geo1", "geo2", "geo3", "geo4", "clim1", "clim2", "clim3", "clim4",
                                "pH1", "pH2", "pH3", "pH4"),
                     labels = c("Geo", "Geo | Clim", "Geo | pH", "Geo | Clim+pH",
                                "Clim", "Clim | Geo", "Clim | pH", "Clim | Geo+pH",
                                "pH", "pH | Geo", "pH | Clim", "pH | Geo+Clim"))) %>% 
  filter(id %in% c("Geo", "Clim", "pH")) 


# 1.2.2 EMF vs. abiotic factors
p5 <- df.tibet %>% 
  select(EMF, EMF.plant, EMF.soil, geo, clim, pH) %>% 
  pivot_longer(geo:pH, names_to = "abiotic.variable", values_to = "abiotic.value") %>% 
  pivot_longer(EMF:EMF.soil, names_to = "biotic.variable", values_to = "biotic.value") %>% 
  mutate(biotic.variable = factor(biotic.variable,
                                  levels = c("EMF", "EMF.plant", "EMF.soil"),
                                  labels = c("Overall", "Plant-mediated", "Soil-mediated")),
         abiotic.variable = factor(abiotic.variable,
                                   levels = c("geo", "clim", "pH"),
                                   labels = c("Geographic distance", 
                                              "Climatic distance",
                                              "pH distance"))) %>% 
  ggplot(aes(abiotic.value, biotic.value)) +
  geom_point(size = 0.2, color = "grey70") +
  geom_smooth(method = "lm", se = FALSE, size = 0.3, color = "black") +
  facet_grid(biotic.variable ~ abiotic.variable, scales = "free_x", switch = "both") +
  labs(y = "Turnover of nutrient stocks",
       x = NULL) +
  lims(y = c(0, 15)) +
  scale_x_continuous(breaks = scales::pretty_breaks(3), limits = c(0, NA)) +
  theme_bw(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = 'none')

EMF.part <- partialMantel1(df.tibet[c("EMF", "geo", "clim", "pH")])
EMF.plant.part <- partialMantel1(df.tibet[c("EMF.plant", "geo", "clim", "pH")])
EMF.soil.part <- partialMantel1(df.tibet[c("EMF.soil", "geo", "clim", "pH")])

p6 <- rbind(EMF.part, EMF.plant.part, EMF.soil.part) %>% 
  data.frame() %>% 
  mutate(variable = rep(c("EMF", "EMF.plant", "EMF.soil"), each = 12)) %>% 
  mutate(variable = factor(variable,
                           levels = c("EMF", "EMF.plant", "EMF.soil"),
                           labels = c("Overall", "Plant-mediated", "Soil-mediated")),
         id = factor(id,
                     levels = c("geo1", "geo2", "geo3", "geo4", "clim1", "clim2", "clim3", "clim4",
                                "pH1", "pH2", "pH3", "pH4"),
                     labels = c("Geo", "Geo | Clim", "Geo | pH", "Geo | Clim+pH",
                                "Clim", "Clim | Geo", "Clim | pH", "Clim | Geo+pH",
                                "pH", "pH | Geo", "pH | Clim", "pH | Geo+Clim"))) %>% 
  filter(!id %in% c("Geo", "Clim", "pH")) %>%
  ggplot(aes(forcats::fct_rev(id), mantelr)) +
  geom_errorbar(aes(ymin = llim.2.5., ymax = ulim.97.5.), 
                color = "black", alpha = 0.66,
                width = 0.15) +
  geom_hline(yintercept = 0, color = "gray", lty = 2) +
  geom_vline(xintercept = 3.5, color = "gray") +
  geom_vline(xintercept = 6.5, color = "gray") +
  geom_point(size = 3, color = "grey70") +
  facet_grid(~ variable, switch = "x") +
  scale_y_continuous(breaks = scales::pretty_breaks(3), limits = c(-0.1, NA)) +
  labs(x = NULL, y = "Mantel coefficient (rho)") +
  coord_flip() +
  theme_bw(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

plot_grid(p5, p6, ncol = 1, labels = c("a)", "b)"), rel_heights = c(1/2, 1/2))
ggsave("./outs/tibet_beta-EMF.pdf", width = 6.5, height = 9.5)  # Figure 2

rbind(EMF.part, EMF.plant.part, EMF.soil.part) %>% 
  data.frame() %>% 
  mutate(variable = rep(c("EMF", "EMF.plant", "EMF.soil"), each = 12)) %>% 
  mutate(variable = factor(variable,
                           levels = c("EMF", "EMF.plant", "EMF.soil"),
                           labels = c("Overall", "Plant-mediated", "Soil-mediated")),
         id = factor(id,
                     levels = c("geo1", "geo2", "geo3", "geo4", "clim1", "clim2", "clim3", "clim4",
                                "pH1", "pH2", "pH3", "pH4"),
                     labels = c("Geo", "Geo | Clim", "Geo | pH", "Geo | Clim+pH",
                                "Clim", "Clim | Geo", "Clim | pH", "Clim | Geo+pH",
                                "pH", "pH | Geo", "pH | Clim", "pH | Geo+Clim"))) %>% 
  filter(id %in% c("Geo", "Clim", "pH"))


# 1.2.3 beta-diversity --> beta-multifunctionality
# EMF vs. Sorensen
partialMantel2 <- function(data) {
  temp <- data
  names(temp) <- c("y", "x", "geo", "clim", "pH")
  nul <- ecodist::mantel(y ~ x, data = temp, nperm = 9999, mrank = TRUE)
  geo <- ecodist::mantel(y ~ x + geo, data = temp, nperm = 9999, mrank = TRUE)
  clim <- ecodist::mantel(y ~ x + clim, data = temp, nperm = 9999, mrank = TRUE)
  pH <- ecodist::mantel(y ~ x + pH, data = temp, nperm = 9999, mrank = TRUE)
  geo.clim <- ecodist::mantel(y ~ x + geo + clim, data = temp, nperm = 9999, mrank = TRUE)
  geo.pH <- ecodist::mantel(y ~ x + geo + pH, data = temp, nperm = 9999, mrank = TRUE)
  clim.pH <- ecodist::mantel(y ~ x + clim + pH, data = temp, nperm = 9999, mrank = TRUE)
  geo.clim.pH <- ecodist::mantel(y ~ x + geo + clim + pH, data = temp, nperm = 9999, mrank = TRUE)
  res <- rbind(nul, geo, clim, pH, geo.clim, geo.pH, clim.pH, geo.clim.pH)
  res <- data.frame(res)
  res$id <- c("nul", "geo", "clim", "pH", 
              "geo.clim", "geo.pH", "clim.pH", "geo.clim.pH")
  return(res)
}

plant.part <- partialMantel2(df.tibet[c("EMF", "plant.sor", "geo", "clim", "pH")])
bac.part <- partialMantel2(df.tibet[c("EMF", "bac.BC", "geo", "clim", "pH")])
fun.part <- partialMantel2(df.tibet[c("EMF", "fun.BC", "geo", "clim", "pH")])

plant.part2 <- partialMantel2(df.tibet[c("EMF.plant", "plant.sor", "geo", "clim", "pH")])
bac.part2 <- partialMantel2(df.tibet[c("EMF.plant", "bac.BC", "geo", "clim", "pH")])
fun.part2 <- partialMantel2(df.tibet[c("EMF.plant", "fun.BC", "geo", "clim", "pH")])

plant.part3 <- partialMantel2(df.tibet[c("EMF.soil", "plant.sor", "geo", "clim", "pH")])
bac.part3 <- partialMantel2(df.tibet[c("EMF.soil", "bac.BC", "geo", "clim", "pH")])
fun.part3 <- partialMantel2(df.tibet[c("EMF.soil", "fun.BC", "geo", "clim", "pH")])

p3 <- df.tibet %>% 
  select(EMF, EMF.plant, EMF.soil, plant.sor, bac.BC, fun.BC) %>% 
  pivot_longer(plant.sor:fun.BC, names_to = "div_variable", 
               values_to = "div_value") %>% 
  pivot_longer(EMF:EMF.soil, names_to = "EMF.variable", 
               values_to = "EMF.value") %>% 
  mutate(div_variable = factor(div_variable,
                               levels = c("plant.sor", "bac.BC", "fun.BC"),
                               labels = c("Plant", "Bacteria", "Fungi"))) %>% 
  mutate(EMF.variable = factor(EMF.variable,
                               levels = c("EMF", "EMF.plant", "EMF.soil"),
                               labels = c("Overall", "Plant-mediated", "Soil-mediated"))) %>% 
  ggplot(aes(div_value, EMF.value)) +
  geom_point(size = 0.2, aes(color = div_variable)) +
  geom_smooth(method = "lm", se = FALSE, size = 0.3, color = "black") +
  facet_grid(EMF.variable ~ div_variable, switch = "both") +
  scale_color_manual(values = c("#8FC93E", "#EF7318", "#105E7F")) +
  lims(y = c(0, 15)) +
  scale_x_continuous(breaks = scales::pretty_breaks(3), limits = c(0, NA)) +
  labs(x = expression(paste(beta, "-diversity")),
       y = "Turnover of nutrient stocks") +
  theme_bw(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = 'none')

p4 <- rbind(plant.part, bac.part, fun.part,
            plant.part2, bac.part2, fun.part2,
            plant.part3, bac.part3, fun.part3) %>% 
  data.frame() %>% 
  mutate(div_variable = rep(rep(c("Plant", "Bacteria", "Fungi"), each = 8), 3)) %>% 
  mutate(EMF_variable = rep(c("Overall", "Plant-mediated", "Soil-mediated"), each = 24)) %>% 
  mutate(div_variable = factor(div_variable,
                               levels = c("Plant", "Bacteria", "Fungi")),
         EMF_variable = factor(EMF_variable,
                               levels = c("Overall", "Plant-mediated", "Soil-mediated")),
         id = factor(id,
                     levels = c("nul", "geo", "clim", "pH", 
                                "geo.clim", "geo.pH", "clim.pH", "geo.clim.pH"),
                     labels = c("X", "X | Geo", "X | Clim", "X | pH", 
                                "X | Geo+Clim", "X | Geo+pH", "X | Clim+pH", 
                                "X | Geo+Clim+pH"))) %>% 
  filter(!id %in% c("X")) %>% 
  ggplot(aes(forcats::fct_rev(id), mantelr)) +
  geom_errorbar(aes(ymin = llim.2.5., ymax = ulim.97.5.), 
                color = "black", alpha = 0.66,
                width = 0.15) +
  geom_hline(yintercept = 0, color = "gray", lty = 2) +
  geom_point(size = 3.0, aes(color = div_variable)) +
  scale_color_manual(values = c("#8FC93E", "#EF7318", "#105E7F")) +
  scale_y_continuous(breaks = scales::pretty_breaks(3), limits = c(-0.1, NA)) +
  facet_grid(EMF_variable ~ div_variable, switch = "both") +
  labs(x = NULL, y = "Mantel coefficient (rho)") +
  coord_flip() +
  theme_bw(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")

plot_grid(p3, p4, ncol = 1, labels = c("a)", "b)"), rel_heights = c(4.5/10, 5.5/10))
ggsave("./outs/tibet_beta-BEMF.pdf", width = 6.5, height = 9.5)  # Figure 3
# Simple Mantel tests
rbind(plant.part, bac.part, fun.part,
      plant.part2, bac.part2, fun.part2,
      plant.part3, bac.part3, fun.part3) %>% 
  data.frame() %>% 
  mutate(div_variable = rep(rep(c("Plant", "Bacteria", "Fungi"), each = 8), 3)) %>% 
  mutate(EMF_variable = rep(c("Overall", "Plant-mediated", "Soil-mediated"), each = 24)) %>% 
  mutate(div_variable = factor(div_variable,
                               levels = c("Plant", "Bacteria", "Fungi")),
         EMF_variable = factor(EMF_variable,
                               levels = c("Overall", "Plant-mediated", "Soil-mediated")),
         id = factor(id,
                     levels = c("nul", "geo", "clim", "pH", 
                                "geo.clim", "geo.pH", "clim.pH", "geo.clim.pH"),
                     labels = c("X", "X | Geo", "X | Clim", "X | pH", 
                                "X | Geo+Clim", "X | Geo+pH", "X | Clim+pH", 
                                "X | Geo+Clim+pH"))) %>% 
  filter(id %in% c("X"))


###########################################################
# Part 2: Structural equation model
###########################################################

library(lavaan)
library(lavaanPlot)

sem.data <- df.tibet %>%
  select(EMF, EMF.plant, EMF.soil, geo, clim, pH, plant.sor, bac.sor, bac.BC, fun.sor, fun.BC) %>%
  scale() %>%
  data.frame()

mod_formula_A <- '
# regressions
EMF.plant ~ clim + bac.BC + plant.sor + fun.BC
EMF.soil ~ clim + pH + bac.BC + plant.sor
plant.sor ~ geo + clim
bac.BC ~ clim + pH + plant.sor
fun.BC ~ geo + clim + plant.sor
pH ~ geo + clim
clim ~ geo
# residual correlations
EMF.plant ~~ EMF.soil
fun.BC ~~ bac.BC
'
fitA <- sem(mod_formula_A, data = sem.data, se = "boot", bootstrap = 9999)
summary(fitA, standardized = T, rsq = T, fit.measures = F)
fitMeasures(fitA, c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr"))
lavaanPlot(model = fitA,
           node_options = list(shape = "box", fontname = "Arial"),
           edge_options = list(color = "grey"),
           coefs = TRUE, stars = TRUE)

# An alternative model
mod_formula_B <- '
# regressions
EMF.plant ~ clim
EMF.soil ~ clim + pH
plant.sor ~ geo + clim + EMF.plant + EMF.soil
bac.BC ~ clim + pH + plant.sor + EMF.plant + EMF.soil
fun.BC ~ geo + clim + plant.sor
pH ~ geo + clim
clim ~ geo
# residual correlations
EMF.plant ~~ EMF.soil
fun.BC ~~ bac.BC
'
fitB <- sem(mod_formula_B, data = sem.data, se = "boot", bootstrap = 9999)
summary(fitB, standardized = T, rsq = T, fit.measures = F)
fitMeasures(fitB, c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr"))
lavaanPlot(model = fitB,
           node_options = list(shape = "box", fontname = "Arial"),
           edge_options = list(color = "grey"),
           coefs = TRUE, stars = TRUE)
AIC(fitA, fitB)
BIC(fitA, fitB)

###########################################################
# Part 3: Linear models with permutation tests
###########################################################
library(lmPerm)

# "EMF", "EMF.plant", "EMF.soil"
res.lmp2.repl <- NULL
for (j in c("EMF", "EMF.plant", "EMF.soil")) {
  print(j)
  res.lmp1.repl <- NULL
  for (i in seq(40, 1000, 20)) {
    print(i)
    res.repl <- NULL
    temp <- df.tibet %>%
      select(j, geo, clim, pH,
             plant.repl, plant.rich,
             bac.repl, bac.rich, fun.repl, fun.rich) %>%
      filter(geo <= i) %>%
      data.frame()
    temp <- apply(temp, 2, scale)
    temp <- data.frame(temp)
    names(temp)[1] <- "y"
    lmp.fit <- lmp(y ~ geo + clim + pH +
                     plant.repl + 
                     bac.repl + 
                     fun.repl,
                   center = FALSE,
                   data = temp)
    res.repl <- summary(lmp.fit)$coef %>% 
      data.frame() %>%
      slice(2L:10L) %>% 
      mutate(geoD = i,
             yVar = j,
             xVar = rownames(.)) %>% 
      select(Estimate, Iter, "Pr.Prob.", geoD, yVar, xVar)
    
    res.lmp1.repl <- rbind(res.lmp1.repl, res.repl)
  }
  res.lmp2.repl <- rbind(res.lmp2.repl, res.lmp1.repl)
}

res.lmp2.rich <- NULL
for (j in c("EMF", "EMF.plant", "EMF.soil")) {
  print(j)
  res.lmp1.rich <- NULL
  for (i in seq(40, 1000, 20)) {
    print(i)
    res.rich <- NULL
    temp <- df.tibet %>%
      select(j, geo, clim, pH,
             plant.rich,
             bac.rich, fun.rich) %>%
      filter(geo <= i) %>%
      data.frame()
    temp <- apply(temp, 2, scale)
    temp <- data.frame(temp)
    names(temp)[1] <- "y"
    lmp.fit <- lmp(y ~ geo + clim + pH +
                     plant.rich + 
                     bac.rich + 
                     fun.rich,
                   center = FALSE,
                   data = temp)
    res.rich <- summary(lmp.fit)$coef %>% 
      data.frame() %>%
      slice(2L:10L) %>% 
      mutate(geoD = i,
             yVar = j,
             xVar = rownames(.)) %>% 
      select(Estimate, Iter, "Pr.Prob.", geoD, yVar, xVar)
    
    res.lmp1.rich <- rbind(res.lmp1.rich, res.rich)
  }
  res.lmp2.rich <- rbind(res.lmp2.rich, res.lmp1.rich)
}

res.lmp2.repl <- res.lmp2.repl %>% 
  mutate(betas = "Replacement")

res.lmp2.rich <- res.lmp2.rich %>% 
  mutate(betas = "Richness difference")


rbind(res.lmp2.repl, res.lmp2.rich) %>% 
  filter(!xVar %in% c("geo", "clim", "pH")) %>% 
  droplevels() %>% 
  mutate(yVar = factor(yVar, 
                       levels = c("EMF", "EMF.plant", "EMF.soil"),
                       labels = c("Overall", "Plant-mediated", "Soil-mediated")),
         xVar = factor(xVar,
                       levels = c("plant.repl", "bac.repl", "fun.repl",
                                  "plant.rich", "bac.rich", "fun.rich"),
                       labels = c("Plant", "Bacteria", "Fungi",
                                  "Plant", "Bacteria", "Fungi")),
         betas = factor(betas,
                        levels = c("Replacement", "Richness difference"))) %>% 
  ggplot(aes(geoD, Estimate, color = xVar)) +
  geom_line() +
  facet_grid(yVar~.) +
  scale_color_manual(values = c("#8FC93E", "#EF7318", "#105E7F")) +
  scale_x_continuous(breaks = seq(0, 1000, 250)) +
  geom_hline(yintercept = 0, lty = 2, color = "gray") +
  facet_grid(yVar ~betas, switch = "y") +
  labs(y = "Standardized regression coefficients",
       x = "Geographic distance (km)") +
  ylim(c(-0.2, 0.6)) +
  theme_bw(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")
ggsave("./outs/tibet_lmp_betadiversity.pdf",
       width = 182, height = 157, units = "mm")


rbind(res.lmp2.repl, res.lmp2.rich) %>% 
  filter(xVar %in% c("geo", "clim", "pH")) %>% 
  droplevels() %>% 
  mutate(yVar = factor(yVar, 
                       levels = c("EMF", "EMF.plant", "EMF.soil"),
                       labels = c("Overall", "Plant-mediated", "Soil-mediated")),
         xVar = factor(xVar,
                       levels = c("geo", "clim", "pH"),
                       labels = c("Geographic distance",
                                  "Climatic distance",
                                  "pH distance"))) %>% 
  ggplot(aes(geoD, Estimate, color = xVar)) +
  geom_line() +
  facet_grid(yVar~.) +
  scale_color_manual(values = c("#BFEB4E", "#61A6FF", "#E23726")) +
  scale_x_continuous(breaks = seq(0, 1000, 250)) +
  geom_hline(yintercept = 0, lty = 2, color = "gray") +
  facet_grid(yVar ~betas, switch = "y") +
  labs(y = "Standardized regression coefficients",
       x = "Geographic distance (km)") +
  ylim(c(-0.2, 0.6)) +
  theme_bw(base_size = 13.5) +
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")

ggsave("./outs/tibet_lmp_abiotic_factors.pdf",
       width = 182, height = 157, units = "mm")

res.lmp2.repl %>% select(xVar, geoD, yVar, Estimate, Pr.Prob.) %>% 
  mutate(xVar = factor(xVar,
                       levels = c("plant.repl", "bac.repl", "fun.repl", 
                                  "geo", "clim", "pH"),
                       labels = c("Plant", "Bacteria", "Fungi", 
                                  "Geographic distance", "Climatic distance", "pH distance"))) %>% 
  write.csv(., "./outs/tibet_lmp_replacement.csv")
res.lmp2.rich %>% select(xVar, geoD, yVar, Estimate, Pr.Prob.) %>% 
  mutate(xVar = factor(xVar,
                       levels = c("plant.rich", "bac.rich", "fun.rich", 
                                  "geo", "clim", "pH"),
                       labels = c("Plant", "Bacteria", "Fungi", 
                                  "Geographic distance", "Climatic distance", "pH distance"))) %>% 
  write.csv(., "./outs/tibet_lmp_richness_diff.csv")


###########################################################
# Part 4: Inspect data distribution of replacement and richness difference
###########################################################
library(ggdist)

# slabinterval
# https://mjskay.github.io/ggdist/articles/slabinterval.html
p.replacement <- df.tibet %>% 
  select(plant.repl, bac.repl, fun.repl) %>% 
  pivot_longer(cols = plant.repl:fun.repl,
               names_to = "metric",
               values_to = "value") %>% 
  mutate(metric = factor(metric,
                         levels = c("plant.repl", "bac.repl", "fun.repl"),
                         labels = c("Plant", "Bacteria", "Fungi"))) %>% 
  ggplot(aes(value, fill = metric, color = metric)) +
  lims(x = c(0, 1)) +
  theme_bw(base_size = 13.5) +
  labs(x = "Replacement",
       y = "Probability density") +
  scale_fill_manual(values = c("#8FC93E", "#EF7318", "#105E7F")) +
  scale_color_manual(values = c("#8FC93E", "#EF7318", "#105E7F")) +
  stat_slab(alpha = .6, linetype = 0, slab_type = "pdf") +
  stat_pointinterval(.width = c(0.66, .95),
                     point_interval = median_qi,
                     position = position_dodge(width = .4, preserve = "single")) +
  theme(legend.position = "none",
        panel.grid = element_blank())

p.richdiff <- df.tibet %>% 
  select(plant.rich, bac.rich, fun.rich) %>% 
  pivot_longer(cols = plant.rich:fun.rich,
               names_to = "metric",
               values_to = "value") %>% 
  mutate(metric = factor(metric,
                         levels = c("plant.rich", "bac.rich", "fun.rich"),
                         labels = c("Plant", "Bacteria", "Fungi"))) %>% 
  ggplot(aes(value, fill = metric, color = metric)) +
  lims(x = c(0, 1)) +
  theme_bw(base_size = 13.5) +
  labs(x = "Richness difference",
       y = "Probability density") +
  scale_fill_manual(values = c("#8FC93E", "#EF7318", "#105E7F")) +
  scale_color_manual(values = c("#8FC93E", "#EF7318", "#105E7F")) +
  stat_slab(alpha = .6, linetype = 0, slab_type = "pdf") +
  stat_pointinterval(.width = c(0.66, .95),
                     point_interval = median_qi,
                     position = position_dodge(width = .4, preserve = "single")) +
  theme(legend.title = element_blank(),
        legend.position = c(0.7, 0.75),
        panel.grid = element_blank())

plot_grid(p.replacement, p.richdiff,
          nrow = 1, labels = c("a)", "b)"))
ggsave("./outs/tibet_beta_diversity_density_plot.pdf",
       width = 182, height = 90, units = "mm")

df.tibet %>% 
  select(plant.rich, bac.rich, fun.rich, plant.repl, bac.repl, fun.repl) %>% 
  pivot_longer(cols = plant.rich:fun.repl,
               names_to = "metric",
               values_to = "value") %>% 
  group_by(metric) %>% 
  summarise(median_qi(value, .width = c(0.66, 0.95)), .groups = "drop")


###########################################################
#                End of the Script                        #
###########################################################