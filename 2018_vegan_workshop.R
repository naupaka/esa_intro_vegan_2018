# Example code for vegan workshop at ESA 2018 in New Orleans
# August 5, 2018
# Naupaka Zimmerman
# nzimmerman@usfca.edu

# load libraries'
library("vegan")

# load in broken data
BCI_broken <- read.csv("esa_intro_vegan_2018-master/01-intro-basics/data/BCI_small_broken.csv")

summary(BCI_broken)

MLM_otus <- read.csv("esa_intro_vegan_2018-master/01-intro-basics/data/MLM_data_otus.csv",
                     row.names = 1,
                     header = TRUE)
# check names
row.names(MLM_otus)
names(MLM_otus)

head(MLM_otus[1:3])

# sum of the rows
sum_of_rows <- apply(MLM_otus, 1, sum)
sort(sum_of_rows)[1:8]

# sum of the columns
sum_of_columns <- apply(MLM_otus, MARGIN = 2, sum)
sort(sum_of_columns, decreasing = TRUE)[1:8]

# count the number of trees each species occurs in
spec_pres <- apply(MLM_otus > 0, 2, sum)
sort(spec_pres, decreasing = TRUE)[1:8]

# square root transformation
head(MLM_otus[1:3 , 1:3])
spec_sqrt <- sqrt(MLM_otus)
head(spec_sqrt[1:3 , 1:3])

# using decostand in vegan

## total standardization
head(MLM_otus[1:3 , 1:3])
site_total <- decostand(MLM_otus,
                        method = "total",
                        MARGIN = 1) # by row (sites)
head(site_total[1:3 , 1:3])

## total standardization
head(MLM_otus[1:3 , 1:3])
spec_total <- decostand(MLM_otus,
                        method = "total",
                        MARGIN = 2) # by cols (species)
head(spec_total[1:3 , 1:3])

## max standardization
head(MLM_otus[1:3 , 1:3])
spec_max <- decostand(MLM_otus,
                        method = "max",
                        MARGIN = 2) # by cols (species)
head(spec_max[1:3 , 1:3])

## presence absence standardization
head(MLM_otus[1:3 , 1:3])
spec_pa <- decostand(MLM_otus,
                      method = "pa",
                      MARGIN = 2) # by cols (species)
head(spec_pa[1:3 , 1:3])

## Hellinger standardization
# Legendre and Gallagher 2001
head(MLM_otus[1:3 , 1:3])
spec_hellinger <- decostand(MLM_otus,
                     method = "hellinger",
                     MARGIN = 2) # by cols (species)
head(spec_hellinger[1:3 , 1:3])

# Wisconsin standardization
# species to maximum and then sites by totals
spec_wisc <- wisconsin(MLM_otus)
head(spec_wisc[1:3 , 1:3])

# calculating distances
spec_jaccpa <- vegdist(MLM_otus,
                       method = "jaccard",
                       binary = TRUE)

# read in the environmental data
MLM_env <- read.csv("esa_intro_vegan_2018-master/01-intro-basics/data/MLM_data_env.csv",
                    header = TRUE)

rank_elev <- rankindex(MLM_env$elevation_m,
                       MLM_otus,
                       indices = c("bray",
                                   "euclidean",
                                   "manhattan",
                                   "horn"),
                       method = "spearman")
rank_elev

rank_elev <- rankindex(MLM_env$elevation_m,
                       wisconsin(MLM_otus),
                       indices = c("bray",
                                   "euclidean",
                                   "manhattan",
                                   "horn"),
                       method = "spearman")
rank_elev


## Diversity metrics
site_richness <- apply(MLM_otus > 0,
                       1,
                       sum)
site_richness

site_fisher <- fisher.alpha(MLM_otus)
site_fisher

site_shannon <- diversity(MLM_otus,
                          index = "shannon",
                          MARGIN = 1)
site_shannon
sort(site_shannon, decreasing = TRUE)[1:5]

# Rarefaction
MLM_S <- specnumber(MLM_otus)
MLM_S

MLM_raremax <- min(apply(MLM_otus, 1, sum))

MLM_Srare <- rarefy(MLM_otus, MLM_raremax)
plot(x = MLM_S,
     y = MLM_Srare,
     xlab = "Observed number of species",
     ylab = "Rarefied number of species")
abline(0, 1)

# rarecurve
rarecurve(MLM_otus,
          step = 20,
          sample = MLM_raremax,
          col = "blue",
          cex = 0.6)

## Beta diversity
MLM_bray <- vegdist(MLM_otus,
                    method = "euclidean")
MLM_bray_bdisp <- betadisper(MLM_bray,
                             group = as.factor(MLM_env$site_ID))
MLM_bray_bdisp
permutest(MLM_bray_bdisp)
plot(MLM_bray_bdisp)



boxplot(MLM_bray_bdisp, las = 3)

# PERMANOVA
adonis(MLM_otus ~ MLM_env$flow_age * MLM_env$elevation_m)


## Ordinations!
set.seed(42)
MLM_bry_ord <- metaMDS(MLM_otus,
                       distance = "bray",
                       k = 2,
                       trymax = 50)
plot(MLM_bry_ord)

plot(MLM_bry_ord,
     display = "sites")

plot(MLM_bry_ord,
     display = "sites",
     type = "t")

plot(MLM_bry_ord,
     display = "sites")
set.seed(50)
ordipointlabel(MLM_bry_ord,
               display = "sites",
               scaling = 3,
               add = TRUE)

plot(MLM_bry_ord,
     display = "sites",
     cex = 0.5)

colors_vector <- c("green", "blue")
plot(MLM_bry_ord,
     display = "sites",
     type = "n")
points(MLM_bry_ord,
       display = "sites",
       cex = 2,
       pch = 21,
       col = colors_vector[MLM_env$side_of_island],
       bg = colors_vector[MLM_env$side_of_island])
legend("topright",
       legend = levels(MLM_env$side_of_island),
       bty = "n",
       col = colors_vector,
       pch = 21,
       pt.bg = colors_vector)

plot(MLM_bry_ord,
     display = "sites",
     cex = 2)
ordihull(MLM_bry_ord,
         groups = MLM_env$site_ID,
         label = FALSE,
         col = "blue")
ordispider(MLM_bry_ord,
           groups = MLM_env$site_ID,
           label = TRUE)

plot(MLM_bry_ord,
     type = "n",
     display = "sites")
points(MLM_bry_ord,
       display = "sites",
       cex = 2)
ordiellipse(MLM_bry_ord,
            groups = MLM_env$site_ID,
            label = FALSE,
            col = "blue",
            lwd = 3)

plot(MLM_bry_ord,
     type = "n",
     display = "sites")
points(MLM_bry_ord,
       display = "sites",
       cex = 2)
ordisurf(MLM_bry_ord,
         MLM_env$elevation_m,
         add = TRUE)

?ordisurf

MLM_bray_ord_elev_fit <- envfit(MLM_bry_ord ~ elevation_m,
                                data = MLM_env,
                                permutations = 999)
MLM_bray_ord_elev_fit

MLM_bray_ord_rain_fit <- envfit(MLM_bry_ord ~ approx_annual_rainfall_mm,
                                data = MLM_env,
                                permutations = 999)
MLM_bray_ord_rain_fit


png("esa_intro_vegan_2018-master/ordination.png")
plot(MLM_bry_ord,
     type = "n",
     display = "sites")
points(MLM_bry_ord,
       display = "sites",
       cex = 2)
ordisurf(MLM_bry_ord,
         MLM_env$elevation_m,
         add = TRUE)
plot(MLM_bray_ord_elev_fit,
     add = TRUE,
     lwd = 3)
plot(MLM_bray_ord_rain_fit,
     add = TRUE,
     lwd = 3)
dev.off()


