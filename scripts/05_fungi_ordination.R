# Prepare environment ####
rm(list = ls()) # Reset R`s brain

library(vegan)
library(tidyverse)

# Load data, using species-sample matrix (in wide format) 
data_spe <- read.csv("./data/sh_abundance.csv") %>% 
  as.data.frame()

# Assign the first column as row names
rownames(data_spe) <- data_spe[, 1]
# Remove the first column
data_spe <- data_spe[, -1]

# Transpose data to make it suitable for {vegan}
data_spe <- t(data_spe)
data_spe

# Identify row names where no SH in a sample
rows_to_remove <- rownames(data_spe)[rowSums(data_spe) == 0]

# Remove these rows from data_spe
data_spe <- data_spe[!(rownames(data_spe) %in% rows_to_remove), , drop = FALSE]


# Load environmental data
env_data <- read_csv("./data/sampled_soil_envir.csv") %>% 
  select(Site_nm, Ecosystem = Ecsystm, clay, organic, pH, precipitation, temperature) %>% 
  tibble::column_to_rownames(var = "Site_nm") %>% 
  # Replace NAs with the mean values by columns
  mutate(across(c(clay, organic, pH, precipitation, temperature), 
                ~replace_na(., mean(., na.rm=TRUE))
                )
         )

# Remove the same rows from env_data that we removed from data_spe earlier
env_data <- env_data[!(rownames(env_data) %in% rows_to_remove), , drop = FALSE]
env_data


# Basic stats for communities ####
## Analysis of Similarities ####
?vegan::anosim()

str(env_data)
spe.dist <- vegdist(data_spe, method = "jaccard")
spe.ano <- with(env_data, anosim(spe.dist, pH))
summary(spe.ano)
plot(spe.ano)



## PERMANOVA ####
# Read more: https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/permanova/
?vegan::adonis2

adonis2(
  data_spe ~ precipitation,  # formula
  data = env_data     # environmental data
)



# Transformations and distances ####

# Identify row names where less than 10 reads in a sample
rows_to_remove <- rownames(data_spe)[rowSums(data_spe) <= 10]

# Remove these rows from data_spe
data_spe <- data_spe[!(rownames(data_spe) %in% rows_to_remove), , drop = FALSE]

# Remove the same rows from env_data that we removed from data_spe earlier
env_data <- env_data[!(rownames(env_data) %in% rows_to_remove), , drop = FALSE]

reads <- rowSums(data_spe)
View(reads)

hist(reads, breaks = 3000)
summary(reads)

hist(log1p(reads), breaks = 30)
hist(sqrt(reads), breaks = 30)


## Raw data transformations ####
data_spe_log <- log1p(data_spe) # log transformation
# log1p(x) is the same as log(x+1)
head(data_spe_log)

reads_log <- rowSums(data_spe_log)
hist(reads_log)


data_spe_sqrt <- sqrt(data_spe) # square root transformation
head(data_spe_sqrt)


## Dissimilarity ####
#Jaccard
d.jac <- vegdist(data_spe, method = "jaccard")
d.jac

#Bray
d.bry <- vegdist(data_spe_log, method = "bray")
d.bry

#Horn-Morisita
d.hor <- vegdist(data_spe, method = "horn")
d.hor




## NMDS (Non-metric Multi-Dimensional Scaling) ####
set.seed(123)
nmds.nuts <- metaMDS(data_spe_log,             # Our community-by-species matrix
                     distance = "bray",    # Ecological distance (respects abundance)
                     # distance = "jaccard", # Species presence-absence
                     k = 2,                # The number of reduced dimensions
                     trymax = 100)         # The number of default iterations

# metaMDS has automatically applied a square root transformation and calculated the 
# Bray-Curtis distances for our community-by-site matrix.

# Shepard plot, which shows scatter around the regression between the interpoint distances 
# in the final configuration (i.e., the distances between each pair of communities) against 
# their original dissimilarities.

stressplot(nmds.nuts) # Large scatter around the line suggests that original dissimilarities are not well preserved in the reduced number of dimensions

# Basic plot function
plot(nmds.nuts)
plot(nmds.nuts, type = "t")


# dbzp02 <- as.data.frame(data_spe["dbzp02", ])
# summary(dbzp02)


# an10 <- as.data.frame(data_spe["an_10", ])


# # Remove these rows from data_spe
# data_spe_log <- log1p(
#   data_spe[!(rownames(data_spe) %in% "dbzp02"), 
#            , 
#            drop = FALSE]
#   )

# # Remove the same rows from env_data that we removed from data_spe earlier
# env_data <- env_data[!(rownames(env_data) %in% "dbzp02"), , drop = FALSE]



# Filter and fit environmental data to the ordination
envfit_vars <- env_data %>%
  select(clay, organic, pH, precipitation, temperature)

fit <- envfit(nmds.nuts, envfit_vars, permutations = 999)

# Plot environmental data vectors
plot(fit)

# Step-by-step plotting
ordiplot(nmds.nuts, type="n")                                   # coordinate plot
orditorp(nmds.nuts, display="species",col="red", air=0.01)      # species with names
orditorp(nmds.nuts, display="sites",cex=1,air=0.01)             # groups





## Environmental factors as contour lines# Use the function ordisurf to plot contour lines
# png('nmds_pf_altitude.png', units = "cm", width = 22, height = 18, res = 300)
ordisurf(nmds.nuts ~ temperature, envfit_vars)
ordisurf(nmds.nuts ~ temperature, envfit_vars, main = "", cex = 0, col = "blue")
orditorp(nmds.nuts, display = "sites", col="grey30", air=0.1, cex=0.8)
# dev.off()# Use the function ordisurf to plot contour lines



ordiplot(nmds.nuts, type = "n")                        # basic plot
ordisurf(nmds.nuts ~ precipitation, envfit_vars, main = "", cex = 0, col = "blue") # Altitude contour lines
ordisurf(nmds.nuts ~ temperature, envfit_vars, main = "", cex = 1, col = "red")  # species + deadwood contour lines






# STEP 6: Plot with ggplot2
library(ggplot2)
library(ggrepel)

# Extract NMDS scores
site_scores <- as.data.frame(scores(nmds.nuts, display = "sites"))
site_scores$sample_id <- rownames(site_scores)

# species_scores <- as.data.frame(scores(nmds.nuts, display = "species"))
# species_scores$species <- rownames(species_scores)

vectors <- as.data.frame(scores(fit, display = "vectors"))
vectors$variable <- rownames(vectors)

# Combine with environmental metadata
# Assumes rows are in same order!
nmds_df <- bind_cols(site_scores, env_data)

# Plot ordination with samples as points and environmental metadata as vectors
ggplot(nmds_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(
    colour = Ecosystem, 
    shape = Ecosystem
    ), 
    alpha = 0.5,
    size = 4) +
  # scale_size_continuous(range = c(2, 6)) +
  theme_minimal() +
  labs(title = "NMDS Ordination of EcM Fungi data",
       colour = "Ecosystem",
       shape = "Ecosystem"
       ) + 
  geom_segment(data = vectors,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  geom_text(data = vectors,
            aes(x = NMDS1, y = NMDS2, label = variable),
            hjust = 1.1, vjust = 1.1)


ggsave("./figures/fungi_nmds.png", width = 8.3, height = 5, units = "in", dpi = 300)




# We can also add species scores:
ggplot(nmds_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(colour = region, shape = origin), size = 5) +
  labs(title = "NMDS Ordination of Nut & Fruit mix data",
       colour = "Region",
       shape = "Origin") +
  geom_text_repel(data = species_scores,
                  aes(x = NMDS1, y = NMDS2, label = species),
                  colour = "darkred", size = 3.5) +
  geom_segment(data = vectors,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")),
               colour = "black") +
  geom_text(data = vectors,
            aes(x = NMDS1, y = NMDS2, label = variable),
            hjust = 1.1, vjust = 1.1, size = 4.5) +
  theme_bw() 

ggsave("./figures/nuts_nmds.png", width = 8.3, height = 5, units = "in", dpi = 300)
