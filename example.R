# R-code of "Structural changes within trophic levels are constrained by within-family assembly rules at lower trophic levels" by:
# Chuliang Song, Florian Altermatt, Ian Pearse, Serguei Saavedra
# published in: Ecology Letters

# load the toolboxs
source("toolbox.R")

# load necessary R packages
library(tidyverse)
library(readxl)

# read (binary) meta herbivore-plant interaction matrix
matrix <- read.table("matrix.txt", header = TRUE)
matrix <- t(matrix[-1])
matrix <- (matrix > 0) * 1

# read the arrival order of plants
arrivals <- read_csv("Plant_arrival.csv") %>%
  .[, c(2:5)] %>%
  as_tibble() %>%
  filter(!is.na(time)) %>% #native plants are filtered
  filter(time <= 2000) %>%
  arrange(desc(time)) %>%
  filter(species %in% rownames(matrix)) %>%
  mutate(
    order = 1:nrow(.),
    degree = rowSums(matrix[rownames(matrix) %in% .$species, ])
  )

#extract the time-dependent matrix of herbivore competition and
#compute the structural stability of community persistence (for each arrival time) 
times <- unique(arrivals$time) # the arrival time is # of years ago, counting from 2015
year <- 2015 - times[-1] # all the arrival years
structural_stability <- 2:length(times) %>% 
  map_dbl( ~ arrivals[arrivals$time >= times[.x], ] %>%
    extract_binary_Inte(matrix) %>%
    get_interaction_matrix() %>%
    Omega())

#Plot Figure 2 of main text 
tibble(year, structural_stability) %>% 
  ggplot(aes(x=year, y =structural_stability)) +
  geom_point(size=1, col='grey') +
  geom_smooth(method = "lm", se = FALSE, linetype='dotted') +
  xlab("Year") +
  ylab("Structural Stability") +
  theme(
    strip.background = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    text = element_text(size=15),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5, linetype="solid"))