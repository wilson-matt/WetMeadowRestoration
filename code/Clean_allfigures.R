####Read in final csv for quads and start here####
setwd("C:/GitHub/WetMeadowRestoration/data")
quads_complete <- read.csv("Quads_seeding_meta_fnl.csv")


#Functional groups:
library(dplyr)
library(ggplot2)

# Subset -- Filter was removing all zero values, even though that's not what your code said to do.
  #This is why tidyverse is garbage and I'll never use it :)
wetland_data <- subset(quads_complete, !quads_complete$NSC == "Contr" & !quads_complete$Wetland.Indicator == "NA")

# Set logical order on x-axis
wetland_data$Wetland.Indicator <- factor(
  wetland_data$Wetland.Indicator,
  levels = c("OBL", "FACW", "FAC", "FACU", "UPL")
)

#Wetland Functional Groups - pretty figure
ggplot(wetland_data,        ###two samples for biennial.beeblossom (FACU, 20% and 7.5% success) omitted to reduce y-axis for clarity.
       aes(x = Wetland.Indicator,
           y = seeding_success)) +
  geom_boxplot(fill = "grey75",
               width = 0.6,
               outlier.shape = 16,
               alpha = 0.8) +
  geom_jitter(width = 0.15,
              alpha = 0.4,
              size = 1) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Wetland Functional Group",
    y = "Proportional Seeding Success",
#    title = "Quad ZOOMED"
  ) +
  scale_y_continuous(
    limits = c(0, 0.035),
    breaks = seq(0, 0.08, by = 0.01),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black")
  )

#use next line to resave the plot if there are changes to be made:
#ggsave(filename = "C:/GitHub/WetMeadowRestoration/figures/FAC_success.jpg", width = 8, height = 10, device='jpeg', dpi=700)


#Wetland functional group stats - Kruskal-Wallis and Dunn's test to compare differences between groups:
library(rstatix)

kruskal.test(seeding_success ~ Wetland.Indicator, data = wetland_data)
dunn_test(wetland_data, seeding_success ~ Wetland.Indicator, p.adjust.method = "bonferroni")
pairwise.wilcox.test(wetland_data$seeding_success,
                     wetland_data$Wetland.Indicator,
                     p.adjust.method = "BH")

####Success by species ####
setwd("C:/GitHub/WetMeadowRestoration/data")
quads_complete <- read.csv("Quads_seeding_meta_fnl.csv")
no_circs <- subset(quads_complete, !quads_complete$NSC == "Contr")

#Species success with SCIENTIFIC names (instead of common names)#
sci.name <- read.csv("WetMeadowScientific Names.csv")
colnames(sci.name)[1:2] <- c("CommonName", "ScientificName")

seeding_merged <- merge(
  no_circs,
  sci.name,
  by.x = "Species.sitecircrat",
  by.y = "CommonName",
  all.x = TRUE
)
seeding_merged$seeding_success <- seeding_merged$Count / seeding_merged$Wet.Meadow
seeding_merged$seeding_success[is.infinite(seeding_merged$seeding_success)] <- NA
seeding.short <- subset(seeding_merged, !is.na(seeding_success))
seeding.shorter <- subset(seeding.short, seeding_success > 0)
scaleFactor <- 250  # 50 / 0.20


ggplot(
  seeding.shorter %>%
    dplyr::filter(
      !is.na(ScientificName),
      !ScientificName %in% 
        c("Schizachyrium scoparium",
          "Juncus spp.",
          "Penstemon digitalis",
          "Asclepias spp.",
          "Baptisia australis",
          "Lycopus americanus",
          "Chasmanthium latifolium",
          "Monarda fistulosa",
          "Panicum rigidulum")),
  aes(x = reorder(ScientificName, seeding_success))
) +
  geom_boxplot(aes(y = Count), 
               fill = "grey55",
               position = position_nudge(x = -.21), width = 0.4) +
  geom_boxplot(aes(y = seeding_success * scaleFactor), 
               fill = "grey90",
               position = position_nudge(x = .21), width = 0.4) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "italic")) +  # scientific formatting
  scale_y_continuous(
    name = "Stems/m^2 (dark grey)",
    limits = c(0, 50),
    breaks = seq(0, 50, by = 10),
    sec.axis = sec_axis(
      transform = ~ . / scaleFactor,   # updated (no deprecation warning)
      name = "Success rate (stems per seeds planted; light grey)",
      breaks = seq(0, 0.20, by = 0.05),
      labels = scales::number_format(accuracy = 0.01)
    )
  ) +
  labs(title = "", x = "Species", y = NULL)
#use next line to resave the plot if there are changes to be made:
#ggsave(filename = "C:/GitHub/WetMeadowRestoration/figures/Species_success.jpg", width = 6, height = 8, device='jpeg', dpi=700)

#Success by species stats:
#This needs some work - let's think about what we're really trying to say statistically in association with this figure.
#Just added this as a starting point.
library(mgcv)
mod <- glm(seeding_success ~ Species.sitecircrat, data = seeding.shorter)
summary(mod)


####Ordination ####
setwd("C:/GitHub/WetMeadowRestoration/data")
quad <- read.csv("Quadrats.csv")
quad$m2 <- quad$Count/1


library(vegan)
quad$Sample.ID.event <- paste(quad$Sample.ID,quad$Sample.event)
ordtable <- as.data.frame.matrix(xtabs(m2~Sample.ID.event+Spp, quad))
ordtable$Sample.ID <- rownames(ordtable)
ordtable$NSC <- substr(ordtable$Sample.ID, 1, 5)
ordtable$Season <- substr(ordtable$Sample.ID, 8, 30)

#Variance Partitioning models:
mod_site <- cca(ordtable[,1:33]~ordtable$NSC)
summary(mod_site)
anova(mod_site)

mod_season <- cca(ordtable[,1:33]~ordtable$Season)
summary(mod_season)
anova(mod_season)

mod_both <- cca(ordtable[,1:33]~ordtable$NSC+ordtable$Season)
summary(mod_both)
anova.cca(mod_both)

site_scores <- scores(mod_both, display="sites") 
species_scores <- scores(mod_both, display="species")


#Combined a/b CCA figure
jpeg(filename = "C:/GitHub/WetMeadowRestoration/figures/CCA_combined.jpg",
     width = 1600, height = 800, units = "px", pointsize = 12, quality = 500)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
#panel A (NSC)
plot(site_scores, type = "n",
     xlab = "Axis 1 (13.7%)",
     ylab = "Axis 2 (6.9%)",
     xlim = c(-1.3, 2),
     cex.lab = 1.5)
abline(h = 0, col = "white")
abline(v = 0, col = "white")
box()
cols <- gray.colors(3, start = 0.05, end = 0.6)
points(mod_both,
       col = "black",
       pch = c(8,11,15)[as.numeric(as.factor(ordtable$NSC))],
       cex = 1.8)
ordiellipse(site_scores, ordtable$NSC,
            kind = "se", conf = 0.95,
            lwd = 0.75,
            draw = "polygon",
            col = cols, border = cols, alpha = 30)

text(-.75,-.9, "South", cex = 1.3)
text(-.8,.6, "North", cex = 1.3)
text(1,.4, "Reference", cex = 1.3)

legend("topright",
       legend=c("Reference", "North","South"),
       col = "black",
       pch = c(8,11,15),
       bty = "n",
       cex = 1.2)
mtext("a", side = 3, adj = 0, line = 0.5, cex = 3.5, font = 2)

#panel b (season)
plot(site_scores, type = "n",
     xlab = "Axis 1 (13.7%)",
     ylab = "Axis 2 (6.9%)",
     xlim = c(-1.3, 2),
     cex.lab = 1.5)

abline(h = 0, col = "white")
abline(v = 0, col = "white")
box()

points(mod_both,
       col = "black",
       pch = c(8,11,15)[as.numeric(as.factor(ordtable$Season))],
       cex = 1.8)

ordiellipse(site_scores, ordtable$Season,
            kind = "se", conf = 0.95,
            lwd = 0.75,
            draw = "polygon",
            col = cols, border = cols, alpha = 30)
text(-.35,1, "September", cex = 1.3)
text(-.2,-0.1, "October", cex = 1.3)
text(0,-1.1, "July", cex = 1.3)

legend("topright",
       legend = levels(as.factor(ordtable$Season)),
       col = "black",
       pch = c(8,11,15),
       bty = "n",
       cex = 1.2)
mtext("b", side = 3, adj = 0, line = 0.5, cex = 3.5, font = 2)

dev.off()

#Functional Donut####
####Create df (species, success rate, wetland indicator, grass/forbes, per/annual/bi)
library(webr)
meta <- read.csv("C:/GitHub/WetMeadowRestoration/Data/Meta.csv")
quads_meta <- read.csv("C:/GitHub/WetMeadowRestoration/Data/Quads_seeding_meta_fnl.csv")
success <- quads_meta[,c(2,13)]
# merged to create table
merge_meta <- merge(meta, success, by=1)
colnames(merge_meta)

#### Rate calcs 
#wetland funct
wetfunct <- aggregate(seeding_success ~ Wetland.Indicator, data= merge_meta, 
                      FUN= function(x) sum(x, na.rm = TRUE)/sum(merge_meta$seeding_success, na.rm= TRUE))
#grasses/forbs/sedge/rush/fern
collective.class <-aggregate(seeding_success ~ Wetland.Indicator + Forbes.grasses.sedge.rush.fern, data= merge_meta, 
                             FUN= function(x) sum(x, na.rm = TRUE)/sum(merge_meta$seeding_success, na.rm= TRUE))
#ALL (ADDING perennial/annual/biennial)
ALL.class <- aggregate(seeding_success ~ Wetland.Indicator + Forbes.grasses.sedge.rush.fern + Annual.perennial.biennial, 
                       data= merge_meta, FUN= function(x) sum(x, na.rm = TRUE)/sum(merge_meta$seeding_success, na.rm= TRUE))

### Donut pie chart 

#### Rate calcs 
#wetland funct
wetfunct <- aggregate(seeding_success ~ Wetland.Indicator, data= merge_meta, 
                      FUN= function(x) sum(x, na.rm = TRUE)/sum(merge_meta$seeding_success, na.rm= TRUE))
#grasses/forbs/sedge/rush/fern
collective.class <-aggregate(seeding_success ~ Wetland.Indicator + Forbes.grasses.sedge.rush.fern, data= merge_meta, 
                             FUN= function(x) sum(x, na.rm = TRUE)/sum(merge_meta$seeding_success, na.rm= TRUE))
#ALL (ADDING perennial/annual/biennial)
ALL.class <- aggregate(seeding_success ~ Wetland.Indicator + Forbes.grasses.sedge.rush.fern + Annual.perennial.biennial, 
                       data= merge_meta, FUN= function(x) sum(x, na.rm = TRUE)/sum(merge_meta$seeding_success, na.rm= TRUE))


###
#Can't handle the "ALL.class" because there are three levels:
collective.class <-aggregate(seeding_success ~ Wetland.Indicator + Forbes.grasses.sedge.rush.fern, data= merge_meta, 
                             FUN= mean)#I think the plotting function does what we were imagining with the fun function? (nice work, btw!)


PieDonut(collective.class, aes(Wetland.Indicator, Forbes.grasses.sedge.rush.fern, count=seeding_success),  title = "Functional Classification Breakdown", 
         showRatioPie = FALSE,
         showRatioDonut = FALSE)

