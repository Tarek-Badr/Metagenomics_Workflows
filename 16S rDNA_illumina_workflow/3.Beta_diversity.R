
# Beta diversity #

#### Package setup ####

library(phyloseq)
library(ggplot2)
library(ggordiplots)
library(plyr)
library(gridExtra)
library(openxlsx)

#### End ####

#### Load RData with filtered and normalized phyloseq object ####

load("Results/RData/2.phyloseq.filtered.RData")

# Create directories for results
dir.create("Results/5.Beta_plots")
dir.create("Results/4.Beta_stats")

#### End ####

#### Ordination ####

# Ordinate with method "NMDS" and "PCoA" on distance "Bray-Curtis" 

ps.prop=physeq_re

ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
ord.pcoa.bray <- ordinate(ps.prop, method="PCoA", distance="bray")

#### End ####

#### Plot ordination with ggordiplots and ggplot2 ####

sampdf = data.frame(sample_data(ps.prop))

ordiplot <- gg_ordiplot(ord.nmds.bray, groups=sampdf$Plate, hull = TRUE, spiders = TRUE, 
                        ellipse = FALSE, plot = FALSE)


# Take the data from the gg_ordiplot object to plot with our preferred parameters

ordiplot$df_spiders$patient <- gsub("-.*", "", rownames(ordiplot$df_spiders))
ordiplot$df_hull$patient <- rep(NA, nrow(ordiplot$df_hull))

ggplot(data = as.data.frame(ordiplot$df_spiders), aes(x = x, y = y, color = patient, shape = Group)) +
  geom_point(size=2) +
  geom_point(size=7, alpha=0.5) +
  geom_segment(aes(yend=cntr.y,xend=cntr.x,group=Group), color = "gray", size=1, alpha=0.5) +
  geom_polygon(data=ordiplot$df_hull,aes(x=x,y=y,fill=Group,group=Group), color="gray", alpha=0.2) +
  #scale_fill_brewer(type="qual", palette="Set2", name = "Plate") +
  guides(fill = guide_legend(title="Plate"), shape=guide_legend(title="Plate"), color =guide_legend(title = "Sample")) +
  xlab("NDMS 1") + ylab("NDMS 2") +
  ggtitle("NMDS ordination using Bray–Curtis dissimilarity index") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=20)) 

#### End ####

#### Plot ordination with plot_ordination from Phyloseq ####

plot_ps_nmds = plot_ordination(ps.prop, ord.nmds.bray, color="Plate", title="NMDS ordination of Plate samples using Bray–Curtis dissimilarity") 

NMDS = plot_ps_nmds + geom_point(size=8, alpha=0.5)  + stat_ellipse(type = "norm", linetype = 2) + stat_ellipse(type = "t") + theme_bw() #default "t" assumes a multivariate t-distribution while norm (fine line) assumes a multivariate normal distribution

plot_ps_pcoa = plot_ordination(ps.prop, ord.pcoa.bray, color="Plate", title="PCoA ordination of Plate samples using Bray–Curtis dissimilarity") 

PCOA = plot_ps_pcoa + geom_point(size=8, alpha=0.5) + stat_ellipse(type = "norm", linetype = 2) + stat_ellipse(type = "t") + theme_bw() #default "t" assumes a multivariate t-distribution while norm (fine line) assumes a multivariate normal distribution

NP <- grid.arrange(NMDS + ggtitle("NMDS") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
                                                  text = element_text(size=15)),
                   PCOA + ggtitle("PCoA") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
                                                  text = element_text(size=15)),
                   nrow=1, top=textGrob("NMDS and PCoA ordination of Plate samples using Bray–Curtis dissimilarity",  gp = gpar(fontsize = 20)))

grid.arrange(NP, get_legend(NMDS + theme(legend.position = "bottom", text = element_text(size=15))), nrow= 2, heights = c(10,1))

#### End ####

#### Plot ordination with only top 5 Abundant Phyla ####

# Keep only the most abundant five phyla to Plot
phylum.sum = tapply(taxa_sums(ps.prop), tax_table(ps.prop)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
phy_top5 = prune_taxa((tax_table(ps.prop)[, "Phylum"] %in% top5phyla), ps.prop)
phy_top5.ord <- ordinate(phy_top5, method="NMDS", distance="bray")

p1 = plot_ordination(phy_top5, phy_top5.ord, type="Taxa", color="Phylum", title="taxa")
print(p1)
p1 + facet_wrap(~Phylum, 3)

p2 = plot_ordination(phy_top5, phy_top5.ord, type="Biplot", color="Phylum", shape="Plate", title="biplot")
print(p2)

p3 = plot_ordination(phy_top5, phy_top5.ord, type="Split", color="Phylum", shape="Plate", title="split") 
print(p3)

#### End ####

#### Plot Beta Diversity with All Methods c("DCA", "CCA", "RDA", "CAP", "NMDS", "MDS", "PCoA") ####

dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "NMDS", "MDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, physeq, dist){
  ordi = ordinate(physeq, method=i, distance=dist)
  plot_ordination(physeq, ordi, "samples", color="Plate")
}, ps.prop, dist)
names(plist) <- ord_meths

pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})

names(pdataframe)[1] = "method"
p5 = ggplot(pdataframe, aes(Axis_1, Axis_2, color=Plate, shape=Plate, fill=Plate))
p5 = p5 + geom_point(size=4) #+ geom_polygon()
p5 = p5 + facet_wrap(~method, scales="free")
p5 = p5 #+ scale_colour_brewer(type="qual", palette="Set2")
print(p5)

#### End ####

####ADONISimplementation of a non-parametric permutation based MANOVA (PERMANOVA)####
#https://rstudio-pubs-static.s3.amazonaws.com/343284_cbadd2f3b7cd42f3aced2d3f42dc6d33.html

dist = phyloseq::distance(ps.prop, method="bray")

set.seed(123456789)
test.adonis <- adonis2(dist ~ Plate, data = sampdf, permutations = 99999)
test.adonis
write.xlsx(test.adonis, file = "Results/4.Beta_stats/test.adonis.xlsx", rowNames = T)

# Beta dispersion
dispersion<-betadisper(dist, group=sampdf$Plate)
dispersion
permutest(dispersion)
anova(dispersion)

#### End ####
