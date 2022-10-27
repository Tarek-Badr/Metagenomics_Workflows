library("reshape2")
library("ggplot2")
library("openxlsx")
library("tidyverse")


setwd("D:/Nano16S")

Matrixsamp <- read.xlsx("00-Nano_clinical_sample-Results_11okt2022.xlsx", sheet = 1, rowNames = T)
Matrixsamp <- Matrixsamp[,-c(1:3,36)]
colnames(Matrixsamp) <- sub("_2_", "2_", colnames(Matrixsamp)) %>% sub("_", ".", .) %>% sub("_.*", "", .) %>% sub("\\.", "_", .)



tot <- summarise_all(Matrixsamp, ~if(is.numeric(.)) sum(.) else "Total")
tot <- tot[rep(seq_len(nrow(tot)), each = 138), ]

Matrixpercent <- Matrixsamp / tot*100

Matrixpercent$name <- rownames(Matrixpercent)
Matrixpercent <- Matrixpercent[,c("name", colnames(Matrixpercent)[-length(colnames(Matrixpercent))])]
Matrixreshaped <- melt(Matrixpercent, id.vars = "name", measure.vars = colnames(Matrixpercent)[-1])
Matrixreshaped <- Matrixreshaped[Matrixreshaped$value != 0,]


#Matrixplot <- Matrixreshaped[Matrixreshaped$value > 5,]

for (i in 1:nrow(Matrixreshaped)) {
  if(Matrixreshaped$value[i] < 5){
    Matrixreshaped$name[i] <- "Other"
  }
}



# Stacked bar plot
ggplot(Matrixreshaped, aes(fill=name, x=variable, y = value)) + 
  geom_bar(position="fill", stat="identity") +
  #ggtitle("Title") +
  labs(fill='Bacteria', x = "Sample") +
  theme(panel.background = element_rect(fill="#F9FCFB", colour="LightGrey"),
        panel.grid = element_line(colour="LightGrey"),
        text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title = element_text(size=10),
        #plot.title = element_text(size = 10),
        legend.title = element_text(size = 10), #angle = 90, vjust = 0.5, hjust=1), 
        legend.text  = element_text(size = 9),
        legend.spacing = unit(0.3, "cm"),
        legend.key.size = unit(0.2, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = 1) 

write.xlsx(Matrixreshaped, "Stacked_bar_plot_noMock_Other_5.xlsx")
ggsave("Plot.pdf")

