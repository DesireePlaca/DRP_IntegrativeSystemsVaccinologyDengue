library(ggpubr)
library(lemon)
library(ggplot2)
library(openxlsx)
library(stringr)
library(reshape2)
library(dplyr)

tabl <- "input.data"
tablebox <- melt(tabl)
names(tablebox)[names(tablebox) == "Groups"] <- "Group"
names(tablebox)[names(tablebox) == "variable"] <- "gene"
names(tablebox)[names(tablebox) == "value"] <- "counts"

#Severity
my_comparisons <- list(c("OFI", "DF"), c("OFI", "DHF"),c("OFI", "DSS"), c("DF", "DHF"), c("DF", "DSS"), c("DHF", "DSS") )#define comparisons between groups
#c("Control","DF"), c("Control", "DHF"), c("DF", "DHF")
#c("OFI", "DF"), c("OFI", "DHF"), c("DF", "DHF"), c("OFI", "DSS"), c("DF", "DSS"), c("DHF", "DSS")
symnum.args <- list(
  cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.01), #1.01 Because in some rare cases, 1 is shown as NS for no reason whatsoever 
  symbols = c("****", "***", "**", "*", waiver()))

tablebox$counts <- as.numeric(tablebox$counts)

tablebox$Group <- factor(tablebox$Group, levels= c("OFI", "DF", "DHF", "DSS")) #define groups
tablebox$gene <- factor(tablebox$gene, levels=str_sort(unique(tablebox$gene)))
#"Control","DF", "DHF"
#png(file = "Fig01-IgA.png", bg = "transparent", width = 4500, height = 2500, units = "px", res = 300)
ggboxplot(tablebox, x = "Group", y = "counts",
          color = "Group", palette = c("sandybrown", "palegreen3","cornflowerblue", "darkblue"),
          add = c("jitter", "median_iqr"), #"mean_sd"
          add.params = list(size = 0.5))+
  scale_y_continuous(breaks = pretty(c(4,20), n = 5), #Values showed from 0 to 8 with 4 breaks will show 0, 2, 4, 6 and 8
                     limits = c(4,20), name="log2 counts")+ # y scale + axis' name
  facet_rep_wrap(~gene, strip.position="bottom", ncol=4)+ #if you have more than 2 lines, use ( facet_rep_wrap(~gene, strip.position="bottom") ) instead of facet_Wrap
  stat_compare_means(method= "wilcox.test", label = "p.format", 
                     label.y = c(16, 17, 18, 19, 20, 21, 22), bracket.size = 0.4, tip.length = 0.006, #location and brackets customization
                     hide.ns = T, aes(Group = Group), comparisons = my_comparisons, 
                     vjust = 0.5, size = 5, # * customization: vjust == 0 -> in the middle of the bracket; accepts negative values
                     symnum.args = symnum.args)+
  theme(strip.text.x = element_text(size = 12))+
  theme(axis.line.x=element_line(), axis.ticks.x=element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(),
        strip.background = element_rect(color="white", fill="white"), strip.placement = "outside",
        panel.border = element_blank(), panel.background = element_blank(), panel.grid = element_blank(), 
        panel.spacing.x = unit(0,"line"), panel.spacing.y = unit(0,"line"))

dev.off()

#Vaccine 
my_comparisons <- list( c("Day0", "Day8"), c("Day0", "Day28"), c("Day8", "Day28"))#define comparisons between groups

symnum.args <- list(
  cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.01), #1.01 Because in some rare cases, 1 is shown as NS for no reason whatsoever 
  symbols = c("****", "***", "**", "*", waiver()))

tablebox$counts <- as.numeric(tablebox$counts)

tablebox$Group <- factor(tablebox$Group, levels= c("Day0","Day8", "Day28")) #define groups
tablebox$gene <- factor(tablebox$gene, levels=str_sort(unique(tablebox$gene)))

#png(file = "Fig01-IgA.png", bg = "transparent", width = 4500, height = 2500, units = "px", res = 300)
ggboxplot(tablebox, x = "Group", y = "counts",
          color = "Group", palette = c("sandybrown", "palegreen3","cornflowerblue", "darkblue"),
          add = c("jitter", "median_iqr"), #"mean_sd"
          add.params = list(size = 0.5))+
  scale_y_continuous(breaks = pretty(c(0,22), n = 5), #Values showed from 0 to 8 with 4 breaks will show 0, 2, 4, 6 and 8
                     limits = c(0,22), name="log2 counts")+ # y scale + axis' name
  facet_rep_wrap(~gene, strip.position="bottom", ncol=5)+ #if you have more than 2 lines, use ( facet_rep_wrap(~gene, strip.position="bottom") ) instead of facet_Wrap
  stat_compare_means(method= "wilcox.test", label = "p.format", 
                     label.y = c(17, 18.5, 19), bracket.size = 0.4, tip.length = 0.006, #location and brackets customization
                     hide.ns = F, aes(Group = Group), comparisons = my_comparisons, 
                     vjust = 0.5, size = 5, # * customization: vjust == 0 -> in the middle of the bracket; accepts negative values
                     symnum.args = symnum.args)+
  theme(strip.text.x = element_text(size = 12))+
  theme(axis.line.x=element_line(), axis.ticks.x=element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(),
        strip.background = element_rect(color="white", fill="white"), strip.placement = "outside",
        panel.border = element_blank(), panel.background = element_blank(), panel.grid = element_blank(), 
        panel.spacing.x = unit(0,"line"), panel.spacing.y = unit(0,"line"))
dev.off()


