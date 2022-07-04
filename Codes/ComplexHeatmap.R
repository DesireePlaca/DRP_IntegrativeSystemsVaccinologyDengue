
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
install.packages("circlize")
library(circlize)
library(writexl)

#inputdata
heatEX<-read.delim("overlapFC_ALL.txt", header=T, sep="\t", row.names = 1)

#convert the dataframe in matrix
heatEX <-as.matrix(heatEX)

ha1 = HeatmapAnnotation(
  dist1 = anno_barplot(
    colSums(heatEX), 
    bar_width = 1, 
    gp = gpar(col = "white", fill = "lightblue"), 
    border = FALSE),
  show_annotation_name = FALSE)
ha2 = rowAnnotation(
  dist2 = anno_barplot(
    rowSums(heatEX), 
    bar_width = 1, 
    gp = gpar(col = "white", fill = "lightblue"), 
    border = FALSE
  ), show_annotation_name = FALSE)

col_fun = RColorBrewer::brewer.pal(n =9, name = "RdYlBu") 
#to see more colors: https://www.r-graph-gallery.com/38-rcolorbrewers-palettes.html

ht_list = Heatmap(heatEX, name = "expression",
                  cluster_columns = TRUE, show_row_dend = TRUE, rect_gp = gpar(col= "white"), 
                  show_column_names = TRUE, column_names_gp = gpar(fontsize = 6),
                  row_names_side = "left", row_names_gp = gpar(fontsize = 9),
                  column_title = 'Gene Expression',
                  top_annotation = ha1, col = col_fun) + ha2
                   
draw(ht_list, ht_gap = unit(1, "mm"))


#decorate_heatmap_body("cases", {
 # i = which(colnames(heatEX) == "1961")
  #x = i/ncol(heatEX)
  # grid.lines(c(x, x), c(0, 1), gp = gpar(lwd = 2, lty = 2))
  #grid.text("Vaccine introduced", x, unit(1, "npc") + unit(5, "mm"))
#})


#Export table
heatEX <-as.data.frame(heatEX)
write_xlsx(heatEX,file.path(outDir, cond,"heatEX.xlsx"))


#heatEX <-t(heatEX) #transpose the matrix

#View(heatEX)
#write_xlsx(heatEX,"/Users/otavio.cmarques/Documents/Otavio/DISCIPLINA POS/R - Atualizado/Complex Heatmap with bars/heatEX.xlsx")
#transform the matrix in dataframe

#heatEX[1]= NULL
#boxplot(heatEX)  