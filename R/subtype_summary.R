#' Define the row colours for the heatmap.
#'
#' @param m A data frame with gene symbols as row names.
plot_summary=function(m){
  print('Plotting summary...')
  print(head(m))
  mm=melt(m,id.vars=1,measure.var=c(5:9))
  #print(head(mm))
  #print(dim(mm))
  pdf(paste0(inputDir,"subtype_summary.pdf"))
  g = ggplot(data = mm, aes(x = sub("_.*","",Sample), fill = value)) + geom_bar(position="fill")
  g = g + labs(title = "Classification counts", y = "Classification Percentage", x = "Sample", fill = "PAM50 Subtype")
  g = g + theme(text = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, size=4))
  print(g)
  dev.off()
}

m=read.delim('ccle_data/subtype_summary.tsv')
mm=plot_summary(m)
