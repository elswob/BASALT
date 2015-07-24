#' Define the row colours for the heatmap.
#'
#' @param m The final subtype summary dataframe
plot_summary=function(m,outDir){
  print('Plotting summary...')
  print(outDir)
  print(head(m))
  mm<-melt(m,id.vars=1,measure.var=c(5:9))
  mm$value=factor(mm$value,levels=sort(levels(mm$value)))
  pdf(paste0(outDir,"/subtype_summary.pdf"))
  g = ggplot(data = mm, aes(x = sub("_.*","",Sample), fill = value)) + geom_bar(position="fill")
  g = g + labs(title = "Classification counts", y = "Classification Percentage", x = "Sample", fill = "PAM50 Subtype")
  g = g + theme(text = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, size=4))
  print(g)
  dev.off()
}

#m=read.delim('subtype_summary.tsv')
#plot_summary(m,"...")
