#' Function to combine plots
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' Generate the PAM50 correlation coefficient plots
#'
#' @param p50 The PAM50 results text file
#' @param id The ID of the samples to extract
#' @param name The name to represent the samples
#' @param outDir The output directory
cor_plot = function(p,id,name,outDir){
  print(paste0("Generating correlation plots - ",id," - ",name))
  p$X=rownames(p)
  print(head(p))
  #get colors
  #myPal <- brewer.pal(6, "Spectral")
  myPal = c("red","green","blue","orange","purple","black")
  names(myPal) <- c("LumB","LumA","Basal","Normal","Her2","Confidence")

  subtypes=c("Basal","Her2","LumA","LumB","Normal","Confidence")
  plotter=function(orderType){
    g = ggplot(fm, aes(x=X,y=value,group=variable,color=variable))
    g = g + geom_line()
    g = g + scale_color_manual(values=myPal[as.character(fm$variable)],name = "")
    g = g + ggtitle(paste0(name," - sorted by ",orderType))

    #g = g + ylab("Correlation coefficient") + xlab("Single Cells")
    #g = g + theme(axis.text.x = element_blank(), axis.ticks.x = element_line(size = 3, colour = colorOrder))
    g = g + ylab("Correlation coefficient") + xlab("Samples")
    g = g + theme(axis.ticks.x = element_line(size = 3, colour = colorOrder), axis.text.x  = element_text(size=3,angle=90))

    if (orderType == "Call"){
      g = g + theme(legend.position="bottom", axis.text.y = element_blank())
    }else{
      g = g + theme(legend.position="bottom")
    }
    g = g + theme(plot.title=element_text(family="Times", face="bold", size=10), axis.title.x = element_text(size=6), axis.text.y  = element_text(size=6), axis.title.y = element_text(size=6), legend.text = element_text(size=5))
    return(g)
  }
  p=p[grep(id,p$X),]

  #order by Call or confidence
  #for Call use p$Call, confidence use -p$Confidence
  p$X <- factor(p$X, levels = p[order(-p$Confidence), "X"])
  fm=melt(p,id.vars=c("X","Call"),measure.vars=subtypes)
  #set the colours for the ticks
  colorOrder = myPal[as.character(p[order(-p$Confidence),]$Call)]
  p1=plotter("Confidence")

  p$X <- factor(p$X, levels = p[order(p$Call), "X"])
  fm=melt(p,id.vars=c("X","Call"),measure.vars=subtypes)
  print(head(fm))
  #set the colours for the ticks
  colorOrder = myPal[as.character(p[order(p$Call),]$Call)]
  p2=plotter("Call")

  pdf(paste0(outDir,"/",short,"_cor_coef_plot.pdf"))
  multiplot(p1,p2,cols=2)
  dev.off()
}
