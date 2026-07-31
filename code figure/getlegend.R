library(ggplot2)
library(gridExtra)

# Function to extract legend
get_legend <- function(myplot) {
  tmp <- ggplot_gtable(ggplot_build(myplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Extract legend
legend_only <- get_legend(p)

# Save legend as PDF
pdf(file = file.path(save_path,"figureS2-2_legend.pdf"), width = 1, height = 2)  # adjust size as needed
grid.arrange(legend_only)
dev.off()
