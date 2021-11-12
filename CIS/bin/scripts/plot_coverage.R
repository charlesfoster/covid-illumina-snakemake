#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)

### main df
df <- fread(args[1])
names(df) <- c("Genome", "Site", "Coverage")
df$cov_range <- df$Coverage
df$cov_range[df$Coverage < 10] <- "low"
df$cov_range[df$Coverage >= 10 & df$Coverage < 100] <- "high"
df$cov_range[df$Coverage >= 100 & df$Coverage < 1000] <- "very high"
df$cov_range[df$Coverage >= 1000] <- "highest"
df$cov_range <- factor(df$cov_range, levels = c("low", "high", "very high", "highest"))

df$Coverage_Capped <- df$Coverage
df$Coverage_Capped[df$Coverage >= 50] <- 50

### amplicon wrangling

primers <- fread(args[2])
primers <- primers[,1:6]
names(primers) <- c("Genome", "Start", "End", "Primer", "Num", "Strand")
primers$Amplicon <- gsub("_R|_F", "", primers$Primer)

amplicons <- bind_rows(lapply(unique(primers$Amplicon), function(x){
  start <- primers %>% 
    filter(Amplicon == x, Strand == "+") %>% 
    dplyr::select(Start) %>% 
    as.numeric
  end <- primers %>% 
    filter(Amplicon == x, Strand == "-") %>% 
    dplyr::select(End) %>% 
    as.numeric
  if(is.na(start) || is.na(end)){
    return(NULL)
  }
  res <- data.frame(Amplicon=x, Start=start, End=end)
  res$Coverage <- mean(df$Coverage[res$Start[1]:res$End[1]])
  return(res)
}))

### plotting

y_max <- (mean(df$Coverage)/10)
y_min <- 0

pos1 <- -y_max*3
plot1 <- ggplot(df, aes(x=Site, y=Coverage)) + geom_line() + labs(title = paste0("Coverage plot: ", args[6]), x = "Site in SARS-CoV-2 Genome", y = "Coverage") 
plot1 <- plot1 + annotate("rect", xmin = 0,   xmax = 265, ymin = y_min, ymax = y_max, alpha = .2, fill = "blue")
plot1 <- plot1 + annotate("rect", xmin = 266,   xmax = 21555, ymin = y_min, ymax = y_max, alpha = .2, fill = "red")
plot1 <- plot1 + annotate("rect", xmin = 21563, xmax = 25384, ymin = y_min, ymax = y_max, alpha = .2, fill = "blue")
plot1 <- plot1 + annotate("rect", xmin = 25393, xmax = 26220, ymin = y_min, ymax = y_max, alpha = .2, fill = "red")
plot1 <- plot1 + annotate("rect", xmin = 26245, xmax = 26472, ymin = y_min, ymax = y_max, alpha = .2, fill = "blue")
plot1 <- plot1 + annotate("rect", xmin = 26523, xmax = 27191, ymin = y_min, ymax = y_max, alpha = .2, fill = "red")
plot1 <- plot1 + annotate("rect", xmin = 27202, xmax = 27387, ymin = y_min, ymax = y_max, alpha = .2, fill = "blue")
plot1 <- plot1 + annotate("rect", xmin = 27394, xmax = 27759, ymin = y_min, ymax = y_max, alpha = .2, fill = "red")
plot1 <- plot1 + annotate("rect", xmin = 27894, xmax = 28259, ymin = y_min, ymax = y_max, alpha = .2, fill = "blue")
plot1 <- plot1 + annotate("rect", xmin = 28274, xmax = 29533, ymin = y_min, ymax = y_max, alpha = .2, fill = "red")
plot1 <- plot1 + annotate("rect", xmin = 29558, xmax = 29674, ymin = y_min, ymax = y_max, alpha = .2, fill = "blue")
plot1 <- plot1 + annotate("rect", xmin = 29675, xmax = 29903, ymin = y_min, ymax = y_max, alpha = .2, fill = "red")

plot1 <- plot1 + annotate("text", angle=90, hjust =0, x = (1+265)/2, y = pos1, label = "5' UTR")
plot1 <- plot1 + annotate("text", angle=90, hjust =0, x = (266+21555)/2, y = pos1, label = "ORF1ab")
plot1 <- plot1 + annotate("text", angle=90, hjust =0, x = (21563+25384)/2, y=pos1, label = "S")
plot1 <- plot1 + annotate("text", angle=90, hjust =0, x = (25393+26220)/2, y=pos1, label = "ORF3a")
plot1 <- plot1 + annotate("text", angle=90, hjust =0, x = (26245+26472)/2, y=pos1, label = "E")
plot1 <- plot1 + annotate("text", angle=90, hjust =0, x = (26523+27191)/2, y=pos1, label = "M")
plot1 <- plot1 + annotate("text", angle=90, hjust =0, x = (27202+27387)/2, y=pos1, label = "ORF6")
plot1 <- plot1 + annotate("text", angle=90, hjust =0, x = (27394+27759)/2, y=pos1, label = "ORF7a")
plot1 <- plot1 + annotate("text", angle=90, hjust =0, x = (27894+28259)/2, y=pos1, label = "ORF8")
plot1 <- plot1 + annotate("text", angle=90, hjust =0, x = (28274+29533)/2, y=pos1, label = "N")
plot1 <- plot1 + annotate("text", angle=90, hjust =0, x = (29558+29674)/2, y=pos1, label = "ORF10")
plot1 <- plot1 + annotate("text", angle=0, hjust =0, x = (29675+29903)/2, y =-y_max, label = "3' UTR")

# Note: I'm disabling the coverage capped plot. Not useful.
plot2 <- ggplot(df, aes(x=Site, y=Coverage_Capped)) + 
  geom_line() + labs(x = "Site in SARS-CoV-2 Genome", y = "Coverage (capped at 50)") 

plot3 <- ggplot(amplicons, aes(x=Amplicon, y=Coverage)) + 
  geom_bar(stat="identity", color="blue") + labs(title = paste("Mean coverage per amplicon for: ", args[6], " (", gsub(".bed","",args[5]), " protocol)", sep=""), x = "Amplicon", y = "Mean Coverage") +
  coord_flip() + scale_x_discrete(limits = rev(levels(amplicons$Amplicon)))

coverage_plot <- ggarrange(plot1, plot3, labels = c("A", "B"), ncol = 1, nrow = 2, align = "v")

ggsave(filename = args[3], plot = coverage_plot, width=16, height=12)

#### Plot dodgy amplicons ####

a2 <- amplicons
a2$Amplicon <- as.character(a2$Amplicon)

my_plots <- lapply(a2$Amplicon, function(x){
  coords <- amplicons %>%
    filter(Amplicon==x)
  df2 <- df %>% filter(Site %in% coords$Start:coords$End)
  if(any(df2$Coverage < 10)){
    #print(x)
    p <- ggplot(df2, aes(x=Site, y=Coverage)) + 
      geom_line() + labs(x = "Site in SARS-CoV-2 Genome", y = "Coverage") 
    p<- p + ggtitle(x)
    return(p)
  }
})

my_plots <- my_plots[which(lapply(my_plots,is.null) == FALSE)]

if(args[5] != "swift.bed"){
  if(length(my_plots) > 0){
   bad_amplicons <- ggarrange(plotlist=my_plots, 
                               labels = LETTERS[1:length(my_plots)], 
                               ncol = 2, nrow = ceiling(length(my_plots)/2), 
                               align = "v")
   ggsave(filename = args[4], plot = bad_amplicons)
  }
}
#### EXPLORATORY
#plot1 <- ggplot(df, aes(x=Site, y=Coverage, fill=cov_range)) + 
#  geom_area() + labs(x = "Site in SARS-CoV-2 Genome", y = "Coverage") +
#  scale_fill_manual(name = "Coverage", labels = c("<10", "10 to 99", "100 to 999", "1000+"), values = c("#D55E00", "#F0E442", "#009E73", "#56B4E9"))
#
#plot2 <- ggplot(df, aes(x=Site, y=Coverage_Capped, fill=cov_range)) + 
#  geom_area() + labs(x = "Site in SARS-CoV-2 Genome", y = "Coverage (capped at 100)") +
#  scale_fill_manual(name = "Coverage", labels = c("<10", "10 to 99", "100 to 999", "1000+"), values = c("#D55E00", "#F0E442", "#009E73", "#56B4E9"))
#
#plot3 <- ggplot(amplicons, aes(x=Amplicon, y=Coverage)) + 
#  geom_bar(stat="identity", color="blue") + labs(x = "Amplicon", y = "Mean Coverage") +
#  coord_flip()
#
