library(readr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(plyr)
library(cowplot)
library(grid)
library(lattice)
library(gridExtra)
library(utils)

My_Theme = theme(text = element_text(size=rel(3)))
gpPlain = gpar(fontsize = 9)
gp = gpar(fontsize = 9, fontface='bold')
tagsize = 9
rates = c("lambda", "mu", "psi", "nLTT") #"Reff" #"pulledDR"
rate_symbols = c(expression(lambda), expression(mu), expression(psi),  "nLTT") #expression(tilde(r)),
individual = c("nLTT") #"Reff", 
text_labels = c("speciation", "extinction", "sampling",  "normalised LTT") #"pulled diversification rate", #"effective\nreproductive ratio"
gridExp = NULL
gridLinear = NULL
# Code for replotting FBDsim graphs to adjust visualisation.
for (type in c("linear")){
  max_Reff = 0
  base_plots = list()
  kappa_plots = list()
  graphed_below = 0
  base_count = 0
  #base_plots[[base_count]] = textGrob("Rates", gp = gp, just="center")
  kappa_plot_count = 0
  files = list.files(path=type, pattern="*.tsv", all.files=FALSE, full.names=FALSE)
  for (f in files){
    plot_label = substr(f, 1, nchar(f)-4)
    kappas = c()
    congruence = c()
    con <- file(paste(type,"/",f, sep=""), "r") 
    lines <- c()
    count = 1
    while(TRUE) {
      line = readLines(con, 1)
      if (length(line) == 0){
        break
      }
      if (substr(line,1,1) == "#" & grepl("[[:alnum:]]", substr(line,str_length(line), str_length(line)))){
        simLoc = str_locate(line, "run")
        if(!(is.na(simLoc[2]))){
          kappaLoc = str_locate(line, "kappa_")
          afterKappa = str_split(substr(line,kappaLoc[2]+1,str_length(line)), "_")[[1]]
          kappaVal = as.numeric(afterKappa[1])
          congruentVal = as.numeric(substr(afterKappa[2],2, str_length(afterKappa[2])))
          congruentCategory = TRUE
          if(congruentVal == 1){
            congruentCategory = FALSE
          }
          kappas[count] = kappaVal
          congruence[count] = congruentCategory
          count = count + 1
        }
      }
    }
    close(con)
    
    df <- as.data.frame(read_tsv(paste(type, "/",f, sep=""), comment="#", col_names=TRUE))
    table_splits <- c(0,which(df[1] == "age"))
    allDfs = list()
    start = 1
    for (s in 1:length(table_splits)){
      if (s < length(table_splits)){
        end = table_splits[s + 1] -1 
      }
      else {
        end = nrow(df) 
      }
      # cat("Slicing ", start, end, "\n")
      dfSplit = as.data.frame(df[start:end,])
      for (i in 1:ncol(dfSplit)){
        dfSplit[[i]]<-as.numeric(dfSplit[[i]])
      }
      allDfs[[s]] =dfSplit
      start = table_splits[s+1]+1
    }
    #str(allDfs)
    if(plot_label == "Reff"){
      for(i in 1:length(allDfs)){
        max_val = max(allDfs[[i]]$Reff)
        if (max_val > max_Reff || max_Reff == 0){
          max_Reff = max_val
        }
      }
      cat("max_Reff", max_Reff)
    }
    
    orderedKappas = sort(unique(kappas))
    base_rate = FALSE
    include_individual = FALSE
    if(plot_label %in% individual){
      include_individual = TRUE
    }
    for (k in 1:length(orderedKappas)){
      o_data = allDfs[[which(kappas == orderedKappas[k] & congruence == FALSE)]]
      c_data = allDfs[[which(kappas == orderedKappas[k] & congruence == TRUE)]]
      combine_data = o_data
      colnames(c_data) = c("age", "var")
      combine_data$ageC = c_data$age
      combine_data$varC = c_data$var
      colnames(combine_data) = c("age", "var", "ageC", "varC")
      if (plot_label == rates[1]){
        if(type == "exp"){
          coords = c(0.25, 0.68)
        }
        else {
          coords = c(0.25, 0.37)
        }
        p <- ggplot(data=combine_data) +
          geom_line(aes(x=age, y=var, color="Original", linetype="Original")) +
          geom_line(aes(x=ageC, y=varC, color="Congruent", linetype="Congruent"), linetype=4, lwd=0.9) +
          theme_classic() +
          scale_color_manual(values = c("Original" = "darkblue", "Congruent" = "black"), 
                             guide = guide_legend(override.aes = list(lwd=c(0.9,0.5), linetype=c(4,1)))) +
          guides(linetype ="none") +
          labs(color="Scenario") +
          scale_x_reverse() +
          theme(legend.position=coords, legend.title=element_text(size=9), legend.text=element_text(size=9),
                legend.background = element_blank(), legend.key.height=unit(0.7,"line")) + My_Theme
      }
      else {
        p <- ggplot(data=combine_data) +
          geom_line(aes(x=age, y=var, linetype="Original"), color="darkblue") +
          geom_line(aes(x=ageC, y=varC,  linetype="Congruent"), color="black", linetype=4, lwd=0.9) +
          theme_classic() + 
          guides(linetype ="none") + #linetype
          scale_x_reverse() + My_Theme
      }
      if (include_individual){
        plot_count = length(setdiff(rates, individual)) + graphed_below + 1
      }
      else {
        plot_count = which(rates == plot_label)
      }
      if (plot_label == "nLTT"){
        labP <- p + labs(x = "Time before present", y=plot_label, tag=toupper(letters[plot_count])) + theme(plot.tag=element_text(size=tagsize))
      }
      else{
        #if (plot_label == "lambda"){
        #  labP <- p + labs(x = "Time before present", y=rate_symbols[rates==plot_label], tag=toupper(letters[plot_count]))
        #}
        #else {
        labP <- p + labs(y=rate_symbols[rates==plot_label], x="", tag=toupper(letters[plot_count])) + theme(plot.tag=element_text(size=tagsize))
        #}
      }
      if (include_individual) {
        #if (orderedKappas[k] == 0){
        #  kappa_plot_count = kappa_plot_count + 1
        #  kappa_plots[[kappa_plot_count]] = textGrob(text_labels[rates==plot_label], gp = gp, just="center")
        #}
        kappa_plot_count = kappa_plot_count + 1
        graphed_below = graphed_below + 1
        kappa_plots[[kappa_plot_count]] = labP
      }
      else {
        if (plot_label %in% rates && !(include_individual)){
          base_count = base_count + 1
          base_plots[[base_count]] = labP
          break
        }
      }
    }
  }
  rate_labels = list()
  rate_labels[[1]] = textGrob("")
  for (i in 1:length(rates)){
    if(!(rates[i] %in% individual)){
      rate_labels[[i]] = textGrob(text_labels[i], gp = gp, just="center")
    }
  }
  if(type == "exp"){
    typeLabel = "Exponential parameters"
  }
  else {
    typeLabel = "Linear parameters"
  }
  gtop = grid.arrange(grobs=c(rate_labels, base_plots), ncol = 3, nrow=2, heights=c(0.2, 1))
  gbottom = grid.arrange(grobs=c(list(textGrob(""), textGrob("normalised LTT", gp = gp, just="center"), textGrob(""), textGrob("r = 1", gp = gpPlain, just="center"), textGrob("r = 0.5", gp = gpPlain), textGrob("r = 0", gp = gpPlain)), kappa_plots), nrow=3, ncol = 3, heights=c(0.1,0.1,1))
  if(type == "exp"){
    gridExp = grid.arrange(gtop, gbottom, heights=c(1,1), top=textGrob(typeLabel, just="center", gp = gpar(fontsize = 11, fontface = 'bold')))  
  }
  else {
    gridLinear = grid.arrange(gtop, gbottom, heights=c(1,1), top=textGrob(typeLabel, just="center", gp = gpar(fontsize = 11, fontface = 'bold')))  
  }
}
barwidth = 0.9
bargap = (1-barwidth)/2
gBar = grid.arrange(textGrob(""), textGrob(""), textGrob(""), widths=c(bargap,barwidth,bargap))
gAll = grid.arrange(textGrob(""), gridExp, textGrob(""), gBar,  textGrob(""), gridLinear, textGrob(""), heights = c(0.01,1,0.02,0.001,0.02,1, 0.01))




