library(readr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(plyr)
library(cowplot)
library(grid)
library(lattice)
library(gridExtra)


gp = gpar(fontsize = 13, fontface = 'bold')
rates = c("lambda", "mu", "psi")
rate_symbols = c(utf8_print("\u03BB"), utf8_print("\u03BC"),  utf8_print("\u03C8"))
individual = c("PDR", "Reff", "nLTT")
# Code for replotting FBDsim graphs to adjust visualisation.
for (type in c("linear", "exp")){
  max_Reff = 0
  base_plots = list()
  kappa_plots = list()
  base_count = 1
  base_plots[[base_count]] = textGrob("Rates", gp = gp, just="center")
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
        simLoc = str_locate(line, "sim")
        if(!(is.na(simLoc[2]))){
          kappaLoc = str_locate(line, "kappa_")
          afterKappa = str_split(substr(line,kappaLoc[2]+1,str_length(line)), "_")[[1]]
          kappaVal = as.numeric(afterKappa[1])
          congruentVal = as.numeric(substr(afterKappa[2],2, str_length(afterKappa[2])))
          cat(congruentVal, "\n")
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
      cat("Slicing ", start, end, "\n")
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
    if(plot_label %in% rates){
      base_rate = TRUE
    }
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
        p <- ggplot(data=combine_data) +
          geom_line(aes(x=age, y=var, color="Original", linetype="Original")) +
          geom_line(aes(x=ageC, y=varC, color="Congruent", linetype="Congruent"), linetype=4, lwd=0.9) +
          theme_classic() +
          scale_color_manual(values = c("Original" = "darkblue", "Congruent" = "black"), 
                             guide = guide_legend(override.aes = list(lwd=c(0.9,0.5), linetype=c(4,1)))) +
          guides(linetype ="none") +
          labs(color="Scenario") +
          scale_x_reverse() +
          theme(legend.position=c(0.7,0.5))
      }
      else {
      p <- ggplot(data=combine_data) +
        geom_line(aes(x=age, y=var, linetype="Original"), color="darkblue") +
        geom_line(aes(x=ageC, y=varC,  linetype="Congruent"), color="black", linetype=4, lwd=0.9) +
        theme_classic() +             ) +
        guides(linetype ="none") + #linetype
        scale_x_reverse()
    }
      labP <- p + labs(x = "Time before present", y="Density") #title=paste(plot_label, ", kappa =", orderedKappas[k]) 
      if (base_rate){
        base_count = base_count + 1
        base_plots[[base_count]] = labP
        break
      }
      else if (include_individual) {
        if (orderedKappas[k] == 0){
          kappa_plot_count = kappa_plot_count + 1
          kappa_plots[[kappa_plot_count]] = textGrob(plot_label, gp = gp, just="center")
        }
        kappa_plot_count = kappa_plot_count + 1
        kappa_plots[[kappa_plot_count]] = labP
      }
    }
  }
   rate_labels = list()
  rate_labels[[1]] = textGrob("")
  for (i in 1:length(rates)){
    rate_labels[[i + 1]] = textGrob(rate_symbols[i], gp = gp, just="center")
  }
  if(type == "exp"){
    typeLabel = "Exponential"
  }
  else {
    typeLabel = "Linear"
  }
  gtop = grid.arrange(grobs=c(rate_labels, base_plots), ncol = 4, heights=c(0.1, 0.3), widths=c(0.5,1,1,1))
  gbottom = grid.arrange(grobs=c(list(textGrob(""), textGrob("kappa = 0", gp = gp, just="center"), textGrob("kappa = 0.5", gp = gp), textGrob("kappa = 1", gp = gp)), kappa_plots), ncol = 4,  widths=c(0.5,1,1,1), heights=c(0.1, rep(0.3, 3)))
  grid.arrange(gtop, gbottom, heights=c(1.3,3.3), top=textGrob(typeLabel, just="center", gp = gpar(fontsize = 15, fontface = 'bold')))  
}




