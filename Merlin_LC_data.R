###
#Analyzes Drosophila melanogaster live/dead data from a 96-well bioassay 
#author James Klimavicz
#TyraTech
#Version 2.1.05
#Jan 19 2021
###

library(tidyverse)
library(drc)
library(ggplot2)
library(binom)
library(readxl)
library(hash)
library(RColorBrewer)
# library("png")
library(grid)
library(gridExtra)
library(pdftools)

#Input/output
# Merlin.path <- "C:\\Users\\jklimavicz\\OneDrive - TyraTech, Inc\\Merlin\\Data_analysis"

input.file <- "C:/Users/jklimavicz/Downloads/Merlin_assay_data_stats.csv"
merlin.key.file <- "C:\\Users\\jklimavicz\\OneDrive - TyraTech, Inc\\Merlin\\Data_analysis\\Merlin_key.csv"
out.path <- "C:\\Users\\jklimavicz\\OneDrive - TyraTech, Inc\\Merlin\\Data_analysis"
LC.out.file <- "Merlin_LC50_data.csv"

earliest.date <- as.Date("10/30/2020", "%m/%d/%Y")

# input.file <- "geraniol_10162020.xlsx"
# LC.out.file <- "geraniol_LC50_data.csv"


#Function for plotting the worm data
plot_worm_data <- function(worm.model, cmpd.title = ""){
   conf.level = 0.8 #for binomial CI error bars
   
   #Determine log-2 values of the min and max concentrations
   max.val <- log2(max(worm.model$data$ppm,na.rm =T))
   min.val <- log2(min(worm.model$data$ppm,na.rm =T))
   
   #Make a vector of concentration values on log scale for plotting the dose-response curve
   demo.fits <- expand.grid(conc=exp(seq(log(2^(min.val-1.5)), 
                                         log(2^(max.val+1.5)), length=181)))
   #Use the above set of concentrations to make an array of estimated values and CI bounds
   pm <- predict(worm.model, newdata=demo.fits, interval="confidence") 
   demo.fits$p <- pm[,1] #estimate dose-response curve
   demo.fits$pmin <- pm[,2] #LB
   demo.fits$pmax <- pm[,3] #UB
   
   demo.fits$pmin[demo.fits$pmin < -0.1] <- -0.1
   demo.fits$pmax[demo.fits$pmax > 1.1] <- 1.1
   
   #Jittering is useful for multiple data points overlapping
   # jitter.factor <- 0.05    
   # nem.melt.func<- nem.melt.func %>% 
   #    mutate(Conc.jitter = runif(Conc, min = (1-jitter.factor)*Conc, 
   #                               max = (1+jitter.factor)*Conc))
   
   #Extract experimental concentrations, counts, and mortalities
   worm.data <- worm.model$data[c(1,2,5)]
   
   #Add binomial error bar data to the experimental mortalities 
   worm.data <- worm.data %>% mutate(CI = 
                                        binom.confint(ratio.dead*weights, weights, 
                                                      conf.level = conf.level, methods="profile")) 
   worm.data<- worm.data %>% mutate(LB = CI$lower, UB = CI$upper)
   
   # breaks <- c(2^c(seq(max.val,min.val,-1)))
   # labels <- ifelse(breaks>=1,breaks, parse(text =(paste((1/breaks),'^-{1}',sep=""))))
   # breaks <- c(breaks, 2^(min.val-2))
   # labels <- c(labels, 0)
   
   #Make break/labels for the graph. Note that 0 is included on the graph at 
   #two log conc below the lowest data value
   breaks <- c(2^seq(max.val,min.val,-1))
   labels <- signif(breaks,4)
   breaks <- c(breaks, 2^(min.val-2))
   labels <- c(labels, 0)
   
   p <- worm.data %>% ggplot(aes(x=ppm, y=ratio.dead)) + 
      geom_point() + #experimental data points
      coord_trans(x="log2") +
      scale_x_continuous(name = "Concentration (ppm)", 
                         breaks=breaks,
                         minor_breaks = NULL,
                         labels=labels )+
      geom_errorbar(data=worm.data, aes(x= ppm, ymin=LB, ymax=UB), width= 0) + #experimental error bars
      geom_ribbon(data=demo.fits, aes(x=conc, y=p, ymin=pmin, ymax=pmax), alpha=0.2, fill = "grey10") + #dose-response curve confidence intervals
      geom_line(data=demo.fits, aes(x=conc, y=p)) + #dose response curve
      theme_minimal() +
      scale_y_continuous(breaks=c(seq(0,3,.2)),
                         labels=as.character(c(seq(0,3,.2))),
                         name = "Portion Dead",
                         minor_breaks = c(seq(-0.1,2.9,.2)), 
                         limits = c(-0.1,1.1)) + 
      theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5)) +
      ggtitle(cmpd.title)
   # time.string <- strptime(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "EST5EDT")
   # p <- p + annotate("text", label = time.string, size = 4, alpha = 0.3, y=0.01, x=-0.2)
   return(p)
}


#data entry and cleaning
worm_input <- read.csv(input.file) %>% filter(!is.na(Date))
# worm_input <- read_excel(input.file)
worm_input$Row <- factor(worm_input$Row)
worm.names <- colnames(worm_input)
worm.names[1] <- "Ref.ID"
colnames(worm_input) <- worm.names
worm_input$Date <- as.Date(worm_input$Date, "%m/%d/%Y")

worm_input <- worm_input %>% filter(Date > earliest.date)

worm_input$ratio.dead <- worm_input$Dead/worm_input$Count


merlin.key <- read.csv(merlin.key.file)
worm_input <- worm_input %>%left_join(merlin.key, by="Ref.ID")

worm_input$Compound[is.na(worm_input$Compound)] <- worm_input$Ref.ID[is.na(worm_input$Compound)]
worm_input$Compound <- as.factor(worm_input$Compound)
worm_input$Parent <- as.factor(worm_input$Parent)
#remove controls
worm_input <- worm_input %>% filter(!(ppm == 0)) %>% droplevels()
#remove bad datapoint
worm_input <- worm_input %>% filter(!grepl("OVERCROWDED",Notes))
worm_input <- worm_input %>% filter(Date >= earliest.date)%>% droplevels()
worm_input <- worm_input %>% mutate(ID = ifelse(is.na(Plate), 
                                       paste(Date, ".id=", substr(Ref.ID,1,nchar(Ref.ID)-1), sep=""), 
                                       paste(Date, ".plate=", Plate, ".id=", substr(Ref.ID,1,nchar(Ref.ID)-1),sep="")))

##QUALITY CONTROL: BASED ON MALATHION CONTROL
mal.data <- worm_input %>% filter(Compound=="malathion") %>% droplevels()


jitter.factor <- 0.05
mal.data <- mal.data %>%
   mutate(ppm.jitter = runif(ppm, min = (1-jitter.factor)*ppm, max = (1+jitter.factor)*ppm))

mal.curve.data <- data.frame()
methods = c( "BFGS", "Nelder-Mead","CG", "SANN")
for (m in levels(as.factor(mal.data$ID))) {
   curr.data <- mal.data %>% filter(ID == m)
   for (n in methods){
      #fit the experimental data to a dose-response curve. For noisy data, curve-
      #fitting may be difficult; therefore, we can try several different nonlinear
      #optimization methods if one or more fails to converge.
      fit <- try(drm(ratio.dead~ppm, 
                     data = curr.data, 
                     fct= LL.2(), 
                     weights = curr.data$Count,
                     start = c(-1,0.1), lowerl = c(-10,0.00001), upperl = c(-0.05,100),
                     control=drmc(method = n, noMessage=T, maxIt = 5000)),silent = T)
      if(!inherits(fit, "try-error")) {break} #fit was found
   }
   if(inherits(fit, "try-error")) {next} #no fit found
   # fit <- try(drm(ratio.dead~ppm, 
   #                data = curr.data, 
   #                fct= LL.2(), 
   #                weights = curr.data$Count,
   #                start = c(-0.5,0.1), lowerl = c(-10,0.00001), upperl = c(-0.05,100),
   #                control=drmc(method = n, noMessage=T, maxIt = 5000)),silent = T)
   mal.check <- ED(fit, c(50), interval = "tfls", level = 0.95, display=F)
   slope <- signif(fit$coefficients[["b:(Intercept)"]],3)
   mal.EC <- signif(mal.check[1,1],4)
   mal.str <- paste("EC50: ", mal.EC, ", slope: ", slope, ", ", m, sep="")
   print(mal.str)

   max.val <- log2(max(curr.data$ppm,na.rm =T))
   min.val <- log2(min(curr.data$ppm,na.rm =T))
   demo.fits <- expand.grid(conc=exp(seq(log(2^(min.val-2.5)),
                                         log(2^(max.val+2.5)), length=201)))
   #Use the above set of concentrations to make an array of estimated values and CI bounds
   pm <- predict(fit, newdata=demo.fits, interval="confidence")
   demo.fits$p <- pm[,1] #estimate dose-response curve
   demo.fits$pmin <- pm[,2] #LB
   demo.fits$pmax <- pm[,3] #UB
   demo.fits$ID <- m
   mal.curve.data <- rbind(mal.curve.data,demo.fits)

}

max.val <- log2(max(mal.data$ppm,na.rm =T))
min.val <- log2(min(mal.data$ppm,na.rm =T))
breaks <- c(2^seq(max.val,min.val,-1))
labels <- signif(breaks,4)
breaks <- c(breaks, 2^(min.val-2))
labels <- c(labels, 0)

n.lines <- nlevels(as.factor(mal.data$ID))
n.mal.colors <- 8
colors <- rep(brewer.pal(n = n.mal.colors, name = "Set1"), times = n.lines%/%n.mal.colors + 1)
lines <- rep(c("solid", "dashed", "twodash", "dotted", "longdash"), each = n.mal.colors)
shapes <- rep(c(16, 17, 15, 8, 14, 4), each = n.mal.colors)
# mal.curve.data$factorID <- as.numeric(as.factor(mal.curve.data$ID))
# mal.curve.data <- mal.curve.data %>% mutate(colorID = colors[(factorID-1)%%n.mal.colors+1], 
#                                             lineID=lines[(factorID-1)%/%n.mal.colors+1])
# mal.data <- mal.data %>% mutate(colorID = colors[(as.numeric(as.factor(ID))-1)%%n.mal.colors+1])

q <- mal.curve.data %>% ggplot(aes(x=conc, y=p, fill=as.factor(ID)))+
   coord_trans(x="log2") +
   scale_x_continuous(name = "Concentration (ppm)",
                      breaks=breaks,
                      minor_breaks = NULL,
                      labels=labels,
                      limits = c(2^(min.val-1), 16))+
   theme_minimal() +
   #geom_ribbon(aes(x=conc, y=p, ymin=pmin, ymax=pmax), alpha=0.2) + #dose-response curve confidence intervals
   geom_line(data = mal.curve.data, 
             aes(x=conc, y=p, color =as.factor(ID), linetype = as.factor(ID)), size = 0.4) +  #dose response curve
   scale_y_continuous(breaks=c(seq(0,3,.2)),
                      labels=as.character(c(seq(0,3,.2))),
                      name = "Portion Dead",
                      minor_breaks = c(seq(-0.1,2.9,.2)),
                      limits = c(-0.1,1.1)) +
   theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5),
         legend.text=element_text(size=5)) +
   geom_point(data = mal.data, aes(x=ppm.jitter, y=ratio.dead, color=as.factor(ID), shape = as.factor(ID)),  size=1) + 
   scale_fill_manual(values=colors) +
   scale_color_manual(values=colors) +
   scale_shape_manual(values=shapes) +
   scale_linetype_manual(values=lines) +
   theme(legend.position = "none") + 
   # labs(color = "ID", linetype = "ID", shape = "ID", fill ="ID")
   ggtitle("malathion")

ggsave(path = out.path, units = "in", width = 6.4, height = 4.134,
       filename = paste("malathion_curves.png"), 
       plot=q, device = "png")



#available non-linear optimization solvers for drm

model <- LL.2() #model selection
#determine number of compounds for preallocation of sufficient dataframe memory
n.cmpd <- length(levels(worm_input$Ref.ID)) 
worm.LC.data <- data.frame(Compound  = levels(worm_input$Compound), 
                           LC50 = NaN,
                           LC50.CI = NA,
                           LC75 = NaN,
                           LC75.CI = NA,
                           LC90 = NaN,
                           LC90.CI = NA,
                           slope = NaN,
                           slope.CI = NA,
                           n.trials = NA,
                           curve.qual = 0,
                           malathion.potency.rel.to.cmpd = NA,
                           potency.rel.to.geraniol = NA,
                           potency.rel.to.parent = NA,
                           fit = NA,
                           class = NA,
                           parent = NA,
                           PAINS = NA,
                           modifications = "",
                           min.conc = NaN)


for (m in levels(worm_input$Compound)) {
   #find index of worm.LC.data array that matches the correct treatment compound
   ind <- which(worm.LC.data$Compound==m) 
   old.ind <- which(worm_input$Compound==m)[1]
   curr.class <- worm_input$Class[old.ind]
   curr.subclass <- worm_input$Subclass[old.ind]
   curr.gen <- worm_input$Generation[old.ind]
   curr.parent <- as.character(worm_input$Parent[old.ind])
   worm.LC.data$class[ind] <- ifelse(curr.class %in% c("Conventional","Natural"), paste(curr.class, " (", curr.subclass, ")", sep=""), curr.class)
   worm.LC.data$parent[ind] <- curr.parent
   worm.LC.data$PAINS[ind] <- worm_input$PAINS[old.ind]
   worm.LC.data$modifications[ind] <- worm_input$Modifications[old.ind]
   #filter the data to that matching the compound
   cmpd.dat<- worm_input %>% filter(Compound==m) %>% droplevels()
   #count number of unique trials
   # worm.LC.data$n.trials[ind] <- cmpd.dat %>% group_by(Date, Row) %>% 
   #    summarise(sum(Dead), .groups="keep") %>% nrow()
   worm.LC.data$n.trials[ind] <- cmpd.dat %>% group_by(Date) %>% 
      summarise(sum(Dead), .groups="keep") %>% nrow()
   max_ppm  <- max(cmpd.dat$ppm)
   min_ppm  <- min(cmpd.dat$ppm)
   worm.LC.data$min.conc[ind] <- 2^(log(min_ppm,2)-2)
   
   
   for (n in methods){
      #fit the experimental data to a dose-response curve. For noisy data, curve-
      #fitting may be difficult; therefore, we can try several different nonlinear
      #optimization methods if one or more fails to converge.
      fit <- try(drm(ratio.dead~ppm, 
                     data = cmpd.dat, 
                     fct= model, 
                     type="continuous",
                     weights = cmpd.dat$Count,
                     start = c(-2.0,(max_ppm+min_ppm)/2),
                     lowerl = c(-3.5,min_ppm/2),
                     upperl = c(-0.05,max_ppm*2),
                     control=drmc(method = n, noMessage=T, maxIt = 5000)),silent = T)
      if(!inherits(fit, "try-error")) {break} #fit was found
      else{
         #If no optimized converged, give up and keep NAs/NaNs in dataframe, and 
         #ensure the fit is NA instead of error class.
         fit <- NA
      }
   }
   if(is.na(fit)){
      if(sum(cmpd.dat$Live)==0) {
         worm.LC.data$LC50[ind] <- paste("<", signif(min_ppm,3), sep = "")
         worm.LC.data$LC75[ind] <- paste("<", signif(min_ppm,3), sep = "")
         worm.LC.data$LC90[ind] <- paste("<", signif(min_ppm,3), sep = "")
      next
      }
   }
   worm.LC.data$fit[ind] <- list(fit) #Store the fit as a list object
   #extract the slope of the dose-response curve
   slope <- worm.LC.data$fit[ind][[1]]$coefficients[["b:(Intercept)"]]
   if(slope>0){ 
      #Here, the data was so noisy that our dose-response curve goes the wrong way, e.g.
      #we predict that more larvae die at low concentrations than at high concentrations.
      #Thus, the estimated LC50/75/90 are not meaningful, so we will not calculate 
      #them because this could be misleading later.
      next 
   } else { 
      #Calculate LC50, LC75, and LC90 values and CIs. Note that we're on a log
      #scale, but we want non-log intervals, so we use 'interval = "tfls"'.
      #We are calculating 95% confidence intervals.
      worm.EC <- ED(worm.LC.data$fit[ind][[1]], c(50,75,90), interval = "tfls", 
                    level = 0.95, display=F)
      if(worm.EC[1,1] > 0.995*sqrt(2)*max_ppm) { #curve maxed out without fitting data.
         worm.LC.data$LC50[ind] <- paste(">", signif(sqrt(2)*max_ppm,3), sep = "")
         worm.LC.data$LC75[ind] <- paste(">", signif(sqrt(2)*max_ppm,3), sep = "")
         worm.LC.data$LC90[ind] <- paste(">", signif(sqrt(2)*max_ppm,3), sep = "")
         worm.LC.data$slope[ind] <- "NA" #store slope in the dataframe
      } else {
         worm.LC.data$LC50[ind] <- signif(worm.EC[1,1],4) #LC50 estimate
         worm.LC.data$LC75[ind] <- signif(worm.EC[2,1],4) #LC75 estimate
         worm.LC.data$LC90[ind] <- signif(worm.EC[3,1],4) #LC90 estimate
         #LC50 95% confidence interval
         worm.LC.data$LC50.CI[ind] <- paste("[", signif(worm.EC[1,3],3), ", ", 
                                            signif(worm.EC[1,4],3), "]",sep = "")
         #LC75 95% confidence interval
         worm.LC.data$LC75.CI[ind] <- paste("[", signif(worm.EC[2,3],3), ", ", 
                                            signif(worm.EC[2,4],3), "]",sep = "")
         #LC90 95% confidence interval
         worm.LC.data$LC90.CI[ind] <- paste("[", signif(worm.EC[3,3],3), ", ", 
                                            signif(worm.EC[3,4],3), "]",sep = "")
         worm.LC.data$slope[ind] <- signif(slope,3) #store slope in the dataframe
         worm.LC.data$slope.CI[ind] <- paste("[", signif(confint(fit,"b")[1],3), ", ", 
                                             signif(confint(fit,"b")[2],3), "]",sep = "")
         worm.LC.data$curve.qual[ind] <- round(100/(1+((confint(fit,"b")[2]-confint(fit,"b")[1])/2)^2)* 
                                                   ifelse(worm.LC.data$LC50[ind]<min_ppm*2 | 
                                                             worm.LC.data$LC50[ind]>max_ppm/2, 0.5, 1))
         
      }
      
      
   }
   
   if(m=="malathion"){malathion.ref <- worm.EC[1,1]}
   if(m=="geraniol"){geraniol.ref <- worm.EC[1,1]}
}
options(scipen = 50)
worm.LC.data$malathion.potency.rel.to.cmpd <- formatC(signif(as.numeric(worm.LC.data$LC50)/malathion.ref,3), format="fg")
worm.LC.data$potency.rel.to.geraniol<- formatC(signif(geraniol.ref/as.numeric(worm.LC.data$LC50),3), format="fg")


#determine potency relative to parent compound. Stores LD50 values with a hash with compound name as a key. 
h <- hash()
parents.levels <- levels(as.factor((worm.LC.data %>% filter(parent != ""))$parent))
for(m in parents.levels){
   cmpd <- tolower(m)
   ind <- which(tolower(worm.LC.data$Compound)==cmpd)[1]
   h[[cmpd]] <- worm.LC.data$LC50[ind]
}

#now use hash to determine relative potency 
for(ind in seq(1,nrow(worm.LC.data))){
   if(is.na(worm.LC.data$parent[ind])){next}
   if(worm.LC.data$parent[ind] == ""){next}
   if(!(tolower(worm.LC.data$parent[ind]) %in% tolower(parents.levels))){next}
   if(is.na(as.numeric(worm.LC.data$LC50[ind]))){next}
   if(is.na(as.numeric(h[[tolower(worm.LC.data$parent[ind])]]))){next}
   worm.LC.data$potency.rel.to.parent[ind] <- with(worm.LC.data, ifelse(!is.na(as.numeric(LC50[ind])),
                                                     formatC(signif(as.numeric(h[[tolower(parent[ind])]])/as.numeric(LC50[ind]),3), format="fg"), NA))
}

#determine if valid curve found
worm.LC.data <- worm.LC.data %>% mutate(curve.found = ifelse(is.na(fit), "No", ""))

#save data as a csv file
worm.LC.data %>% dplyr::select(-c(fit, curve.qual, LC75, LC75.CI, LC90, LC90.CI, min.conc)) %>% 
   write.csv(paste(out.path, LC.out.file, sep="\\"), row.names=F)

plots <- vector(mode = "list", length = 1+length(levels(as.factor(worm.LC.data$Compound))))
q$layers[[1]]$geom$default_aes$size <- 0.3
q$layers[[2]]$geom$default_aes$width <- 0.1
q$layers[[2]]$geom$default_aes$size <- 0.3
new.data <- q$data %>% mutate(ppm = runif(conc, min = (1-jitter.factor)*conc, max = (1+jitter.factor)*conc))
q$data$conc <- new.data$conc
q<- q+ theme(panel.grid.minor = element_line(size = 0.1, colour="grey85"), 
             panel.grid.major = element_line(size = 0.3, colour="grey85"),
             legend.key.height = unit(0.2, "cm"), 
             legend.position = c(0.8, 0.4), 
             legend.key.size = grid::unit(2, "lines"),
             legend.title = element_blank(),
             plot.title = element_text(hjust = 0.5))
q$theme$plot.title$size <- 8
q$theme$text$size <- 6

plots[[1]]<-q
i<-2
jitter.factor <- 0.1
for (m in levels(as.factor(worm.LC.data$Compound))) {
   save.name <-  str_replace_all(m, "[.]", "-") %>% str_replace_all("[%]", "")
   ind <- which(tolower(worm.LC.data$Compound)==tolower(m))
   

   if(is.na(worm.LC.data$fit[ind])){
      print(paste("Graph for ", m, " not printed; no fit found.", sep=""))
      next
   } else{
      print(paste("Making graph for ", m, ".",sep=""))
   }
   fit <- worm.LC.data$fit[ind]
   
   p<-plot_worm_data(fit[[1]],m)
   curr.plot.name <- paste(save.name, ".png", sep="")
   # ggsave(path = out.path, units = "in", width = 6.4, height = 4.134,
   #        filename = curr.plot.name, 
   #        plot=p, device = "png")
   p$layers[[1]]$geom$default_aes$size <- 0.3
   p$layers[[2]]$geom$default_aes$width <- 0.1
   p$layers[[2]]$geom$default_aes$size <- 0.3
   if(worm.LC.data$n.trials[ind]>1){ #jitter data if there are multiple trials
      new.data <- p$data %>% mutate(ppm = runif(ppm, min = (1-jitter.factor)*ppm, max = (1+jitter.factor)*ppm))
      p$data$ppm <- new.data$ppm
      p$layers[[2]]$data$ppm <- new.data$ppm
   }
   p$theme$plot.title$size <- 8
   p$theme$text$size <- 6
   p <- p +  theme(panel.grid.minor = element_line(size = 0.1, colour="grey85"), 
                   panel.grid.major = element_line(size = 0.3, colour="grey85"),
                   plot.title = element_text(hjust = 0.5)) +
      annotate(geom = 'text', label = paste("LD50:", worm.LC.data$LC50[ind], "ppm"), 
               x = worm.LC.data$min.conc[ind], y = 0.95, hjust = -0.05, vjust = 0.5, size = 1.75)
   plots[[i]] <- p
   i <- i+1
}
plots<-plots[lengths(plots) != 0] #remove empty info from array

#number of pages determined by integer division of number of plots plus three for large malathion image.
n.pages <- (length(plots)+3)%/%8 + if_else((length(plots)+3)%%8==0, 0, 1) 
grobs1 <- grobTree(gp = gpar(fontsize = 11), 
                   textGrob(label= "CONFIDENTIAL: Merlin Data", name="title1",
                            x = unit(1, "lines"), y = unit(0, "lines"), 
                            just = "left", hjust = 0, vjust = .1))

make.time <-format(Sys.time(), "%d %b %Y %H:%M:%S")

# ggsave(paste(out.path,"dose-response_curves1.pdf", sep="\\"),
#        marrangeGrob(grobs=plots[1:5], nrow=4, ncol=2, units = "in", 
#                     layout_matrix = matrix(c(1,1,1,1:5), 4, 2, TRUE), top= NULL, bottom=arrangeGrob(grobs1, grobs2, ncol = 2)),
#        width=7.5, height = 10)

ggsave(paste(tempdir(),"dose-response_curves1.pdf", sep="\\"),
       marrangeGrob(grobs=plots[1:5], nrow=4, ncol=2, units = "in", 
                    layout_matrix = matrix(c(1,1,1,1:5), 4, 2, TRUE), 
                    bottom=textGrob(label=paste(make.time, "\nPage 1 of ", n.pages, sep=""), name="title3",
                                                x = unit(2.5, "lines"), y = unit(0, "lines"),
                                                just = c("right", "bottom"), hjust = 0, vjust = -.6,
                                                gp = gpar(fontsize = 7)),
                    top=arrangeGrob(grobs1)),
              width=7.5, height = 10, colormodel = "cmyk")

p1.length <- pdf_length(paste(tempdir(),"dose-response_curves1.pdf", sep="\\"))
if(p1.length > 1){
   pdf_subset(paste(tempdir(),"dose-response_curves1.pdf", sep="\\"),
                             pages = 2, output = paste(tempdir(),"dose-response_curves1_temp.pdf", sep="\\"))
   file.rename(paste(tempdir(),"dose-response_curves1_temp.pdf", sep="\\"),
               paste(tempdir(),"dose-response_curves1.pdf", sep="\\"))
}


ggsave(paste(tempdir(),"dose-response_curves2.pdf", sep="\\"),
       marrangeGrob(grobs=plots[6:length(plots)], nrow=4, ncol=2, units = "in",
                    layout_matrix = matrix(c(1:8), 4, 2, TRUE), 
                    bottom=quote(grid::textGrob(label=paste(make.time, "\nPage ", g+1, " of ", npages+1, sep=""), name="title3",
                                                x = unit(2.5, "lines"), y = unit(0, "lines"),
                                                just = c("right", "bottom"), hjust = 0, vjust = -.6,
                                                gp = gpar(fontsize = 7))),
                    top=arrangeGrob(grobs1)),
       width=7.5, height = 10, colormodel = "cmyk")


pdf_combine(c(paste(tempdir(),"dose-response_curves1.pdf", sep="\\"), 
              paste(tempdir(),"dose-response_curves2.pdf", sep="\\")), 
            output = paste(out.path,"Merlin_DR-curves.pdf", sep="\\"))

do.call(file.remove, list(paste(tempdir(),"dose-response_curves1.pdf", sep="\\"),paste(tempdir(),"dose-response_curves2.pdf", sep="\\")))

