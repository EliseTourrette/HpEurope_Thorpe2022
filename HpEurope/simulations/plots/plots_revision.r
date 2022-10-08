## last version of the plots
## figure 4 of the HpEurope paper
## put all panels together, in one figure


setwd("~/Desktop/HpEurope/article_repeatedOutAfricaHPylori")

library(reshape2)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(ggpubr)

selCoeff <- -c(0, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4)
bottleneck <- c(5000, 5500)
admixture <- 8000
rec <- c('0bp', '500bp', '5000bp', '50000bp')


########### FIGURE 4 ###########

############# 4.A ###########

irec <- 3

load("ANALYSIS/DATA/pairwiseDiff.RData")
data <- data[data$generation == 8000 & data$rec == rec[irec] & data$rep == 1, ]

selCoeff <- c(0, -0.005, -0.002, -0.001, -0.0005, -0.0002, -0.0001, 0, -0.005, -0.002, -0.001, -0.0005, -0.0002, -0.0001)
mutType <- paste0("m", 1:length(selCoeff))
for (i in 1:length(mutType)) {
    data$s[data$mutType == mutType[i]] <- selCoeff[i]
}

data1 <- data[data$mutType %in% mutType[1:7], -c(1:3, 5)]
data2 <- data[data$mutType %in% mutType[8:14], -c(1:3, 5)]

sub <- data.frame(pop = data1$pop, s = data1$s, data1[, 2:4951] + data2[, 2:4951])

data <- NULL
for (i in 2:7) {
    data <- rbind(data, data.frame(ratio = unlist(sub[i, 3:ncol(sub)]) / unlist(sub[1, 3:ncol(sub)]), neutral = unlist(sub[1, 3:ncol(sub)]), s = sub[i, 2], pop = sub[i, 1]))
    data <- rbind(data, data.frame(ratio = unlist(sub[i + 7, 3:ncol(sub)]) / unlist(sub[8, 3:ncol(sub)]), neutral = unlist(sub[8, 3:ncol(sub)]), s = sub[i + 7, 2], pop = sub[i + 7, 1]))
}

pathwd <- paste0("20201105_", irec, "_1")

data <- NULL
data <- rbind(data, data.frame(ratio = unlist(apply(sub[2:7, 3:ncol(sub)], 2, sum)) / unlist(sub[1, 3:ncol(sub)]), neutral = unlist(sub[1, 3:ncol(sub)]), pop = 1))
data <- rbind(data, data.frame(ratio = unlist(apply(sub[9:14, 3:ncol(sub)], 2, sum)) / unlist(sub[8, 3:ncol(sub)]), neutral = unlist(sub[8, 3:ncol(sub)]), pop = 2))

data2 <- aggregate(data[c("ratio", "neutral")], by = list(pop = data$pop), mean, na.rm = TRUE)

p <- ggplot(data) +
    geom_point(mapping = aes(x = neutral, y = ratio, col = as.factor(pop)), alpha = 0.1, shape = 1, size = 0.5) +
    geom_point(data2,mapping = aes(x = neutral, y = ratio, fill = as.factor(pop)),alpha = 1,shape = 21,size = 0.75) +
    labs(x = "dS",y = "dN/dS",color = "Population",fill = "Population",linetype = " ", tag = "a") +
    scale_colour_manual(values = c("forestgreen", "Darkgoldenrod1"),labels = c("non-bottlenecked \npopulation \n","bottlenecked \npopulation")) +
    scale_fill_manual(values = c("forestgreen", "Darkgoldenrod1"),labels = c("non-bottlenecked \npopulation \n","bottlenecked \npopulation"))
pa <- p + theme_cowplot()  + theme(legend.position = "none",text = element_text(size = 7),axis.text = element_text(size = 7))
#ggsave(paste0("4B.jpg"),width = 8.5,height = 8.5,units = "in",dpi = 500,scale = 0.75)



############# 4.B ###########

i <- 3

pathwd <- paste0('20201105_', i, '_1')
load(file = paste0('ANALYSIS/DATA/dNdS/dNdSancestorSumSelCoeff_', pathwd, '_g8000.RData'))

data2 <- aggregate(data[c("ratio", "neutral")], by = list(pop = data$pop), mean)

p <- ggplot(data) +
    geom_point(mapping = aes(x = neutral, y = ratio, col = as.factor(pop)),alpha = 0.1,shape = 1,size = 0.5) +
    geom_point(data2,mapping = aes(x = neutral, y = ratio, fill = as.factor(pop)),alpha = 1,shape = 21,size = 0.75) +
    labs(x = "dS",y = "dN/dS",color = "Population", fill = "Population",linetype = " ", tag = "b") +
    scale_colour_manual(values = c("forestgreen", "Darkgoldenrod1"),labels = c("non-bottleneck \npopulation", "bottleneck \npopulation")) +
    scale_fill_manual(values = c("forestgreen", "Darkgoldenrod1"),labels = c("non-bottleneck \npopulation", "bottleneck \npopulation"))
pb <- p + theme_cowplot() + theme(legend.position = "none", text = element_text(size = 7), axis.text = element_text(size = 7))
#ggsave(paste0("4C.jpg"),width = 8.5,height = 8.5,units = "in",dpi = 500,scale = 0.75)




############# 4.F ###########

g <- 8300
irow <- 1
rep <- 1
cutOff <- 1000

irec <- 3

pathwd <- paste0("20201105_", irec, "_", rep)
load(paste0("ANALYSIS/DATA/ancestry_g15000/",pathwd,"/sampleData.RData"))

load(paste0("ANALYSIS/DATA/20211021-ancestry/dataAncestry",cutOff,"_",pathwd,".RData"))
anc <- data

anc[4:5, -(1:2)] <- anc[4:5, -(1:2)] / 10
## to get the ancestry in percent
## in generation 15000, corresponding to lines 5 to 7, the sample size is 1000

anc <- anc[, -(1:2)]

pos <- unlist(anc[irow + 2, ])
x1 <- unlist(anc[irow, ])
x2 <- unlist(anc[irow + 1, ])

data <- cbind(sub,ancestry1 = unlist(lapply(sub$pos, function(x) x1[pos == x])),ancestry2 = unlist(lapply(sub$pos, function(x) x2[pos == x])),pos2 = unlist(lapply(sub$pos, function(x) pos[pos == x])))

data$ratio1 <- data$ancestry1 / (data$ancestry1 + data$ancestry2)

## average over the delta fr bins
bins = c(-1.01,-0.99,-0.75,-0.5,-0.25,-0.01,0.01,0.25,0.5,0.75,0.99,1.01)
#bins = seq(from = -1.01, to = 1.01, by = 0.1)
data$bins <- NA
for (i in 2:length(bins)) {
    data$bins[data$fr > bins[i - 1] &data$fr <= bins[i]] <- mean(c(bins[i - 1], bins[i]))
}

p <-ggplot(data, aes(x = -bins,y = 1 - ratio1,col = as.factor(mtype))) +
    stat_summary(fun = mean, geom = "point", size = 0.5) +
    stat_summary(fun = mean, geom = "line") +
    stat_summary(fun.data = mean_se,fun.args = list(mult = 1), geom = "errorbar", width = 0.05) +
    geom_point(mapping = aes(x = -fr,y = 1 - ratio1,col = factor(mtype)),alpha = 0.1, size = 0.25) +
    #expand_limits(y = c(0, 1.1)) +
    labs(x = "Mutation score", y = " Bottleneck population \nancestry", color = "Mutation type", tag = "f") #+
    #geom_segment(aes(x = -0.01,y = 1.2,xend = -0.1,yend = 1.2),
    #             col = "black",
    #             arrow = arrow(length = unit(0.05, "in"))) +
    #geom_text(x=-0.6, y=1.2, label="Excess \nbottleneck mutations", col = "black", size = 1.5) +
    #geom_segment(aes(x = 0.01,y = 1.2,xend = 0.1,yend = 1.2),
    #             col = "black",
    #             arrow = arrow(length = unit(0.05, "in"))) +
    #geom_text(x=0.58, y=1.2, label="Excess \nnon-bottleneck mutations", col = "black", size = 1.5)
pf <- p + theme_cowplot() + theme(text = element_text(size = 7),axis.text = element_text(size = 7))
#ggsave(paste0("4D.jpg"),width = 10,height = 8.5,units = "in",dpi = 500,scale = 0.75)




############# 4.C ###########

load("./ANALYSIS/DATA/fitnessOverall.RData")
colnames(data) <- c("rec", "rep", "pop", "generation", "mean", "var")

data$varnorm <- sqrt(data$var) / data$mean

#data <- data[-which(data$rec == "8e+05bp"),]

data$rec <- as.numeric(unlist(strsplit(as.character(data$rec), split = "bp")))
data$rec <- paste("import size per generation: ", data$rec, "bp")

p <- ggplot(data[data$rec == unique(data$rec)[3],]) +
    stat_summary(mapping = aes(x = generation, y = mean, col = as.factor(pop)),fun.data = mean_se,geom = "line") +
    #stat_summary(fun.data = mean_se, fun.args = list(mult = 1.96), geom = "ribbon", alpha = 0.2) +
    #facet_wrap(~ rec) +
    labs(x = "generation", y = "fitness", color = "Population", tag = "c") +
    scale_colour_manual(values = c("forestgreen", "Darkgoldenrod1"),labels = c("non-bottleneck \npopulation \n","bottleneck \npopulation")) +
    geom_segment(data = data.frame(admixture = admixture), aes(x = bottleneck[1], y = 0.01 , xend = bottleneck[2], yend = 0.01)) +
    geom_segment(data = data.frame(admixture = admixture), aes(x = admixture, y = 0.01 , xend = admixture, yend = 0.2),
                 arrow = arrow(length = unit(0.05, "in")))
pc <- p + theme_cowplot() + theme(legend.justification = c(0.9, 0),legend.position = c(1, 0.6),text = element_text(size = 7),axis.text = element_text(size = 7))
#ggsave("4.A.jpg",width = 8.5,height = 8.5,units = "in",dpi = 500,scale = 0.75)

############# 4.D ###########

ancestryString <- c("bottleneck", "unknown", "non_bottleneck")
ancestryColor <- c("Darkgoldenrod1", "grey", "forestgreen")
ancestryLab <- c("bottleneck", "unknown", "non-bottleneck")

load("REVISION/simu10/analysesAncestry/ANALYSIS/proportionInterRec.RData")


irec <- "interRec"

p <- ggplot(proportionTot[!(proportionTot$population %in% c("pre_split")) & proportionTot$rec == irec,]) +
    geom_point(aes(x = G, y = proportion, color = population),size = 0.5) +
    scale_color_manual(name = "Ancestry", breaks = ancestryString, values = ancestryColor, labels = ancestryLab) +
    geom_segment(data = data.frame(admixture = seq(from = 8000, to = 15000, by = 500)), aes(x = admixture, y = 0 , xend = admixture, yend = 5), arrow = arrow(length = unit(0.05, "in"))) +
    labs(x = "generation", y = "Ancestry proportion", tag = "d") +
    guides(color=guide_legend(nrow=1, byrow=TRUE))
pd <- p + theme_cowplot() + theme(legend.spacing = unit(.001, 'mm'), legend.justification = c(0,0), legend.direction = "horizontal",legend.position = c(0,1), panel.background = element_rect(fill = 'white', colour = 'white'), text = element_text(size = 7),axis.text = element_text(size = 7)) 
    #ggsave(paste0("proportion_", irec, ".jpeg"),width = 8.5,height = 8.5,units = "in",dpi = 500,scale = 0.75)


############# 4.E ###########

plotList <- vector("list",4)
iplot <- 1

for (ig in c(8000, 8100, 8200, 8300)) {
    
    
    print(ig)
    
    nameFile <- paste0("genome2_", ig, ".txt") ## name of the file you will look at (number = generation)
    parameters <- read.table("ANALYSIS/painting/parameters.csv", stringsAsFactor = FALSE, header = TRUE, sep = ",")
    
    nrow <- as.numeric(parameters$nrow[parameters$nameFile == nameFile])## number of rows of the file; get it from parameters.csv; can also get it with the "system" function 
    nameFile <- paste0("ANALYSIS/painting/genome2_", ig, ".txt")
    ncol <- as.numeric(system(paste("awk '{print NF}'", nameFile, "| sort -nu | tail -n 1"), intern = TRUE)) ## maximum number of column in the file

    
    ## load the population data
    pop <- read.table(nameFile, stringsAsFactor = FALSE, skip = nrow - 100, col.names = 1:ncol, fill = TRUE)
    pop <- as.matrix(pop[,-(1:2)])
    ## pop: matrix, rows = strains, columns = mutations contained by each strains
    
    ## load the mutation data
    mut <- read.table(nameFile,stringsAsFactor = FALSE, skip = 2,nrows = nrow - 100 - 3, fill = TRUE )
    mut <- mut[order(mut$V4),]
    ## mut: data.frame, rows = mutations, column = different information about the given mutations
    
    ##if want to look at a portion of the genome (or the whole genome), independently on whether a site has a mutation or not
    nbSite <- 1600000 ## if want to look at the whole genome, nbSite = 1600000
    ancestry <- matrix(0, nrow = nrow(pop), ncol = nbSite) ## if want to look at a portion of the chromosomes of nbSite sites
    posBegin <- 1 ## first position, in bp of the interval looked at
    
    ##fill the ancestry matrix
    for (i in 1:nbSite) {
        ## for-loop on each site
        j <- i + posBegin - 1 ## get the physical position of each site = index of the site in the interval plus the first position of the interval (minus 1, else would add 1 to every position)
        id <- mut$V1[mut$V4 == j] ## get the (within-file) id of the mutation(s) at this position (there can be more than one or none at all)
        if (length(id) > 0) {
            ## do the following part only if there is at least one mutation at this position
            #print("mut!")
            for (k in id) {
                ## as there can be more than mutation at this position, loop on the possible mutation
                ind <- which(pop == k) %% 100 ## get the index of the strains with the mutation (use the modulo 100 as there are 100 strains in the sample; would need another way if the number of strains in the sample was different)
                ind[ind == 0] <- 100 ## for the element that are a multiple of 100, the modulo will result in 0; thus need to change the strain index back to 100 (could also do ind <- ind + 1, as we need an index that begin at 1)
                # fitness[ind,i] <- 1 + mut$V5[mut$V1 == k]
                ##fill in the ancestry matrix depending on the mutation
                if (mut$V7[mut$V1 == k] == "p1" & mut$V8[mut$V1 == k] < 5000)
                    ancestry[ind, i] <- 1
                if (mut$V7[mut$V1 == k] == "p1" & mut$V8[mut$V1 == k] >= 5000 & mut$V8[mut$V1 == k] < 8000)
                    ancestry[ind, i] <- 2
                if (mut$V7[mut$V1 == k] == "p1" & mut$V8[mut$V1 == k] >= 5000 & mut$V8[mut$V1 == k] >= 8000)
                    ancestry[ind, i] <- 3
                if (mut$V7[mut$V1 == k] == "p2" & mut$V8[mut$V1 == k] >= 5000 & mut$V8[mut$V1 == k] < 8000)
                    ancestry[ind, i] <- 4
                if (mut$V7[mut$V1 == k] == "p2" & mut$V8[mut$V1 == k] >= 5000 & mut$V8[mut$V1 == k] >= 8000)
                    ancestry[ind, i] <- 5
            }
        }
    }
    positionSite <- posBegin:(posBegin + nbSite - 1) ## vector of position, in bp, of continuous sites
    
    ## initialize the selcoeff matrix with 0, which is the default for a site in a strain without mutation
    selcoeff <- matrix(0, nrow = nrow(pop), ncol = length(unique(mut$V4)))
    ## put the selective coefficient in selcoeff matrix, rows = strains and columns = site
    for (i in 1:length(unique(mut$V4))) {
        ## for site which has a mutation; for-loop on a index
        j <- unique(mut$V4)[i] ## get the physical position
        id <- mut$V1[mut$V4 == j] ## get the (within-file) id of the mutation(s) at this position (there can be more than one)
        for (k in id) {
            ## as there can be more than mutation at this position, loop on the possible mutation
            ind <- which(pop == k) %% 100 ## get the index of the strains with the mutation (use the modulo 100 as there are 100 strains in the sample; would need another way if the number of strains in the sample was different)
            ind[ind == 0] <- 100 ## for the element that are a multiple of 100, the modulo will result in 0; thus need to change the strain index back to 100 (could also do ind <- ind + 1, as we need an index that begin at 1)
            ## fill in the selcoeff matrix depending on the mutation type(mut$3)
            selcoeff[ind, i] <- mut$V5[mut$V1 == k]
        }
    }
    
    ## matrix of fitness; fitness=1+selective coefficient
    fitness <- selcoeff + 1
    
    ## if want to look at the fitness for the different classes of mutations
    # components <- apply(fitness, 1, table) ## first, for each strain, calculate the number of fitness of each class of mutations
    # effect <- as.numeric(attributes(components)$dimnames[[1]]) ## the get the fitness effect (1+s) of each class of mutations present
    # fitnessPerClass <- effect^components ## then calculate the fitness per strain per clas of mutations (which is simply effect^(nbr mutations of this class) as the fitness components are multiplied between the different sites)
    
    ##fitness_i = mult{from j = 1 to N} (1 + s_ij), with i = strain, j = site (N = number of site) and s_ij the selection coeff of strain i at site j
    fitnessTotal <- apply(fitness, 1, prod)## here, calculate the total fitness of each strain
    
    ## for the test, sample 10 strains among the 100 strains
    Nstrains <- c(1, 11, 21, 31, 41, 51, 61, 71, 81, 91)
    
    ancestry10 <- ancestry[Nstrains,]
    fitness10 <- fitnessTotal[Nstrains]
    
    widthBlock <- apply(ancestry10, 2, sum) ## if the value > 0, then there is at least one mutation at this site for this subset of strain
    widthBlock[widthBlock > 0] <- 1 ## only use 0 and 1
    
    anc10 <- melt(ancestry10) ## need to reshape the matrix to get a data.frame with col1 = strain, col2 = site, col3 = ancestry (note that you could directly put it in this format instead of the ancestry matrix
    ## add another variable to anc10 saying if there is a mutation at the site or not
    anc10$Var3 <- rep(widthBlock, each = length(Nstrains)) ## use each as we first look at every strains for a position; would use times if we were looking first at every site for one strain
    ## replace the site number by their physical positions
    anc10$Var2 <- rep(positionSite, each = length(Nstrains)) ## use each as we first look at every strains for a position; would use times if we were looking first at every site for one strain
    
    ## want to look at the global proportion mutations non bottleneck/bottleneck population for the mutations post split / pre admixture
    origin10 <- NULL
    for (i in 1:length(Nstrains)) {
        origin <- anc10$value[anc10$Var1 == i]
        origin10 <- c(origin10, length(origin[origin == 4]) / (length(origin[origin == 2])+length(origin[origin == 4])))
    }
    print(max(origin10))
    
    anc10 <- anc10[anc10$Var2 <= 20000,] ## only look at the first 20000 sites
    
    lab <- c( "ancestral population \npre split \n",
              "ancestral population \npost split \npre admixture \n",
              "ancestral population \npost split \npost admixture \n",
              "bottleneck population \npost split \npre admixture \n",
              "bottleneck population \npost split \npost admixture \n") ## possible labels for the ancestry
    colLab <- c("1" = "green4",  "2" = "red", "3" = "orange", "4" = "blue", "5" = "purple") ## color to be used for the different labels; same color scheme
    

    if(iplot == 1) {
        
        plot1 <- ggplot(data = anc10[anc10$value != 0,], aes(x = Var2, y = as.factor(Var1))) +
            scale_fill_manual(labels = lab, values = colLab, drop = FALSE) +
            labs(title = "At admixture") +
            geom_tile(aes(fill = as.factor(value), height = 0.7, width = Var3)) +  ## add the plus if want to add one of the geom_point ## with width = Var3, you only put a bar when there is a mutation
            # geom_raster(aes(fill = as.factor(value))) +
            #geom_point(data = s10, aes(x = Var2, y = as.factor(Var1), pch = as.factor(round(value,4)))) ## for each element of the ancestry matrix, will add a point whose shape will depend on the value of the selection coefficient
            # geom_point(data = s10, aes(x = Var2, y = as.factor(Var1), size = as.factor(abs(round(value,4)))))  ## for each element of the ancestry matrix, will add a point whose size will depend on the value of the selection coefficient
            theme_set(theme_bw()) + theme(
                panel.grid = element_line(colour = NA),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text = element_blank(),
                #axis.text.x = element_text(size = rel(0.5), angle = 45),
                axis.line.x = element_line(),
                panel.border = element_blank(),
                axis.ticks.y = element_blank(),
                text = element_text(size = 7)) + ##get rid of the grey background and grid lines
            guides(fill = "none")
        
        plot2 <-
            ggplot(data = data.frame(strain = 1:length(fitness10), fitness = fitness10), aes(x = fitness, y = as.factor(strain))) +
            geom_point(size = 0.5) +
            theme_set(theme_bw()) + theme(
                panel.grid = element_line(colour = NA),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text = element_blank(),
                #axis.text.x = element_text(size = rel(0.5), angle = 45),
                panel.border = element_blank(),
                axis.line = element_line(),
                text = element_text(size = 7)) + ##get rid of the grey background and grid lines
            xlim(0.31, 0.37)
        
        plot3 <- ggplot(data = data.frame(strain = 1:length(origin10), ratio = origin10), aes(x = ratio, y = as.factor(strain))) +
            geom_point(size = 0.5) +
            theme_set(theme_bw()) + theme(
                panel.grid = element_line(colour = NA),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text = element_blank(),
                #axis.text.x = element_text(size = rel(0.5), angle = 45),
                panel.border = element_blank(),
                axis.line = element_line(),
                text = element_text(size = 7)) + ##get rid of the grey background and grid lines
            xlim(0,1)
        
    }
    
    if(iplot == 2) {
        
        plot1 <- ggplot(data = anc10[anc10$value != 0,], aes(x = Var2, y = as.factor(Var1))) +
            labs(title = "After 100 generations", fill = "population and \nperiod of origin") +
            scale_fill_manual(labels = lab, values = colLab, drop = FALSE) +
            geom_tile(aes(fill = as.factor(value), height = 0.7, width = Var3)) +  ## add the plus if want to add one of the geom_point ## with width = Var3, you only put a bar when there is a mutation
            # geom_raster(aes(fill = as.factor(value))) +
            #geom_point(data = s10, aes(x = Var2, y = as.factor(Var1), pch = as.factor(round(value,4)))) ## for each element of the ancestry matrix, will add a point whose shape will depend on the value of the selection coefficient
            # geom_point(data = s10, aes(x = Var2, y = as.factor(Var1), size = as.factor(abs(round(value,4)))))  ## for each element of the ancestry matrix, will add a point whose size will depend on the value of the selection coefficient
            theme_set(theme_bw()) + theme(
                panel.grid = element_line(colour = NA),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text = element_blank(),
                #axis.text.x = element_text(size = rel(0.5), angle = 45),
                axis.line.x = element_line(),
                panel.border = element_blank(),
                axis.ticks.y = element_blank(),
                text = element_text(size = 7))  ##get rid of the grey background and grid lines
            
        
        plot2 <-
            ggplot(data = data.frame(strain = 1:length(fitness10), fitness = fitness10), aes(x = fitness, y = as.factor(strain))) +
            geom_point(size = 0.5) +
            theme_set(theme_bw()) + theme(
                panel.grid = element_line(colour = NA),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text = element_blank(),
                #axis.text.x = element_text(size = rel(0.5), angle = 45),
                panel.border = element_blank(),
                axis.line = element_line(),
                text = element_text(size = 7)) + ##get rid of the grey background and grid lines
            xlim(0.31, 0.37)
        
        plot3 <- ggplot(data = data.frame(strain = 1:length(origin10), ratio = origin10), aes(x = ratio, y = as.factor(strain))) +
            geom_point(size = 0.5) +
            theme_set(theme_bw()) + theme(
                panel.grid = element_line(colour = NA),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text = element_blank(),
                #axis.text.x = element_text(size = rel(0.5), angle = 45),
                panel.border = element_blank(),
                axis.line = element_line(),
                text = element_text(size = 7)) + ##get rid of the grey background and grid lines
            xlim(0,1)
        
    }
    
    if(iplot == 3) {
        
        plot1 <- ggplot(data = anc10[anc10$value != 0,], aes(x = Var2, y = as.factor(Var1))) +
            labs(title = "After 200 generations", fill = "population and \nperiod of origin") +
            scale_fill_manual(labels = lab, values = colLab, drop = FALSE) +
            geom_tile(aes(fill = as.factor(value), height = 0.7, width = Var3)) +  ## add the plus if want to add one of the geom_point ## with width = Var3, you only put a bar when there is a mutation
            # geom_raster(aes(fill = as.factor(value))) +
            #geom_point(data = s10, aes(x = Var2, y = as.factor(Var1), pch = as.factor(round(value,4)))) ## for each element of the ancestry matrix, will add a point whose shape will depend on the value of the selection coefficient
            # geom_point(data = s10, aes(x = Var2, y = as.factor(Var1), size = as.factor(abs(round(value,4)))))  ## for each element of the ancestry matrix, will add a point whose size will depend on the value of the selection coefficient
            theme_set(theme_bw()) + theme(
                panel.grid = element_line(colour = NA),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text = element_blank(),
                #axis.text.x = element_text(size = rel(0.5), angle = 45),
                axis.line.x = element_line(),
                panel.border = element_blank(),
                axis.ticks.y = element_blank(),
                text = element_text(size = 7))  ##get rid of the grey background and grid lines
        
        
        plot2 <-
            ggplot(data = data.frame(strain = 1:length(fitness10), fitness = fitness10), aes(x = fitness, y = as.factor(strain))) +
            geom_point(size = 0.5) +
            theme_set(theme_bw()) + theme(
                panel.grid = element_line(colour = NA),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text = element_blank(),
                #axis.text.x = element_text(size = rel(0.5), angle = 45),
                panel.border = element_blank(),
                axis.line = element_line(),
                text = element_text(size = 7)) + ##get rid of the grey background and grid lines
            xlim(0.31, 0.37)
        
        plot3 <- ggplot(data = data.frame(strain = 1:length(origin10), ratio = origin10), aes(x = ratio, y = as.factor(strain))) +
            geom_point(size = 0.5) +
            theme_set(theme_bw()) + theme(
                panel.grid = element_line(colour = NA),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text = element_blank(),
                #axis.text.x = element_text(size = rel(0.5), angle = 45),
                panel.border = element_blank(),
                axis.line = element_line(),
                text = element_text(size = 7)) + ##get rid of the grey background and grid lines
            xlim(0,1)
        
    }
    
    if(iplot == 4) {
        
        plot1 <- ggplot(data = anc10[anc10$value != 0,], aes(x = Var2, y = as.factor(Var1))) +
            labs(title = "After 300 generations", fill = "population and \nperiod of origin") +
            xlab("Physical position (bp)") +
            scale_fill_manual(labels = lab, values = colLab, drop = FALSE) +
            geom_tile(aes(fill = as.factor(value), height = 0.7, width = Var3)) +  ## add the plus if want to add one of the geom_point ## with width = Var3, you only put a bar when there is a mutation
            # geom_raster(aes(fill = as.factor(value))) +
            #geom_point(data = s10, aes(x = Var2, y = as.factor(Var1), pch = as.factor(round(value,4)))) ## for each element of the ancestry matrix, will add a point whose shape will depend on the value of the selection coefficient
            # geom_point(data = s10, aes(x = Var2, y = as.factor(Var1), size = as.factor(abs(round(value,4)))))  ## for each element of the ancestry matrix, will add a point whose size will depend on the value of the selection coefficient
            theme_set(theme_bw()) + theme(
                panel.grid = element_line(colour = NA),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.line.x = element_line(),
                panel.border = element_blank(),
                axis.ticks.y = element_blank(),
                text = element_text(size = 7))   ##get rid of the grey background and grid lines
           
        
        plot2 <-
            ggplot(data = data.frame(strain = 1:length(fitness10), fitness = fitness10), aes(x = fitness, y = as.factor(strain))) +
            geom_point(size = 0.5) +
            xlab("Fitness") +
            theme_set(theme_bw()) + theme(
                panel.grid = element_line(colour = NA),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.text.x = element_text(size = 5),
                panel.border = element_blank(),
                axis.line = element_line(),
                text = element_text(size = 7)) + ##get rid of the grey background and grid lines
            xlim(0.31, 0.37)
        
        plot3 <- ggplot(data = data.frame(strain = 1:length(origin10), ratio = origin10), aes(x = ratio, y = as.factor(strain))) +
            geom_point(size = 0.5) +
            theme_set(theme_bw()) + theme(
                panel.grid = element_line(colour = NA),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.text.x = element_text(size = 2, angle = 45),
                panel.border = element_blank(),
                axis.line = element_line(),
                text = element_text(size = 7)) + ##get rid of the grey background and grid lines
            xlim(0,1)
        
    }
    
   
    
    #plotTot <- plot_grid( plot2, NULL,  plot3, NULL, plot1,  nrow = 1, rel_widths = c(1, -0.01, 1, -0.01, 4),  align = "hv")
    
    # now add the title
    #title <- ggdraw() + draw_label(paste0("g = ", ig), size = 12)
    #plotList[[iplot]] <- plot_grid(title, plotTot,  ncol = 1, rel_heights = c(1, 9))
    plotList[[iplot]] <- vector("list",3)
    plotList[[iplot]][[1]] <- plot2
    plotList[[iplot]][[2]] <- plot3
    plotList[[iplot]][[3]] <- plot1
    iplot <- iplot + 1
    
    # save_plot(paste0("generation8000.pdf"), p, base_height = 1.8)
    
    #save_plot(paste0("generation", ig, ".pdf"), p, base_height = 2.5,  dpi = 2000)
    
}

pe1 <- ggarrange(plotList[[1]][[3]], plotList[[2]][[3]], plotList[[3]][[3]], plotList[[4]][[3]], ncol = 1, common.legend = TRUE, legend = "right")
pe3 <- ggarrange(plotList[[1]][[1]], plotList[[2]][[1]], plotList[[3]][[1]], plotList[[4]][[1]], ncol = 1)
pe <- plot_grid(pe3, NULL, pe1, nrow = 1, rel_widths = c(1, -0.01, 4),  align = "hv", labels = "e", label_size = 7)

############# 4 ###########
## put all the subplots together

ptot <- ggdraw() +
    draw_plot(pd, 0.01, 0.01, .43, .23, scale = 1) +
    draw_plot(pc, 0.01, .26, .43, .23, scale = 1) +
    draw_plot(pb, 0.01, .51, .43, .23, scale = 1) +
    draw_plot(pa, 0.01, .76, .43, .23, scale = 1) +
    draw_plot(pf, .46, 0.01, .53, .23, scale = 1) +
    draw_plot(pe, .46, .26, .53, .73, scale = 1)
ggsave(ptot, file = "Fig4.pdf", width = 180, height = 185, units = "mm")


