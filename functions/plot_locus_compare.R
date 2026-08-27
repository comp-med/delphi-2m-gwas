###########################################
## function to plot stacked locus zoom plot
## colouring for different variants 

plot.locus.compare <- function(sum.stat, sum.coloc, embQTL, ld_mat){
  
  ## 'sum.stat'  -- summary statistics for plotting
  ## 'sum.coloc' -- summary from coloc
  ## 'embQTL'    -- SNP to be plotted
  ## 'ld_mat'    -- LD matrix as imported with the help of REGENIE
  
  ## package for faster data reading
  require(data.table)
  
  #-----------------------#
  ##--  perpare input  --##
  #-----------------------#
  
  ## convert to data frame if needed
  sum.stat                <- as.data.frame(sum.stat)
  ## dummy
  a.vars                  <- sum.coloc[ type.snp == "embQTL"]$ID
  print(a.vars)
  
  ## compute estimates LOG10P value for GWAS traits to avoid underflow
  sum.stat$log10p         <- -pchisq((as.numeric(sum.stat$BETA.drug)/sum.stat$SE.drug)^2, df=1, lower.tail=F, log.p=T)/log(10)
  
  ## create a couple of colour gradients
  # cls <- c("#5F4690", "#1D6996", "#38A6A5", "#0F8554", "#73AF48", "#EDAD08", "#E17C05", "#CC503E", "#94346E", "#6F4070", "#994E95")
  cls <- c("#E58606", "#5D69B1", "#52BCA3", "#99C945", "#CC61B0", "#24796C", "#DAA51B", "#2F8AC4", "#764E9F", "#ED645A", "#CC3A8E", "#A5AA99", "#00A4CC", "#F95700", rainbow(10))
  # cls <- c("#7F3C8D", "#11A579", "#3969AC", "#F2B701", "#E73F74", "#80BA5A", "#E68310", "#008695", "#CF1C90", "#f97b72", "#4b4b8f", "#A5AA99")
  
  ## create colouts list
  col.list <- lapply(1:20, function(x){
    colorRampPalette(c( "white", cls[x]))(101)
  })
  
  ## --> add LD column for SNPs to be added <-- ##
  
  for(j in embQTL){
    ## obtain from LD matrix
    sum.stat[, paste0("R2.", j)] <- ld_mat[sum.stat$ID, j]
  }
  
  print(head(sum.stat))
  
  #------------------------------------#
  ##--        oppose p-values       --##
  #------------------------------------#
  
  cat("\n plot oppsing p-value \n")
  
  plot(log10p ~ LOG10P, 
       xlim=c(0,max(sum.stat$LOG10P)*1.05),
       data=sum.stat,
       cex=.4,
       col="grey90", 
       xlab=bquote(.(sum.coloc$embedding[1])~"GWAS"~~-log[10]("p-value")),
       ylab=bquote(.(sum.coloc$disease[1])~"GWAS"~~-log[10]("p-value"))
  )
  
  ## add variants to be highlighted
  for(j in 1:length(a.vars)){
    
    cat("\n plot LD with ", a.vars[j], "\n")
    
    ## rsidentify all up to R2>.3
    tmp <- sum.stat[which(sum.stat[, paste0("R2.",embQTL)] >= .3),]
    print(summary(sum.stat[, paste0("R2.",embQTL)]))
    ## add points
    points(tmp$LOG10P, tmp$log10p, cex=.5, pch=21, col="grey10", 
           bg=col.list[[j]][ceiling(tmp[, paste0("R2.",embQTL)]*100)],
           lwd=.2)
    
    ## add lead variant as diamond
    tmp <- subset(sum.stat, ID == a.vars[j])
    print(tmp)
    # print(tmp)
    points(tmp$LOG10P, tmp$log10p, pch=23, lwd=.8, cex=.7, 
           col=col.list[[j]][ceiling(tmp[, paste0("R2.",embQTL)]*100)], bg="white", 
           type="p")
    
    ## annotate
    text(tmp$LOG10P, tmp$log10p*1.05, 
         labels=a.vars[j], cex=.4, pos=4, xpd=NA, offset=.15)
    # print(tmp)
    ## add arrow to combine both
    arrows(tmp$LOG10P, tmp$log10p, 
           tmp$LOG10P, tmp$log10p*1.05,
           lwd=.5, length = 0, xpd=NA)
  }
  
  ## add legend for Coloc
  legend(ifelse(sum.coloc[1, "PP.H4.abf"] > .5, "topleft", "topright"), lty=0, pch=NA, cex=.5,
         legend=paste(paste0("H", 1:4), "=", sprintf("%.1f", sum.coloc[1, paste0("PP.H",1:4, ".abf")]*100)),
         title="Posterior prob. [%]")
  
  #------------------------------------#
  ##--  regional association - OUT  --##
  #------------------------------------#
  
  ## plot the Seer
  par(mar=c(.1,1.5,1,.5))
  
  ## add all variants
  plot(sum.stat[, "BP"], sum.stat[, "log10p"], xlab="", ylab=expression(-log[10]("p-value")),
       ylim=c(0, max(sum.stat[, "log10p"], na.rm=T)*1.1), cex=.3, col="grey90",
       xaxt="n", yaxt="n")
  axis(2, lwd=.5)
  
  ## rsidentify all up to R2>.3
  tmp <- sum.stat[which(sum.stat[, paste0("R2.",embQTL)] >= .3),]
  ## add points
  points(tmp[, "BP"], tmp[, "log10p"], cex=.5, pch=21, col="grey10", 
         bg=col.list[[j]][ceiling(tmp[, paste0("R2.",embQTL)]*100)],
         lwd=.2)
  
  ## add legend
  legend("topleft", lty=0, pch=NA, cex=.5, bty="n", legend = sum.coloc$disease[1])
  
  ## plotting coordinates
  pm <- par("usr")
  
  ## add colour gradient legend
  l  <- seq(pm[1]+(pm[2]-pm[1])*.75, pm[1]+(pm[2]-pm[1])*.95, length.out = 100)
  ## rectangle for the colours
  rect(l-(l[2]-l[1])/2, pm[3]+(pm[4]-pm[3])*(.93-(j*.08)), l+(l[2]-l[1])/2, pm[3]+(pm[4]-pm[3])*(1-(j*.08)), border=NA, col=col.list[[1]])
  ## box
  rect(l[1]-(l[2]-l[1])/2, pm[3]+(pm[4]-pm[3])*(.93-(j*.08)), l[100]+(l[2]-l[1])/2, pm[3]+(pm[4]-pm[3])*(1-(j*.08)), border="black", col=NA, lwd=.3)
  ## add header
  text(pm[1]+(pm[2]-pm[1])*.75, pm[3]+(pm[4]-pm[3])*(.97-(1*.08)), cex=.4, labels = a.vars[1], pos=4,
       offset = .2)
  ## genetic variant
  text(pm[1]+(pm[2]-pm[1])*.75, pm[3]+(pm[4]-pm[3])*.96, cex=.5, labels = "r2 with", pos=4,
       offset = .2)
  ## simple axis
  text(l[c(1,20,40,60,80,100)], pm[3]+(pm[4]-pm[3])*(.93-(1*.08)), labels=c(0,.2,.4,.6,.8,1), pos=1, cex=.4, offset = .1)
  
  #------------------------------------#
  ##--  regional association - EXP  --##
  #------------------------------------#
  
  
  ## add all variants
  plot(sum.stat[, "BP"], sum.stat[, "LOG10P"], xlab="", ylab=expression(-log[10]("p-value")),
       ylim=c(0,max(sum.stat[, "LOG10P"], na.rm=T)*1.1), cex=.3, col="grey90",
       xaxt="n", yaxt="n")
  axis(2, lwd=.5)
  
  ## rsidentify all up to R2>.3
  tmp <- sum.stat[which(sum.stat[, paste0("R2.",embQTL)] >= .3),]
  ## add points
  points(tmp[, "BP"], tmp[, "LOG10P"], cex=.5, pch=21, col="grey10", 
         bg=col.list[[j]][ceiling(tmp[, paste0("R2.",embQTL)]*100)],
         lwd=.2)
  
  ## add legend
  legend("topleft", lty=0, pch=NA, cex=.5, bty="n", legend = sum.coloc$embedding[1])
  
  #-----------------------#
  ##--    plot genes   --##
  #-----------------------#
  
  cat("\n plot genes ", a.vars[j], "\n")
  
  #------------------------------------#
  ##--        gene assignment       --##
  #------------------------------------#
  
  ## import entire list, given that biomart acts weirdly
  tmp.genes <- rtracklayer::readGFF("<path_to_file>",
                                    filter = list(type=c("gene")))
  ## subset to what is needed
  tmp.genes <- subset(tmp.genes, seqid == ifelse(sum.stat$CHR_ENSEMBL[1] == 23, "X", sum.stat$CHR_ENSEMBL[1]) & start >= min(sum.stat$BP, na.rm=T)-3e6 & end <= max(sum.stat$BP, na.rm=T)+3e6)
  
  ## restrict to protein encoding genes for now
  if(nrow(tmp.genes) > 20){
    tmp.genes <- subset(tmp.genes, gene_biotype %in% c("protein_coding", "processed_transcript"))
  }
  ## sort by start
  tmp.genes <- tmp.genes[order(tmp.genes$start),]
  
  print(tmp.genes)
  
  ## dummy for the line in the plot
  tmp.genes$line    <- NA
  tmp.genes$line[1] <- 1
  
  ## start sorting
  l <- 0
  
  ## loop over everything, collect genes row-wise
  while(sum(is.na(tmp.genes$line)) > 0){
    ## increase line 
    l <- l+1
    e <- 0
    for(j in 1:nrow(tmp.genes)){
      ## test whether new start is larger than the current end
      if(tmp.genes$start[j] > e+4e4 & is.na(tmp.genes$line[j])){
        ## assign line to be drawn
        tmp.genes$line[j] <- l
        ## assign new end
        e                 <- tmp.genes$end[j]
      }
    }
  }
  
  print(head(tmp.genes))
  
  ## now plot it just below
  par(mar=c(1.5,1.5,.1,.5), tck=-.02, bty="o", lwd=.5)
  plot(range(sum.stat$BP), c(0,max(tmp.genes$line)+.5), yaxt="n", xaxt="n", ylab="", 
       xlab=paste("Genomic position on chrosome", sum.stat$CHROM[1]), type="n",
       ylim=rev(c(0,max(tmp.genes$line)+.5)))
  axis(1, lwd=.5)
  ## position of the lead variant
  for(j in 1:length(a.vars)){
    ii <- which(sum.stat$ID == a.vars[j])
    abline(v=sum.stat$BP[ii], lwd=.3, col=cls[j])
  }
  ## add genes
  arrows(tmp.genes$start, tmp.genes$line, tmp.genes$end, tmp.genes$line, lwd=.5, length = 0)
  ## add Gene Names on top
  text(tmp.genes$start+(tmp.genes$end - tmp.genes$start)/2, tmp.genes$line, cex=.3, font=3,
       labels=tmp.genes$gene_name, pos=3, offset=.1)
}