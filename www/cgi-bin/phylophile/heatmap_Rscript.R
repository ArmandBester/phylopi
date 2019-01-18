
        library(seqinr)
        library(ape)
        library(d3heatmap)
        library(reshape2)
        library(RColorBrewer)
        library(htmlwidgets)
        library(jsonlite)
        library(yaml)
        library(plotly)
        
	aln <- read.alignment('output/replicate_3_repeat_22.02.2018_Iqra/replicate_3_repeat_22.02.2018_Iqra-mafft_trimal.cleaned.aln', format = 'fasta')

        alnBin <- as.DNAbin(aln)
        alnDist <- dist.dna(alnBin)
        alnDist <- as.matrix(dist.dna(alnBin))
        alnDistLong <- melt(alnDist)
        alnDistLong$sameSample <- alnDistLong$Var1 == alnDistLong$Var2
        alnDistLong <- subset(alnDistLong, sameSample %in% FALSE)
        alnDistLong <- alnDistLong[,-4]
        colnames(alnDistLong) <- c('sample1', 'sample2', 'distance')
        alnDistLong$combined <- paste(alnDistLong$sample1, alnDistLong$sample2, sep = ",")
        alnDistLong$sorted <- apply(alnDistLong,
                                    1,
                                    function(x) {paste(sort(unlist(strsplit(x[4], split=","))), collapse = ",")})
        alnDistDedup <- alnDistLong[!duplicated(alnDistLong$sorted),]
        alnDistDedup <- alnDistDedup[, c(-4,-5)]
        alnDistDedup <- alnDistDedup[order(alnDistDedup$distance),]

        alnDist <- data.frame(alnDist,check.names = FALSE) #new
    
	qSample <- names(read.fasta("input/replicate_3_repeat_22.02.2018_Iqra.fas"))

        rowNames <- rownames(alnDist)
        colNames <- colnames(alnDist)

        t <- list(
        family = "sans serif",
        size = 10,
        color = toRGB("black"))
        anno.df <- data.frame(q = character(),
                      rowNames = character())

        for (q in qSample){
            tmp <- cbind(q, rowNames)
            anno.df <- rbind(anno.df, tmp)
        }

        focus.df <- alnDistLong
        focus.df$prim <- focus.df$sample1 %in% qSample

        focus.df <- focus.df %>% 
          filter(prim == TRUE | sample1 == sample2)

        p <- plot_ly(x = focus.df$sample2,
                    y = focus.df$sample1,
                    z = focus.df$distance,
                    type = "heatmap", colors = brewer.pal(11, "RdYlGn"), 
                    zmin = 0.0, zmax = 0.03,  xgap = 2, ygap = 1) %>% 
            layout(margin = list(l = 100, r = 10, b = 100, t = 10, pad = 4), 
                             yaxis = list(tickfont = list(size = 10), showspikes = TRUE),
                             xaxis = list(tickfont = list(size = 10), showspikes = TRUE))

write.csv(alnDist,"/var/www/cgi-bin/phylophile/output/replicate_3_repeat_22.02.2018_Iqra/replicate_3_repeat_22.02.2018_Iqra_heatmap.csv")
write.csv(alnDistDedup,"/var/www/cgi-bin/phylophile/output/replicate_3_repeat_22.02.2018_Iqra/replicate_3_repeat_22.02.2018_Iqra_dists.csv", row.names = FALSE)
htmlwidgets::saveWidget(as_widget(p), '/var/www/cgi-bin/phylophile/output/replicate_3_repeat_22.02.2018_Iqra/replicate_3_repeat_22.02.2018_Iqra_heatmap.html')