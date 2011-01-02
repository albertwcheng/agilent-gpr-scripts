# Run script from directory with agilent files, in other words, working directory should be set to where your agilent files live
library(limma)

setwd("..")

setwd("raw")

# Read GRP and GAL files
gprFiles = dir(pattern = "*.gpr$")
RG <- read.maimages(gprFiles, source="genepix")
RG$genes <- readGAL()
RG$printer <- getLayout(RG$genes)

# Preprocess data
RG.b <-backgroundCorrect(RG,method="none",offset=0) #no background subtract. RG.b <-backgroundCorrect(RG,method="minimum")
MA.p <-normalizeWithinArrays(RG.b, method="loess")
MA.pq <- normalizeBetweenArrays(MA.p, method="quantile")
RG.pq <- RG.MA(MA.pq)


setwd("..")
dir.create("normalized")
setwd("normalized")

# Write data
plotDensities(RG.pq)
data <- cbind(RG$genes, log2(RG.pq$R), log2(RG.pq$G))
write.table(data, file = "G.bgc0.norm.gn.txt", sep="\t")
data.nolog <- cbind(RG$genes, RG.pq$R, RG.pq$G)
write.table(data.nolog, file = "G.bgc0.norm.gn.nolog.txt", sep="\t")
#plotDensities(RG.pq)

setwd("..")
