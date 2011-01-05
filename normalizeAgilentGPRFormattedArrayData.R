library(limma)

setwd("..")

setwd("raw")

# Read GRP and GAL files
gprFiles = dir(pattern = "*.gpr$")
RG <- read.maimages(gprFiles, source="genepix")
RG$genes <- readGAL()
RG$printer <- getLayout(RG$genes)

# Preprocess data
RG.b <-backgroundCorrect(RG,method="minimum") #use TOBI method Bell Method: RG.b <-backgroundCorrect(RG,method="none",offset=0) #no background subtract.
MA.p <-normalizeWithinArrays(RG.b, method="loess")
MA.pq <- normalizeBetweenArrays(MA.p, method="quantile")
RG.pq <- RG.MA(MA.pq)


setwd("..")
dir.create("normalized")
setwd("normalized")

# Write data
plotDensities(RG.pq)

numChnCols <- dim(RG.pq$R)[2]
numInfoCols <- dim(RG$genes)[2]

data <- cbind(RG$genes, log2(RG.pq$R), log2(RG.pq$G))

colnam=colnames(data)

for(i in (numInfoCols+1):(numInfoCols+numChnCols))
{
	#cat("settingR",i)
	colnam[i]=paste(colnames(data)[i],".","R",sep="")
}
for(i in (numInfoCols+numChnCols+1):(numInfoCols+numChnCols*2))
{
	#cat("settingG",i)
	colnam[i]=paste(colnames(data)[i],".","G",sep="")	
}
colnames(data)=colnam
write.table(data, file = "G.bgc0.norm.gn.txt", sep="\t",quote=FALSE,row.names=FALSE)
data.nolog <- cbind(RG$genes, RG.pq$R, RG.pq$G)
colnames(data.nolog)=colnam
write.table(data.nolog, file = "G.bgc0.norm.gn.nolog.txt", sep="\t",quote=FALSE,row.names=FALSE)
#plotDensities(RG.pq)

setwd("..")
