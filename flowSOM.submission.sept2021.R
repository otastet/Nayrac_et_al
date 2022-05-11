library(flowCore)
library(FlowSOM)
library(M3C)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(ggfortify)
library(grid)
library(ggplotify)
library(ggpubr)

################################
#     Convert CSV to FCS       #
################################    


## Get CSV files
outDir <- 'new_FCS'
PrimaryDirectory <- getwd()
FileNames <- list.files(path=PrimaryDirectory, pattern = ".csv")
DataList=list()

for (File in FileNames) {       
    tempdata <- fread(File, check.names = FALSE)
    File <- gsub(".csv", "", File)
    DataList[[File]] <- tempdata
}

rm(tempdata)
AllSampleNames <- names(DataList)

## Set Time feature for FCS
x <- Sys.time()
x <- gsub(":", "-", x)
x <- gsub(" ", "_", x)

newdir <- paste(PrimaryDirectory,outDir,sep='/')

setwd(PrimaryDirectory)
dir.create(paste0(PrimaryDirectory,newdir), showWarnings = FALSE)
setwd(paste0(PrimaryDirectory,newdir))

## Iterate on samples and save as FCS 
for(i in c(1:length(AllSampleNames))){
		data_subset <- DataList[i]
		data_subset <- rbindlist(as.list(data_subset))
		colnames(data_subset)=cols
		dim(data_subset)
		a <- names(DataList)[i]

		metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))
		
		## Create FCS file metadata - ranges, min, and max settings
		# metadata$range <- apply(apply(data_subset,2,range),2,diff)
		metadata$minRange <- apply(data_subset,2,min)
		metadata$maxRange <- apply(data_subset,2,max)
		
		data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
		print(colnames(data_subset.ff))
		write.FCS(data_subset.ff, paste0(a, ".fcs"))
}

#### Analysis on newly created FCS files ####

## Set working directories 
aim.cd4.dir = ...
ics.cd4.dir = ...

## Extract FCS files 
aim.files = paste(aim.cd4.dir,list.files(aim.cd4.dir),sep='/')
ics.files = paste(ics.cd4.dir,list.files(ics.cd4.dir),sep='/')

## Read all AIM files as a FlowSET
fs.aim = read.flowSet(fcs_files,ignore.text.offset=T,truncate_max_range = FALSE)

## Channels of interest for AIM analysis 
ids=c(9,11,12,13,16,18,19)

## Print total number of gated events 
k=fsApply(fs.aim, function(ff) {
	print(nrow(ff))
})

## Establish ideal downsample size and subset each FCS files by that number 
sampling.ceiling.cd4 <- 1000
fs.ds <- fsApply(fs.aim, function(ff) {
	idx <- sample.int(nrow(ff), min(sampling.ceiling.cd4, nrow(ff)))
	ff[idx,]  
})

## Print channels corresponding to ids, to validate markers
print(paste(c('Channels = ', paste(colnames(fs[[1]])[ids],collapse=', ')),collapse=''))

## Use FlowSOM to integrated the FCS files, apply comppensation, scaling and transformation 
fs = FlowSOM(fs.ds, compensate=FALSE, transform=TRUE, toTransform=ids, scale=TRUE,colsToUse=ids,nClus=10)

## Parse file names 
names(fs$metaData)[!grepl('V0|V1|V3', names(fs$metaData))]=gsub('_S',' V2_S',names(fs$metaData)[!grepl('V0|V1|V3', names(fs$metaData))])
names(fs$metaData)=gsub('export V2','export',names(fs$metaData))
names(fs$metaData)=gsub('CF021','CF21',names(fs$metaData))
names(fs$metaData)=gsub('CF023','CF23',names(fs$metaData))
names(fs$metaData)=gsub('V2 V2','V2',names(fs$metaData))
names(fs$metaData)=gsub('CF02 ','CF2 ',names(fs$metaData))

## Run analysis function 

## Extract the transformed MFI values from the flow set 
df = data.frame(fs$data)

## Perform UMAP only on the channels of interest
um = umap(t(df[,ids]))

## Parse sample-wise metadata info
out = c()
for(j in names(fs$metaData)){
	idx = fs$metaData[j]
	sID = gsub(' ','_',unlist(strsplit(j,'_'))[3])
	r = (idx[[1]][2]-idx[[1]][1])+1
	out=c(out,rep(sID,r))
}
	
## Create output data frame with UMAP coordinates 
dat = um$data[,1:2]
colnames(dat)=c('UMAP_1','UMAP_2')

pheno_cluster = Rphenograph(df[,ids],k=150)
dat$pheno_cluster=paste0('C',pheno.cluster[[2]]$membership)

## Add sample annotation based on the out variable 
dat$group = out
dat$condition = substr(dat$group,1,1)
dat$time = unlist(strsplit(dat$group,'_'))[seq(2,nrow(dat)*2,2)]
dat$group2 = paste(dat$condition,dat$time,sep='.')
	
## Add transformed MFI
a = data.frame(fs$data)[,ids]
dat.num = data.frame(cbind(dat,a)) 
	
## Plot UMAPs (p2) splitted by condition
p2=ggplot(dat.num,aes(UMAP_1,UMAP_2,col=cluster.k9))+geom_point()+facet_grid(condition~time)+theme_bw()
	
## Subset only channel MFIs and aggregate with mean across clusters 
sub = dat.num[,colnames(df)[ids]]
sub$cluster = dat.num[,]
sub = melt(sub)
sub = aggregate(sub$value,by=list(sub$cluster,sub$variable),mean)
d = dcast(sub,formula=Group.1~Group.2)
rownames(d)=d[,1]

## Plot mean MFI across clusters and scale to see relative intensity  
ph = as.ggplot(pheatmap(d[,2:ncol(d)],scale='column',border_col='white',cluster_col=F,cluster_row=F))


## Make boxplots accross timepoints 
t=table(dat.num[,c('group','pheno_cluster')])
for(i in 1:nrow(t)){t[i,]=t[i,]/sum(t[i,])}
t=data.frame(t)
t$condition  = substr(t[,1],1,1)
t$timepoint= unlist(strsplit(as.character(t[,1]),'_'))[seq(2,nrow(t)*2,2)]
t$con = paste(t$condition,t$timepoint,sep='_')
t$sample = unlist(strsplit(as.character(t[,1]),'_'))[seq(1,nrow(t)*2,2)]
pl=ggplot(t,aes(con,Freq*100,fill=con))+geom_boxplot(outlier.shape=NA)+facet_wrap(~t[,'cluster.k9'],scales='free_y',ncol=5)+theme_minimal()+	stat_compare_means(comparisons=list(c('C_V0','C_V1'),c('S_V0','S_V1'),c('C_V1','C_V2'),c('S_V1','S_V2'),c('C_V2','C_V3'),c('S_V2','S_V3'))	,size=1.5,color='red')+geom_point(size=.1)+ylab('Percentage (%)')+xlab('')+theme(axis.text.x=element_text(angle=60,hjust=1))+geom_line(aes	(group=sample),size=.1)
pl2=ggplot(t,aes(timepoint,Freq*100,fill=condition))+geom_boxplot(outlier.shape=NA)+facet_wrap(~t[,'cluster.k9'],scales='free_y',ncol=5)+theme_minimal()+stat_compare_means(comparisons=list(c('C_V0','C_V1'),c('S_V0','S_V1'),c('C_V1','C_V2'),c('S_V1','S_V2'),c('C_V2','C_V3'),c('S_V2','S_V3')),size=1.5,color='red')+geom_point(position=position_dodge(width=.85),size=.1)+ylab('Percentage (%)')+xlab('')+theme(axis.text.x=element_text(angle=60,hjust=1))+geom_line(aes(col=condition,group=sample),size=.1)

## Save sample-wise cluster proportion table 
props = dcast(t,group~cluster.k9,value.var='Freq')

output.aim=list(
	plots = list(umap=p2, heatmap=ph, boxplots=list(b1=pl, b2=pl2)),
	data = list(facs=dat.num, props=props)
)	



## Read FCS files corresponding to ICS experiment
fs.ics = read.flowSet(files.ics,ignore.text.offset=T)

## Channels of interest for ICS analysis
ids=c(8, 18, 19, 20, 21, 23)

## Print number of events per FCS file
fsApply(fs.ics, function(ff) {
	print(nrow(ff))
})

## Downsample to appropriate value 
sampling.ceiling.cd4 <- 100
fs.ds <- fsApply(fs.ics, function(ff) {
	idx <- sample.int(nrow(ff), min(sampling.ceiling.cd4, nrow(ff)))
	ff[idx,]  
})

ids=1:13
## Print channels corresponding to ids, to validate markers
print(paste(c('Channels = ', paste(colnames(fs.ics[[1]])[ids],collapse=', ')),collapse=''))

## Use FlowSOM to integrated the FCS files, apply comppensation, scaling and transformation 
fs.ics = FlowSOM(fs.ds, compensate=TRUE,transform=TRUE,  scale=TRUE,colsToUse=ids,nClus=10)

## Parse encoded file names 
names(fs$metaData)=gsub('_S.fcs|_Spike.fcs','.fcs',names(fs$metaData))
names(fs$metaData)=gsub('_Spike_','_S_',names(fs$metaData))
names(fs$metaData)[!grepl('V0|V1|V3', names(fs$metaData))]=gsub('_S_',' V2_S_',names(fs$metaData)[!grepl('V0|V1|V3', names(fs$metaData))])
names(fs$metaData)=gsub('V2 V2','V2',names(fs$metaData))
names(fs$metaData)[grepl(' ',names(fs$metaData))]
names(fs$metaData)[!grepl(' ',names(fs$metaData))]=gsub('.fcs',' V2.fcs',names(fs$metaData)[!grepl(' ',names(fs$metaData))])
bool=grepl('^export',names(fs$metaData))
names(fs$metaData)[bool]=paste0(unlist(strsplit(names(fs$metaData)[bool],'_'))[seq(2,sum(bool)*4,4)],'.fcs')

## Run analysis function 
## Extract the transformed MFI values from the flow set 
df = data.frame(fs$data)

## Perform UMAP only on the channels of interest
um = umap(t(df[,ids]))

## Parse sample-wise metadata info
out = c()
for(j in names(fs$metaData)){
	idx = fs$metaData[j]
	sID = gsub(' ','_',unlist(strsplit(j,'_'))[3])
	r = (idx[[1]][2]-idx[[1]][1])+1
	out=c(out,rep(sID,r))
}
	
## Create output data frame with UMAP coordinates 
dat = um$data[,1:2]
colnames(dat)=c('UMAP_1','UMAP_2')

pheno_cluster = Rphenograph(df[,ids],k=150)
dat$pheno_cluster=paste0('C',pheno.cluster[[2]]$membership)

## Add sample annotation based on the out variable 
dat$group = out
dat$condition = substr(dat$group,1,1)
dat$time = unlist(strsplit(dat$group,'_'))[seq(2,nrow(dat)*2,2)]
dat$group2 = paste(dat$condition,dat$time,sep='.')
	
## Add transformed MFI
a = data.frame(fs$data)[,ids]
dat.num = data.frame(cbind(dat,a)) 
	
## Plot UMAPs (p2) splitted by condition
p2=ggplot(dat.num,aes(UMAP_1,UMAP_2,col=cluster.k9))+geom_point()+facet_grid(condition~time)+theme_bw()
	
## Subset only channel MFIs and aggregate with mean across clusters 
sub = dat.num[,colnames(df)[ids]]
sub$cluster = dat.num[,]
sub = melt(sub)
sub = aggregate(sub$value,by=list(sub$cluster,sub$variable),mean)
d = dcast(sub,formula=Group.1~Group.2)
rownames(d)=d[,1]

## Plot mean MFI across clusters and scale to see relative intensity  
ph = as.ggplot(pheatmap(d[,2:ncol(d)],scale='column',border_col='white',cluster_col=F,cluster_row=F))

## Make boxplots accross timepoints 
t=table(dat.num[,c('group','pheno_cluster')])
for(i in 1:nrow(t)){t[i,]=t[i,]/sum(t[i,])}
t=data.frame(t)
t$condition  = substr(t[,1],1,1)
t$timepoint= unlist(strsplit(as.character(t[,1]),'_'))[seq(2,nrow(t)*2,2)]
t$con = paste(t$condition,t$timepoint,sep='_')
t$sample = unlist(strsplit(as.character(t[,1]),'_'))[seq(1,nrow(t)*2,2)]
pl=ggplot(t,aes(con,Freq*100,fill=con))+geom_boxplot(outlier.shape=NA)+facet_wrap(~t[,'cluster.k9'],scales='free_y',ncol=5)+theme_minimal()+	stat_compare_means(comparisons=list(c('C_V0','C_V1'),c('S_V0','S_V1'),c('C_V1','C_V2'),c('S_V1','S_V2'),c('C_V2','C_V3'),c('S_V2','S_V3'))	,size=1.5,color='red')+geom_point(size=.1)+ylab('Percentage (%)')+xlab('')+theme(axis.text.x=element_text(angle=60,hjust=1))+geom_line(aes	(group=sample),size=.1)
pl2=ggplot(t,aes(timepoint,Freq*100,fill=condition))+geom_boxplot(outlier.shape=NA)+facet_wrap(~t[,'cluster.k9'],scales='free_y',ncol=5)+theme_minimal()+stat_compare_means(comparisons=list(c('C_V0','C_V1'),c('S_V0','S_V1'),c('C_V1','C_V2'),c('S_V1','S_V2'),c('C_V2','C_V3'),c('S_V2','S_V3')),size=1.5,color='red')+geom_point(position=position_dodge(width=.85),size=.1)+ylab('Percentage (%)')+xlab('')+theme(axis.text.x=element_text(angle=60,hjust=1))+geom_line(aes(col=condition,group=sample),size=.1)

## Save sample-wise cluster proportion table 
props = dcast(t,group~cluster.k9,value.var='Freq')

output.ics=list(
	plots = list(umap=p2, heatmap=ph, boxplots=list(b1=pl, b2=pl2)),
	data = list(facs=dat.num, props=props)
)	

## Save results of analysis as a R data file 
saveRDS(output.ics,'ICS.CD4.analysis.rds')
saveRDS(output.aim,'AIM.CD4.analysis.rds')


