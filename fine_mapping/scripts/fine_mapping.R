#   This file is part of the analysis for the publication 'CRISPR-directed mitotic recombination 
#   enables genetic mapping without crosses'. Code written by Joshua S Bloom.
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

# 05/01/15 LOH  Image processing and mapping pipeline
# note: saves and loads are commented. uncomment as necessary

require(EBImage)
require(gdata)
library(locfit)
options(browser='chromium-browser')
library(foreach)
library(doMC)
registerDoMC(cores=6)
library(Matrix)
library(rrBLUP)

# a couple of global variables -------------------------------------------
plate.types  <- list('96'=c(8,12), '384'=c(16,24), '1536'=c(32,48) )
dimensions = c(5184,3456)
corners = data.frame(X=c(525,4618,525,4624),
                     Y=c(388,400,3060,3064))
#-----------------------------------------------------------------------

fx.file1='/data/CRISPR_LOH/04/phenotyping/plategrowthSPFx.R'
source(fx.file1)

image.path =  '/data/CRISPR_LOH/04/phenotyping/Images/'
layout.path = '/data/CRISPR_LOH/04/phenotyping/keys/'
key.file =    '/data/CRISPR_LOH/04/phenotyping/Key.xls'
Results.file = '/data/CRISPR_LOH/04/phenotyping/Results.RData'
Results.processed.file = '/data/CRISPR_LOH/04/phenotyping/ResultsProcessed.RData'
outDir='/data/CRISPR_LOH/04/phenotyping/out/'
dir.create(outDir)

key=read.xls(key.file,stringsAsFactors=F)
# formatting
# key$Concentration[key$Concentration=='0']=''
# Run image processing
Results=foreach(i=1:nrow(key)) %dopar% {  processImages(i, key, image.path, layout.path, outDir, plate.types[['384']], corners) }
plate.names=apply(key, 1, paste, collapse='::' )
names(Results)=plate.names

save(Results, file=Results.file)


# check what 10% quantile of radius size looks like for each plate to check for plates with alignment issues or no growth
bg=sapply(Results, function(x) {quantile(x$s.radius.mean, na.rm=T,.10)})

# Filtering out edges (this is conservative)-----------------------------------------------------------------------------
Results=filterEdges(Results)
# run plate normalization but ultimately not using the normalized values
Results.processed=normalizePlate(Results)  
save(Results.processed, file=Results.processed.file)

#SEG_Results.df= getResults.df(Results, key, grep('^B', key$StrainLayout))
LOH_Results.df= getResults.df(Results, key, grep('^M', key$StrainLayout))

#!!!! pick feature of interest!!!!!!!
foi='s.radius.mean'
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# look for plate effects and effect of condition ------------------------------------------------------------
# in LOH panel
par(oma=c(20,1,1,1))
x=(split(LOH_Results.df[,foi], paste(LOH_Results.df$plate, LOH_Results.df$condition, sep='::')))
#x=x[order(do.call('rbind', strsplit(names(x), '::'))[,2])]
#names(x)=paste(unique(LOH_Results.df$plate),  unique(LOH_Results.df$condition))
boxplot(x, las=2)
#-------------------------------------------------------------------------------------------------------------

# extract radius and strain names
LOHpheno=split(LOH_Results.df[,foi], LOH_Results.df$condition)
LOHpheno.names=split(LOH_Results.df$strain, LOH_Results.df$condition)

Lpheno=mapply(function(x,n){split(x,n)}, LOHpheno, LOHpheno.names, SIMPLIFY=F)

pheno.file='/data/CRISPR_LOH/04/phenotyping/Lphenos.RData'
#save(Lpheno, file=pheno.file)
load(pheno.file)

# global variables
all.strain.names.L=names(Lpheno[[1]])


# Load genotype data for fine-mapping panel that incorporates Sanger data 
load('/home/jbloom/Dropbox/Public/CRISPR_LOH/LOH.384.BYxRM.fm.w.Sanger.RData')
#str(LOH.384.BYxRM.fm.filt)
LOH.384.BYxRM.fm.filt=LOH.384.BYxRM.fm.w.Sanger
LOH.384.BYxRM=t(LOH.384.BYxRM.fm.filt)

#rownames(LOH.384.BYxRM)=toupper(rownames(LOH.384.BYxRM))
la=do.call('rbind', strsplit(colnames(LOH.384.BYxRM), ':'))[,2]
L.pos=as.numeric(do.call('rbind', strsplit(la, '_'))[,1])
L.chr=do.call('rbind', strsplit(colnames(LOH.384.BYxRM), ':'))[,1]
LOH.384.BYxRM=LOH.384.BYxRM[rownames(LOH.384.BYxRM) %in% all.strain.names.L,]
all.strain.names.L=all.strain.names.L[all.strain.names.L %in% rownames(LOH.384.BYxRM)]
print(identical(rownames(LOH.384.BYxRM), all.strain.names.L))
LOH.384.BYxRM.sdx=apply(LOH.384.BYxRM,2,sd, na.rm=T)

BYBY.cnt=apply(LOH.384.BYxRM, 2, function(x) {sum(x==0)})
BYRM.cnt=apply(LOH.384.BYxRM, 2, function(x) {sum(x==1)})
RMRM.cnt=apply(LOH.384.BYxRM, 2, function(x) {sum(x==2)})

pdf(file='~/Dropbox/Public/CRISPR_LOH/LOH_af.pdf' ,width=8, height=8)
par(mfrow=c(3,1))
plot(BYBY.cnt, ylim=c(0,384))
        abline(v=cumsum(rle(L.chr)$lengths), col='blue')
plot(BYRM.cnt, ylim=c(0,384))
        abline(v=cumsum(rle(L.chr)$lengths), col='blue')
plot(RMRM.cnt, ylim=c(0,384))
        abline(v=cumsum(rle(L.chr)$lengths), col='blue')
dev.off()
#-------------------------------------------------------------------------------------


# for now only look at cases where conditions were tested in both LOH and in segregants
# for mgso4
#1190
search.space=which(L.chr=='chrVII')[950:1300]
#search.space=which(L.chr=='chrVII')[1:2000]

# Manganese 10nM 
pp=6
#stats.galore=list()   
#for(pp in 1:length(Lpheno) ) {
    # for debugging
    # pp=7
pheno=names(Lpheno)[pp]

ylist.L=processPhenos_for_MM(Lpheno[[pheno]], all.strain.names.L)
    
l.geno=LOH.384.BYxRM[ylist.L$strains.with.phenos,]
    # code dominance as 0,1,0
l.geno.d=-1*(abs(l.geno-1))+1

# extract avg for replicates
y.avg.l=avg.nvec(ylist.L$y )
L.genome.LOD=get.LOD.by.COR(ylist.L$n.strains, y.avg.l,l.geno, betas=TRUE, sdx=as.vector(LOH.384.BYxRM.sdx))

plot(L.pos[search.space], L.genome.LOD$LOD[1,search.space], xlab='chr VII position', ylab='LOD')
save.image('~/Dropbox/Public/CRISPR_LOH/010516_fine_mapping_workspace.RData')

