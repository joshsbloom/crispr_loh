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
source('~/Dropbox/Public/CRISPR_LOH/calcMM.R')


# a couple of global variables -------------------------------------------
plate.types  <- list('96'=c(8,12), '384'=c(16,24), '1536'=c(32,48) )
dimensions = c(5184,3456)
corners = data.frame(X=c(525,4618,525,4624),
                     Y=c(388,400,3060,3064))
#-----------------------------------------------------------------------

fx.file1='/data/CRISPR_LOH/04/phenotyping/plategrowthfx.R'
source(fx.file1)

image.path =  '/data/CRISPR_LOH/phenotyping/Images/'
layout.path = '/data/CRISPR_LOH/phenotyping/keys/'
key.file =    '/data/CRISPR_LOH/phenotyping/Key.xls'
Results.file = '/data/CRISPR_LOH/phenotyping/Results.RData'
Results.processed.file = '/data/CRISPR_LOH/phenotyping/ResultsProcessed.RData'
outDir='/data/CRISPR_LOH/phenotyping/out/'
dir.create(outDir)

key=read.xls(key.file,stringsAsFactors=F)
# formatting
# key$Concentration[key$Concentration=='0']=''
# Run image processing
Results=foreach(i=1:nrow(key)) %dopar% {  processImages(i, key, image.path, layout.path, outDir, plate.types[['384']], corners) }
plate.names=apply(key, 1, paste, collapse='::' )
names(Results)=plate.names
#save(Results, file=Results.file)

# check what 10% quantile of radius size looks like for each plate to check for plates with alignment issues or no growth
bg=sapply(Results, function(x) {quantile(x$s.radius.mean, na.rm=T,.10)})
# BYxRM_C.txt::1::MS_0066_15-04-19_14-56-28.JPG::SDS::0.10%::0.10% , #66 agar detatched from plate
# BYxRM_D.txt::1::MS_0085_15-04-19_15-23-13.JPG::SDS::0.05%::0.10%   #85 agar detatched from plate
# MS_1.txt::1::SDS_01_MS_0001_15-04-21_10-49-06.JPG::SDS::0.10%::0.10% #93  agar detatched from plate
# typical filter is  bad.plates=which(bg<5) but here we hand picked
bad.plates=c(66,85,93)
Results=Results[-bad.plates]
key=key[-bad.plates,]

# visualize different image processing outputs
glob.filt=sapply(Results, function(x) {(x[,'s.perimeter'])})
# perimeter filter 
hist(glob.filt, breaks=10000, xlim=c(0,quantile(glob.filt,.9995)), ylim=c(0,200))
# eccentricity filter
glob.filt=sapply(Results, function(x) {(x[,'s.area'])})
hist(log(glob.filt), breaks=100, xlim=c(0,15))
abline(v=12)
glob.filt=sapply(Results, function(pa) { (log(pa[,'s.radius.max']-pa[,'s.radius.min']))  } )

# Post processing 
# Filtering out irregular colonies -----------------------------------------------------------------------
# turned off some filters for MS
for( i in 1:length(Results) ){
    pa =Results[[i]]
    # previous filters
    # (pa[,'m.eccentricity']>.8) |  (log(pa[,'s.radius.max']-pa[,'s.radius.min'])>3.5)
    badstuff = (pa[,'s.perimeter']>550) | (log(pa[,'s.area']>12)) 
    Results[[i]][badstuff,-c(1,2)] = NA
}
#-----------------------------------------------------------------------------------------------------------------------

# Filtering out edges ( this is conservative)-----------------------------------------------------------------------------
Results=filterEdges(Results)
# run plate normalization but ultimately not using the normalized values
Results.processed=normalizePlate(Results)  
save(Results.processed, file=Results.processed.file)

SEG_Results.df= getResults.df(Results, key, grep('^B', key$StrainLayout))
LOH_Results.df= getResults.df(Results, key, grep('^M', key$StrainLayout))

#!!!! pick feature of interest!!!!!!!
foi='s.radius.mean'
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# look for plate effects and effect of condition ------------------------------------------------------------
# in LOH panel
par(oma=c(20,1,1,1))
x=(split(LOH_Results.df[,foi], paste(LOH_Results.df$plate, LOH_Results.df$condition, sep='::')))
x=x[order(do.call('rbind', strsplit(names(x), '::'))[,2])]
names(x)=paste(unique(LOH_Results.df$plate),  unique(LOH_Results.df$condition))
boxplot(x, las=2)

#in segregants
par(oma=c(20,1,1,1))
x=(split(SEG_Results.df[,foi], paste(SEG_Results.df$plate, SEG_Results.df$condition, sep='::')))
x=x[order(do.call('rbind', strsplit(names(x), '::'))[,2])]
boxplot(x, las=2)
names(x)=paste(unique(SEG_Results.df$plate),  unique(SEG_Results.df$condition))
boxplot(x, las=2)
#-------------------------------------------------------------------------------------------------------------

# extract radius and strain names
LOHpheno=split(LOH_Results.df[,foi], LOH_Results.df$condition)
LOHpheno.names=split(LOH_Results.df$strain, LOH_Results.df$condition)

Lpheno=mapply(function(x,n){split(x,n)}, LOHpheno, LOHpheno.names, SIMPLIFY=F)

SEGpheno=split(SEG_Results.df[,foi], SEG_Results.df$condition)
SEGpheno.names=split(SEG_Results.df$strain, SEG_Results.df$condition)

Spheno=mapply(function(x,n){split(x,n)}, SEGpheno, SEGpheno.names, SIMPLIFY=F)

phenos=list(Spheno=Spheno, Lpheno=Lpheno)
pheno.file='/data/CRISPR_LOH/phenotyping/phenos.RData'
load(pheno.file)
attach(phenos)
#save(phenos, file=pheno.file)

# global variables
all.strain.names.S=names(Spheno[[1]])
all.strain.names.L=names(Lpheno[[1]])

# Load segregant genotype data---------------------------------------------------------
# substitute with dropbox version
load('/data/CRISPR_LOH/1000BYxRM_with_names.RData')
BYxRM_orig=BYxRM_orig[order(rownames(BYxRM_orig)),]
#str(BYxRM_orig)
BYxRM_orig=BYxRM_orig+1
rownames(BYxRM_orig)=do.call('rbind', strsplit(rownames(BYxRM_orig), ':'))[,1]
S.pos=as.numeric(do.call('rbind', strsplit(colnames(BYxRM_orig), '_'))[,3])
S.chr=(do.call('rbind', strsplit(colnames(BYxRM_orig), '_'))[,2])

# might as well ignore segregants that never have a phenotype
BYxRM_orig=BYxRM_orig[rownames(BYxRM_orig) %in% all.strain.names.S,]
all.strain.names.S=all.strain.names.S[all.strain.names.S %in% rownames(BYxRM_orig)]
print(identical(rownames(BYxRM_orig), all.strain.names.S))
#woot (precompute genotype SDs)
BYxRM_orig.sdx=apply(BYxRM_orig,2,sd, na.rm=T)
#-------------------------------------------------------------------------------------

#Load LOH genotype data -------------------------------------------------------------
# substitute with dropbox version
load('/data/CRISPR_LOH/BYxRM_LOH_segs.RData')
str(LOH.384.BYxRM)
LOH.384.BYxRM=t(LOH.384.BYxRM)
rownames(LOH.384.BYxRM)=toupper(rownames(LOH.384.BYxRM))
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



# power calculation ----------------------------------------------------------------------------
search.space=which(L.chr=='chrVII')[1:2000]

fasterLOD=function(n.pheno, pheno.s,gdata.s, betas=FALSE, sdx=1, pheno=NULL){
   r=crossprod(pheno.s, gdata.s)/(n.pheno-1)
   LOD=(-n.pheno*log(1-r^2))/(2*log(10))
   if(betas==FALSE) {
       return(LOD)
   } else {
       beta=r*apply(cbind(pheno),2, sd,na.rm=T)/sdx
       return(list(r=r, LOD=LOD, beta=beta))
   }
}

#h2=.1
lg=LOH.384.BYxRM[,search.space]
lg.s=scale(lg)
powerLOD=list()
powerLOD[[1]]=matrix(0,2000,100)
powerLOD[[2]]=matrix(0,2000,100)
powerLOD[[3]]=matrix(0,2000,100)
powerLOD[[4]]=matrix(0,2000,100)
powerLOD[[5]]=matrix(0,2000,100)
powerLOD[[6]]=matrix(0,2000,100)
powerLOD[[7]]=matrix(0,2000,100)
powerLOD[[8]]=matrix(0,2000,100)

for(n in 1:100) {
    print(n)
    eff.size=c(.1, .25, .5, 1, 1.5, 3, 5)
    h2mat=matrix(0, 2000, length(eff.size))
    ylist=list()
    for(j in 1:length(eff.size)){
        #print(j)
        ymat=matrix(0,384,2000)
        for(i in 1:2000) {
        g=eff.size[j]*lg[,i]
        varg=var(g)
        e=rnorm(384,mean=0, sd=1)
        vare=var(e)
        h2=var(g)/(var(g)+var(e))
        h2mat[i,j]=h2
        y=g+e
        ymat[,i]=y
        }
       ylist[[j]]=ymat
    }
    ynorm=lapply(ylist, scale)
    yLOD=lapply(ynorm, function(x) {fasterLOD(384, x, lg.s)})
    mLOD=lapply(yLOD, function(x) {apply(x,1,max)})
    for(j in 1:length(eff.size)){
        powerLOD[[j]][,n]=mLOD[[j]]
    }
}
powerOUT=lapply(powerLOD, function(x) { apply(x,1, function(y) sum(y>3)/100 )})
pdf(file='/data/CRISPR_LOH/power_calc.pdf', width=8.5, height=11)
    par(mfrow=c(2,1))
    par(yaxs='i')
    plot(L.pos[search.space], h2mat[,7], ylim=c(0,1), xlab='QTL position, chrVII', ylab=expression(h^2), col='black', cex=.2)
    points(lowess(L.pos[search.space], h2mat[,7], f=.05), type='l', lwd=2)
    points(L.pos[search.space], h2mat[,6], ylim=c(0,1), col='lightblue', cex=.2)
    points(lowess(L.pos[search.space], h2mat[,6], f=.05),col='lightblue', type='l', lwd=2)

    points(L.pos[search.space], h2mat[,5], ylim=c(0,1), col='pink', cex=.2)
    points(lowess(L.pos[search.space], h2mat[,5], f=.05),col='pink', type='l', lwd=2)

    points(L.pos[search.space], h2mat[,4], ylim=c(0,1), col='lightgreen', cex=.2)
    points(lowess(L.pos[search.space], h2mat[,4], f=.05),col='lightgreen', type='l', lwd=2)

    points(L.pos[search.space], h2mat[,3], ylim=c(0,1), col='yellow', cex=.2)
    points(lowess(L.pos[search.space], h2mat[,3], f=.05),col='yellow', type='l', lwd=2)

    points(L.pos[search.space], h2mat[,2], ylim=c(0,1), col='orange', cex=.2)
    points(lowess(L.pos[search.space], h2mat[,2], f=.05),col='orange', type='l', lwd=2)
    par(yaxs='i')
    plot(L.pos[search.space], powerOUT[[7]], ylim=c(0,1.1), xlab='QTL position, chrVII', ylab='power', col='black', cex=.2)
    points(lowess(L.pos[search.space],  powerOUT[[7]], f=.05), type='l', lwd=2)

    points(L.pos[search.space],powerOUT[[6]], ylim=c(0,1), col='lightblue' , cex=.2  )
    points(lowess(L.pos[search.space],  powerOUT[[6]], f=.05), col='lightblue' ,type='l', lwd=2)

    points(L.pos[search.space],powerOUT[[5]], ylim=c(0,1), col='pink'  , cex=.2  )
    points(lowess(L.pos[search.space],  powerOUT[[5]], f=.05),col='pink'  , type='l', lwd=2)

    points(L.pos[search.space],powerOUT[[4]], ylim=c(0,1), col='lightgreen' , cex=.2 ) 
    points(lowess(L.pos[search.space],  powerOUT[[4]], f=.05),col='lightgreen' , type='l', lwd=2)

    points(L.pos[search.space],powerOUT[[3]], ylim=c(0,1), col='yellow', cex=.2 )
    points(lowess(L.pos[search.space],  powerOUT[[3]], f=.05), col='yellow', type='l', lwd=2)

    points(L.pos[search.space],powerOUT[[2]], ylim=c(0,1), col='orange', cex=.2 )
    points(lowess(L.pos[search.space],  powerOUT[[2]], f=.1),col='orange', type='l', lwd=2)

    legend('bottomleft', title=expression(beta), legend=rev(eff.size)[-1][-7], fill=c('black','lightblue',
                                                              'pink',   
                                                              'lightgreen', 
                                                              'yellow',
                                                              'orange'))
dev.off()
save.image('/data/CRISPR_LOH/power.RData')

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# for now only look at cases where conditions were tested in both LOH and in segregants
search.space=which(L.chr=='chrVII')[1:2000]

stats.galore=list()   
for(pp in 1:length(Lpheno) ) {
    # for debugging
    # pp=7
    
    # mapping on segregants--------------------------------------------------------------------
    pheno=names(Lpheno)[pp]
    print(pheno)
    # phenotype to genotype matching and constructs other useful things 
    ylist.S=processPhenos_for_MM(Spheno[[pheno]], all.strain.names.S)
    # extract avg for replicates
    y.avg.s=avg.nvec(ylist.S$y)
    seg.geno=BYxRM_orig[ylist.S$strains.with.phenos,]
    # construct relatedness matrices 
    # Additive, genome
    sA       = A.mat(seg.geno-1,shrink=F)/2
    # Additive (without chr7)
    sA.loco7 = A.mat(seg.geno[,-which(S.chr=='chrVII')]-1,shrink=F)/2
    # Additive (chr7 only)
    sA.7     = A.mat(seg.geno[,which(S.chr=='chrVII')]-1,shrink=F)/2
    # Addititive x Additive
    sAA=sA*sA

    # calculate repeated measures mixed model y=XB+Z(aloco_7)+Z(a7)+Z(AA)+Z(Repeatibility)+e
    s.mm=calcMM(ylist.S$y, B=list(A.loco7=sA.loco7, A.7=sA.7,AA=sAA),alg='ai', reps=TRUE)
    s.mm=extractVarCompResults(s.mm)
   
    # subtract out the additive effect for the rest of the genome and the epistatic effect 
    G.f=s.mm$sigma[1]*sA.loco7 + s.mm$sigma[3]*sAA
    X=rep(1,length(ylist.S$y))
    blup.residuals=ylist.S$y-as.vector(ylist.S$Z%*%(calc.BLUPS(G.f, ylist.S$Z, s.mm$W, ylist.S$y, X, t(s.mm$Bhat))))
    # big and annoying to store this
    s.mm$W=NULL
    #and average the residuals for each strain for mapping
    bresid.avg.s=avg.nvec(blup.residuals)

    # LOD scores, genome-wide for trait averages
    S.genome.LOD=get.LOD.by.COR(ylist.S$n.strains, y.avg.s, seg.geno,
                                betas=TRUE, sdx=as.vector(BYxRM_orig.sdx))
    # LOD scores, genome-wide for blup residuals
    S.br7.LOD=get.LOD.by.COR(ylist.S$n.strains, bresid.avg.s, seg.geno,
                             betas=TRUE, sdx=as.vector(BYxRM_orig.sdx))

    #LOD scores, genome-wide for blup residuals (chrVII only)
    S.br7.LOD.c7=get.LOD.by.COR(ylist.S$n.strains, bresid.avg.s,seg.geno[,which(S.chr=='chrVII')])
    # emprical null
    print('permutations for blup residuals')
    S.br7.LOD.c7.null=get.LOD.by.COR(ylist.S$n.strains, replicate(1000, sample(bresid.avg.s)), seg.geno[,which(S.chr=='chrVII')])
    # attach FWER<.05 threshold
    attr(S.br7.LOD.c7, 'thresh')=quantile(  apply(S.br7.LOD.c7.null$LOD,1,max), .95)
                           
   
    ylist.L=processPhenos_for_MM(Lpheno[[pheno]], all.strain.names.L)
    
    l.geno=LOH.384.BYxRM[ylist.L$strains.with.phenos,]
    # code dominance as 0,1,0
    l.geno.d=-1*(abs(l.geno-1))+1
    # extract avg for replicates
    y.avg.l=avg.nvec(ylist.L$y )
    
    bck.split=grepl('^B', names(y.avg.l))
    #plot(y.avg.l, bck.split)
    #y.avg.l.r=y.avg.l-predict(lm(y.avg.l~as.numeric(bck.split)))

    # scanone for all of LOH panel   
    L.genome.LOD=get.LOD.by.COR(ylist.L$n.strains, y.avg.l,l.geno, betas=TRUE, sdx=as.vector(LOH.384.BYxRM.sdx))

    # Test on chromosome VII first 2000 markers
    l.geno.s=l.geno[,search.space]
    l.geno.s.d=l.geno.d[,search.space]
    L.LOD.c7=get.LOD.by.COR(ylist.L$n.strains, y.avg.l, l.geno.s)
    # empirical null, permutations structured by strain background
    print('permutations for LOH panel .. stuctured')
    L.LOD.c7.structured_null=get.LOD.by.COR(ylist.L$n.strains, replicate(1000, structured.perm(y.avg.l, bck.split)),  l.geno.s)
    attr(L.LOD.c7, 'structured_thresh')=quantile(  apply(L.LOD.c7.structured_null$LOD,1,max), .95 )
    rm(L.LOD.c7.structured_null)
   
    # empirical null, full permutations
    print('permutations for LOH panel .. unstuctured')
    L.LOD.c7.null=get.LOD.by.COR(ylist.L$n.strains, replicate(1000,sample(y.avg.l)),l.geno.s)
    # attach FWER thresholds
    attr(L.LOD.c7, 'thresh')=quantile(  apply(L.LOD.c7.null$LOD,1,max), .95 )
    rm(L.LOD.c7.null)

    # calculate effects in both strain backgrounds separately                       
    L.genome.LOD.B=get.LOD.by.COR(sum(bck.split), y.avg.l[bck.split], LOH.384.BYxRM[ylist.L$strains.with.phenos,][bck.split,], 
                                betas=TRUE, sdx=as.vector(LOH.384.BYxRM.sdx))
    L.genome.LOD.R=get.LOD.by.COR(sum(!bck.split), y.avg.l[!bck.split], LOH.384.BYxRM[ylist.L$strains.with.phenos,][!bck.split,], 
                                betas=TRUE, sdx=as.vector(LOH.384.BYxRM.sdx))
    # Scan for dominance 
    domList=dominanceScan(ylist.L, l.geno.s, l.geno.s.d, search.space)
    
    stats.galore[[pheno]]=list(ylist.S=ylist.S,
                               ylist.L=ylist.L,
                               s.mm   =s.mm,
                               S.genome.LOD= S.genome.LOD,
                               S.br7.LOD=    S.br7.LOD,
                               S.br7.LOD.c7= S.br7.LOD.c7,
                               L.genome.LOD=        L.genome.LOD,
                               L.LOD.c7= L.LOD.c7,
                               L.genome.LOD.B=       L.genome.LOD.B,    
                               L.genome.LOD.R=        L.genome.LOD.R,
                               domList=domList
                               )    
}

#save(stats.galore, file='/data/CRISPR_LOH/LOH_stats.RData')
load('/home/jbloom/Dropbox/Public/CRISPR_LOH/LOH_stats.RData')

#stats.galore
#for(i in 1:length(stats.galore)){
#    stats.galore[[i]]$s.mm$W=NULL
#    stats.galore[[i]]$domList$add.model=NULL
#    stats.galore[[i]]$domList$dom.model=NULL
#}

makePlots.1('~/Dropbox/Public/CRISPR_LOH/LOH_plots_v3.pdf',stats.galore)

vc.results= lapply(stats.galore, function(x) {
    s.mm=x$s.mm
    # scale variance components
    nf=sum(s.mm$sigma)
    vcs=t(s.mm$sigma/nf)
    vcs.se =t(sqrt(diag(s.mm$sigma.cov))/nf)
    colnames(vcs)=c('A.loco7',  'A.7', 'AA', 'Strain', 'E')
    return(rbind(vcs,vcs.se))})

vcs=sapply(vc.results, function(x) x[1,])
vc.se=sapply(vc.results, function(x) x[2,])

pdf(file='~/Dropbox/Public/CRISPR_LOH/seg_var_comp.pdf', width=11, height=11)
par(oma=c(8,1,1,1))
vc.cum=apply(vcs,2, cumsum)
bp=barplot(vcs[1:4,], las=2, col=c('lightblue', 'blue', 'lightgreen', 'pink' ), ylim=c(0,1),
           ylab='fraction of phenotypic variance',
            legend=c(rownames(vcs)[1:4]),
args.legend=list(x='topleft') )
segments(bp-.2,  vc.cum[1,]-vc.se[1,], bp-.2, vc.cum[1,]+vc.se[1,], lwd=1.5, col='black')
segments(bp-.1,  vc.cum[2,]-vc.se[2,], bp-.1, vc.cum[2,]+vc.se[2,], lwd=1.5, col='black')
segments(bp+.2,  vc.cum[3,]-vc.se[3,], bp+.2, vc.cum[3,]+vc.se[3,], lwd=1.5, col='black')
segments(bp+.3,  vc.cum[4,]-vc.se[4,], bp+.3, vc.cum[4,]+vc.se[4,], lwd=1.5, col='black')
dev.off()


#test for sig QTL 
traits.with.sig.LOH.QTL=sapply(names(stats.galore), function(pp) {
   #print(pp)
   # print(    attr(stats.galore[[pp]]$L.LOD.c7, 'structured_thresh'))
   # print(max(stats.galore[[pp]]$L.LOD.c7$LOD[1,]))    
    max(stats.galore[[pp]]$L.LOD.c7$LOD[1,]) >  attr(stats.galore[[pp]]$L.LOD.c7, 'structured_thresh')  
    }  ) 
sig.traits=which(traits.with.sig.LOH.QTL)[c(1,2,3,4,5,8:11)]

# if we collapse to unique compounds 
#9 traits have a LOH linkage out of 12

pdf(file='~/Dropbox/Public/CRISPR_LOH/full_panel_LOD_plots_segs_blue_new_thresh_new_yaxs.pdf', width=18,height=21)
    par(mfrow=c(4,3))
    par(xaxs='i', yaxs='i')
    for(i in 1:9) {
        pp=sig.traits[i]

        #str(stats.galore[[pp]])
        sL=stats.galore[[pp]]$S.br7.LOD.c7$LOD
        sL.pos= as.numeric(sapply(strsplit(colnames(sL), '_'), function(x) x[3]))
        lL=stats.galore[[pp]]$L.LOD.c7$LOD[1,] 

        #blue is segregants
        plot(sL.pos, sL[1,], col='blue', main=gsub('.95%', '', names(sig.traits))[i], 
             ylim=c(0, max(sL[1,], lL,10) ), xlim=c(0,470000), type='l', ylab='LOD', xlab='pos' ,lwd=1.5)
        lL.pos = as.numeric(sapply(strsplit(sapply(strsplit(names(lL), '_'), function(x) x[1]), ':'), function(x) x[2]))
        points(lL.pos, lL, type='l',lwd=1.5)
        abline(h=attr(stats.galore[[pp]]$S.br7.LOD.c7, 'thresh'), lty=3, col='blue')
        abline(h=attr(stats.galore[[pp]]$L.LOD.c7, 'thresh'), lty=2, col='black')
        #readline()
    }

    nsig.traits=which(!traits.with.sig.LOH.QTL)[c(1,2,5)]
    for(i in 1:3) {
        pp=as.vector(nsig.traits[i])

        #str(stats.galore[[pp]])
        sL=stats.galore[[pp]]$S.br7.LOD.c7$LOD
        sL.pos= as.numeric(sapply(strsplit(colnames(sL), '_'), function(x) x[3]))
        lL=stats.galore[[pp]]$L.LOD.c7$LOD[1,] 

        #blue is segregants
        plot(sL.pos, sL[1,], col='blue', main=gsub('.95%', '', names(nsig.traits))[i], 
             ylim=c(0, max(sL[1,], lL,10) ), xlim=c(0,470000), type='l', xlab='pos', ylab='LOD',lwd=1.5)
        lL.pos = as.numeric(sapply(strsplit(sapply(strsplit(names(lL), '_'), function(x) x[1]), ':'), function(x) x[2]))
        points(lL.pos, lL, type='l',lwd=1.5)

        abline(h=attr(stats.galore[[pp]]$S.br7.LOD.c7, 'thresh'), lty=3, col='blue')
        abline(h=attr(stats.galore[[pp]]$L.LOD.c7, 'thresh'), lty=2, col='black')

        #readline()
    }
dev.off()


pp=names(Lpheno)[7]
sL=stats.galore[[pp]]$S.br7.LOD.c7$LOD
sL.pos= as.numeric(sapply(strsplit(colnames(sL), '_'), function(x) x[3]))
lL=stats.galore[[pp]]$L.LOD.c7$LOD[1,] 
lL.pos = as.numeric(sapply(strsplit(sapply(strsplit(names(lL), '_'), function(x) x[1]), ':'), function(x) x[2]))
mgso4.LODs=list()
mgso4.LODs[['segs']]=data.frame(pos=sL.pos,
                                LOD=sL[1,]
                                )
mgso4.LODs[['LOH_full']]=data.frame(pos=lL.pos,
                                LOD=lL)

save(mgso4.LODs, file='/home/jbloom/Dropbox/Public/CRISPR_LOH/LODs_mgso4.RData')
L.genome.LOD

fm.LODs.mgso4=L.genome.LOD$LOD[1,]
fm.LODs.mgso4=fm.LODs.mgso4[grep('chrVII:', names(fm.LODs.mgso4))]
fm.LODs.mgso4.pos=as.vector

fm.M.pos=as.numeric(sapply(strsplit(sapply((strsplit(names(fm.LODs.mgso4), ':')), function(x)x[2]), '_'), function(x) x[1]))
mgso4.LODs[['LOH_fm']]=data.frame(pos=fm.M.pos,
                                  LOD=as.vector(fm.LODs.mgso4))
save(mgso4.LODs, file='/home/jbloom/Dropbox/Public/CRISPR_LOH/LODs_mgso4.RData')

save.image('/data/CRISPR_LOH/for_publication_020816.RData')
