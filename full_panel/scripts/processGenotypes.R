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

source('/data/rr/R_hmm_code/hmm_funcs.R')
library(VariantAnnotation)

# read table of cut sites
cut.sites = read.delim('/data/CRISPR_LOH/03/cut_sites.csv', header=F, sep=',', stringsAsFactors=F)
cut.sites[,2]=toupper(cut.sites[,2])

vcf=readVcf('/data/CRISPR_LOH/03/MS_chrVII_1234.recode.vcf', "sacCer3")
# extract CHROM, POS, ID, and REF, alt and QUAL
vcf.rd=rowData(vcf)
vcf.info=info(vcf)
# hist(vcf.rd$QUAL)
# filter on QUALITY
# rle(vcf.rd$FILTER)
# QUAL was set to 222
# histogram shows bimodal distribution so made this filter more aggressive
# hist(vcf.rd$QUAL)
vcf.filtered=vcf[(vcf.rd$QUAL>2500) & (vcf.rd$QUAL<18000) & (vcf.info$comp==TRUE)]

# for memory concerns
rm(vcf)
gc()

# seqlevels(vcf.filtered)=sortSeqlevels(seqlevels(vcf.filtered))
# re-extract rowData
vcf.filtered.rd=rowData(vcf.filtered)
# extract all geno info (equivalent to getSNPInfo)
gall=geno(vcf.filtered)
# get matrix of genotype calls (using 0,1,2 coding (ref/ref, ref/alt, alt/alt)
res=genotypeToSnpMatrix(vcf.filtered)
gdata=t(as(res$genotype, "numeric"))

# there is probably a cleaner way to get this marker annotation information
#extract marker names
mname= rownames(gdata)
mname.split=strsplit(rownames(gdata), ':')
chr=sapply(mname.split, function(x)x[1])
pos.by.chr=split(start(vcf.filtered.rd), seqnames(vcf.filtered))

# extract depth (annoying/slow that this is stored by VariantAnnotation as an array of lists)
REF.d=apply(gall$AD, c(1,2), function(x) {  x[[1]][1] })
ALT.d=apply(gall$AD, c(1,2), function(x) {  x[[1]][2] })

# diagnostics
#alt.sum=apply(ALT.d, 1, sum, na.rm=T)
#ref.sum=apply(REF.d, 1, sum, na.rm=T)
#na.count.A=apply(ALT.d, 1, function(x) {sum(is.na(x)) })
#na.count.R=apply(REF.d, 1, function(x) {sum(is.na(x)) })

# set NA counts to 0
REF.d[is.na(REF.d)]=0
ALT.d[is.na(ALT.d)]=0


# extract genotype likelihoods
REF.l=apply(gall$PL, c(1,2), function(x) { 10^( x[[1]][1]/-10 )})
HET.l=apply(gall$PL, c(1,2), function(x) { 10^( x[[1]][2]/-10) })
ALT.l=apply(gall$PL, c(1,2), function(x) { 10^( x[[1]][3]/-10) })

# more memory cleanup
rm(gall)
gc()

# might want to play with this ... flat prior for NA, was 0.5
# if ref and het and alt are all NA
REF.l[is.na(REF.l)]=.333
HET.l[is.na(HET.l)]=.333
ALT.l[is.na(ALT.l)]=.333

#typically this is something like this
#tmatrix=matrix(.01,3,3)
#diag(tmatrix)=.98
tmatrix=matrix((1/5e20)/2,3,3)
diag(tmatrix)=1-((1/5e20))

# Go through by chromosome, and then by segregant
# fit HMM and run baum welch for each segregant
doHMM = function(REF.l, HET.l, ALT.l, chr, tmatrix) {
    unique.chrs=unique(chr)
    REFliks.c  = split(data.frame(REF.l), chr)
    HETliks.c  = split(data.frame(HET.l), chr)
    ALTliks.c  = split(data.frame(ALT.l), chr)
    REF.post=list()
    ALT.post=list()
    HET.post=list()
    CALL.hmm=list()

    #split by chromosome 
    # usually use lapply for cases like this but for returning multiple lists this is more straightforward
    # iterate through each chromosome
    for (cchrom in unique.chrs ) {
            print(cchrom)
            REF.c = REFliks.c[[cchrom]]
            HET.c = HETliks.c[[cchrom]]
            ALT.c = ALTliks.c[[cchrom]]
            REFp  = matrix(NA, ncol=ncol(REF.c), nrow=nrow(REF.c))
            ALTp  = REFp
            HETp  = REFp
            CALL.HMM = REFp
            colnames(REFp) = colnames(REF.c[[1]])
            colnames(ALTp) = colnames(REFp)
            colnames(HETp) = colnames(REFp)
            colnames(CALL.HMM) = colnames(REFp)
            # for each individual
            for (ccol in 1:ncol(REF.c) ) {
               hmmOUT=  BWtrans(rbind(REF.c[,ccol], HET.c[,ccol], ALT.c[,ccol]), tmatrix, 1)
               REFp[,ccol]=hmmOUT$pos[1,]
               HETp[,ccol]=hmmOUT$pos[2,]
               ALTp[,ccol]=hmmOUT$pos[3,]
               # ref/het/alt as 0,1,2
               # pick the most likely state
               # could add filter here, e.g. only ouput if posterior prob > threshold
               # or posterior prob > next most likely state by certain amount
               CALL.HMM[,ccol]=apply(hmmOUT$post, 2, function(x) { which.max(x) })-1
            }
            REF.post[[cchrom]] = REFp
            HET.post[[cchrom]] = HETp
            ALT.post[[cchrom]] = ALTp
            CALL.hmm[[cchrom]] = CALL.HMM
    }
        return( CALL.hmm)
}
CALL.hmm=doHMM(REF.l, HET.l, ALT.l, chr, tmatrix)

# map indices up to start of chrVII 
c7.anchor=sum( seqnames(vcf.filtered)@lengths[1:6])

# store the distance between the closest breakpoint and the cut
cut.breakpoint.dists=rep(NA,384)
names(cut.breakpoint.dists)=colnames(gdata)
#store the closest breakpoint to a cut
cut.breakpoint=rep(NA,384)
names(cut.breakpoint)=colnames(gdata)

# blech, break this up into separate bits 
# plot chrVII only and find breakpoint dist to targeted cut
pdf(file='/home/jbloom/Dropbox/Public/CRISPR_LOH/CHR_VII_LOH_chrVII_only.pdf', width=15, height=8)
for(i in 1:384) {
    # which strain?
    strain.name=colnames(gdata)[i]
    strain.name.split=strsplit(strain.name, '-')[[1]]
    strain.name.base=paste(toupper(strain.name.split[1:2]), collapse='-')
    # find target cut site from table
    cut.pos= cut.sites[which(cut.sites[,2]==strain.name.base),1]
    # find index of target cut site
    c7.ind=findInterval(cut.pos, pos.by.chr[['chrVII']])
    # extract hmm calls
    hmm.calls.chr=(CALL.hmm[['chrVII']][,i])-1
    
    # breakpoint identification
    r=rle(hmm.calls.chr)
    breakpoints=cumsum(r$lengths)
    breakpoints=breakpoints[r$values==1 | r$values==-1]
    closest.breakpoint=breakpoints[which.min(abs(breakpoints-c7.ind))]
    # positive = right of cut, negative = left of cut
    breakpoint.dist.to.cut=pos.by.chr[['chrVII']][closest.breakpoint]-cut.pos
    cut.breakpoint.dists[i]= ifelse(!length(breakpoint.dist.to.cut), NA, breakpoint.dist.to.cut)
    cut.breakpoint[i]      = ifelse(!length(breakpoint.dist.to.cut), NA, pos.by.chr[['chrVII']][closest.breakpoint])
    #plot REF reads 
    plot(pos.by.chr[['chrVII']],  
         -REF.d[which(seqnames(vcf.filtered.rd)=='chrVII'),i],
         main=colnames(gdata)[i], 
         type='h', ylim=c(-10,10), 
         ylab='RM reads +, BY reads -', xlab='pos',col='#11111133',
         sub=paste('dist cut to breakpoint = ',  breakpoint.dist.to.cut, 'bp')
         )
    # plot ALT reads
    points(pos.by.chr[['chrVII']],  ALT.d[which(seqnames(vcf.filtered.rd)=='chrVII'),i], type='h',col='#11111133' )
    # plot HMM calls
    points(pos.by.chr[['chrVII']], (hmm.calls.chr)*5, col='orange', type='l', lwd=3)
    # indicate cut site
    abline(v=cut.pos, col='red', lwd=2)
    abline(h=0, lty=2, col='grey')
    #readline() # for step through plots in interactive session
}
dev.off()

#aggregate all calls 
CALL.hmm.genome=do.call('rbind', CALL.hmm)
png(file='~/Desktop/breakpoints_vs_cuts_chrVII_hist.png', width=1024, height=1024)
hist(cut.breakpoint.dists, breaks=100, main='',
     xlab='bp distance from gRNA target to observed breakpoint',
     sub='telomere ------- centromere')
dev.off()
png(file='~/Desktop/breakpoints_vs_cuts_chrVII_hist_zoom.png', width=1024, height=1024)

hist(cut.breakpoint.dists, breaks=1000, main='', xlim=c(-4e4, 4e4),
     xlab='bp distance from gRNA target to observed breakpoint',
     sub='telomere ------- centromere')
dev.off()
cpi=rep(0,384)
for(i in 1:384) {
    strain.name=colnames(gdata)[i]
    strain.name.split=strsplit(strain.name, '-')[[1]]
    strain.name.base=paste(toupper(strain.name.split[1:2]), collapse='-')
    cpi[i]= cut.sites[which(cut.sites[,2]==strain.name.base),1]
}


tcov=REF.d+ALT.d
tcov=apply(tcov,2,sum, na.rm=T)

png(file='~/Desktop/breakpoints_vs_cuts_chrVII.png', width=1024, height=1024)
par(xaxs='i', yaxs='i')
plot(cpi, cut.breakpoint, xlab='gRNA target position', ylab='observed breakpoint', 
     main='(targetted cuts) vs (observed LOH breakpoint) in 384 BYxRM ChrVII hybrid LOH panel',
     xlim=c(0, max( cut.breakpoint, na.rm=T)),
     ylim=c(0, max( cut.breakpoint, na.rm=T)), col=ifelse(tcov>25000,'black', 'red')
     #cut.breakpoint
         )
abline(0,1, lty=2)
dev.off()

png(file='~/Desktop/breakpoints_vs_cuts_chrVII_L_of_centromere.png', width=1024, height=1024)
par(xaxs='i', yaxs='i')
plot(cpi, cut.breakpoint, xlab='gRNA target position', ylab='observed breakpoint', 
     main='(targetted cuts) vs (observed LOH breakpoint) in 384 BYxRM ChrVII hybrid LOH panel',
     xlim=c(0, max( cpi, na.rm=T)+100000),
     ylim=c(0, max( cpi, na.rm=T)+100000),col=ifelse(tcov>25000,'black', 'red')


         )
abline(0,1, lty=2)
dev.off()

cut.breakpoint


# plot whole genome 
pdf(file='/home/jbloom/Dropbox/Public/CRISPR_LOH/CHR_VII_LOH_genome.pdf', width=15, height=8)
for( i in 1:384){
    strain.name=colnames(gdata)[i]
    strain.name.split=strsplit(strain.name, '-')[[1]]
    strain.name.base=paste(toupper(strain.name.split[1:2]), collapse='-')
    cut.pos= cut.sites[which(cut.sites[,2]==strain.name.base),1]
    c7.ind=findInterval(cut.pos, pos.by.chr[['chrVII']])
    plot( -REF.d[,i],main=colnames(gdata)[i], type='h', ylim=c(-10,10),ylab='RM reads +, BY reads -',
         xlab='pos', col='#11111122')
    points(ALT.d[,i], type='h',col='#11111122' )
    points(((CALL.hmm.genome[,i])-1)*5, col='orange', type='l', lwd=3)
    abline(v=c7.anchor+c7.ind, col='red', lwd=2)
    abline(v=cumsum(rle(chr)$lengths), col='blue', lwd=2)
    abline(h=0, lty=2, col='grey')
}
dev.off()

save.image(file='~/Dropbox/Public/CRISPR_LOH/041415.RData')

colnames(CALL.hmm.genome)=colnames(gdata)
rownames(CALL.hmm.genome)=rownames(gdata)
LOH.384.BYxRM=CALL.hmm.genome
save(LOH.384.BYxRM, file='/data/CRISPR_LOH/BYxRM_LOH_segs.RData')


## chr copy number analysis
total.coverage=REF.d+ALT.d
total.cov.sum=apply(total.coverage,2,sum)
total.ref.sum=apply(REF.d,2,sum)
total.alt.sum=apply(ALT.d,2,sum)

cov.by.chr=split(data.frame(total.coverage),chr) 
ref.by.chr=split(data.frame(REF.d),chr) 
alt.by.chr=split(data.frame(ALT.d),chr) 

chr.norm=lapply(cov.by.chr, function(x) apply(x,2,sum)/total.cov.sum )
ref.norm=lapply(ref.by.chr, function(x) apply(x,2,sum)/total.ref.sum )
alt.norm=lapply(alt.by.chr, function(x) apply(x,2,sum)/total.alt.sum )
