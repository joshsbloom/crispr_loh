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


#Supplementary functions for image processing
# Part 1 functions start ##################################################################################
# make grid for expected colony locations
makeGrid = function(plate.type, corners) {
    x.c=sort(corners$X)
    y.c=sort(corners$Y)
    av=list( xmin=round(mean(x.c[c(1,2)])), xmax=round(mean(x.c[c(3,4)])),
                     ymin=round(mean(y.c[c(1,2)])), ymax=round(mean(y.c[c(3,4)])))
    gridme = expand.grid( 
                    round(seq(av$xmin, av$xmax, length.out=max(plate.type)) ), 
                    round(seq(av$ymin, av$ymax, length.out=min(plate.type))))
    colnames(gridme)=c('X','Y') 
    return(gridme)
}     

# segment plate into features and non-features using kmeans clustering
segmentPlate=function(img, spot.seeds){
    print('Segmenting plate')
    tplate=img

#    km = kmeans(as.vector(img), centers=c(.3906,.7812))
    km = kmeans(as.vector(img), centers=c(.13,.7812))
    kmin = which.min(km$centers)
    kmax = which.max(km$centers)
    if(kmin==2) {km$cluster = abs(km$cluster-3)}
    tplate@.Data=matrix(km$cluster-1,nrow(tplate),ncol(tplate))

   # EBImage::display(matrix(km$cluster-1, 5184,3456))

    seeds = matrix(0, nrow(img), ncol(img))
    seeds[is.numeric(seeds)]=0
   # ycirc    <- makeBrush(9, 'disc', step=TRUE)
    # brush was 9
    ycirc    <- makeBrush(49, 'disc', step=TRUE)

    for(i in 1:nrow(spot.seeds) ){
            x = spot.seeds[i,'X']
            y = spot.seeds[i, 'Y']
            #seeds[(x-4):(x+4), (y-4):(y+4) ] =ycirc*i
            seeds[(x-24):(x+24), (y-24):(y+24) ] =ycirc*i
    }

    # find voroni region for each seed
    plate.mask = propagate(img, seeds, mask=tplate)
    #EBImage::display(plate.mask)
 
    return(plate.mask)
}

# calculate colony features and attach names
getColonyFeatures=function(plate.mask, img, layout.file){
    print('Getting features')
    strain.layout=read.delim(layout.file, sep='\t', header=F, stringsAsFactors=F)
    mfeatures=computeFeatures.moment(plate.mask, img )
    sfeatures=computeFeatures.shape(plate.mask)
    nfeatures=data.frame(mfeatures, sfeatures)
    rownames(nfeatures)=paste(seq(1,384), as.vector(t(strain.layout)), sep=':')
    #as.vector(t(strain.layout))
    return(nfeatures)
}

# spline fit for position effects
# calculate colony features and attach names
localPfit = function(vecin, x,y) {
   m = locfit(vecin~lp(x,y, scale=FALSE),lfproc=locfit.robust)
   predx = predict(m, data.frame(vecin, x,y))
   return(vecin-predx) }

# do everything
processImages = function(i, key, image.path, layout.path, outdir, plate.type, corners){
    print(i)
    image.file   = paste(image.path, key$Filename[i], sep='')
    layout.file  = paste(layout.path, key$StrainLayout[i], sep='')
    raw.img      = rotate( readImage(image.file), 180)
    #EBImage::display(raw.img)
    img          = EBImage::channel(raw.img, 'grey')
    spot.seeds   = makeGrid(plate.type, corners)
    
    plate.mask   = segmentPlate(img,spot.seeds)

    features     = getColonyFeatures(plate.mask, img, layout.file)
    #display(paintObjects(plate.mask, raw.img ), title=key$Filename[i] )
    writeImage(paintObjects(plate.mask, raw.img ), file=paste(outdir, i, '.jpeg', sep=''))
    return(features)
}
#  text files ... 
writeTXTfiles = function(outD, Results.norm) {
    dir.create(outDir)
    for(n in names(Results.norm)){
        n2 =n
        n2=gsub('/', '__', n2)
        n2=gsub(':', '_', n2)
        write.table(Results.norm[[n]], file=paste(outDir,n2,'.txt',sep=''), 
                sep='\t', row.names=T,col.names=NA,quote=F) }
}
################### Part 1 functions end #####################################################################


# Part 2 functions start ####################################################################################
# reconstruct key from 
reconstructKey=function(Results) {
    ndf=c("StrainLayout" ,    "PermutationGroup", "Filename"   ,      "Condition"     ,   "Concentration"  ,  "Control")
    df=data.frame(do.call('rbind', strsplit(names(Results), '::'))) 
    names(df)=ndf
    return(ndf)
}

getResults.df=function(Results, key, ind, plate.type=384) {
    R_list=Results.processed[ind]
    key_sub=key[ind,]
    #key_sub$Concentration[key_sub$Concentration=='0']=''
    plate.factor=as.factor(rep(key_sub$Filename, each=plate.type))
    cond_conc.factor=rep(paste(key_sub$Condition, key_sub$Concentration, sep=';'), each=plate.type)
    R_list.df=do.call('rbind', R_list)
    strain.names = as.vector(unlist(lapply(R_list, function(x) rownames(x) )))
    strain.names=do.call('rbind', strsplit(strain.names, ':') )[,2] 
    R_list.df=data.frame(strain=strain.names, R_list.df, plate=plate.factor, condition=cond_conc.factor, stringsAsFactors=F)
    return(R_list.df)
}

filterEdges=function(Results, plate.type=c(16,24)){
    rows.n=plate.type[1]
    cols.n=plate.type[2]
    col.row = expand.grid(seq(1,cols.n), seq(1,rows.n))
    edges = list(    
        top.row = col.row[,2]==1,
       bottom.row = col.row[,2]==rows.n,
       left.col = col.row[,1]==1,
       right.col =col.row[,1]==cols.n)                                                      

    for( i in 1:length(Results) ){
        pa = Results[[i]]
        pagn = pa[,'s.radius.mean']
        e.ps= sapply(edges, function(x) {
                     tryCatch( {t.test(pagn[x], pagn)$p.value}, error=function(e) {return(1) } )  
         } )
     for(j in 1:4){   if(e.ps[j]<.05) {Results[[i]][edges[[j]], -c(1,2)]=NA;   } } }
    return(Results)
}

# within plate normalization --------------------------------------------------------------------------------------------
# this would be a good spot to adjust for plate effects, if desired
normalizePlate=function(Results){
   Results.processed=list()
   for( i in 1:length(Results) ){
        print(i)
        pa = Results[[i]]
        normalized =  apply(pa[,c('s.radius.mean', 's.area')], 2, function(x){ 
                    tryCatch( { return(localPfit(x, pa[,'m.cx'], pa[,'m.cy'])+mean(x,na.rm=T))}, 
                       error = function(e){return(x)} ) })
        colnames(normalized)=paste(colnames(normalized),'norm', sep='.')
        r = cbind(pa, normalized)
        Results.processed[[i]]=r
    }
    names(Results.processed)=names(Results)
    return(Results.processed)
}   
#-------------------------------------------------------------------------------------------------------------------------



#################### Part 2 functions end ######################################################################


# Part 3 Mapping functions ######################################################################################
processPhenos_for_MM=function(s, all.strain.names){
    #s=(pheno_raw[[phenotype]])
    # fix this -----------------------------------
    s=s[match(all.strain.names, names(s))]
    srle=sapply(s, length)
    sr=list()
    sr$lengths=as.vector(srle)
    sr$values=names(srle)
    attr(sr, 'class')='rle'
    ny=as.vector(unlist(s))
    names(ny)=inverse.rle(sr)
    #-----------------------------------------------
    y=ny[!is.na(ny)]
    strain.names=(names(y))
    unique.sn=unique(strain.names)
    n.to.m=match(strain.names, unique.sn)
    strain.ind  = seq_along(strain.names)
    strain.cnt  = length(unique.sn)
   
    #for constructing Strain Variance component
    Strain      = diag(strain.cnt)
    Z=Matrix(0, length(y), strain.cnt,sparse=T);   Z[cbind(strain.ind, n.to.m)]=1 
    strains.with.phenos=match(unique.sn, all.strain.names)
    n.strains=length(strains.with.phenos)
    return(list(y=y, Z=Z, Strain=Strain, strains.with.phenos=strains.with.phenos, n.strains=n.strains))
}  


# MAPPING
extractScaledPhenotype=function(impcross,scaleVar=FALSE){apply(impcross$pheno, 2, scale, scale=scaleVar)}
extractGenotype=function(impcross){ (do.call('cbind', sapply(impcross$geno, function(x) { x$data }))*2)-3 }
extractGenotype.argmax=function(impcross){ (do.call('cbind', sapply(impcross$geno, function(x) { x$argmax }))*2)-3 }
countStrainsPerTrait = function(pheno) {apply(pheno, 2, function(x){sum(!is.na(x))})}

get.LOD.by.COR = function(n.pheno, pheno, gdata,betas=FALSE, sdx=1) {
   # Lynch and Walsh p. 454
   r=cor(pheno, gdata, use='pairwise.complete.obs')
   LOD=(-n.pheno*log(1-r^2))/(2*log(10))
   if(betas==TRUE) {
       # convert pearson R to B from y=u+Bx+e model
       beta=r * apply(cbind(pheno),2, sd,na.rm=T)/sdx
       return(list(r=r, LOD=LOD, beta=beta))
   }
   else {   return( list(r=r, LOD=LOD) ) } 
}

calc.pval.LRtest=function(null,full) { pchisq(-2*(null-full),1, lower.tail=FALSE) }
calc.BLUPS= function(G,Z,Vinv,y,X,B ){    G%*%t(Z)%*%Vinv%*%(y- X%*%B)     }

extractVarCompResults = function(r) {list(sigma=r$Var, sigma.cov=r$invI, W=r$W,
                                          Bhat=as.vector(r$Bhat), llik=as.vector(r$llik)) }

# can generalize for more than two groups if desired
structured.perm=function(y, bck.split){
    yperm=y
    yperm[bck.split]=sample(y[bck.split])
    yperm[!bck.split]=sample(y[!bck.split])
    return(yperm)
}

# calculate averages for a named vector with replicates
avg.nvec=function(y, un=NULL){
    if(is.null(un)) {un=names(y) }
    y.avg.l=(by(y, un, mean, na.rm=T))
    nvec=names(y.avg.l)
    y.avg.l=as.vector(y.avg.l)
    names(y.avg.l)=nvec
    return(y.avg.l)
}

dominanceScan=function(ylist, l.geno.s, l.geno.s.d, search.space) {
    print("scanning for dominance QTL")
    X=matrix(NA, nrow=ylist$n.strains, ncol=2)
    dom.model=list()
    add.model=list()
    dom.v.add.F.nlps=rep(NA,2000)
    #meh, this is fast enough, could be more fancy here
    for(i in 1:2000) {
        if(i%%100==0) print(i)
        X[,1]=l.geno.s[,i]
        X[,2]=l.geno.s.d[,i]
        dom.model[[i]]=lm(y.avg.l~X)
        add.model[[i]]=lm(y.avg.l~X[,1])
        dom.v.add.F.nlps[i]= -log10(anova(add.model[[i]], dom.model[[i]])$'Pr(>F)'[2])
    }
    dominance.scan=list(search.space=search.space,
                        add.model=add.model,
                        dom.model=dom.model,
                        dom.v.add.F.nlps=dom.v.add.F.nlps)

}

# X=matrix(NA, nrow=ylist.L$n.strains, ncol=2)
#    dom.model=list()
#    add.model=list()
#    dom.v.add.F.nlps=rep(NA,2000)
#    #meh, this is fast enough, could be more fancy here
#    for(i in 1:2000) {
#        print(i)
#        X[,1]=l.geno.s[,i]
#        X[,2]=l.geno.s.d[,i]
#        dom.model[[i]]=lm(y.avg.l~X)
#        add.model[[i]]=lm(y.avg.l~X[,1])
#        dom.v.add.F.nlps[i]= -log10(anova(add.model[[i]], dom.model[[i]])$'Pr(>F)'[2])
#    }
    # uncomment to retain results of additive and dominance scans
    # these are stored in a grossly inefficient way... advise recomputing if you see something interesting
    # from -log10(p)
#    dominance.scan=list(search.space=search.space,
                        #add.model=add.model,
                        #dom.model=dom.model,
#                        dom.v.add.F.nlps=dom.v.add.F.nlps)


makePlots.1=function(pdf.file, stats.galore, chr='chrVII') {
       
    pdf(file=pdf.file, width=16, height=14)
    for(pp in 1:length(stats.galore)){
    
    pheno=names(stats.galore)[pp]
    ylist.S=stats.galore[[pp]]$ylist.S
    ylist.L=stats.galore[[pp]]$ylist.L
    s.mm =stats.galore[[pp]]$s.mm
    S.genome.LOD = stats.galore[[pp]]$S.genome.LOD
    S.br7.LOD    = stats.galore[[pp]]$ S.br7.LOD
    L.genome.LOD = stats.galore[[pp]]$ L.genome.LOD
    L.genome.LOD.B=stats.galore[[pp]]$ L.genome.LOD.B   
    L.genome.LOD.R=stats.galore[[pp]]$L.genome.LOD.R
    domList=stats.galore[[pp]]$domList
    search.space=domList$search.space

        #x11()
        # genome wide LOD plot
        par(mfrow=c(3,1))
        plot(S.genome.LOD$LOD[1,],
             main=paste(pheno, ylist.S$n.strains, "BYxRM Segregants", "trait avg"),
             ylab='LOD')
        abline(v=cumsum(rle(S.chr)$lengths), col='blue')
        plot(S.br7.LOD$LOD[1,],
             main=paste(pheno, ylist.S$n.strains, "BYxRM Segregants", "(chr7 loco) avg blup residuals") , 
             ylab='LOD' )
        abline(v=cumsum(rle(S.chr)$lengths), col='blue')
        #par(mfrow=c(2,1))
        plot(L.genome.LOD$LOD[1,], main=paste(pheno, ylist.L$n.strains, "BYxRM LOH Segregants"),ylab='LOD')
        abline(v=cumsum(rle(L.chr)$lengths), col='blue')

        # ChrVII LOD plot
        #par(mfrow=c(3,1))
        #  overlapping  chrVII LOD scores 
        par(mfrow=c(1,1))
       ymax=max( c(S.genome.LOD$LOD[1,which(S.chr=='chrVII')], 
                   S.br7.LOD$LOD[1,which(S.chr=='chrVII')] , 
                   L.genome.LOD$LOD[1,which(L.chr=='chrVII')]  
                   ), na.rm=T)
        plot(S.pos[which(S.chr=='chrVII')],
             S.genome.LOD$LOD[1,which(S.chr=='chrVII')],
             main=paste(pheno, 'chrVII'),
             ylab='LOD',
             xlab='pos' , ylim=c(0,ymax), col='grey' ,type='l', lwd=1.5)

        #abline(v=cumsum(rle(S.chr)$lengths), col='blue')
       points(S.pos[which(S.chr=='chrVII')],
             S.br7.LOD$LOD[1,which(S.chr=='chrVII')], 
             #main=paste(pheno, ylist.S$n.strains, "BYxRM Segregants", "chr7 loco avg blup residuals") ,
             ylab='LOD',xlab='pos' , col='black', type='l', lwd=1.5)
       # fine to use threshold from other object ... bit on chr VII is identical 
       abline(h=attr(S.br7.LOD.c7, 'thresh'), col='black')

        #abline(v=cumsum(rle(S.chr)$lengths), col='blue',ylim=c(-2,0))
        points(L.pos[which(L.chr=='chrVII')], 
             L.genome.LOD$LOD[1,which(L.chr=='chrVII')], 
            # main=paste(pheno, ylist.L$n.strains, "BYxRM LOH Segregants", 'chrVII'),
             ylab='LOD',xlab='pos', col='red', type='l', lwd=1.5)
       # fine to use threshold from other object ... bit on chr VII is identical 
        abline(h=attr(L.LOD.c7, 'structured_thresh'), col='red')
        abline(h=attr(L.LOD.c7, 'thresh'), col='pink')

         legend('topright', legend=c('BYxRM segregants (blup residuals)', 
                                    'BYxRM segregants (scanone on strain averages)',
                                    'LOH panel' ), 
                                    col=c('black','grey','red'), cex=1.1, pch=20,
                                    text.col=c('black','grey','red'))

       # ChrVII effect size plot 
       par(mfrow=c(1,1))

       beta.range=range(c(S.genome.LOD$beta[1,which(S.chr=='chrVII')],
                          S.br7.LOD$beta[1,which(S.chr=='chrVII')],
                          L.genome.LOD$beta[1,which(L.chr=='chrVII')]))
       
       plot(S.pos[which(S.chr=='chrVII')],
             S.genome.LOD$beta[1,which(S.chr=='chrVII')],
             main=paste(pheno, 'chrVII'),
             ylab='effect size B',
             xlab='pos' ,ylim=beta.range,type='l', lwd=1.5, col='grey')
       points(S.pos[which(S.chr=='chrVII')],
             S.br7.LOD$beta[1,which(S.chr=='chrVII')], 
             ylab='effect size B',type='l', lwd=1.5, col='black' )
       points(L.pos[which(L.chr=='chrVII')], 
             L.genome.LOD$beta[1,which(L.chr=='chrVII')],
             ylab='effect size B',type='l', lwd=1.5, col='red')
        legend('topright', legend=c('BYxRM segregants (blup residuals)', 
                                    'BYxRM segregants (scanone on strain averages)',
                                    'LOH panel' ), 
                                    col=c('black','grey','red'), cex=1.1, pch=20,
                                    text.col=c('black','grey','red'))

       # effect size split by strain background
       par(mfrow=c(1,1))
       plot(L.pos[which(L.chr=='chrVII')], 
             L.genome.LOD$beta[1,which(L.chr=='chrVII')],
             ylab='effect size B',type='l', lwd=1.5, col='red',ylim=beta.range, xlab='pos',
             main=paste(pheno, 'chrVII', 'LOH effect sizes')
             )

        points(L.pos[which(L.chr=='chrVII')],
            L.genome.LOD.B$beta[1,which(L.chr=='chrVII')],type='l', lwd=1.5, col='orange')

        points(L.pos[which(L.chr=='chrVII')],
               L.genome.LOD.R$beta[1,which(L.chr=='chrVII')],type='l', lwd=1.5, col='purple')
        legend('topright', legend=c('LOH panel',
                                    'LOH panel BY_GFP',
                                    'LOH panel RM_GFP'
                                    ), 
                                    col=c('red', 'orange', 'purple'), cex=1.1, pch=20,
                                    text.col=c('red', 'orange', 'purple'))


       par(mfrow=c(1,1))
       plot(L.pos[search.space],  domList$dom.v.add.F.nlps, main='Dominance vs Additive Model -log10(p)',
            xlab='pos (right end of chrVII only)',
            ylab='-log10(p)', type='l'
            )

    }
    dev.off()
}



