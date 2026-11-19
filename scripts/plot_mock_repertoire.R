# Purpose: plot final generated-repertoire summaries without rerunning simulation.
# Context: descriptive raw-generator truth; no empirical light-chain comparison.
a <- commandArgs(trailingOnly=TRUE)
if(length(a)!=2L) stop('Usage: Rscript plot_mock_repertoire.R summary_directory fresh_output_directory')
stopifnot(requireNamespace('jsonlite',quietly=TRUE),requireNamespace('digest',quietly=TRUE),requireNamespace('ragg',quietly=TRUE))
src <- a[1];out <- a[2];if(file.exists(out)) stop('Fresh output directory required')
m <- jsonlite::fromJSON(file.path(src,'manifest.json'));stopifnot(m$status=='complete')
for(n in names(m$artifact_sha256)) stopifnot(digest::digest(file.path(src,n),algo='sha256',file=TRUE)==m$artifact_sha256[[n]])
read <- function(n) read.delim(file.path(src,paste0(n,'.tsv')),stringsAsFactors=FALSE)
pop <- read('population_counts');ids <- read('identity_counts');lens <- read('length_distributions');genes <- read('gene_usage');shm <- read('shm_distributions');families <- read('family_sizes')
loc <- c('IGH','IGK','IGL');lab <- c('Heavy','Kappa','Lambda');cols <- c('#7898B8','#B0A0BF','#A2B6AC');names(cols)<-loc
stopifnot(identical(ids$locus,loc),sum(ids$chain_records)==m$observed_chains,all(ids$chain_records>0))
ns <- setNames(ids$chain_records,loc)
for(l in loc)stopifnot(sum(pop$observed_chain_records[pop$locus==l])==ns[l])
len <- subset(lens,quantity=='anchored_junction_aa_length');g <- subset(genes,gene_namespace=='native_GenAIRR_gene_calls' & segment %in% c('V','J'))
s <- subset(shm,quantity=='root_substitutions' & component=='related_shm_family')
for(l in loc){stopifnot(sum(len$count[len$locus==l])==ns[l]);for(k in c('V','J'))stopifnot(sum(g$numerator[g$locus==l & g$segment==k])==ns[l])}
dir.create(out,recursive=TRUE)
for(n in names(m$artifact_sha256))file.copy(file.path(src,n),file.path(out,n),overwrite=FALSE)
write.table(data.frame(locus=loc,records=ns,junction_AA_classes=ids$anchored_junction_AA_classes,variable_IL_classes=ids$variable_AA_IL_classes),file.path(out,'panel_b_source.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
ranked <- do.call(rbind,lapply(c('V','J'),function(k)do.call(rbind,lapply(loc,function(l){z<-g[g$locus==l & g$segment==k,];z<-z[order(-z$numerator,z$gene_call),];z$rank<-seq_len(nrow(z));z}))))
write.table(ranked,file.path(out,'gene_rank_source.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
fmt <- function(x)format(x,big.mark=',',scientific=FALSE,trim=TRUE)
panel <- function(letter,title){title(main=title,adj=0,font.main=2,cex.main=.93,line=1.0);mtext(letter,side=3,line=1.0,adj=-.16,font=2,cex=1)}
axes <- function(){axis(1);axis(2,las=1);box(bty='l')}
curve_panel <- function(z,x,y,letter,title,xlab,ylab){
 plot(range(z[[x]]),c(0,max(z[[y]])*1.13),type='n',xlab=xlab,ylab=ylab,axes=FALSE);axes();panel(letter,title)
 for(l in loc){q<-z[z$locus==l,];q<-q[order(q[[x]]),];lines(q[[x]],q[[y]],col=cols[l],lwd=1.5,lty=match(l,loc))}
 legend('topright',lab,col=cols,lty=1:3,lwd=1.5,bty='n',cex=.73)
}
draw <- function(){
 par(mfrow=c(2,3),mar=c(4.1,4.3,2.7,.8),oma=c(1.5,0,0,0),mgp=c(2.7,.65,0),tcl=-.22,cex=.83,family='sans',las=1)
 mat <- rbind(sapply(loc,function(l)pop$observed_chain_records[pop$locus==l & pop$component=='unmutated_singleton']),sapply(loc,function(l)pop$observed_chain_records[pop$locus==l & pop$component=='related_shm_family']))
 bp<-barplot(mat,col=c('#D6D9DD','#7898B8'),border=NA,names.arg=lab,ylim=c(0,max(ns)*1.32),ylab='Generated chain records',axes=FALSE);axis(2,at=seq(0,40000,10000),labels=c('0','10k','20k','30k','40k'),las=1);box(bty='l');panel('a','Observed model population')
 text(bp,ns+max(ns)*.04,fmt(ns),cex=.78)
 legend('topright',c('Independent singleton','SHM family descendant'),fill=c('#D6D9DD','#7898B8'),bty='n',cex=.69)
 fractions<-rbind(ids$anchored_junction_AA_classes/ns,ids$variable_AA_IL_classes/ns)
 bp<-barplot(fractions,beside=TRUE,col=c('#A69BB6','#7898B8'),border=NA,names.arg=lab,ylim=c(0,1.27),ylab='Unique classes / chain records',axes=FALSE);axis(2,at=seq(0,1,.25),las=1);box(bty='l');panel('b','Sequence diversity')
 text(bp,fractions+.035,fmt(rbind(ids$anchored_junction_AA_classes,ids$variable_AA_IL_classes)),srt=90,adj=0,cex=.60)
 legend('topleft',c('Exact junction AA','Full variable AA (I/L)'),fill=c('#A69BB6','#7898B8'),bty='n',cex=.68)
 len$percent<-100*len$count/len$denominator
 curve_panel(len,'value','percent','c','Junction lengths','Anchored junction length (AA)','Chain records (%)')
 for(k in c('V','J')){z<-ranked[ranked$segment==k,];z$percent<-100*z$numerator/z$denominator;curve_panel(z,'rank','percent',if(k=='V')'d' else 'e',paste(k,'gene usage'),paste(k,'gene rank (within locus)'),'Chain records (%)')}
 s$percent<-100*s$count/s$denominator
 curve_panel(s,'value','percent','f','Within-family mutation burden','Substitutions from latent root (nt)','Family descendant records (%)')
 mtext('Generated model truth; latent ancestors excluded. Junction includes conserved anchors; CDR3 excludes them.',side=1,outer=TRUE,cex=.71,line=.25)
}
# PDF retains text; native Cairo SVG embeds outlined glyphs, as in existing figures.
grDevices::cairo_pdf(file.path(out,'mock_repertoire.pdf'),width=183/25.4,height=136/25.4,pointsize=8);draw();dev.off()
grDevices::svg(file.path(out,'mock_repertoire.svg'),width=183/25.4,height=136/25.4,pointsize=8);draw();dev.off()
ragg::agg_png(file.path(out,'mock_repertoire.png'),width=183,height=136,units='mm',res=300,pointsize=8);draw();dev.off()
files <- list.files(out,full.names=TRUE)
receipt<-list(status='complete',backend='R',dimensions_mm=c(183,136),input_manifest_sha256=digest::digest(file.path(src,'manifest.json'),algo='sha256',file=TRUE),script_sha256=digest::digest('scripts/plot_mock_repertoire.R',algo='sha256',file=TRUE),artifact_sha256=setNames(lapply(files,digest::digest,algo='sha256',file=TRUE),basename(files)),session=capture.output(sessionInfo()),SVG_text='outlined glyphs; PDF text retained')
jsonlite::write_json(receipt,file.path(out,'manifest.json'),pretty=TRUE,auto_unbox=TRUE)
