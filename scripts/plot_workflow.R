# Purpose: draw the declared AIRR/proteomics architecture from documented contracts.
# Context: architecture only; no simulation, calculated result or completion claim.
a<-commandArgs(trailingOnly=TRUE)
if(length(a)!=1L || file.exists(a[1]))stop('Usage: Rscript scripts/plot_workflow.R fresh_output_directory')
stopifnot(requireNamespace('ragg',quietly=TRUE),requireNamespace('jsonlite',quietly=TRUE),requireNamespace('digest',quietly=TRUE))
out<-a[1];dir.create(out,recursive=TRUE)
nodes<-data.frame(id=c('airr','database_a','physical_a','generation','reads','arda','database_b','proteins','truth','digest','search','exports','fit','evaluation'),x=c(.16,.55,.55,.13,.39,.65,.88,.36,.76,.20,.50,.81,.28,.76),y=c(.87,.92,.82,.65,.65,.65,.65,.45,.45,.25,.25,.25,.08,.08),w=c(.26,.38,.38,.23,.23,.23,.20,.36,.35,.29,.25,.28,.45,.43),h=c(.10,.085,.085,.13,.13,.13,.13,.09,.07,.11,.11,.11,.10,.10),label=c('Experimental AIRR\nSequences as supplied','D: fixed candidate database\n33,052 supplied records','P: physical protein substrate\n12 H references OR all 33,052','GenAIRR 2.4.0\nH / kappa / lambda\nSHM families + singletons','InSilicoSeq 2.0.1\nAmplicon OR\ntranscript-aware RNA','fastp + ARDA 2.36.0\nRetain unmerged reads\nNative correction','D: observed H/K/L DB\nObserved segments\n+ public constants','P: physical H/K/L proteins\nConstants + declared protein amounts','T: generating truth retained\nEvaluation only','P → separate protease acquisitions\nTrypsin | chymotrypsin\nMSSimulator (OpenMS 3.0.0)','Internal spectra + D\nSimpleSearchEngine\n(OpenMS 3.0.0)','Spectrum-linked peptide exports\nValidate contributor linkage\nThen discard spectra','Yield + per-reference region coverage\nSelect assay settings / conditional reconstruction\nH/K/L capture: descriptive compatibility','Heavy-chain scoring; synthetic evaluation uses T\nPeptide accuracy; H identity recovery\nPrecision/recall, AP/ROC-AUC\nRecovery excluded from assay selection'),stringsAsFactors=FALSE)
edges<-data.frame(from=c('airr','airr','generation','reads','arda','generation','generation','digest','search','exports','exports'),to=c('database_a','physical_a','reads','arda','database_b','proteins','truth','search','exports','fit','evaluation'),type=c(rep('data',6),'truth_only',rep('data',4)))
write.table(nodes,file.path(out,'nodes.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(edges,file.path(out,'edges.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
arrow<-function(x,y,col='#68737D',lty=1){lines(x,y,col=col,lwd=.85,lty=lty);n<-length(x);arrows(x[n-1],y[n-1],x[n],y[n],length=.055,angle=23,col=col,lwd=.85,lty=lty)}
draw<-function(){
 par(mar=c(.1,.1,.1,.1),xaxs='i',yaxs='i',family='sans');plot.new();plot.window(c(0,1),c(0,1))
 rect(.005,.765,.995,.995,col='#F7F8F9',border=NA);rect(.005,.365,.995,.755,col='#F7F8F9',border=NA);rect(.005,.005,.995,.355,col='#F7F8F9',border=NA)
 text(.025,.975,'a  Primary: experimental AIRR-conditioned MS',adj=0,font=2,cex=.9)
 text(.025,.735,'b  Complementary: known-truth synthetic repertoires',adj=0,font=2,cex=.9)
 text(.025,.335,'c  Shared proteomics and evaluation architecture',adj=0,font=2,cex=.9)
 arrow(c(.29,.33,.33,.36),c(.89,.89,.92,.92));arrow(c(.29,.33,.33,.36),c(.85,.85,.82,.82))
 arrow(c(.245,.275),c(.65,.65));arrow(c(.505,.535),c(.65,.65));arrow(c(.765,.78),c(.65,.65))
 arrow(c(.13,.13,.18),c(.585,.45,.45));arrow(c(.08,.08,.76,.76),c(.585,.39,.39,.415),col='#A69BB6',lty=2)
 arrow(c(.345,.375),c(.25,.25));arrow(c(.625,.67),c(.25,.25))
 arrow(c(.81,.81,.28,.28),c(.195,.16,.16,.13));arrow(c(.81,.81,.76,.76),c(.195,.16,.16,.13))
 for(i in seq_len(nrow(nodes))){n<-nodes[i,];fill<-if(n$id=='truth')'#EAE5EF' else if(n$id %in% c('fit','evaluation'))'#E6EDE9' else '#E4EAF0';rect(n$x-n$w/2,n$y-n$h/2,n$x+n$w/2,n$y+n$h/2,col=fill,border='#99A6B1',lwd=.6);text(n$x,n$y,n$label,cex=.70)}
 text(.76,.543,'Matched reconstruction only: append 12 H queries + reported VL auxiliaries.\nGeneric synthetic backgrounds remain entirely FASTQ-derived.',cex=.64)
 text(.76,.383,'Dashed branch: truth never supplies observed sequences.',cex=.62,col='#675A75')
 text(.50,.178,'Ports D / P select the matching branch database / physical substrate; no truth-to-database link.',cex=.62)
}
grDevices::cairo_pdf(file.path(out,'workflow.pdf'),width=183/25.4,height=166/25.4,pointsize=9);draw();dev.off()
grDevices::svg(file.path(out,'workflow.svg'),width=183/25.4,height=166/25.4,pointsize=9);draw();dev.off()
ragg::agg_png(file.path(out,'workflow.png'),width=183,height=166,units='mm',res=300,pointsize=9);draw();dev.off()
f<-list.files(out,full.names=TRUE);jsonlite::write_json(list(status='complete',type='declared_architecture_not_results',backend='R',dimensions_mm=c(183,166),source_script_sha256=digest::digest('scripts/plot_workflow.R',algo='sha256',file=TRUE),artifact_sha256=setNames(lapply(f,digest::digest,algo='sha256',file=TRUE),basename(f)),session=capture.output(sessionInfo()),SVG_text='outlined glyphs; PDF retains text'),file.path(out,'manifest.json'),pretty=TRUE,auto_unbox=TRUE)
