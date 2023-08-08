# graph blast bit scores of kogs with self and non-self hits

#generated data with:
#blastplustable.py -q Core_genes_Parra.fta -d ~/ncbi-blast-2.2.29+/db/Core_genes_Parra.fta -b blastp -M -m 16 > Core_genes_Parra_self_blasttable_with_errors.txt
btxt = "~/cogdata/Core_genes_Parra_self_blasttable.txt"
btxt2 = "~/cogdata/Core_genes_Parra_self_blasttable_with_errors.txt"

# 2022-05-19
# this would probably be better redone with blastp -outfmt 6

#pdf(file="~/cogdata/parra_blast_v_bitscore_chart_v3.pdf", width=6, height=6)

alldat = read.table(btxt2,header=TRUE,sep='\t',stringsAsFactors=FALSE)

selfhits = alldat[alldat$Resulttitle==alldat$Query,]
otherhits = alldat[alldat$Resulttitle!=alldat$Query,]
#
#get the n-th items from a list of vectors, use sapply
koghits = otherhits[sapply(strsplit(otherhits$Resulttitle,"_"),"[[",4)==sapply(strsplit(otherhits$Query,"_"),"[[",4),]
nonkoghits = otherhits[sapply(strsplit(otherhits$Resulttitle,"_"),"[[",4)!=sapply(strsplit(otherhits$Query,"_"),"[[",4),]

xmax=500
ymax=2*xmax
#ymax=3000

#plot query length vs bitscore
plot(selfhits[,12], selfhits[,3], xlab="query length", ylab="bitscore", xlim=c(0,xmax), ylim=c(0,ymax), type='p', col='#00000080', frame.plot=FALSE, main="KOG all-v-all BLAST")
points(koghits[,12], koghits[,3], col=rainbow(7438,alpha=0.5), pch=20)
points(nonkoghits[,12], nonkoghits[,3], col="#00000044", pch=20)
lines(c(0,2500),c(200,200),lty=2)
lines(c(100,100),c(0,1000),lty=2)
lines(c(0,2500),c(0,1000),lty=3)
legend(100,ymax*0.9, legend=c("self-hits","same-kog hits","off-kog hits"), col=c("black","red","darkgray"), pch=c(1,20,20))

#plot(selfhits[,12], selfhits[,3]/(selfhits[,12]), xlab="query length", ylab="bitscore", xlim=c(0,xmax), ylim=c(0,2), type='p', col='#00000080', frame.plot=FALSE, main="KOG all-v-all BLAST")
#points(koghits[,12], koghits[,3]/(koghits[,12]), col=rainbow(7438,alpha=0.5), pch=20)
#points(nonkoghits[,12], nonkoghits[,3]/(nonkoghits[,12]), col="#00000044", pch=20)
#lines(c(100,100),c(0,1000),lty=2)
#lines(c(0,2500),c(0.4,0.4),lty=2)

#
#dev.off()