library(phytools)
library(igraph)
library(dplyr)
library(cowplot)
library(data.table)
library(ggtree)
library(ggraph)
library(getopt)

#Get call input
spec <- matrix(c('path', 'p', 1, "character",
                 'phylogeny', 't', 1, "character",
                 'gene_pa', 'g', 1, "character",
                 'output', 'o', 1, "character"), byrow=TRUE, ncol=4)
opt <- getopt(spec)
setwd(opt$path)

#Read in
nodstr <- paste(opt$output, "_nodes.tsv", sep="")
nodes  <- read.table(nodstr, header=TRUE, sep="\t", quote=NULL)
colnames(nodes) <- c("alphas", "D")
nodes$alphas <- make.names(nodes$alphas)
edgstr <- paste(opt$output, "_edges.tsv", sep="")
edges <- read.table(edgstr, header=TRUE, sep="\t", quote=NULL)
colnames(edges) <- c("alpha1", "alpha2", "p")
edges$alpha1 <- make.names(edges$alpha1)
edges$alpha2 <- make.names(edges$alpha2)
edges$p <- as.numeric(as.character(edges$p))
#genepa <- read.csv(opt$gene_pa, header=T, row.names=1)
tree <- read.newick(opt$phylogeny)
#Remove any nodes that do not have a value for D and edges which don't have a p-value
nodes <- nodes[complete.cases(nodes),]
edges <- edges[complete.cases(edges),]
edges <- subset(edges, (edges$alpha1 %in% nodes$alphas == TRUE) & (edges$alpha2 %in% nodes$alphas == TRUE))

#Define CCs and order by D-value
g <- graph_from_data_frame(d = edges, directed = FALSE)
CCs <- data.frame(components(g)$membership)
colnames(CCs) <- c("CC")
CCs$alphas <- rownames(CCs)
#Remove any spaces in the alpha names
CCs$alphas <- gsub(" ", "", CCs$alphas)
#Add CC information to nodes table
nodes <- merge(nodes,CCs, by="alphas")
#Reorder by D, then by CC
ord    <- nodes %>% arrange(D) %>% mutate(ord=order(D, decreasing=TRUE))
ord2   <- ord %>% group_by(CC) %>% summarize(min_ord=min(ord))
ord.CC <- merge(ord,ord2)
ord.CC <- ord.CC %>% arrange(min_ord,ord)
#Colour by connectedness
ord.CC.cnts <- count(ord.CC, CC)
tmp <- merge(ord.CC, ord.CC.cnts)
tmp <- tmp %>% arrange(min_ord)
tmp2 <- setDT(tmp)[, unique(n), by = CC]
CC.size <- tmp2$V1
#Define colour array based on the size of the components
colour.array <- c(
"#e6932e", # orange #normally skipped
"#37a2bf", # blue
"#a30b49", # pink
"#d4af37", # gold
"#663ac7", # purple
"#e6932e", # orange
"#0a7a6d", # teal
"#a30b49", # pink
"#37a2bf", # blue
"#e6932e", "#109e57","#a30b49","#34accf","#cc9e08","#006842","#984db8","#3243c7","#7c2fa8","#0245ec","#e3db00","#671a00","#00bf72","#ff9dd6","#f386ff","#d57600","#20000a","#ace967","#840019","#0094f7","#006695","#ff5bf7","#b9007e","#ff89c7","#360093","#001733","#005a69","#d5dca2","#61cb00","#00338b","#26a600","#fdd268","#ffbdca","#b9e590","#ff711b","#8ce9c7","#331300","#a40053","#02d6cf","#640015","#c90039","#fecdb7","#230025","#008c51","#ea21dd","#7a63ff","#b43b00","#0068a7","#01cf41","#f6d633","#93a200","#61f0b9","#ffbf91","#990015","#c4a8ff","#02d9b8","#f2007e","#ac8b00","#d0e072","#ffa368","#02b096","#02c1cc","#ff8de1","#fc0072","#dea3ff","#c37d00","#03c4f4","#ff5dad","#8990ff","#d2d9d9","#b744f6","#380006","#5a3700","#002a23","#392ad1","#238100","#ff673c","#9bddff","#b2e3c0","#89ef68","#00726c","#004376","#5d9fff","#210052","#009e40","#ffac48","#9d1fd7","#d7d5ef","#7cf22c","#b7ddef","#a8a7ff","#008634","#8c3300","#8400a7","#783900","#0099b3","#ff3628","#386700","#003c45","#7f7100","#001f5c","#b4a200","#272200","#ff9ea3","#ff325c","#e19100","#c80024","#ff6a61","#3c5300","#b9d400","#fc2109","#3a3b00","#024fd6","#5d7b00","#a90009","#006027","#6ab200","#016989","#4f005e","#016dbf","#29000b","#002700","#659700","#ffc545","#d7d0ff","#910085","#0064db","#ffa9c8","#dc0034","#a30b49","#ff5664","#cc9e08","#441500","#ff7588","#9ce5da","#ff55cd","#c27fff","#004f31","#e2da91","#7bb4ff","#7c0051","#34accf","#001106","#01ca10","#e56800","#0005a7","#0a3700","#0245ec","#e3db00","#671a00","#00bf72","#ff9dd6","#f386ff","#d57600","#20000a","#ace967","#840019","#0094f7","#006695","#ff5bf7","#b9007e","#ff89c7","#360093","#001733","#005a69","#d5dca2","#61cb00","#00338b","#26a600","#fdd268","#ffbdca","#b9e590","#ff711b","#8ce9c7","#331300","#a40053","#02d6cf","#640015","#c90039","#fecdb7","#230025","#008c51","#ea21dd","#7a63ff","#b43b00","#0068a7","#01cf41","#f6d633","#93a200","#61f0b9","#ffbf91","#990015","#c4a8ff","#02d9b8","#f2007e","#ac8b00","#d0e072","#ffa368","#02b096","#02c1cc","#ff8de1","#fc0072","#dea3ff","#c37d00","#03c4f4","#ff5dad","#8990ff","#d2d9d9","#b744f6","#380006","#5a3700","#002a23","#392ad1","#238100","#ff673c","#9bddff","#b2e3c0","#89ef68","#00726c","#004376","#5d9fff","#210052","#009e40","#ffac48","#9d1fd7","#d7d5ef","#7cf22c","#b7ddef","#a8a7ff","#008634","#8c3300","#8400a7","#783900","#0099b3","#ff3628","#386700","#003c45","#7f7100","#001f5c","#b4a200","#272200","#ff9ea3","#ff325c","#e19100","#c80024","#ff6a61","#3c5300","#b9d400","#fc2109","#3a3b00","#024fd6","#5d7b00","#a90009","#006027","#6ab200","#016989","#4f005e","#016dbf","#29000b","#002700","#659700","#ffc545","#d7d0ff","#910085","#0064db","#ffa9c8","#dc0034","#a30b49","#ff5664","#cc9e08","#441500","#ff7588","#9ce5da","#ff55cd","#c27fff","#004f31","#e2da91");#,"#7bb4ff");
node.colour <- NULL
for(i in c(1:length(CC.size))){
  #if(i > length(colour.array)) {
  #  print("Error: Add more colours to colour.array")
  #  quit()
  #}
  node.colour <- c(node.colour,rep(colour.array[(i%%length(colour.array))+1],times=CC.size[i]))
}
node.order <- unique(ord.CC$alphas)
names(node.colour) <- node.order
node.order <- factor(node.order, levels=node.order)

#Write CCs to coincident_components.csv
CC_out <- data.frame(id = character(0), stringsAsFactors = FALSE)
#CCs$alpha <- make.names(CCs$alpha)
#for(i in c(1:nrow(CCs))) {
#  if (is.na(CC_out[CCs$CC[i],1])) {
#    CC_out[CCs$CC[i],1] <- CCs$alphas[i]
#  } else {
#    CC_out[CCs$CC[i],1] <- paste(CC_out[CCs$CC[i],1], ",", CCs$alphas[i], sep="")
#  }
#}
ord.CC$alphas <- make.names(ord.CC$alphas)
for(i in c(1:nrow(ord.CC))) {
	  curCC <- which(tmp2$CC == ord.CC$CC[i])
  if (is.na(CC_out[curCC,1])) {
	      CC_out[curCC,1] <- ord.CC$alphas[i]
    } else {
	        CC_out[curCC,1] <- paste(CC_out[curCC,1], ",", ord.CC$alphas[i], sep="")
      }
}
outstr <- paste(opt$output, "_components.tsv", sep="")
write.table(CC_out, file=outstr, sep="\t", quote=FALSE, row.names=TRUE, col.names=FALSE)

#Create annot
#genepa[,1:14] <- NULL
#rownames(genepa) <- gsub(" ", "", rownames(genepa))
#annot <- genepa[(rownames(genepa) %in% as.character(names(node.colour))),]
#genepa <- NULL
con = file(opt$gene_pa, "r")
header=readLines(con, n=1) #read in line
header.sp = strsplit(header,"\",\"") #split line into list
header.sp[[1]][1] = gsub("\"", "", header.sp[[1]][1]) #remove first \" from line
header.sp[[1]][length(header.sp[[1]])] = gsub("\"", "", header.sp[[1]][length(header.sp[[1]])]) #and last
header.sp[[1]] <- make.names(header.sp[[1]])
annot <- matrix(ncol=length(header.sp[[1]]))
colnames(annot) <- header.sp[[1]]
flag <- 1
while(TRUE) {
  line=readLines(con, n=1) #read in line
  if (length(line)==0) {
    break #file done
  }
  line.sp = strsplit(line,"\",\"") #split line into list
  line.sp[[1]][1] = gsub("\"", "", line.sp[[1]][1]) #remove first \" from line
  line.sp[[1]][length(line.sp[[1]])] = gsub("\"", "", line.sp[[1]][length(line.sp[[1]])]) #and last
  line.sp[[1]] <- make.names(line.sp[[1]])
  if (line.sp[[1]][1] %in% names(node.colour)) { #its a gene of interest, keep around
    if (flag == 1) {
      annot[1,] <- as.character(line.sp[[1]])
      flag <- 0
    } else {
      annot <- rbind(annot, as.character(line.sp[[1]]))#, stringsAsFactors=FALSE)
    }
  }
}
close(con)
annot <- as.data.frame(annot)
rownames(annot) <- annot[,1]
annot[,1:14] <- NULL

annot <- as.data.frame(t(annot), stringsAsFactors = FALSE)
for(a in 1:ncol(annot)) {
  annot[,a][annot[,a]!="X"] <- as.character(colnames(annot)[a])
}
setcolorder(annot, as.character(names(node.colour)))
#setcolorder(annot, as.character(node.order))

#Draw heatmap
heatmap.breaks <- colnames(annot)
heatmap.breaks <- factor(heatmap.breaks, levels=node.order)
#node.colour[length(node.colour)+1] <- paste0("name", "black")
node.colour  <- setNames(c(node.colour, "white"), c(names(node.colour), "name"))

#
#
#
# GENERATE TREE
#
#
#
tree$tip.label <- make.names(tree$tip.label)
p.tree <- ggtree(tree, layout='circular') +
geom_tiplab(size=5, hjust=0) + #hjust=0.1 
# geom_text(aes(label=node), hjust=-.3) + # add node numbers to tree - used to change highlighting pattern
# shade the tree to highlight different species/genera
# PECTOBACTERIUM
geom_hilight(node=513, fill="#a6004d", extendto=0.55) + # P.brasillience
geom_hilight(node=605, fill="#e3026a", extendto=0.55) + # P.aqua
geom_hilight(node=599, fill="#ff479d", extendto=0.55) + # P.quasi
geom_hilight(node=591, fill="#ff6eb2", extendto=0.55) + # p.polaris
geom_hilight(node=587, fill="#ff82de", extendto=0.55) + # P.parvum
geom_hilight(node=421, fill="#d466b7", extendto=0.55) + # P.versatile
geom_hilight(node=472, fill="#b54597", extendto=0.55) + # P.carv
geom_hilight(node=496, fill="#962678", extendto=0.55) + # P.odor
geom_hilight(node=417, fill="#9926a3", extendto=0.55) + # P.actin
geom_hilight(node=392, fill="#9812a3", extendto=0.55) + #P.parmmentieri
geom_hilight(node=415, fill="#b808c7", extendto=0.55) + # P.punjabense
geom_hilight(node=414, fill="#da17eb", extendto=0.55) + # p.polonicun
geom_hilight(node=39, fill="#f047ff", extendto=0.55) + # P.wasab
geom_hilight(node=382, fill="#f473ff", extendto=0.55) + # P.astrospectum
geom_hilight(node=381, fill="#e68eed", extendto=0.55) + # P.peruo
geom_hilight(node=378, fill="#ce97de", extendto=0.55) + # P.zant
geom_hilight(node=375, fill="#bb81cc", extendto=0.55) + # p.colocasium
# DICKEYA
geom_hilight(node=690, fill="#0f49a6", extendto=0.55) + # d.solani
geom_hilight(node=727, fill="#145fd9", extendto=0.55) + # d.fang
geom_hilight(node=719, fill="#256ee6", extendto=0.55) + # d.diadantii
geom_hilight(node=688, fill="#4388fa", extendto=0.55) + # d.undi
geom_hilight(node=656, fill="#619bfa", extendto=0.55) + # d.dian
geom_hilight(node=649, fill="#79bef7", extendto=0.55) + # d.chry
geom_hilight(node=648, fill="#539edb", extendto=0.55) + # d.poaceiphilia
geom_hilight(node=635, fill="#459ee6", extendto=0.55) + # d.aquatica
geom_hilight(node=261, fill="#3b87c4", extendto=0.55) + # d.lacutsns
geom_hilight(node=618, fill="#2773b0", extendto=0.55) + # d.zeae
geom_hilight(node=613, fill="#1897d6", extendto=0.55) + # d.oryzae
# MUSICOLA
geom_hilight(node=632, fill="#694bd6", extendto=0.55) +
# BRENNERIA
geom_hilight(node=643, fill="#c99a20", extendto=0.55) +
geom_hilight(node=638, fill="#c99a20", extendto=0.55) +
# LONSDALEA
geom_hilight(node=642, fill="#339615", extendto=0.55)

outstr <- paste(opt$output, "_tree", a, ".pdf", sep="")
pdf(outstr,height=58,width=55)
#print(p.out)
print(p.tree)
dev.off()

# 433




arc.breaks = as.numeric(c(0, 0.005, 1))
#Split the heatmap into sections that equal roughly 500 genes each;
#however, must be sure not to split a cluster into parts so must look
#at the size of the clusters and split the input at a location close
#to multiples of 500 based on the cluster size.
CC.cumsum <- cumsum(CC.size)
countie <- 1
a <- 0
while (countie < length(colnames(annot))) {
	#for(i in seq(0,length(colnames(annot)), by=500)) {
	#j <- min(i+500, length(colnames(annot)))
	i <- countie
  	j.loc <- which.min(abs(CC.cumsum - (countie+500)))
  	j <- CC.cumsum[j.loc]
	#However, if there is an element with >500 elements heatmap creation will spin. Inforce j > i.
	if (j<i) {
		j <- CC.cumsum[j.loc+1]
	}
	#However, however, if the distance between i and j is > 1000, we won't be able to see anything
	#if (j-i > 1000) {
	#	j <- i+500
	#}
  	countie <- j+1
	#
	# Edit this code to alter the relative placement of the heatmap to the tree
	#
	p.heat <- gheatmap(p.tree, annot[i:j], offset=0.1, width=2.4, font.size=10, colnames_angle=-90, hjust=0) +
	  guides(fill=FALSE) +
	  scale_fill_manual(breaks=heatmap.breaks, values=node.colour) +
	  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
	  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

	outstr <- paste(opt$output, "_heatmap", a, ".pdf", sep="")
	a <- a + 1
	pdf(outstr,height=200,width=60)
	#print(p.out)
	print(p.heat)
	dev.off()
}
