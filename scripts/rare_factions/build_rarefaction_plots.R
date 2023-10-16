#!/usr/bin/env Rscript
#
# (c) University of St Andrews 2023
# (c) University of Strathclyde 2023
# (c) James Hutton Institute 2023
#
# Author:
# Emma E. M. Hobbs
#
# Contact
# eemh1@st-andrews.ac.uk
#
# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK
#
# The MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

library("rjson")
library("vegan")

json_file <- "../results/cazy_families/famcounts.json"
json_data <- fromJSON(file=json_file)

fam.df <- read.csv("../results/cazy_families/cazy_fam_freqs.csv")
all.fams <- colnames(fam.df)[5:length(colnames(fam.df))]
row.names(fam.df) <- fam.df$Genome
fam.only.df <- fam.df[all.fams]

get.colours <- function(genera){
  colours <- c()
  for(genus in fam.df$Genus){
    if (genus == 'Pectobacterium'){
      colours <- append(colours, "#940113")
    }
    else if (genus == 'Dickeya'){
      colours <- append(colours, "#0a7a6d")
    }
    else if (genus == 'Lonsdalea'){
      colours <- append(colours, "#d4af37")
    }
    else if (genus == 'Brenneria'){
      colours <- append(colours, "#7844b8")
    }
    else if (genus == 'Musicola'){
      colours <- append(colours, "#13851e")
    }
    else if (genus == 'Samsonia'){
      colours <- append(colours, "#c22176")
    }
    else if (genus == 'Affinibrenneria'){
      colours <- append(colours, "#3888e0")
    }
    else { # Acerihabitans
      colours <- append(colours, "#e07e38")
    }
  }
  
  return(colours)
}

famAbund <- rowSums(fam.only.df)  #gives the number of individuals found in each plot
# view observations per plot 

raremin <- min(rowSums(fam.only.df))  #rarefaction uses the smallest number of observations per sample to extrapolate the expected number if all other samples only had that number of observations # view smallest # of obs (site 17)

sRare <- rarefy(fam.only.df, raremin) # now use function rarefy
#gives an "expected"rarefied" number of species (not obs) if only 15 individuals were present


pecto.colours <- get.colours(fam.df$Genus)
pdf("../results/cazy_families/pectobacteriaceae_rareCurve_genus_large.pdf", width = 20, height = 30)
rarecurve(
  fam.only.df,
  col = pecto.colours,
  sample=raremin,
  xlab='Number of genomes',
  ylab='Number of CAZy families',
)
legend(
  x = 0, y = 70,
  legend = c("Pectobacterium", "Dickeya", 'Lonsdalea', 'Brenneria', 'Musicola', 'Samsonia', 'Affinibrenneria', 'Acerhabitans'),
  col = c(
    "#940113",
    "#0a7a6d",
    "#d4af37",
    "#7844b8",
    "#13851e",
    "#c22176",
    "#3888e0",
    "#e07e38"
  ),
  lwd = 2,
  cex = 0.8,
  title = "Genus"
)
dev.off()

pdf("../results/cazy_families/pectobacteriaceae_rareCurve_genus.pdf", width = 8, height = 12)
rarecurve(
  fam.only.df,
  col = pecto.colours,
  sample=raremin,
  xlab='Number of genomes',
  ylab='Number of CAZy families',
  label=FALSE,
)
legend(
  x = 0, y = 70,
  legend = c("Pectobacterium", "Dickeya", 'Lonsdalea', 'Brenneria', 'Musicola', 'Samsonia', 'Affinibrenneria', 'Acerhabitans'),
  col = c(
    "#940113",
    "#0a7a6d",
    "#d4af37",
    "#7844b8",
    "#13851e",
    "#c22176",
    "#3888e0",
    "#e07e38"
  ),
  lwd = 2,
  cex = 0.8,
  title = "Genus"
)
dev.off()



## _Pectobacterium_

p.fam.df <- subset(fam.df, Genus=='Pectobacterium')
p.fam.only.df <- p.fam.df[all.fams]
famAbund <- rowSums(p.fam.only.df)
raremin <- min(rowSums(p.fam.only.df))
sRare <- rarefy(p.fam.only.df, raremin)
pdf("../results/cazy_families/pectobacterium_rareCurve.pdf", width = 12, height = 16)
rarecurve(
  p.fam.only.df,
  col = "#940113",
  sample=raremin,
  xlab='Number of genomes',
  ylab='Number of CAZy families',
)
dev.off()

## _Dickeya_

d.fam.df <- subset(fam.df, Genus=='Dickeya')
d.fam.only.df <- d.fam.df[all.fams]
famAbund <- rowSums(d.fam.only.df)
raremin <- min(rowSums(d.fam.only.df))
sRare <- rarefy(d.fam.only.df, raremin)
pdf("../results/cazy_families/dickeya_rareCurve.pdf", width = 12, height = 16)
rarecurve(
  d.fam.only.df,
  col = "#0a7a6d",
  sample=raremin,
  xlab='Number of genomes',
  ylab='Number of CAZy families',
)
dev.off()


## _Brenneria_

d.fam.df <- subset(fam.df, Genus=='Brenneria')
d.fam.only.df <- d.fam.df[all.fams]
famAbund <- rowSums(d.fam.only.df)
raremin <- min(rowSums(d.fam.only.df))
sRare <- rarefy(d.fam.only.df, raremin)
pdf("../results/cazy_families/brenneria_rareCurve.pdf", width = 12, height = 16)
rarecurve(
  d.fam.only.df,
  col = "#7844b8",
  sample=raremin,
  xlab='Number of genomes',
  ylab='Number of CAZy families',
)
dev.off()


## _Lonsdalea_

d.fam.df <- subset(fam.df, Genus=='Lonsdalea')
d.fam.only.df <- d.fam.df[all.fams]
famAbund <- rowSums(d.fam.only.df)
raremin <- min(rowSums(d.fam.only.df))
sRare <- rarefy(d.fam.only.df, raremin)
pdf("../results/cazy_families/lonsdalea_rareCurve.pdf", width = 12, height = 16)
rarecurve(
  d.fam.only.df,
  col = "#d4af37",
  sample=raremin,
  xlab='Number of genomes',
  ylab='Number of CAZy families',
)
dev.off()
