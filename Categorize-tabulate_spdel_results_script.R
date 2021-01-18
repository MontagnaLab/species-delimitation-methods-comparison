##################################################################################
##### Factors affecting the efficiency of molecular species delimitation in  #####
##### a species-rich insect family ###############################################

## Giulia Magoga, Diego Fontaneto & Matteo Montagna

##################################################################################

## R code to categorize species delimitation results as match, merge, mixture, split and tabulate them using R 3.6.3 (in this example ASAP partition is used as input).
## Last update: Milano, 18 January 2021

##################################################################################

# Loading R packages -------------------------------------------------------------

library('ape')
library('spider')

# Working directory --------------------------------------------------------------

# check the folder

getwd()

# change it if needed
setwd("") # <- Change me

# Loading the files --------------------------------------------------------------

##### files for the analyses #############################
# best_partition.txt is the file with the output of ASAP #
# seq.fas is the file with the alignment in fasta format #
##########################################################

# Read the output of ASAP into a list --------------------------------------------

# each line is an element of the list
x <- scan("best_partition.txt", what="", sep="\n") 

# spaces between names to separate elements of the list
y <- strsplit(x, "[[:space:]]+") 

# Read the dna alignment ---------------------------------------------------------

# read the alignment
dna <- ape::read.dna("seq.fas", format="fasta")

# Categorize results as match, merge, mixture or split -------------------------- 

# obtain species names
strsp <- strsplit(dimnames(dna)[[1]], split="_")
spp <- sapply(strsp, function(x) paste(x[1],x[2],sep="_"))
un <- unique(spp)

# loop to match names and their position in the list
ch <- vector("list", length(un))
for (i in seq(length(un))) {
			ch[[i]]=grep(un[i], 
			y, 
			ignore.case=FALSE, 
			perl=FALSE, 
			value=FALSE, 
			fixed=FALSE, 
			useBytes=FALSE, 
			invert=FALSE)
			} 
ch

# loop to extract from 'ch' elements of length 1  
names(ch) <- un  
li <- list()
for(i in 1:length(ch)) {
			if(lengths(ch[i], use.names=TRUE) == 1) {
					li=c(li,ch[[i]])
					}
			}
li
 
# MATCH
# loop to obtain match, i.e. elements appearing only once in the list 'ch'
match = list()
times = 0
posCh = 0
posFound = 0
for(elementLi in li) {
 	for(line in ch) {
 		posCh = posCh + 1
 		for(elementCh in line) {
 			if (elementLi == elementCh) { 
 				times = times + 1
 				posFound = posCh
 			}
 		}
 	}
 	if (times == 1) {
 		match = c(match, names(ch[posFound]))
 	}
 	posCh = 0
 	posFound = 0
 	times = 0
	}   
match


# MERGE
# obtain numbers that appear more than once in the list 'li'
te <- unlist(li)
o <- table(te)
thr <- 1
ko <- names(o)[o > thr] 

# loop to obtain the list 'limo'composed by elements with length > 1 in 'ch'
limo <- list()
for(i in 1:length(ch)){
			if(lengths(ch[i], use.names=TRUE) > 1) {
				limo[i] = ch[i]
				}
			}
limo

# find elements of 'ko' in 'limo'
lio <- lapply(limo, intersect, ko) 

# if any element of 'ko' is present in 'limo' write merge
if (all(lapply(lio, function(x) identical(x, character(0))))) {
	merge = lapply(ch, intersect, ko)
	} 

# if any element of 'ko' are present in 'limo' do not include it in merge (if an error occurs ignore it, means merges are already saved in previous vector)

lio <- lio[lapply(lio, length) > 0] 
lio <- unlist(lio)
if (lapply(lio, length) > 0) {mergepre = ko[!ko %in% lio]} 
merge <- lapply(ch, intersect, mergepre) 

# substitute numbers in 'merge' with names
merge <- names(merge[lapply(merge, length) > 0]) 

# MIXTURE
# find elements of 'li' in 'limo'
comli <- lapply(limo, intersect, li) 
uncomli <- unlist(comli)

# find elements of 'uncomli' in 'ch' and save in mixture
mixture <- lapply(ch, intersect, uncomli) 

# assign to 'mixture' the species names
mixture <- names(mixture[lapply(mixture, length) > 0])

# obtain numbers that appear more than once in the list 'tunlimo' 
unlimo <- unlist(limo)
tunlimo <- table(unlimo)
mlimo <- names(tunlimo)[tunlimo > thr] 

# find elements of 'mlimo' in 'ch'
mixlo <- lapply(ch, intersect, mlimo) 

# substitute numbers in 'mixlo' with names
mixlo <- names(mixlo[lapply(mixlo, length) > 0])

# add 'mixlo' to mixture and delate multiple elements
mixture <- append(mixture, mixlo) 
mixture <- unique(mixture)

# SPLIT
# find names of elements that appear more than once in the list 'tunlimo'
onelimo <- names(tunlimo)[tunlimo == thr]

# find elements of 'onelimo' in 'li' and remove them from 'onelimo'
mli <- lapply(li, intersect, onelimo)
unmli <- unlist(mli)
onelimor <- onelimo[!onelimo %in% unmli]

# find elements of 'onelimor' in 'ch'
msplit <- lapply(ch, intersect, onelimor) 

# substitute numbers in 'msplit' with names
splitpre <- names(msplit[lapply(msplit, length) > 0])

# find elements of 'unmli' in 'ch', find their name and remove them from 'splitpre' by name
cut <- lapply(ch, intersect, unmli)
cut <- names(cut[lapply(cut, length) > 0]) 
split <- splitpre[!splitpre %in% cut] 

# find elements of 'mlimo' in 'ch', find their name and remove them from 'split' by name
cut2 <- lapply(ch, intersect, mlimo)
cut2 <- names(cut2[lapply(cut2, length) > 0])
split <- split[!split %in% cut2] 

# Bind match, merge, mixture and split in a data frame --------------------------- 

# find length of 'match', 'merge', 'mixture' and 'split' as vectors
match <- unlist(match)
merge <- unlist(merge)
max.len <- max(length(match), length(merge), length(mixture), length(split))

# concatenate NA to 'match', 'merge', 'mixture' and 'split' according to length 'max.len'
match <- c(match, rep(NA, max.len - length(match)))
merge <- c(merge, rep(NA, max.len - length(merge)))
mixture <- c(mixture, rep(NA, max.len - length(mixture)))
split <- c(split, rep(NA, max.len - length(split)))

# bind 'match', 'merge', 'mixture', 'split' in a data frame
table <- data.frame(match, merge, mixture, split)

# Save table in csv --------------------------------------------------------------

write.csv(table, "table_results.csv")
