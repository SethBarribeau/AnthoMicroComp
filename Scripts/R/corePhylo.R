#Core Phylogroups
#11 November 2022

library(edgeR)
library(ggplot)
library(reshape)

#Sample metadata
met <- read.table("input/SRA/MetaData_Edit_Oct22.tsv",
                  header = T, sep ="\t", quote = "")

#Microbial metadata
genkey <- read.table("output/Keys/GenusFamilySpeciesKey.tsv",
                     sep = "\t", header = T)
#Taxonomy
tax <- read.table("input/Phylo_Misc/rankedlin_Edit_Nov22.tsv",
                  header = T, sep = "\t", quote = "")

#core phylogroup lactobacillus firm -5
firm5 <- c("Lactobacillus apis", "Lactobacillus melliventris", "Lactobacillus kimbladii",
           "Lactobacillus kullabergensis", "Lactobacillus panisapium", 
           "Lactobacillus bombicola", "Lactobacillus helsingborgensis")

#load in species level count table
spec.cnt <- read.table("output/Counts/RawReads/All_rawCounts_Nov22.tsv")

#check the length of the IDs to make a new one
length(unique(genkey$GenHitID))
#apply this to the firm5 group
genkey$GenHitID[genkey$MicroHit %in% firm5] <- "GRH761"

#change phylogeny
genkey$PhyloHit[genkey$MicroHit %in% firm5] <- "Lactobacillus: Firm-5"

#add to taxonomy
tax[tax$genus == "Lactobacillus",]
n <- nrow(tax)+1
tax[n,] <- c("Lactobacillus: Firm5", "", "Lactobacillus: Firm-5", "Lactobacillaceae",
             "Lactobacillales", "Bacilli", "Firmicutes", NA, "Bacteria")

#Add the microhit ids to the genkey so I can recount the single species count table
microkey <- read.table("output/Keys/MicrobialSpeciesKey.tsv",
                       sep = "\t", header = T)
for (i in 1:nrow(genkey)){
  genkey$MicroHitID[i] <- microkey$HitID[microkey$MicroHit == genkey$MicroHit[i]]
}

#Now make a new, phylo collapsed count table ....
#compile samples
samples <- names(spec.cnt)
#compile microbial genus hit ids
phyloz <- unique(genkey$GenHitID)
#make a matrix of sufficient size
cnt <- matrix(ncol = length(samples), nrow = length(phyloz))
#run through each of the samples
for (i in 1:length(samples)){
  #print to screen so progress can be checked
  print(paste(i, "/317: ", samples[i], sep = ""))
  #make a vector for the sample read to populate with counts per microbial genus group
  x <- vector(length = length(phyloz))
  #iterate through the tax hits, pulling out the associated species as we go
  for (j in 1:length(x)){
    #extract the microbial IDs within this microbial genus grouping
    y <- unique(genkey$MicroHitID[genkey$GenHitID == phyloz[j]])
    #pull out the rows of the count table that match these microbial ids
    t <- spec.cnt[rownames(spec.cnt) %in% y, i]
    #combine multiple species from a phylogenetic classification into one
    if (length(t) > 1){
      t <- sum(t)
    }
    #populate sample's vector with counts
    x[j] <- t
    #add sample's counts to count matrix
    cnt[,i] <- x
  }
}
#covert to a dataframe
cnt <- as.data.frame(cnt)
#add row and colomn names
rownames(cnt) <- phyloz
names(cnt) <- samples
#write up
write.table(cnt, "output/Counts/RawReads/All_rawCounts_CorePhylo_PhyloCollapsed.tsv",
            sep = "\t", col.names = T, row.names = T, quote = F)

#remove DNA samples
dna <- met$Sample.ID[met$NucleotideType =="DNA"]
cnt <- cnt[,! names(cnt) %in% dna]

#keep only prokaryotic hits
pro <- genkey$GenHitID[genkey$Class == "bacteria" | genkey$Class == "archaea"]
cnt <- cnt[rownames(cnt) %in% pro,]

#remove any samples with < 1000 reads
keep <- colSums(cnt) > 1000
cnt <- cnt[,keep]

#make a TMM table using edgeR
#I cannot have more more entries in the met object that are in the analysis
met2 <- met[met$Sample.ID %in% names(cnt),]
#set analysis parameters (by sociality)
group <- factor(met2$Sociality)
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
#Creates a DGEList object from a table of counts (rows=features, columns=samples) with group indicator for each column
work <- DGEList(counts = cnt, group = group)
#I would normally try to remove "low expression" rows here but it breaks all the later steps so I removed that step and skipped to 
#calculate normalization factors to align column of a count matrix
#each sample will have a normalisation factor to try and combat effective library size effects
#(no idea if this is going to effect the table / make it useless)
work2 <- calcNormFactors(work)
#estimate common, trended and tagwise negative binomial dispersions by weight likelihood empirical Bayes
work2 <- estimateDisp(work2, design)
#fit a quasi-likelihood negative binomial generalised log-linear model to count data
fit <- glmQLFit(work2, design)

#extract tmm
cnt.tmm <- fit$fitted.values

#write up
write.table(cnt.tmm, "output/Counts/EdgeR/AllMicrobials_Filtered_Nov22.tsv",
            sep = "\t", col.names = T, row.names = T, quote = F)

#extract the core microbiota phylotypes
beecore <- c("Lactobacillus: Firm5", "Apilactobacillus", "Bombilactobacillus",
             "Gilliamella", "Snodgrassella", "Bifidobacterium", "Frischella",
             "Bartonella", "Bombiscardovia", "Schmidhempelia")
bcoreID <- unique(genkey$GenHitID[genkey$PhyloHit %in% beecore])
bc.tmm <- cnt.tmm[rownames(cnt.tmm) %in% bcoreID,]

#prepare to plot
tmp <- melt(bc.tmm)
names(tmp) <- c("GenusHitRead", "Sample", "TMM")
for (i in 1:nrow(tmp)){
  tmp$MicroPhylo[i] <- unique(genkey$PhyloHit[genkey$GenHitID == tmp$GenusHitRead[i]])
}
bc.plot <- merge(met2, tmp, by.x="Sample.ID", by.y = "Sample", all.y = T)

ggplot(data = bc.plot[bc.plot$TMM < 5000,], aes(x = MicroPhylo, y = TMM, fill = Sociality)) +
  geom_boxplot() +
  coord_flip()

bc.plot[bc.plot$TMM > 2e+05,]
