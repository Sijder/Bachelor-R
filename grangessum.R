library(GenomicRanges)
my_table = read.table("A1.txt", fill = TRUE)
my_table2 = my_table

aa = read_excel("mart_export (1).xls")

gr2 = GRanges(
  seqnames = aa$ID,
  ranges = IRanges(start = aa$start, end = aa$end, names = aa$ID),
  strand = aa$strand
)

sum(my_table2$V7[my_table2$V3 >= start(gr) & my_table2$V4 <= end(gr) & my_table2$V12 == as.character(strand(gr))])

findOverlaps
split

Nab2_Table = read.table("GSM1442550_Nab2-121211-4-B4.txt", header = TRUE)
gsub("chr","",Nab2_Table$chromosome)
Genes_table = read.csv("BFE.txt")
gr3 = GRanges(
  seqnames = Rle(Genes_table$Chromosome.scaffold.name),
  ranges = IRanges(start = Genes_table$Gene.start..bp., end = Genes_table$Gene.end..bp., names = Genes_table$Gene.stable.ID),
  strand = Rle(Genes_table$Strand),
  Nab2OccupancySumm = NULL,
  HalfLife = NULL
)
gr3$HalfLife = AB$`Half-Life`
# Normal for loop for occupancy count

n = 0
for(unit in names(gr3)){
  n = n + 1
  gr3$Nab2Occupancy[n] = (sum(Nab2_Table$occupancy[Nab2_Table$position >= start(gr3)[n] & Nab2_Table$position <= end(gr3)[n] & Nab2_Table$strand == as.vector(strand(gr3)[n]) & Nab2_Table$chromosome == as.vector(seqnames(gr3)[n])]))
}

#findoverlaps function for occupancy count

gr4 = GRanges(
  seqnames = Rle(Nab2_Table$chromosome),
  ranges = IRanges(start = Nab2_Table$position, end = Nab2_Table$position),
  strand = Rle(Nab2_Table$strand),
  Occupancy = Nab2_Table$occupancy
)


Overlaps = findOverlaps(gr4, gr3, type = "within")
n = 0
for(unit in as.vector(seqnames(gr3))){
  n = n+1
  gr3$Nab2OccupancySumm[n] = sum(gr4$Occupancy[queryHits(Overlaps[subjectHits(Overlaps) == n])])
}


plot(density(GenesDataFrame2$`GSM1442550_Nab2-121211-4-B4.txt`), main = "Density for Nab2 Occupancy", xlab = "Nab2 Occupancy", col = "green")

plot(density(AB$`Half-Life`), col  = "red", xlim = c(0,120),ylim = c(0,0.16), main = "Density samples", xlab = "Half-life in min", sub = "Occupancy number")
lines(density(GenesDataFrame2$`GSM1442550_Nab2-121211-4-B4.txt`), col = "green")
legend(62,0.1, legend = c("Half-life", "Occupancy"), col = c("red","green"), lty=1:1, cex=0.8)
abline(v = 10, col = "blue")
abline(v = 60, col = "red")



#Trying to find correlation

CompareFrame = data.frame(GenesDataFrame2$`GSM1442550_Nab2-121211-4-B4.txt`, gr3$HalfLife, GenesDataFrame2$`GSM1442550_Nab2-121211-4-B4.txt` / (end(gr3)-start(gr3)))
colnames(CompareFrame) = c("Occ","HL","OccDLength")
CompareFrame$Occ[CompareFrame$Occ < quantile(CompareFrame$Occ, 0.90)] = 0
CompareFrame$Occ[CompareFrame$Occ > quantile(CompareFrame$Occ, 0.90)] = 1
cor.test(CompareFrame$HL, CompareFrame$Occ)

# Maybe if i divide sum to length
CompareFrame = data.frame(gr3$Nab2OccupancySumm, gr3$HalfLife, gr3$Nab2OccupancySumm / (end(gr3)-start(gr3)))
colnames(CompareFrame) = c("Occ","HL","OccDLength")
CompareFrame$OccDLength[CompareFrame$OccDLength < quantile(CompareFrame$OccDLength, 0.90)] = 0
CompareFrame$OccDLength[CompareFrame$OccDLength > quantile(CompareFrame$OccDLength, 0.90)] = 1
cor.test(CompareFrame$HL, CompareFrame$OccDLength)


p = ggscatter(Nab2Compare, x = "Occ90", y = "HL90", 
          add = "reg.line", conf.intw = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "OccTrue", ylab = "Numbers")
