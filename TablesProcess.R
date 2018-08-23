library(readxl)
library(xlsx)
AB = read_excel("mmc1.xlsx")
AB$Description = NULL
colnames(AB) = c("chrom" ,"systematic name" ,"common name", "Simple Type", "Half-Life")
hist(AB$`Half-Life`, col = "blue", freq = FALSE, breaks = 50, xlim = c(0,150), main = "Histogram", xlab = "Half-Life in min")
abline(v = quantile(AB$`Half-Life`, probs = c(0.1,0.9), names = FALSE), col = "red")
abline(v = quantile(AB$`Half-Life`, probs = 0.9, names = FALSE), col = "red")
d = density(AB$`Half-Life`)
plot(d, xlim = c(0,150), main = "Density Plot", sub = paste( "Min = 10%:",,"; Max = 90%"))
abline(v = quantile(AB$`Half-Life`, probs = 0.1, names = FALSE), col = "red")
abline(v = quantile(AB$`Half-Life`, probs = 0.9, names = FALSE), col = "red")
polygon(d, col = "blue")

table_sort = round(AB$`Half-Life`)
table_sort = table(table_sort)
nn = as.data.frame(table_sort)

summary(AB$`Half-Life`)
quantile(AB$`Half-Life`)
quantile(AB$`Half-Life`, probs = seq(0,1, 0.0))

probs = 0.1
aa=  quantile(AB$`Half-Life`, probs = probs, names = FALSE)
if(probs == 0.9){
  var = AB$`Half-Life` > aa
  }else{var = AB$`Half-Life` < aa}

n = 0
my_vektor = vector()
my_vektor2 = vector()
for(i in var){n = n+1
  if(i == TRUE){ 
    my_vektor = append(my_vektor, AB$`common name`[n])
    my_vektor2 = append(my_vektor2, AB$`Half-Life`[n])
  }
}

my_vektor = AB$`common name`[var]
my_vektor = AB$`common name`[which(var)]


data = data.frame(my_vektor, my_vektor2)
colnames(data) = c("Name","Half-Life")
write.xlsx(data, "data.xlsx")











my_vektor = AB$`common name`[AB$`Half-Life` < 10]
my_vektor2 = AB$`common name`[AB$`Half-Life` > 60]
write.xlsx(list(peter=my_vektor, my_vektor2), "data.xlsx")
bka = data.frame(my_vektor, my_vektor2)


gr <- GRanges(
  seqnames = "chr6",
  ranges = IRanges(start = 53260, end = 54696, names = "Actin"),
  strand = Rle(strand(c("-"))))


fileName = "S288C_YFL039C_ACT1_genomic.fsa"
my_file = readChar(fileName, file.info(fileName)$size)
my_file = gsub("\n","",ss)






























plot(density(AB$`Half-Life`), col  = "red", xlim = c(0,120), main = "Density samples", xlab = "Half-life in min", sub = "Full size = 4947")
lines(density(sample(AB$`Half-Life`,  500)), col = "green")
lines(density(sample(AB$`Half-Life`, 50)), col = "blue")
lines(density(sample(AB$`Half-Life`,10)), col = "yellow")
legend(62,0.03, legend = c("Full", "Sample 500", "Sample 50", "Sample 10"), col = c("red","green", "blue", "yellow"), lty=1:1, cex=0.8)


Sample_500 = sample(AB$`Half-Life`, 500)
Sample_4000 = sample(AB$`Half-Life`, 4000)
Sample_50 = sample(AB$`Half-Life`, 50)
Sample_10 = sample(AB$`Half-Life`, 10)
list_box = list(S4000 = Sample_4000,S500 = Sample_500,S50 = Sample_50,S10 = Sample_10)
boxplot(list_box,  main = "Boxplot samples", xlab = "Sample size", ylab = "Half-Life in min",ylim = c(0, 100), col = c("red", "green", "blue", "yellow")) 


plot(density(AB$`Half-Life`[AB$`Half-Life` > 20]))

plot(density(sample(AB$`Half-Life`[AB$`Half-Life` > 20], 3000)), col = "green", main = "Density, more then 20 min, sample 3000+500", xlab = "Time in min")
lines(density(sample(AB$`Half-Life`[AB$`Half-Life` > 20], 500)), col = "red", lty = 2)
legend(200,0.03, legend = c("3000", "500"), col = c("green", "red"), lty=1:2, cex=0.8)

plot(density(sample(AB$`Half-Life`[AB$`Half-Life` > 10], 3000)), col = "green")
lines(density(sample(AB$`Half-Life`[AB$`Half-Life` > 10], 500)), col = "red")

plot(density(AB$`Half-Life`),xlim = c(1,120), main = "Density from Half-Life data", xlab = "Half-life in min")
abline(v = 10, col = "blue")
abline(v = 60, col = "red")
legend(62,0.03, legend = c("Least stable - t1/2<10 min", "Most stable - t1/2>60min"), col = c("blue", "red"), lty=1:1, cex=0.8)

liset = AB$`Half-Life`[AB$`Half-Life` > 10 & AB$`Half-Life` < 60]
plot(density(liset), col  = "red", xlim = c(0,70), main = "Density samples between 10 and 60 min", xlab = "Half-life in min", sub = "Full size = 4392")
lines(density(sample(liset,  500)), col = "green")
lines(density(sample(liset, 50)), col = "blue")
lines(density(sample(liset,10)), col = "yellow")
legend(50,0.035, legend = c("Full", "Sample 500", "Sample 50", "Sample 10"), col = c("red","green", "blue", "yellow"), lty=1:1, cex=0.8)
