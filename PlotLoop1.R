#Start for cycle
i = "Adaptors Nab2"
number = 1
for(i in ProtNames){
  
  Nab2Frame = data.frame(AB$`systematic name`,AB$`Half-Life`, ProtData[i])
  names(Nab2Frame) = c("Name","HalfLife","Occ")
  
  
  VariableHigher = quantile(Nab2Frame$Occ, 0.80)
  VariableLower = quantile(Nab2Frame$Occ, 0.30)
  data=data.frame(name=c(paste("Bind, >",names(VariableHigher)),paste("Not Bind, <", names(VariableLower))) ,  value=c(sum(Nab2Frame$HalfLife[Nab2Frame$Occ >= VariableHigher])/length(Nab2Frame$HalfLife[Nab2Frame$Occ >= VariableHigher]),sum(Nab2Frame$HalfLife[Nab2Frame$Occ <= VariableLower])/length(Nab2Frame$HalfLife[Nab2Frame$Occ <= VariableLower])))
  
  
  
  boxplot(list(Bindet = Nab2Frame$HalfLife[Nab2Frame$Occ >= VariableHigher],BindetNicht = Nab2Frame$HalfLife[Nab2Frame$Occ <= VariableLower], GanzeDataset = Nab2Frame$HalfLife), ylim = c(0, 80), col = c("green", "red", "blue"), ylab = "Halbwertszeit in min", main = paste("Vergleich Summe von Occupancy Events\n", i ,"mehr als", names(VariableHigher), "und kleiner als", names(VariableLower)), sub =  paste("Unterschied in Mittelwert zwischen Bindet und BindetNicht =",abs(round(100 - (sum(data[1,2]/sum(data[2,2])*100)),2)),"%","\n Unterschied in Minuten:",abs(round(mean(Nab2Frame$HalfLife[Nab2Frame$Occ <= VariableLower])-mean(Nab2Frame$HalfLife[Nab2Frame$Occ >= VariableHigher]),1)),"\n Wilcox-Test: W Wert =",abs(round(wilcox.test(Nab2Frame$HalfLife[Nab2Frame$Occ >= VariableHigher],Nab2Frame$HalfLife[Nab2Frame$Occ <= VariableLower])[[1]],2)),"; p Wert:",abs(signif(wilcox.test(Nab2Frame$HalfLife[Nab2Frame$Occ >= VariableHigher],Nab2Frame$HalfLife[Nab2Frame$Occ <= VariableLower])[[3]]))))
  means = list(Bind = round(mean(Nab2Frame$HalfLife[Nab2Frame$Occ >= VariableHigher]),1),Not_Bind = round(mean(Nab2Frame$HalfLife[Nab2Frame$Occ <= VariableLower]),1), All_data = round(mean(Nab2Frame$HalfLife),1))
  means = c(means$Bind,means$Not_Bind,means$All_data)
  points(1:3, means, col = "black")
  text(1:3, means + 4,labels = means)
  
  PlotName = paste(number,i,"BoxPlot",".jpg")
  dev.copy(jpeg,filename=PlotName)
  dev.off ()
  
  
  plot(density(Nab2Frame$Occ), main = paste("Ganze Dataset fur Occ."), xlab = "Summe von Occupancy Events", xlim = c(0,0.02))
  abline(v = VariableLower, col = "green")
  abline(v = VariableHigher, col = "red")
  legend(0.012,max(density(Nab2Frame$Occ)[[2]]), legend = c(paste(names(VariableLower), "Quantile"),paste(names(VariableHigher),"Quantile")), col = c("green","red"), lty=1:1, cex=0.8)
  PlotName = paste(number,i,"Density",".jpg")
  dev.copy(jpeg,filename=PlotName)
  dev.off ()
  print("fine")
  
  
  if(sum(Nab2Frame$Occ[Nab2Frame$Occ <= VariableLower]) != 0 && sum(Nab2Frame$HalfLife[Nab2Frame$Occ <= VariableLower]) != 0){
    print(plot_replicate_smoothscatter_jpeg(Nab2Frame$Occ[Nab2Frame$Occ >= VariableHigher],Nab2Frame$HalfLife[Nab2Frame$Occ >= VariableHigher],paste(i,"Bindet")),split = c(1,1,2,1), more = TRUE)
    
    
    
    'PlotName = paste(number,i,"Scatter_Bind",".jpg")
    dev.copy(jpeg,filename=PlotName)
    dev.off ()
    graphics.off()
    
    print("fine2")
    print(i)'
    
    
    print(plot_replicate_smoothscatter_jpeg(Nab2Frame$Occ[Nab2Frame$Occ <= VariableLower],Nab2Frame$HalfLife[Nab2Frame$Occ <= VariableLower],paste(i,"BindetNicht")),split = c(2,1,2,1))
    
    print("fine3")
    
    PlotName = paste(number,"Scatter_Not_Bind",".jpg")
    dev.copy(jpeg,filename=PlotName)
    dev.off ()
    graphics.off()
  }else{
    print(plot_replicate_smoothscatter_jpeg(Nab2Frame$Occ[Nab2Frame$Occ >= VariableHigher],Nab2Frame$HalfLife[Nab2Frame$Occ >= VariableHigher],paste(i,"Bindet")))
    
    
    
    PlotName = paste(number,i,"Scatter_Bind",".jpg")
    dev.copy(jpeg,filename=PlotName)
    dev.off ()
    graphics.off()
  }
  
  number = number + 1
  
}

ProtNames = c("Cbc2","BBPU2AF65_Msl5","BBPU2AF65_Mud2","U1 snRNP_Luc7","U1 snRNP_Mud1","U1 snRNP_Nam8","U1 snRNP_Snp1","U2 snRNP_Ist3","CFIA_Rna15","CPF_Cft2","CPF_Mpe1","CPF_Yth1","Poly(A)_Pab1","Poly(U)_Pub1","THO_Hpr1","THO_Tho2","TREX_Gbp2","TREX_Hrb1","TREX_Mex67","TREX_Sub2","TREX_Yra1","Adaptors_Nab2","Adaptors_Npl3")



plot_replicate_smoothscatter_jpeg <- function(rep1, rep2, rep3){
  xx = log2(rep2)
  yy = log2(rep1)
  sel = !is.infinite(xx) & !is.infinite(yy)
  
  corcoef = cor.test(xx[sel], yy[sel], method="pearson")
  
  df <- data.frame(x1=xx[sel], x2=yy[sel])
  x <- densCols(xx[sel], yy[sel], colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(x)[1,] + 1L
  cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
  df$col <- cols[df$dens]
  df = df[order(df$dens),]
  pl = xyplot(x1 ~ x2, data=df, scales=list(tck=c(0.1,0)),
              pch=20, cex=0.1,
              panel = function(x,y,...){
                panel.xyplot(x,y,..., col=df$col)
                panel.lmline(x,y,..., col="red", lty=3, lwd=0.7)
              },
              xlab=paste("Logarithmierte Occupancy Events\n Person Korrelationskoeffizient:", round(as.numeric(corcoef[4]),2)), ylab="Logarithmierte Halbwertszeit", aspect=1, main = rep3,grid = TRUE)
  return(pl)
}

n = 1
for(a in names(ProtData)){
  for(i in c(1:4947)){
    ProtData[n,a] = ProtData[n,a]/(end(gr3)[n] - start(gr3)[n])
    n = n + 1
  }
  n = 1
  
}

for(a in c(1:5)){
  print("a")
}



getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

TestList = Nab2Frame$Occ[Nab2Frame$Occ < VariableLower]
n = 1
for(i in TestList2){
  TestList2[n] = VariableLower + (VariableLower - TestList2[n])
  n = n+1
}

lines(density(Listtt), col = "blue")

n = 1
for(i in Listtt){
  Listtt[n] = abs((Listtt[n]-mean(Listtt))/sd(Listtt))
  n = n+1
}