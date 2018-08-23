processFile = function(filepath) {
  counter = 0
  first_line = vector()
  file_seq = vector()
  my_file = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(my_file, n = 1, warn = FALSE)
    if ( length(line) == 0 ) {
      break
    }
    if(counter == 1){
      file_seq = c(file_seq, line)
    }else{
      first_line = line
    }
    counter = 1
    
  }
  
  close(my_file)
  

  first_line = gsub(" |, |:|\\..", ";", first_line)
  file_seq = gsub(" |,|c|\"|\n|\\(|\\)", "",file_seq)
  first_line = gsub("SGDID;|>", "", first_line)
  ComName = strsplit(first_line, ";")[[1]][1]
  Name = strsplit(first_line, ";")[[1]][2]
  SGDID = strsplit(first_line, ";")[[1]][3]
  Chromosome = strsplit(first_line, ";")[[1]][4]
  Start = as.numeric(strsplit(first_line, ";")[[1]][5])
  End = as.numeric(strsplit(first_line, ";")[[1]][6])
  
  
  
  return(list(CommonName = ComName,Name = Name,SGDID = SGDID,Chromosome = Chromosome,Start=  Start,End = End,  Sequence = file_seq))
}


x = processFile("Actin.fsa")

library(GenomicRanges)

gr <- GRanges(
  seqnames = x$Name,
  ranges = IRanges(start = x$Start, end = x$End, names = x$CommonName),
  strand = "-",
  SGDID = x$SGDID,
  chromosome = x$Chromosome)



