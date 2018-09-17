#Functions


  
  ####Function to count alleles and populations parameters###
tassel5_to_params=function(x="hap matrix", y="columns to skip", z="population number"){
      geno=x
      #recount allele A and B and het
      alleleA=rowSums(geno[,(y+1):ncol(geno)]!=substring(geno$alleles, 3, 3) & geno[,(y+1):ncol(geno)]!="N") #only counts what is not allele B and missing.  i.e. counts allele A and various calls for heterozygous
      alleleB=rowSums(geno[,(y+1):ncol(geno)]!=substring(geno$alleles, 1, 1) & geno[,(y+1):ncol(geno)]!="N")
      het=rowSums(geno[,(y+1):ncol(geno)] == "M") + rowSums( geno[,(y+1):ncol(geno)] ==   "R") + rowSums(geno[,(y+1):ncol(geno)] ==  "W") + rowSums(geno[,(y+1):ncol(geno)] ==  "K") + rowSums(geno[,(y+1):ncol(geno)] ==  "S") + rowSums(geno[,(y+1):ncol(geno)] ==  "Y")
      present=1-(rowSums(geno[,(y+1):ncol(geno)]=="N")/z)
      MAF=apply(cbind(((alleleA-het)*2+het), (alleleB-het)*2+het), 1, min)/apply(
    cbind(((alleleA-het)*2+het), ((alleleB-het)*2+het)), 1, sum) 
      percentHet=het/apply(cbind(alleleA-het, alleleB-het, het), 1, sum)
      return(cbind.data.frame(geno[,1:y], "alleleA"=alleleA, "alleleB"=alleleB, "het"=het, "present"= present, "MAF"=MAF, "percentHET"=percentHet, geno[,(y+1):ncol(geno)]))
}

##function to convert hap to 0 and 1
hap_to_G=function(x="hap matrix", y="number of columns of information"){
  ##From Prasana, pulls out first allele for a and second for b
  a = substring(x$alleles,1,1)
  #Checks the frequency of the alleles if the second allele is more frequent it is substitued
  a[x$alleleA<x$alleleB] = substring(x$alleles,3,3)[x$alleleA<x$alleleB]
  #Same thing with the second allele
  b = substring(x$alleles,3,3)
  b[x$alleleA<x$alleleB] = substring(x$alleles,1,1)[x$alleleA<x$alleleB]
  #Checks to make sure all alleles are one or the other
  #print(paste("If 0 all alleles are accounted for: ", sum(a == b), sep=""))
  
  ## Turn into letter matrix for mapping
  #makes a copy of the hap matrix
  hap01 = x
  #sets all allele values to NA
  hap01[,(y+1):ncol(hap01)]=NA
  
  ## Turn allele a and allele b into 1 and -1.  Het into 0
  #line by line if a line is a then it places 1 in hap01 for the allele
  hap01[x == a] = 1
  hap01[x == b] = -1
  hap01[x == "M"] = 0
  hap01[x == "Y"] = 0
  hap01[x == "K"] = 0
  hap01[x == "R"] = 0
  hap01[x == "W"] = 0
  hap01[x == "S"] = 0
  hap01[x== "N"]=NA
  
  return(hap01)}
  



