
#   This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

  
########################################################################
###### estimation of variances  out of an hierarchical experiment #####
######################################################################

## function to estimate variances of each level out of a signal matrix
estVarLevData <- function(k, repli, vals, fac){

    ## reformat to factor in R
    facnew <- apply(fac,1,function(x){as.factor(unlist(x))})

    colnames(facnew)[1:k-1] <- letters[1:k-1]

    ## construct a formula for ANOVA method according to the number of levels
    form <- "a"
    m <- 2
    while(m <k) {
      form <- paste(form,letters[m],sep="/")
      m <- m+1
    }
    form<-as.formula(paste("x", form, sep=" ~ "))
    
    ## prepare class data.frame for one signal and estimate variances of each level for one signal
    varsigs<-t(apply(vals,1, function(val){estVarLevSig(k, repli=repli, formAn=form, data.frame(x=val, facnew ))}))
    
    ## each column according to one variance
    colnames(varsigs)<-c(paste("sigQu",k:1,sep="_"),"sigTotalQu")
    
    ## each row according to one signal
    rownames(varsigs)<-rownames(vals)
    return(varsigs)
}

## function to estimate variances of each level out of one signal
estVarLevSig<-function(k, repli, formAn, data) {
  
  #method ANOVA for one signal
  ano<-summary(aov(formAn, data))
  #ano[[1]]
  
  #estimating unbiased estimator for variances of all levels and total variance for one signal 
  #out of the given ANOVA table (Mean Square Deviations and degrees of freedom)
  sig<-ano[[1]][k,"Mean Sq"]
  for(i in 1:(k-1)){sig<-c(sig,(ano[[1]][k-i,"Mean Sq"]-ano[[1]][(k-i+1),"Mean Sq"])/prod(repli[k:(k-i+1)]))}
  sig<-c(sig,sum(sig))
  return(sig)
}

#####################################################################################
###### nested ANOVA to find differentiell signals in hierarchical experiments  #####
###################################################################################

#function to find differentiell signals out of a signal matrix using nested ANOVA 
#where the first level in hierarchical experiment is the different background, and the secound is the biological variation
diffAnovData<-function(k, repli, file1,file2){
  
  #read measured data out of file
  vals<-read.csv(file=file1,  row.names=1)
  
  #read factors of all levels out of file
  fac<-read.csv(file=file2)
  #support to factor in R
  facnew<-apply(fac,1,function(x){as.factor(x)})
  colnames(facnew)[1:k-1] <-letters[1:k-1]
  
  #construct a formula for ANOVA method according to the number of levels
  form<-"a"
  m<-2
  while(m <k){
    form<-paste(form,letters[m],sep="/")
    m<-m+1
  }
  form<-as.formula(paste("x", form, sep=" ~ "))
  
  #prepare class data.frame for one signal and ask for difference of one signal in to backgrounds
  sigs<-t(apply(vals,1, function(val){diffAnovSig(k, formAn=form, data.frame(x=val, facnew))}))
  return(sigs)
}


#function to find differentiell signal using nested ANOVA 
#where the first level in hierarchical experiment is the different background, and the secound is the biological variation
diffAnovSig<-function(k, formAn, data){
  
  #method ANOVA for one signal
  ano<-summary(aov(formAn, data))
  ano[[1]]
  
  #ask for difference by comparing estimated variances of the first two levels 
  signif<-1-pf(q=ano[[1]][1,"Mean Sq"]/ano[[1]][2,"Mean Sq"], df1=ano[[1]][1,"Df"],df2=ano[[1]][2,"Df"])
  return(signif)
  
}


############################################################
###  Plot "Distribution of estimated variances" ###########
##########################################################

plotStripForVarLev<-function(k, levnames, estvar, file){
  
  #only signals without negativ variances
  pos<-(estvar[,2:k]>0)
  VarNeg<-unique(unlist(sapply(1:(k-2), function(i){which(pos[,i]==FALSE)})))       
  SiganzTrue<- nrow(estvar)-length(VarNeg)
  
  values<-unlist(as.vector(estvar))
  names<-rep(1:(k+1),each=nrow(estvar))
  colpl<-rainbow(k+1)
  ymax<-max(estvar)
  
  #plot in file
  
  if(!missing(file)) {
    pdf(file)
  }
  
  par(oma=c(1,1,1,1))
  stripchart(values~names, method="jitter" ,jitter=0.3,vertical=TRUE,
             pch=8,cex=0.2, ylim=c(-0.2,ymax),
             #col= colpl, 
             at=1:(k+1), 
             group.names=levnames,
             xlab="variances", ylab="signal variances",cex.lab=1.5,
             main="distribution of signal variances", cex.main=1.9)
  
  sapply(1:(k+1), function(i){segments(i-0.4,mean(as.numeric(estvar[-VarNeg,i])), i+0.4, mean(as.numeric(estvar[-VarNeg,i])), lwd=1.8)})
  sapply(1:(k+1), function(i){text(i,-0.2,round(1000*mean(as.numeric(estvar[-VarNeg,i])))/1000, cex=1.2
                                   #, col= colpl[i]
                                   )})
  
  abline(v=k+0.5)
     
  if(!missing(file)) {
       dev.off()
  }      
}

################################################################################################################
###  Powercalculations for given biological and technical variances and experiment-Design choice    ###########
##############################################################################################################

#2-level hierarchical experiment design
#5 functions in one
#A: calculate possible number of technical replicates E for given biological replicates N, delta, power and alpha (given the maximal number of technical replicates)
#B: calculate possible number of biological replicates N, if number of technical replicates E for each biological, delta, power and alpha is given
#C: calculate power for given number of replicates (N, E) and delta
#D: calculate delta for given number of replicates (N, E) and power, alpha
#E: calculate alpha for given number of replicates (N, E), delta and power

#given: 4 of 5 parameters a)number of technical replicates E, b)number of biological replicates N, c)delta and d)power, e)alpha
#always given: biological and technical variance
#get: the free parameter

power.hierarch.ttest<-function(N,Emax,p,d,alpha,sigma_bio_q,sigma_tech_q, type){
  
  switch(type,
         
         #calculation of number of technical replicates, given alpha, power, delta, number of biological replicates, maximal number of technical replicates biological and technical variance
         replicates={
           
           E<-1
           pow<-0
           while((pow<p)&& (E<=Emax) ){
             
             pow<-power.t.test(n=N, delta=d, sd=sqrt(sigma_bio_q+sigma_tech_q/E), alpha,alternative="two.sided", strict=TRUE)$power;
             E<-E+1;
           }
           findE<-E-1;
           if(pow<p){res<-NULL}
           else{
             res<-list(c("biol. Rep."=N,"techn./biol. Rep."=findE),Power=round(pow,2), Delta=d, Alpha=alpha);}
           
         },
         #calculation of number of biological replicates, given alpha, power, delta, number of technical replicates, maximal number of technical replicates biological and technical variance
         n={
           
           n<-power.t.test(n=NULL, delta=d, sd=sqrt(sigma_bio_q+sigma_tech_q/Emax), alpha,power=p,alternative="two.sided", strict=TRUE)$n;
           res<-list(c("biol. Rep."=ceiling(n),"techn./biol. Rep."=Emax),Power=p, Delta=d, Alpha=alpha);
           
         },
         #calculation of power, given alpha, delta, number of biological and technical replicates, biological and technical variance
         power={
           
           pow<-power.t.test(n=N, delta=d, sd=sqrt(sigma_bio_q+sigma_tech_q/Emax), alpha,power=NULL,alternative="two.sided", strict=TRUE)$power;
           res<-list(c("biol. Rep."=N,"techn./biol. Rep."=Emax),Power=round(pow,2), Delta=d, Alpha=alpha);
           
         },
         #calculation of true difference in mean, given alpha, power, number of biological and technical replicates, biological and technical variance
         delta={
        
           delt<-power.t.test(n=N,delta=d,power=p, sd=sqrt(sigma_bio_q+sigma_tech_q/Emax), alpha,alternative="two.sided", strict=TRUE)$delta;
           res<-list(c("biol. Rep."=N,"techn./biol. Rep."=Emax), Power=p,Delta=round(delt,2), Alpha=alpha);
         },
         #calculation error Type I, given effect, power, number of biological and technical replicates, biological and technical variance
         alpha={
           
           alph<-power.t.test(n=N, delta=d,power=p, sd=sqrt(sigma_bio_q+sigma_tech_q/Emax),alpha, alternative="two.sided", strict=TRUE)$sig.level;
           res<-list(c("biol. Rep."=N,"techn./biol. Rep."=Emax), Power=p,Delta=delta,Alpha=round(alph,3));
         }
         
         
  )
  return(res)
  
}



#2-level hierarchical experiment design
#function calculates the possibele combinations of biological and technical replicates
#given: alpha, power, delta, biological and technical variance, maximal number of biological and technical repicates
#get: list of given parameter power, delta and matrix of possible biological and technical replicates 

supportMat<-function(alpha, power, sigma_tech_q, sigma_bio_q, d, Emax, Nmax){
  
  mat<-sapply(2:Nmax, function(N){
    
    temp<-power.hierarch.ttest(N,Emax,p=power,d,alpha,sigma_bio_q,sigma_tech_q, type="replicates")[[1]]
    
  })
  
  if(is.list(mat)){
    
    if(is.null(unlist(mat))){
      mat<-unlist(mat)
    }else
    {mat<-matrix(unlist(mat), nrow=2, ncol=length(unlist(mat))/2, dimnames=list(c("biol. Rep.","techn./biol. Rep.")));
    }
  }
  
  #alle raus, die trotz erhöhter biol. Replikatanzahl auch nur ein techn replikat haben
  
  if(is.null(mat)==FALSE){
    n<-1
    over<-c()
    while(n<ncol(mat)){
      overN<-which(mat["techn./biol. Rep.",]==mat["techn./biol. Rep.",n]);
      over<-c(over,overN[-1]);
      n<-overN[length(overN)]+1;  
    }
    mat<-mat[,-over]
  }
  
  return(list(possibilities=mat,"Power>="=power ,Delta=d, Alpha=alpha))
}



#2-level hierarchical experiment design
#function calculates the possibele combinations of biological and technical replicates with minimal costs
#given: alpha, power, delta, biological and technical variance, maximal number of biological and technical repicates, cost relation of biological and technical replicates
#get: list of given parameter power, delta and possible combination of biological and technical replicates with minimal

minCostPoss<-function(alpha, power, sigma_tech_q, sigma_bio_q, d, Emax, Nmax, cost ){
  
  temp<-supportMat(alpha, power, sigma_tech_q, sigma_bio_q, d, Emax, Nmax)
  
  if(is.null(temp$possibilities)){
    res<-temp$possibilities
  }else{
    if(is.matrix(temp$possibilities)){
    rep<-rbind(temp$possibilities["biol. Rep.",], temp$possibilities["biol. Rep.",]*temp$possibilities["techn./biol. Rep.",])    
    w<-which(colSums(rep*cost)==min(colSums(rep*cost)))
    res<-temp$possibilities[,w]
    }else{
      res<-temp$possibilities
    }
    if(is.matrix(res)){
      res<-t(res)
    }
  }
  
  return (list(res, CostRelation=cost,power=round(power,2), Delta=d))
  
}


# tblofCosts<-function(pow, delta, sigmabioq, sigmatotq, anztrep, relation){
# 
#   stud<-sapply(delta,function(d){rep(ceiling(power.t.test( sd=sqrt(sigmatotq),delta=d,power=pow)$n),2)})
#   coststud<-apply(stud,2,function(a){sum(a*relation)})
#   stud<-rbind(stud, coststud)
#   rownames(stud)<-c("No. plants", "No. measurements", "costs")
#   colnames(stud)<-delta
#   anov<-sapply(delta,function(d){c(ceiling(power.t.test( sd=sqrt(sigmabioq),delta=d,power=pow)$n),
#                                    ceiling(power.t.test( sd=sqrt(sigmabioq),delta=d,power=pow)$n)*anztrep)})
#   costanov<-apply(anov,2,function(a){sum(a*relation)})
#   anov<-rbind(anov, costanov)
#   rownames(anov)<-c("No. plants", "No. measurements", "costs")
#   colnames(anov)<-delta
#   string1<-paste("total variance = ", sigmatotq)
#   string1
#   stud
#   string2<-paste("biological variance = ", sigmabioq)
#   string2
#   anov
#   res<-list(list(string1,stud),list(string2,anov))
#   return(res)
# 
# }

