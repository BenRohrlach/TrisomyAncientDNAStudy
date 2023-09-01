#########################################################
# Analysis of chromosome read count for trisomy study   #
# Created by: AB Rohrlach                               #
#########################################################

# Required packages and custom functions ----
require(extraDistr)
require(tidyverse)
require(ggstar)
require(ggplot2)
require(FactoMineR)
require(factoextra)
require(ggrepel)
require(matrixStats)
require(Rfast)
require(ggside)

centreMatrix <- function(X,diffs=NA){
  n <- ncol(X)
  Y <- X
  
  for(i in 1:n){
    if(length(diffs)==1){
      Y[,i] <- X[,i]-mean(X[,i])
    }else{
      Y[,i] <- X[,i]-diffs[i]
    }
  }
  return(Y)
}
beta.mle_filter <- function (y,which,tol = 1e-09) {
  z <- y[!is.na(y)&(y>0)]
  x <- z[((z>(median(z)-10*IQR(z)))&(z<(median(z)+10*IQR(z))))]
  
  n <- length(x)
  sly1 <- sum(log(x))/n
  sly2 <- sum(log(1 - x))/n
  sy <- sum(x)
  sy2 <- sum(x^2)
  iniphi <- (sy - sy2)/(sy2 - sy^2/n) * (n - 1)/n
  a <- sy * iniphi/n
  b <- iniphi - a
  phi <- a + b
  lik1 <- -n * lbeta(a, b) + (a - 1) * n * sly1 + (b - 1) * 
    sly2 * n
  dera <- sly1 - digamma(a) + digamma(phi)
  derb <- sly2 - digamma(b) + digamma(phi)
  derab <- trigamma(phi)
  dera2 <- -trigamma(a) + derab
  derb2 <- -trigamma(b) + derab
  anew <- c(a, b) - c(derb2 * dera - derab * derb, -derab * 
                        dera + dera2 * derb)/(dera2 * derb2 - derab^2)
  a <- anew[1]
  b <- anew[2]
  phi <- a + b
  lik2 <- -n * lbeta(a, b) + (a - 1) * n * sly1 + (b - 1) * 
    sly2 * n
  i <- 2
  while (lik2 - lik1 > tol) {
    i <- i + 1
    lik1 <- lik2
    dera <- sly1 - digamma(a) + digamma(phi)
    derb <- sly2 - digamma(b) + digamma(phi)
    derab <- trigamma(phi)
    dera2 <- -trigamma(a) + derab
    derb2 <- -trigamma(b) + derab
    anew <- anew - c(derb2 * dera - derab * derb, -derab * 
                       dera + dera2 * derb)/(dera2 * derb2 - derab^2)
    a <- anew[1]
    b <- anew[2]
    phi <- a + b
    lik2 <- -n * lbeta(a, b) + (a - 1) * n * sly1 + (b - 
                                                       1) * sly2 * n
  }
  loglik <- lik2
  names(anew) <- c("alpha", "beta")
  return(anew[which])
}
normVector <- function(x){
  return(x/sum(x))
}
normVal <- function(x,v){
  return(x/v)
}
calcTol <- function(a,b){
  return(a*(a^2+2*a*b+a+b^2+b)*(3*a^2+a^2*(3-9*b)-8*a*b+4*b^2*(3*b+1)))
}
lwrUprBB <- function(a,b,N,k,M){
  mu <- N*a/(a+b)
  sigma <- sqrt(N*a*b*(a+b+N)/((a+b)^2*(a+b+1)))
  l <- floor(mu-k*sigma)
  u <- ceiling(mu+k*sigma)
  
  return(c(seq(l,u,by=max(floor((u-l)/M),1)),u)%>%unique()%>%sort())
}
shiftpca <- function(a,b,i){
  A <- 3*a*(-3*a^2+3*a*(b-1)+2*b*(3*b+1))/(8*b*(a+b))
  B <- (3*a^3+a^2*(3-9*b)-8*a*b+4*b^2*(3*b+1))/(8*b*(a+b))
  return(A)
}
shiftpcb <- function(a,b,i){
  A <- 3*a*(-3*a^2+3*a*(b-1)+2*b*(3*b+1))/(8*b*(a+b))
  B <- (3*a^3+a^2*(3-9*b)-8*a*b+4*b^2*(3*b+1))/(8*b*(a+b))
  return(B)
}

dat.sg.raw <- bind_rows(read_delim('~/Dropbox/MPI/Trisomy21/backup_data/result-dedup-join.txt',delim='\t',col_names=F),
                        read_delim('~/Dropbox/MPI/Trisomy21/backup_data/lost-run-dedup.txt',delim='\t',col_names=F),
                        read_delim('~/Dropbox/MPI/Trisomy21/backup_data/new-runs-2022-11-01_newresults_dedup.txt',delim='\t',col_names=F)) %>%
  `colnames<-`(c("sample",'run','SS','protocol',paste0('chr',1:22),'X','Y','mt')) %>%
  unique() %>%
  dplyr::mutate(protocol=gsub('SsLibrary 2.0 2018','ssLibrary 2.0 2018',protocol)) %>%
  dplyr::filter(!grepl('HEL',sample))

dat.sg.raw %>%
  dplyr::filter(grepl('CDM005',sample))

dat.sg <- dat.sg.raw %>%
  dplyr::mutate(total=rowSums(dplyr::select(.,starts_with('chr')))) %>%
  dplyr::rowwise() %>%
  dplyr::mutate_at(vars(contains('chr')),~ . / !! as.symbol("total")) %>%
  #dplyr::mutate(across(-Detrital_Divisor, ~ . / !! Detrital_Divisor)) %>%
  dplyr::ungroup() %>%
  dplyr::filter((total>1e3)|grepl('YUN039|ERE006|HEL002|THE106',sample),!grepl('HEL002',sample))

dat.sg.supp <- dat.sg.raw %>%
  dplyr::mutate(sample2=paste0('Individual_',1:n())) %>%
  dplyr::mutate(protocol=protocol%>%as.factor()%>%as.numeric()%>%paste0('Protocol_',.),
                sample=ifelse(grepl('CRU001|CRU013|CRU024|ERE004|LAZ019|YUN039',sample),sample,sample2)) %>%
  dplyr::select(-X,-Y,-mt,-run,-SS,-sample2)
write_tsv(dat.sg.supp,file='~/Dropbox/Project_Trisomy/Data/Rohrlach_chrCount.tsv')

dat.sg.raw <- read_tsv('~/Dropbox/Project_Trisomy/Data/Rohrlach_chrCount.tsv')

#### Prob calculations ----
Prior.t13 <- 1/7143
Prior.t18 <- 1/3226
Prior.t21 <- 1/705

SEs <- 5
M <- 1e5

posterior.calcutations <- dat.sg.raw %>%
  dplyr::mutate(total=rowSums(dplyr::select(.,starts_with('chr'))),
                CHR13=chr13/total,
                CHR18=chr18/total,
                CHR21=chr21/total) %>%
  dplyr::filter((total>1e3)|grepl('YUN039|MUM001|ERE006',sample)) %>%
  dplyr::filter(!is.na(CHR21),!is.na(CHR18),!is.na(CHR13)) %>%
  dplyr::group_by(protocol) %>%
  dplyr::mutate(CHR13_F=ifelse(total<1e5,NA,CHR13),
                CHR18_F=ifelse(total<1e5,NA,CHR18),
                CHR21_F=ifelse(total<1e5,NA,CHR21)) %>% 
  dplyr::mutate(alpha.13=beta.mle_filter(CHR13_F,which=1),
                beta.13=beta.mle_filter(CHR13_F,which=2),
                alpha.18=beta.mle_filter(CHR18_F,which=1),
                beta.18=beta.mle_filter(CHR18_F,which=2),
                alpha.21=beta.mle_filter(CHR21_F,which=1),
                beta.21=beta.mle_filter(CHR21_F,which=2)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(factor.13=CHR13/median(CHR13_F,na.rm=T),
                factor.18=CHR18/median(CHR18_F,na.rm=T),
                factor.21=CHR21/median(CHR21_F,na.rm=T)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(alpha.13_t=shiftpca(alpha.13,beta.13,1),
                beta.13_t=shiftpcb(alpha.13,beta.13,2),
                l.13=extraDistr::dbbinom(chr13,total,alpha=alpha.13,beta=beta.13,log=T),
                l.t13=extraDistr::dbbinom(chr13,total,alpha=alpha.13_t,beta=beta.13_t,log=T),
                num.13=log(1-Prior.t13)+l.13,
                num.t13=log(Prior.t13)+l.t13,
                denom.13=logSumExp(c(num.13,num.t13)),
                post.13=(num.13-denom.13)%>%exp(),
                post.t13=(num.t13-denom.13)%>%exp(),
                BF.13=num.t13-num.13) %>%
  dplyr::mutate(alpha.18_t=shiftpca(alpha.18,beta.18,1),
                beta.18_t=shiftpcb(alpha.18,beta.18,2),
                l.18=extraDistr::dbbinom(chr18,total,alpha=alpha.18,beta=beta.18,log=T),
                l.t18=extraDistr::dbbinom(chr18,total,alpha=alpha.18_t,beta=beta.18_t,log=T),
                num.18=log(1-Prior.t18)+l.18,
                num.t18=log(Prior.t18)+l.t18,
                denom.18=logSumExp(c(num.18,num.t18)),
                post.18=(num.18-denom.18)%>%exp(),
                post.t18=(num.t18-denom.18)%>%exp(),
                BF.18=num.t18-num.18) %>%
  dplyr::mutate(alpha.21_t=shiftpca(alpha.21,beta.21,1),
                beta.21_t=shiftpcb(alpha.21,beta.21,2),
                l.21=extraDistr::dbbinom(chr21,total,alpha=alpha.21,beta=beta.21,log=T),
                l.t21=extraDistr::dbbinom(chr21,total,alpha=alpha.21_t,beta=beta.21_t,log=T),
                num.21=log(1-Prior.t21)+l.21,
                num.t21=log(Prior.t21)+l.t21,
                denom.21=logSumExp(c(num.21,num.t21)),
                post.21=(num.21-denom.21)%>%exp(),
                post.t21=(num.t21-denom.21)%>%exp(),
                BF.21=num.t21-num.21) %>%
  dplyr::select(sample,total,protocol,
                CHR21,BF.21,factor.21,post.21,post.t21,
                CHR18,BF.18,factor.18,post.18,post.t18,
                CHR13,BF.13,factor.13,post.13,post.t13) %>%
  dplyr::select(sample,total,protocol,CHR13,post.t13,BF.13,factor.13,CHR18,post.t18,BF.18,factor.18,CHR21,post.t21,BF.21,factor.21,post.13,post.18,post.21)

posterior.scan <- posterior.calcutations %>%
  dplyr::filter(grepl('CRU00[12]|YUN039|ERE00[46]|HKI002|CRU013|CRU024|MUM001|TAF|EML|^PN0|LAZ|CDM005',sample)) %>%
  dplyr::filter(any(c(post.t13>post.13,post.t18>post.18,post.t21>post.21))) %>%
  dplyr::arrange(factor.13,factor.18,factor.21) %>%
  dplyr::select(sample,total,post.t13,post.13,factor.13,post.t18,post.18,factor.18,post.t21,post.21,factor.21)

summary.plot.tib.t <- posterior.calcutations %>%
  dplyr::filter(grepl(paste0(c('PN07','YUN039.A','HKI002.A','LAZ019.A','CRU013.A','CRU001.A','ERE004.A','CRU024.A','CDM005.A'),collapse='|'),sample))

####### Diagnostic ---------
repdat.auto.raw <- dat.sg.raw %>%
  dplyr::mutate(total=rowSums(dplyr::select(.,starts_with('chr')))) 

repdat.auto <- repdat.auto.raw %>%
  dplyr::mutate(total=rowSums(dplyr::select(.,starts_with('chr')))) %>%
  dplyr::rowwise() %>%
  dplyr::ungroup() %>%
  gather(chr,proportion,chr1:chr22)

out.repdat.in <- repdat.auto %>%
  dplyr::filter(total>1e3|grepl('MUM|THE106|HEL',sample)) %>% 
  dplyr::mutate(chr=factor(chr,levels=paste0('chr',1:22)),
                fct=substr(sample,1,3)) %>%
  dplyr::group_by(sample,chr) %>%
  dplyr::mutate(p=proportion/total) %>% 
  dplyr::group_by(chr,protocol) %>%
  dplyr::mutate(a=beta.mle_filter(ifelse(total<1e5,NA,proportion/total),which=1),
                b=beta.mle_filter(ifelse(total<1e5,NA,proportion/total),which=2),
                n=total) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(mp=n*a/(a+b)/n,
                sdp=sqrt(n*a*b*(a+b+n)/(((a+b)^2)*(a+b+1)))/n,
                N=n(),
                Z=(p-mp)/(sdp)) %>% 
  dplyr::mutate(fct=gsub('MSI','MUM',fct),
                pos=grepl(paste0('HEL|SCH|',unique(substr(posterior.scan$sample,1,30)),collapse='|'),sample)) %>%
  dplyr::group_by(fct,protocol) %>% 
  dplyr::mutate(keep=any(pos)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(keep)

small.filter <- out.repdat.in %>%
  dplyr::mutate(pos=grepl(paste0(unique(substr(posterior.scan$sample,1,30)),collapse='|'),sample)) %>%
  dplyr::select(sample,pos,fct) %>%
  unique() %>%
  dplyr::group_by(pos,fct) %>%
  dplyr::sample_n(min(10,n())) %>%
  dplyr::ungroup() %>%
  dplyr::pull(sample) %>%
  paste0(collapse='|')

out.repdat <- out.repdat.in %>%
  dplyr::filter(grepl(small.filter,sample),!grepl('SCHasd',sample)) %>%
  dplyr::mutate(shapeCol=case_when(pos ~ 'Important',
                                   abs(Z)>=2 ~ 'Significant',
                                   T ~ 'Control')) %>%
  ggplot(aes(x=chr,y=Z,shape=substr(sample,1,8),fill=substr(sample,1,8),col=shapeCol))+
  theme_bw()+
  geom_star(aes(starshape=substr(sample,1,8)),size=2.5)+
  xlab(NULL)+
  ylab('Z-score')+
  geom_hline(yintercept=c(-2,2),linetype='dashed')+
  scale_color_manual(name=NULL,values=c('Grey','Red','Black'),guide=F)+
  scale_fill_discrete(name=NULL)+
  scale_starshape_manual(name=NULL,values=rep(setdiff(1:30,c(16,22,14,26)),1e4))+
  scale_shape_manual(name=NULL,values=rep(1:22,1e4))+
  theme(legend.position='bottom',
        axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(starshape=guide_legend(ncol=7,byrow=TRUE),
         fill=guide_legend(ncol=7,byrow=TRUE))+
  facet_wrap(facets=vars(fct),ncol=2,scales='free_y')+
  scale_y_continuous(expand=c(0.1,0.1))
out.repdat + theme(legend.position='none')