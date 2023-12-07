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

getPs <- function(p,x=x,n=n){
  binom.test(x=x,n=n,p=p)$p.value %>% return()
}
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


#### Read in data ----
dat.sg.raw <- read_tsv('Rohrlach_chrCount.tsv')

dat.sg <- dat.sg.raw %>%
  dplyr::mutate(total=rowSums(dplyr::select(.,starts_with('chr')))) %>%
  dplyr::rowwise() %>%
  dplyr::mutate_at(vars(contains('chr')),~ . / !! as.symbol("total")) %>%
  dplyr::ungroup() %>%
  dplyr::filter((total>1e3))

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
  dplyr::filter(grepl('CRU001|YUN039|ERE004|HKI002|CRU013|CRU024|LAZ019|PN07',sample)) %>%
  dplyr::filter(any(c(post.t13>post.13,post.t18>post.18,post.t21>post.21))) %>%
  dplyr::arrange(factor.13,factor.18,factor.21) %>%
  dplyr::select(sample,total,post.t13,post.13,factor.13,post.t18,post.18,factor.18,post.t21,post.21,factor.21)

posterior.scan %>%
  print(n=1e4)

####### Diagnostic ---------
repdat.auto.raw <- dat.sg.raw %>%
  dplyr::mutate(total=rowSums(dplyr::select(.,starts_with('chr')))) 

repdat.auto <- repdat.auto.raw %>%
  dplyr::mutate(total=rowSums(dplyr::select(.,starts_with('chr')))) %>%
  dplyr::rowwise() %>%
  dplyr::ungroup() %>%
  gather(chr,proportion,chr1:chr22)

out.repdat.in <- repdat.auto %>%
  dplyr::filter(total>1e3) %>% 
  dplyr::mutate(chr=factor(chr,levels=paste0('chr',1:22))) %>%
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
  dplyr::mutate(pos=grepl(paste0(unique(substr(posterior.scan$sample,1,30)),collapse='|'),sample)) %>%
  dplyr::group_by(fct,protocol) %>% 
  dplyr::mutate(keep=any(pos)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(keep) %>%
  dplyr::group_by(fct,protocol) %>% 
  dplyr::mutate(site_name=grep('Indiv',sample,value=T,invert=T)%>%substr(1,3)%>%unique()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(point_type=ifelse(grepl('Individ',sample),'Reference',sample) %>%
                  factor(levels=c('CRU001','CRU013','CRU024','ERE004','HKI002','LAZ019','YUN039','Reference')))
  
set.seed(12345)
small.filter <- out.repdat.in %>%
  dplyr::mutate(pos=grepl(paste0(unique(substr(posterior.scan$sample,1,30)),collapse='|'),sample)) %>%
  dplyr::select(sample,pos,fct) %>%
  unique() %>%
  dplyr::group_by(pos,fct) %>%
  dplyr::sample_n(min(10,n())) %>%
  dplyr::ungroup() %>%
  dplyr::pull(sample) %>%
  paste0(collapse='|')

figureS1.tibble <- read_tsv(file='Rohrlach_FigureS1.tsv')

figureS1.tibble %>%
  ggplot(aes(x=chr,y=Z,fill=fct,col=shapeCol))+
  theme_bw()+
  geom_star(aes(starshape=point_type),size=2.5)+
  xlab(NULL)+
  ylab('Z-score')+
  geom_hline(yintercept=c(-2,2),linetype='dashed')+
  scale_color_manual(name=NULL,values=c('Grey','Red','Black'),guide=F)+
  scale_fill_discrete(name=NULL,guide=FALSE)+
  scale_starshape_manual(name=NULL,values=rep(setdiff(1:30,c(16,22,14,26)),1e4))+
  scale_shape_manual(name=NULL,values=rep(1:22,1e4))+
  theme(legend.position='bottom',
        axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(facets=vars(site_name),ncol=2,scales='free_y')+
  scale_y_continuous(expand=c(0.1,0.1))

#### Analyse fold increases ----
fromTo <- posterior.calcutations %>%
  dplyr::filter(factor.21>1.375,sample%in%PCA.like.plot.ID) %>%
  dplyr::select(sample,factor.18,factor.21) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(to18=c(1.07,1.08,1.08,0.96,0.9,0.95,0.88),
                to21=c(1.512947,1.472,1.435,1.52,1.46,1.42,1.494543),
                to18_1col=c(1.07,1.08,1.08,0.96,0.9,0.95,0.88),
                to21_1col=c(1.512947,1.472,1.435,1.52,1.46,1.42,1.494543)) %>%
  dplyr::filter(!grepl('VEN',sample))

PCA.like.plot.ID <- repdat.auto %>%
  dplyr::mutate(chr=factor(chr,levels=paste0('chr',1:22)),
                fct=substr(sample,1,3)) %>%
  dplyr::group_by(sample,chr) %>%
  dplyr::mutate(p=proportion/total) %>%
  dplyr::group_by(chr,protocol) %>%
  dplyr::mutate(a=beta.mle_filter(ifelse(total<1e4,NA,proportion/total),which=1),
                b=beta.mle_filter(ifelse(total<1e4,NA,proportion/total),which=2),
                n=total) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(mp=n*a/(a+b)/n,
                sdp=sqrt(n*a*b*(a+b+n)/(((a+b)^2)*(a+b+1)))/n,
                N=n(),
                Z=(p-mp)/(sdp)) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(sigZ=sum(abs(Z)>2)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(sigZ<=1,total>5e3) %>%
  dplyr::select(sample,protocol,chr,Z,total,mp,p) %>%
  dplyr::mutate(fact=p/mp*100) %>%
  dplyr::arrange(fact) %>%
  dplyr::pull(sample) %>%
  unique() %>%
  c('PN07',.)

gg.1821.factor.tibble <- read_tsv(file='Figure1.tsv')

gg.1821.factor.tibble %>%
  ggplot(aes(x=factor.18,y=factor.21,label=substr(sample,1,6)))+
  theme_bw()+
  geom_point(alpha=0.25)+
  # stat_ellipse(level=0.68,alpha=0.5,col='red',size=1)+
  # stat_ellipse(level=0.90,alpha=0.5,col='red',size=1)+
  # stat_ellipse(level=0.95,alpha=0.5,col='red',size=1)+
  # stat_ellipse(level=0.99,alpha=0.5,col='red',size=1)+
  # stat_ellipse(level=0.999,alpha=0.5,col='red',size=1)+
  geom_point(data=fromTo)+
  geom_point(data=dplyr::filter(posterior.calcutations,factor.18>1.4,sample%in%PCA.like.plot.ID))+
  geom_point(data=fromTo)+
  geom_segment(data=fromTo,aes(x=factor.18,xend=to18,y=factor.21,yend=to21))+
  geom_label(data=fromTo,aes(x=to18,y=to21))+
  geom_label_repel(data=dplyr::filter(posterior.calcutations,factor.18>1.4,sample%in%PCA.like.plot.ID),
                   max.overlaps=1e4,box.padding=1,size=4)+
  xlab('\nFold-increase in reads mapping to chromosome 18')+
  ylab('Fold-increase in reads mapping to chromosome 21\n')+
  # geom_xsidedensity(aes(y=stat(density),col=protocol))+
  # geom_ysidedensity(aes(x=stat(density),col=protocol))+
  # scale_color_discrete(name=NULL)+
  theme(legend.position='bottom',
        text = element_text(size = 20))

#### Prevalence calculations ----
tot.n <- 9855

t13.n <- 0
t18.n <- 1
t21.n <- 6

t13.p <- t13.n/tot.n
t18.p <- t18.n/tot.n
t21.p <- t21.n/tot.n

prev <- c(1/7143,1/3226,1/113,1/705,1/1539)

binom.test(t13.n,tot.n,p=prev[1])
binom.test(t18.n,tot.n,p=prev[2])
binom.test(t21.n,tot.n,p=prev[3])
binom.test(t21.n,tot.n,p=prev[4])
binom.test(t21.n,tot.n,p=prev[5])


age.t21.all <- 1/(1e4/c(7.8,6.9,7.9,12.9,38,121.8,14.2))
age.t21 <- 1/(1e4/c(6.5,6.3,6.9,11.5,28.6,88.6,12.6))

prev.plot.tibble %>%
  ggplot(aes(x=x,y=yall))+
  theme_bw()+
  geom_point(pch=19,size=3)+
  geom_line(aes(group=1))+
  geom_hline(yintercept=t21.p,col='blue',linewidth=2,linetype='dashed')+
  geom_hline(yintercept=age.t21[7],col='red',size=2,linetype='dashed')+
  scale_y_log10(breaks=c(0.0006,0.001,0.0025,0.005,0.0075,0.01,0.013),lim=c(0.0005,0.013),expand=c(0,0))+
  xlab('\nAge of Mother (years)')+
  ylab('Rate of Prevalence of Trisomy 21\n')+
  theme(text = element_text(size = 20))


