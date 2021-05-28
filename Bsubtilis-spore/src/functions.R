####Packages####

library(flowCore)
library(tidyverse)
library(viridis) #color palette
library(ggcyto) #visualization
library(cowplot) #multi-panel plots
library(ggridges) #ridgets
library(mclust) #for GMM clustering
library(beepr) #sounds
library(testthat) #testing functions

if(!require(dsHelper)){
  devtools::install_git("https://gitlab.com/drsudo/drsudo_helper.git")  
}
library(dsHelper) #helper functions

####Functions####

#' firstsecondElement
#'
#' @param x 
#' @description If only one element exists in a vector, the respective element is returned. If it contains more than 2, the second element is returned.
#' 
firstsecondElement<-function(x){
  ifelse(is.na(x[2]),x[1],x[2])
}

expect_equal(firstsecondElement(c(1)),1)
expect_equal(firstsecondElement(c(1,2)),2)

#' firstsecondPair
#'
#' @param x 
#' @description If only one element exists in a vector, a vector containing the respective element twice is returned. Otherwise, The vector itself is returned
#' 
firstsecondPair<-function(x){
  if(is.na(x[2])) {c(x[1],x[1])}else {x}
}

testthat::expect_equal(firstsecondPair(c(1,2)),c(1,2))
testthat::expect_equal(firstsecondPair(c(1)),c(1,1))


#' sdnorm
#'
#' @param x numeric
#' @param lambda scaling parameter 
#' @param ... passed into dnorm
#' @description dnorm with additional scaling parameter lambda used as mixing component
#' 
sdnorm <- function(x, mean=0, sd=1, lambda=1){lambda*dnorm(x, mean=mean, sd=sd)}

#' ref.ds
#'
#' @param strn strain
#' @param time time
#' @param tripl triplicate number
#' @param stn stain type
#' @param df data frame
#' @description subsets dataframe extracted from fcs file, used to get reference dataset
#'
ref.ds<-function(strn="Bs02003",time="24h",tripl=NA,stn=NA,df=NA){
  if(is.na(stn)){
    df%>%
      dplyr::filter(strain==strn)%>%
      dplyr::filter(tripl==tripl)%>%
      dplyr::filter(time==time)%>%
      dplyr::select(asinh.FSC.A,asinh.SSC.A,asinh.FL1.A)%>%
      as.matrix()
  }else{
    df%>%
      dplyr::filter(strain==strn)%>%
      dplyr::filter(time==time)%>%
      dplyr::filter(stain==stn)%>%
      dplyr::select(asinh.FSC.A,asinh.SSC.A,asinh.FL1.A)%>%
      as.matrix()
  }
}


#' findCutoffs
#'
#' @param x_mu mean of first normal
#' @param x_sd sd of first normal
#' @param x_pro scaling of first normal
#' @param y_mu mean of second normal
#' @param y_sd sd of second normal
#' @param y_pro scaling of second normal
#' @description Given parameters of two normals determines threshold value
findCutoffs <-function(x_mu,x_sd,x_pro,y_mu,y_sd,y_pro){
      Vectorize(function(x_mu,x_sd,x_pro,y_mu,y_sd,y_pro) {
            f<-function(x) sdnorm(x, m=x_mu, sd=x_sd, lambda=x_pro) - sdnorm(x, m=y_mu, sd=y_sd, lambda=y_pro)
            uniroot(f, interval=c(min(c(x_mu,y_mu)), max(c(x_mu,y_mu))))$root
      })(x_mu,x_sd,x_pro,y_mu,y_sd,y_pro)
}
