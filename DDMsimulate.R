library(patchwork)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(grid)
library(magrittr)


sim_paths <- function(v=0.6,a=3,b=0.5,s=0.1,nsim=100,tmax=6){

  # Required functions
  rep.row <- function(x,n){
    matrix(rep(x,each=n),nrow=n)
  }

  first.upper <- function(x,bnd){
    r <- min(which(x>bnd))
    return(r)
  }

  first.lower <- function(x,bnd){
    r <- min(which(x<bnd))
    return(r)
  }

  # v: drift
  # a: boundary
  # z: bias
  # s: scaling parameter
  time <- seq(0,tmax,0.01)
  start <- b*a

  # X: accumulated evidence
  # Diffusion process: dX_t = drift*dt + s*dW_t
  # dW_t is a wiener process or borwnian motion where there is Gaussian white noise
  X <- matrix(rnorm(n=nsim*(length(time)-1),sd=s),nrow=nsim,ncol=length(time)-1)
  X <- cbind(rep(start,nsim),X)
  X <- t(apply(X,1,cumsum))
  X <- rep.row(v*time, n=nsim) + X

  # To determine first passage time (where the accumulator first hits the decision boundary)
  end <- length(time)
  r_upper <- vector(length=nsim)
  r_lower <- vector(length=nsim)
  for (i in 1:nsim){
    # check for upper boundary
    r_upper[i] <- first.upper(X[i,],a)
    if (is.finite(r_upper[i])){
      #X[i,r_upper[i]] <- a
      X[i,r_upper[i]:end] <- NA
    }

    #check for lower boundary
    r_lower[i] <- first.lower(X[i,],0)
    if (is.finite(r_lower[i])){
      #X[i,r_lower[i]] <- 0
      X[i,r_lower[i]:end] <- NA
    }
  }

  df <- as.data.frame(X) %>%
    mutate(nsim = row_number(),ref=1) %>%
    pivot_longer(cols = contains("V"),names_to="time",values_to="V") %>%
    mutate(time = as.numeric(str_remove_all(time,"V")))

  # only take the first instance where it hits the boundary
  upper_bnd <- r_upper[!(r_lower<r_upper)]
  lower_bnd <- r_lower[(r_lower<r_upper)]

  # ensure a response is made
  upper_bnd <- upper_bnd[is.finite(upper_bnd)]
  lower_bnd <- lower_bnd[is.finite(lower_bnd)]

  ## REFERENCE simulated data ##
  set.seed(2741) #keep this constant for each set of parameters (but varies with different nsim)
  v0 <- 0.7
  a0 <- 3
  b0 <- 0.5
  X0 <- matrix(rnorm(n=nsim*(length(time)-1),sd=s),nrow=nsim,ncol=length(time)-1)
  X0 <- cbind(rep(start,nsim),X0)
  X0 <- t(apply(X0,1,cumsum))
  X0 <- rep.row(v0*time, n=nsim) + X0

  # To determine first passage time (where the accumulator first hits the decision boundary)
  end <- length(time)
  r_upper0 <- vector(length=nsim)
  r_lower0 <- vector(length=nsim)
  for (i in 1:nsim){
    # check for upper boundary
    r_upper0[i] <- first.upper(X0[i,],a0)
    if (is.finite(r_upper0[i])){
      #X0[i,r_upper0[i]] <- a0
      X0[i,r_upper0[i]:end] <- NA
    }

    #check for lower boundary
    r_lower0[i] <- first.lower(X0[i,],0)
    if (is.finite(r_lower0[i])){
      #X0[i,r_lower0[i]] <- 0
      X0[i,r_lower0[i]:end] <- NA
    }
  }

  df0 <- as.data.frame(X0) %>%
    mutate(nsim = row_number(), ref=0) %>%
    pivot_longer(cols = contains("V"),names_to="time",values_to="V") %>%
    mutate(time = as.numeric(str_remove_all(time,"V")))

  # only take the first instance where it hits the boundary
  upper_bnd0 <- r_upper0[!(r_lower0<r_upper0)]
  lower_bnd0 <- r_lower0[(r_lower0<r_upper0)]

  # ensure a response is made
  upper_bnd0 <- upper_bnd0[is.finite(upper_bnd0)]
  lower_bnd0 <- lower_bnd0[is.finite(lower_bnd0)]

  ## END REFERENCE SIM##

  #combine data
  #Add extra group_reference indexing for plotting purposes
  dat <- rbind(df,df0)
  dat %<>% mutate(grp = paste0(ref,"_",nsim))

  upper_bnd %<>% data.frame(.) %>% rename("val"=".") %>% mutate(ref=1)
  upper_bnd0 %<>% data.frame(.) %>% rename("val"=".") %>% mutate(ref=0)
  upper <- rbind(upper_bnd,upper_bnd0)

  lower_bnd %<>% data.frame(.) %>% rename("val"=".")  %>% mutate(ref=1)
  lower_bnd0 %<>% data.frame(.) %>% rename("val"=".") %>% mutate(ref=0)
  lower <- rbind(lower_bnd,lower_bnd0)

  #Combine simulation and reference simulation

  res <- list(df = dat,upper_bnd=upper,
              lower_bnd=lower,bnd=rbind(a,a0),
              start=rbind(a*b,a0*b0),v=rbind(v,v0))
  return(res)
}

plot_ddmsim <- function(res,histogram=TRUE){

  nsim <- max(res$df$nsim)
  pctA <- nrow(res$upper_bnd[res$upper_bnd$ref==1,])/nsim*100
  pctB <- nrow(res$lower_bnd[res$lower_bnd$ref==1,])/nsim*100

  pctA0 <- nrow(res$upper_bnd[res$upper_bnd$ref==0,])/nsim*100
  pctB0 <- nrow(res$lower_bnd[res$lower_bnd$ref==0,])/nsim*100

  col1 <- alpha("darkgray",alpha=.7)
  col1_alt <- ifelse(nsim<120,"black","white")
  col2 <- alpha("#2c3e50",alpha=.7) #dark blue to match shiny flatly theme
  col2_op <- "#2c3e50"
  col2_alt <- ifelse(nsim<120,"black","white") 
  col2_txt <- "#5b80a4"
  col2_line <- "#364d63"

  if (histogram==TRUE){
    p_upper <- res$upper_bnd %>%
      ggplot(aes(x=val,fill=as.factor(ref))) +
      geom_histogram(color="gray",alpha=.6,bins=100,position="identity") +
      scale_fill_manual(values=c(col1,col2)) +
      scale_x_continuous(limits = c(0, 600), oob = scales::oob_keep) +
      theme_void() + theme(legend.position="none")

    #dmatch ylim of lower bound histogram to first
    upper_y <- max(ggplot_build(p_upper)$data[[1]][,2])

    p_lower <- res$lower_bnd %>%
      ggplot(aes(x=val,fill=as.factor(ref))) +
      geom_histogram(color="gray",alpha=.6,bins=100,position="identity") +
      scale_fill_manual(values=c(col1,col2))+
      scale_y_reverse(limits = c(upper_y, 0), oob = scales::oob_keep) +
      scale_x_continuous(limits = c(0, 600), oob = scales::oob_keep) +
      theme_void()+ theme(legend.position="none")

  } else {
    ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)

    p_upper <- res$upper_bnd %>%
      ggplot(aes(x=val,group=ref,fill=as.factor(ref))) +
      geom_density(alpha=.5,colour='darkgray') +
      scale_fill_manual(values=c(col1,col2))+
      #scale_y_continuous(limits = c(0, 0.02), oob = scales::oob_keep) +
      scale_x_continuous(limits = c(0, 600), oob = scales::oob_keep) +
      theme_void() + theme(legend.position="none")

    upper_y <- ceiling_dec(max(ggplot_build(p_upper)$data[[1]][,4]),3)

    p_lower <-res$lower_bnd %>%
      ggplot(aes(x=val,group=ref,fill=as.factor(ref))) +
      geom_density(alpha=.5,colour='darkgray') +
      scale_fill_manual(values=c(col1,col2)) +
      scale_y_reverse(limits = c(upper_y, 0), oob = scales::oob_keep) +
      scale_x_continuous(limits = c(0, 600), oob = scales::oob_keep) +
      theme_void() + theme(legend.position="none")
  }

  #plot
  bnd <- res$bnd
  start <- res$start
  v <- res$v
  p_ddm <- res$df %>%
    ggplot(aes(x=time,y=V,group=grp,col=as.factor(ref))) +
    geom_line(aes(col=factor(ref)), alpha=0.2) + scale_color_manual(values=c(col1,col2))+ theme(legend.position = "none") +
    scale_y_continuous(limits = c(0, 4), oob = scales::oob_keep, breaks=c(bnd[2],start[2]),
                       labels=c(expression(a[0]),expression(a[0]*b[0]))) +
    scale_x_continuous(limits = c(0, 600), oob = scales::oob_keep) +
    geom_hline(yintercept=bnd[2], color = "black", size=0.7) + #reference boundary
    geom_hline(yintercept=bnd[1], color = col2_line,size=0.7) + #boundary
    geom_hline(yintercept=start[2], linetype="dashed", color = "black", size=0.7) + #reference start point
    geom_hline(yintercept=start[1], linetype="dashed", color = col2_line, size=0.7) + #start point
    theme_classic() + theme(axis.ticks.x=element_blank(),
                            #axis.ticks.y=element_blank(),
                            axis.text.x=element_blank(),
                            axis.text.y=element_text(size=14),
                            axis.title.y=element_blank(),
                            axis.line.y=element_line(size=0.2),
                            axis.title.x=element_blank(),
                            legend.position="none") +
    annotate("segment",x=-30,xend=0,y=start[2]-0.15,yend=start[2]-0.15,colour="black", # Add Ter label + bracket
             size=0.8,arrow=arrow(ends="both",length = unit(0.15, "cm")))+
    annotate("text",x=-15,y=start[2]-0.3,colour="black",label="T[er]",parse=TRUE)+
    annotate("text",x=530,y=0,label=paste0("Response B (",round(pctB0,1),"%)")) + # Add response rates text
    annotate("text",x=600,y=0,label=paste0("(",round(pctB,1),"%)"),colour=col2_txt) +
    annotate("text",x=530,y=4-0.15,label=paste0("Response A (",round(pctA0,0),"%)")) +
    annotate("text",x=600,y=4-0.15,label=paste0("(",round(pctA,0),"%)"),colour=col2_txt) +
    #annotate("text",x=-24,y=bnd[2]+0.2,label="a[0]",parse=TRUE) + # Add boundary label
    #annotate("text",x=-20,y=start[2]+.15,label="a[0]*b[0]", parse=TRUE) + # Add start point label
    annotate("segment", x = 0, xend = 50, y = start[1], yend = start[1]+v[1]*0.5, # Add drift rate labels + vectors
             colour = col2_alt, size = 0.8, arrow = arrow(length = unit(0.3, "cm"))) +
    #annotate("text",x=25,y=start[1]+0.15+v[1]*0.25,label="v", parse=TRUE, size=5,colour=col2_alt) +
    annotate("text",x=25,y=start[2]+0.18+v[2]*0.25,label="v[0]", parse=TRUE, size=5,colour=col1_alt) +
    annotate("segment", x = 0, xend = 50, y = start[2], yend = start[2]+v[2]*0.5,
             colour = col1_alt, size = 0.8, arrow = arrow(length = unit(0.3, "cm"))) +
    annotate("segment", x = 520, xend = 600, y = start[2]-0.2, yend = start[2]-0.2, # Add time arrow on x axis
             colour = "black", size = 0.8, arrow = arrow(length = unit(0.3, "cm"))) +
    annotate("text",x=560,y=start[2]-0.35,label="Time", size=4)

  p_upper + p_ddm + p_lower +
    plot_layout(
      ncol = 1,
      nrow = 3,
      widths = c(1),
      heights = c(1,2,1)
  )
}

## Generate a data frame of parameter values and summary stats to display in a table in the web app

gen_table <- function(v,a,b,nsim,res){
  pctA <- as.character(round(nrow(res$upper_bnd[res$upper_bnd$ref==1,])/nsim*100))
  pctB <- as.character(round(nrow(res$lower_bnd[res$lower_bnd$ref==1,])/nsim*100))
  
  pctA0 <- as.character(round(nrow(res$upper_bnd[res$upper_bnd$ref==0,])/nsim*100))
  pctB0 <- as.character(round(nrow(res$lower_bnd[res$lower_bnd$ref==0,])/nsim*100))
  
  medA <- as.character(round(mean(res$upper_bnd[res$upper_bnd$ref==1,1])))
  medB <- as.character(round(mean(res$lower_bnd[res$lower_bnd$ref==1,1])))
  
  medA0 <- as.character(round(mean(res$upper_bnd[res$upper_bnd$ref==0,1])))
  medB0 <- as.character(round(mean(res$lower_bnd[res$lower_bnd$ref==0,1])))
  
  res <- data.frame(ID = c("Reference","User"),
                    v = c(0.7,v),
                    a = c(3,a),
                    b = c(0.5,b),
                    Prob_A = c(pctA0,pctA),
                    Med_RT_A = c(medA0,medA),
                    Prob_B = c(pctB0,pctB),
                    Med_RT_B = c(medB0,medB))

  return(res)
}
