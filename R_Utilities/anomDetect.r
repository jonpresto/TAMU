anomDetect <- function(X, k = c(10,20), nsub = nrow(X), method = c('lof','ldf','rkof','lpdf'),
                  ldf.param = c(h = 1, c = 0.1),
                  rkof.param = c(alpha = 1, C = 1, sig2 = 1),
                  lpdf.param = c(cov.type = 'full',sigma2 = 1e-5, tmax=1, v=1),
                  treetype='kd',searchtype='standard',eps=0.0,
                  scale.data=TRUE)
{

  if(is.null(k)) stop('k is missing')

  if(!is.numeric(k)) stop('k is not numeric')


  # coerce X to class matrix
  X <- as.matrix(X)

  if(!is.numeric(X)) stop('X contains non-numeric data type')

  if(!is.matrix(X)) stop('X must be of class matrix')

  k <- as.integer(k)

  nsub <- as.integer(nsub)

  # check max k less than than nsub-1
  kmax <-  max(k)
  len.k <- length(k)
  if(kmax > nsub - 1 ){ stop('k is greater than or equal to nsub') }

  if(is.null(nsub)) (nsub <- nrow(X))

  nsub  <- as.integer(nsub)
  kmax  <-  max(k)
  len.k <- length(k)


  if(kmax > nsub -1) stop('k is greater than or equal to nsub')

  if(min(k) < 2 ) stop('k must be greater than 1')





  # number of rows of X
  n <- nrow(X)
  # number of columns of X
  p <- ncol(X)

  # check that nsub <= n
  if(nsub>n){ nsub=n }
  if(nsub<10){ nsub=10 }

  # subsample ids without replacement
  sub_sample_ids = sample(1:n,nsub)


  # scale X
  if(scale.data){  X <- scale(X) }

  # subsample X, Y will subsample of X
  Y <-  as.matrix(X[sub_sample_ids,])



  # compute distance (euclidean) matrix between X and Y and returns kNN ids and kNN distances #
  knn <- RANN::nn2(data = Y,query = X,k = kmax+1,treetype=treetype, searchtype=searchtype, eps=eps)

  # kNN id matrix
  knn_ids <- knn$nn.idx[,-1]

  # kNN distance matrix
  knn_dist_matrix <- knn$nn.dists[,-1]



  # define storeage matrix for each outlier score ................
  if('lof'%in%method){
    store_lrd <- matrix(NA,n,len.k)
    store_lof <- matrix(NA,n,len.k)
  }else{
    store_lrd <- NULL
    store_lof <- NULL
  }

  if('ldf'%in%method){
    store_lde <- matrix(NA,n,len.k)
    store_ldf <- matrix(NA,n,len.k)
  }else{
    store_lde <- NULL
    store_ldf <- NULL
  }

  if('rkof'%in%method){
    store_kde <- matrix(NA,n,len.k)
    store_rkof <- matrix(NA,n,len.k)
  }else{
    store_kde <- NULL
    store_rkof <- NULL
  }

  if('lpdf'%in%method){
    store_lpde <- matrix(NA,n,len.k)
    store_lpdf <- matrix(NA,n,len.k)
    store_lpdr <- matrix(NA,n,len.k)
  }else{
    store_lpde <- NULL
    store_lpdf <- NULL
    store_lpdr <- NULL
  }

  # initialize
  ii <- 0

  ### compute lof and lrd for each k value specified by user
  for(kk in k){

    ii <- ii+1

    # distance to kth nearest neighbor for X to Y
    dist_k <- knn_dist_matrix[,kk]


    ## Reachability distance matrix if outlier score lrd, lof, lde, or ldf has been specified
    if(sum(method%in%c('lof','ldf'))>0){

      # arrange dist_k relative to subsample data set
      dist_k_rel_to_test = t(apply(knn_ids[,1:kk],1,function(x)dist_k[sub_sample_ids[x]]))

      # compute reachability distances between X and kNN of X relative to subsample data set
      # returns reachability distance matrix of size n x kk
      reach_dist_matrix_test <- t(apply( cbind(knn_dist_matrix[,1:kk] , dist_k_rel_to_test),1,function(x){

        x1 = x[1:kk]
        x2 = x[(kk+1):(2*kk)]
        reach_dist = apply(rbind(x1,x2),2,max)
        return(reach_dist)

      }))


      ############  compute outlier scores lrd and lof  ############
      ############  compute outlier scores lrd and lof  ############
      ############  compute outlier scores lrd and lof  ############
      if('lof'%in%method){

        # compute local reachability density for each point in X
        lrd <- 1/(apply(reach_dist_matrix_test,1,mean)+1e-198)

        # compute local outlier factor for each point in X
        lof <- apply(knn_ids[,1:kk],1,function(x)mean(lrd[sub_sample_ids[x]]))/lrd

        # store lof and lrd for each k
        store_lrd[,ii] <- lrd
        store_lof[,ii] <- lof

      }# end if statement for lof

      ############  compute outlier scores lde and ldf  ############
      ############  compute outlier scores lde and ldf  ############
      ############  compute outlier scores lde and ldf  ############
      if('ldf'%in%method){

        # paramters
        h <- ldf.param[names(ldf.param)=='h']
        c <- ldf.param[names(ldf.param)=='c']

        lde <- sapply(1:n,function(id){
                  mean(1/((2*pi)^(p/2))*1/(h*dist_k[sub_sample_ids[knn_ids[id,1:kk]]])^p*
                  exp(-(.5*reach_dist_matrix_test[id,]^2)/(h*dist_k[sub_sample_ids[knn_ids[id,1:kk]]])^2))+1e-198
        })

        ldf <- sapply(1:n,function(id){
                  mean.lde = mean(lde[sub_sample_ids[knn_ids[id,1:kk]]])
                  ldf = mean.lde/(lde[id]+c*mean.lde)
                  return(ldf)
        })

        # store lof and lrd for each k
        store_lde[,ii] <- lde
        store_ldf[,ii] <- ldf

      }# end if statement for ldf


    }# end if statement for reachability distance calculations




    ############  compute outlier scores kde and rkof  ############
    ############  compute outlier scores kde and rkof  ############
    ############  compute outlier scores kde and rkof  ############
    if('rkof'%in%method){

      # parameters
      alpha <- rkof.param[names(rkof.param)=='alpha']
      C     <- rkof.param[names(rkof.param)=='C']
      sig2  <- rkof.param[names(rkof.param)=='sig2']

      ## compute kde for each point in X
      kde <-  sapply(1:n, function(id){

        mean(1/(2*pi)^(p/2)*1/(C*dist_k[sub_sample_ids[knn_ids[id,1:kk]]]^alpha)^2*
               exp(-.5*knn_dist_matrix[id,1:kk]^2/(C*dist_k[sub_sample_ids[knn_ids[id,1:kk]]]^alpha)^2))+1e-198

      })

      ## compute wde for each point in X
      wde <-  sapply(1:n,function(id){

        weights = exp(-(dist_k[sub_sample_ids[knn_ids[id,1:kk]]]/min(dist_k[sub_sample_ids[knn_ids[id,1:kk]]])-1)^2/(2*sig2))

        weights = weights/sum(weights)
        wde = sum(weights*kde[sub_sample_ids[knn_ids[id,1:kk]]])

      })


      ## compute rkof
      rkof <- wde/kde

      # store kde and rkof for each k
      store_kde[,ii] <- kde
      store_rkof[,ii] <- rkof

    }# end if statement for rkof




    ############  compute outlier scores lpde, lpdf, and lpdr  ############
    ###########  compute outlier scores lpde, lpdf, and lpdr  ############
    ###########  compute outlier scores lpde, lpdf, and lpdr  ############
    if('lpdf'%in%method){

      # parameters
      cov.type<- lpdf.param[names(lpdf.param)=='cov.type']
      tmax   <- as.numeric(lpdf.param[names(lpdf.param)=='tmax'])
      sigma2 <- as.numeric(lpdf.param[names(lpdf.param)=='sigma2'])
      v      <- as.numeric(lpdf.param[names(lpdf.param)=='v'])

      # max iterations
      tmax <- tmax+1
      # identity matrix
      II <- diag(1,p,p)
      # regularization matrix
      reg <- sigma2*II

      # matrices to store densities and weights
      store_dens <- matrix(NA,n,tmax)
      store_R <- matrix(NA,n,tmax)

      # compute density and weight function, R, for each iteration in 1:tmax
      for(t in 1:tmax){

        #compute multivarite t density for each point in X referencing subsample data set Y
        dens <-  sapply(1:n,function(id){

          # test point
          x = X[id,]

          # k nearest neighborhood to compute weighted location and scatter
          hood = as.matrix(Y[knn_ids[id,1:kk],])

          # define weights
          if(t==1){weights = rep(1/kk,kk)}
          if(t>1){weights  = R[sub_sample_ids[knn_ids[id,1:kk]]] }

          # normalize weights to sum to 1
          weights <-  weights/sum(weights)

          # comptue weighted sample mean and sample covariance matrix  if cov.type='full'
          if(cov.type=='full'){
            covwt   = cov.wt(hood,wt=weights,method='ML')
            center  = covwt$center
            scatter = covwt$cov+reg

          }

          # comptue weighted sample mean and sample covariance matrix if cov.type='diag'
          if(cov.type=='diag'){
            center     = colSums(weights * hood)
            center.mat = matrix(center,kk,p,byrow=T)
            vars       = colSums(weights*(hood-center.mat)^2)
            scatter    = diag(vars,p,p)+reg
          }


          # compute multivaritae t density with degrees of freedom v
          density = mnormt::dmt(x,mean=center,S=scatter,df=v)+1e-198

        })


        # compute updated weight function, R, for each point in X
        R <- sapply(1:n,function(id){

          # compute weights
          weights = 1/(dist_k[sub_sample_ids[knn_ids[id,1:kk]]]^2+1e-20)
          weights = weights/sum(weights)

          # comptue ratio of density to weighted neighborhood density
          log1p(dens[id])/sum(weights*log1p(dens[sub_sample_ids[knn_ids[id,1:kk]]]))

        })


        # store all densities and R
        store_dens[,t] <- dens
        store_R[,t] <- R

        ## update weight rule
        if(t>2){

          R <- apply(store_R[,2:t],1,max)
        }


      }### end t loop

      # compute outlier scores lpde, lpdf, lpdr
      lpde <- apply(store_dens,1,function(x)log(mean(x[-1])))
      lpdf <- apply(store_R,1,function(x)mean(x[-1]))
      lpdr <- lpde - log(store_dens[,1])


      store_lpde[,ii] <- lpde
      store_lpdf[,ii] <- lpdf
      store_lpdr[,ii] <- lpdr


    }# end if statement for lpdf



  }# k loop

  # return all outlier scores
  return(list(lrd = store_lrd,
              lof = store_lof,
              lde = store_lde,
              ldf = store_ldf,
              kde = store_kde,
              rkof = store_rkof,
              lpde = store_lpde,
              lpdf = store_lpdf,
              lpdr = store_lpdr,
			   kid = knn_ids,
              kdm = knn_dist_matrix
  )
  )



} # end anomDetect function
