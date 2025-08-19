
###binary
get_binary <- function( data , group , g1 , g2, permutation , p.adjust.method , threads ){
  #
  message( Sys.time() , ' | ','Performing permutation test.'  )

  #
  new.meta <- data$parameters$meta.data[  , c(  data$parameters$sample  , group   )  ]
  colnames(new.meta) <- c( 'sample'  , 'group'  )
  groups <- dplyr::distinct( new.meta , sample , group ,.keep_all = T   ) %>% arrange( sample, group  )
  rownames(groups) <- NULL

  #
  raw.score <- data$LRscore
  raw.score <- raw.score[  , match(  groups$sample , colnames( raw.score  )    ) ]

  #
  res <- pbmclapply( 1:nrow( raw.score ) , function(x){
    v = raw.score[x,] %>% unlist() %>% as.numeric()
    #
    op = c( g1.mean = 0 , g2.mean  = 0 , log2fc = NA , p = NA )
    if( sum(v) != 0 ){
      #
      g1.index <- which( groups$group  == g1 )
      g2.index <- which( groups$group  == g2 )
      #
      g1.mean = mean(  v[  g1.index  ] )
      g2.mean = mean(  v[  g2.index  ] )
      log2fc <- log2( g1.mean / g2.mean  )
      #
      ms <- g1.mean - g2.mean
      #
      g1.length = length( which( groups$group  == g1 ) )
      g2.length = length( which( groups$group  == g2 ) )

      #
      ps <- lapply( 1:permutation , function(num){
        set.seed(num)
        perm.data <- sample( v )
        mean(  perm.data[  g1.index  ]  ) - mean( perm.data[  g2.index  ] )

      }) %>% as.numeric()
      #
      pvalue <- (length(which(abs(ps) >= abs(ms))) + 1) / (  permutation + 1 )
      #
      op = c( g1.mean = g1.mean , g2.mean  = g2.mean , log2fc = log2fc , p = pvalue )
    }
    return( data.table( t( op )  )    )

  } ,mc.cores = threads ) %>% rbindlist() %>% as.data.frame()

  colnames( res ) <- c(  paste( 'mean',c(g1, g2), sep='.' ) , 'log2FC' , 'p'  )
  rownames(res ) <- rownames( raw.score )
  res$p.adj <- p.adjust(  res$p, method = p.adjust.method )
  res <- cbind(  CCC.ID = rownames(res)    , res  )
  #
  colnames( groups ) <- c(  data$parameters$sample  , group   )
  return( list(
    meta.data = groups,
    result = res ,
    type = 'binary',
    parameters = list( data = data , group = group , g1  = g1 , g2 = g2,
                       permutation  = permutation , p.adjust.method = p.adjust.method ,
                       threads = threads  )
  ) )
}


###anova
get_anova <-  function( data , group , p.adjust.method, threads ){
  #
  message( Sys.time() , ' | ','Performing ANOVA test.'  )

  #
  new.meta <- data$parameters$meta.data[  , c(  data$parameters$sample  , group   )  ]
  colnames(new.meta) <- c( 'sample'  , 'group'  )
  groups <- dplyr::distinct( new.meta , sample , group ,.keep_all = T   ) %>% arrange( sample, group  )
  rownames(groups) <- NULL
  gs = unique( groups$group  ) %>% sort()

  #
  raw.score <- data$LRscore
  raw.score <- raw.score[  , match(  groups$sample , colnames( raw.score  )    ) ]

  #
  res <- pbmclapply( 1:nrow( raw.score ) , function(x){
    v = raw.score[x,] %>% unlist() %>% as.numeric()
    #
    op = c( F.value = NA , p =  NA )
    op  = c( op , rep(NA,length(gs)  ) )

    if( sum(v) != 0 ){
      #
      df <- data.frame( d = v , g = groups$group )
      aov_res <- aov( d~g,data = df )
      aov_res <- summary( aov_res )[[1]] %>% as.data.frame()
      op = c( F.value = aov_res[1,'F value'] , p = aov_res[1,'Pr(>F)']  )
      #
      d <- lapply( gs, function(m){  mean(   v[ which( groups$group == m  )  ]   )   }  ) %>% as.numeric()
      op <- c( op , d )
      #
    }
    return( data.table( t( op )  )    )

  } ,mc.cores = threads ) %>% rbindlist() %>% as.data.frame()

  #
  colnames( res ) <- c(  'F.value' , 'p' , paste( 'mean', gs, sep='.' )  )
  rownames(res ) <- rownames( raw.score )
  res$p.adj <- p.adjust(  res$p, method = p.adjust.method )
  res <- dplyr::select( res ,F.value , p , p.adj  , everything() )
  res <- cbind(  CCC.ID = rownames(res)    , res  )
  #
  colnames( groups ) <- c(  data$parameters$sample  , group   )
  return( list(
    meta.data = groups,
    result = res ,
    type = 'anova',
    parameters = list( data = data , group = group  ,
                       p.adjust.method = p.adjust.method, threads = threads )
  ) )
}


###glm
get_glm <- function(  data , group , covariance , p.adjust.method , threads ){
  #
  message( Sys.time() , ' | ','Fitting generalized linear model.'  )

  #
  new.meta <- data$parameters$meta.data[  , c(  data$parameters$sample  , c(group , covariance )  )  ]
  groups <- dplyr::distinct( new.meta , .keep_all = T   )
  groups[[ group ]] <- as.numeric( as.character( groups[[ group ]] ) )
  rownames(groups) <- NULL

  #
  raw.score <- data$LRscore
  raw.score <- raw.score[  , match(  groups[[ data$parameters$sample ]] , colnames( raw.score  )    ) ]

  #
  model_res <- mclapply( 1:nrow( raw.score ) , function(x){
    v = raw.score[x,] %>% unlist() %>% as.numeric()
    #
    if( sum(v) != 0 ){
      #
      df <- groups
      df <- cbind( lrscore = v , df  )
      if( is.null( covariance  ) ){
        myformula <- paste0( 'lrscore' ,'~' , colnames(df)[3]   )
      }else{
        myformula <- paste0( 'lrscore' ,'~' , colnames(df)[3] ,'+' , paste( covariance , collapse = '+' )   )
      }
      model_glm <- glm(  as.formula( myformula  ), data = df , family = gaussian  )
      return( list( name = rownames(raw.score)[x] ,  model =  model_glm )  )
    }
  } ,mc.cores = threads )
  model_res <- model_res[ ! unlist( lapply( model_res , is.null   ) )    ]

  #
  names <- lapply( model_res, function(x) x[['name']] ) %>% unlist() %>% as.character()
  models <- lapply( model_res, function(x) x[['model']] )
  names(models) <- names

  #
  res <- pbmclapply( 1 : nrow( raw.score )  , function(x){
    name = row.names( raw.score )[x]
    op = c( coef = NA , p =  NA )
    if ( name %in% names(models)  ){
      model <- models[[name]]
      ds <- summary(  model )$coefficients  %>% as.data.frame()
      op = c( coef = ds[ group , 'Estimate' ],
              p = ds[ group , 'Pr(>|t|)' ]
      )
    }
    return(   data.table( t(op) )   )
  } , mc.cores = threads ) %>% rbindlist() %>% as.data.frame()

  colnames( res ) <- c(  'coef' , 'p' )
  rownames(res ) <- rownames( raw.score )
  res$p.adj <- p.adjust(  res$p, method = p.adjust.method )
  res <- cbind(  CCC.ID = rownames(res)    , res  )
  #
  return( list(
    meta.data = groups,
    result = res ,
    model = models ,
    type = 'glm',
    parameters = list( data = data, group = group, covariance = covariance ,
                       p.adjust.method = p.adjust.method, threads = threads
    )
  ) )
}


###time
get_time <- function(  data , time , replicate , covariance , p.adjust.method , threads ){
  #
  message( Sys.time() , ' | ','Fitting linear mixed-effects model.'  )

  #
  new.meta <- data$parameters$meta.data[  , c(  data$parameters$sample  , c(time , replicate , covariance )  )  ]
  groups <- dplyr::distinct( new.meta , .keep_all = T   )
  groups[[ time ]] <- as.numeric( as.character( groups[[ time ]] ) )
  groups[[ replicate ]] <- as.numeric( as.character( groups[[ replicate ]] ) )
  rownames(groups) <- NULL

  #
  raw.score <- data$LRscore
  raw.score <- raw.score[  , match(  groups[[ data$parameters$sample ]] , colnames( raw.score  )    ) ]

  #
  model_res <- mclapply( 1:nrow( raw.score ) , function(x){
    v = raw.score[x,] %>% unlist() %>% as.numeric()
    #
    if( sum(v) != 0 ){
      #
      df <- groups
      df <- cbind( lrscore = v , df  )
      if( is.null( covariance  ) ){
        myformula <- paste0( 'lrscore' ,'~' , colnames(df)[3] , '+ ( 1 | ' , colnames(df)[4]  ,')'    )
      }else{
        myformula <- paste0( 'lrscore' ,'~' , colnames(df)[3] , '+ ( 1 | ' , colnames(df)[4]  ,')',
                             '+' , paste( covariance , collapse = '+' )   )
      }
      model_lmer <- lmerTest::lmer(  as.formula( myformula  ), data = df  )
      return( list( name = rownames(raw.score)[x] ,  model =  model_lmer )  )
    }
  } ,mc.cores = threads )
  model_res <- model_res[ ! unlist( lapply( model_res , is.null   ) )    ]

  #
  names <- lapply( model_res, function(x) x[['name']] ) %>% unlist() %>% as.character()
  models <- lapply( model_res, function(x) x[['model']] )
  names(models) <- names

  #
  res <- pbmclapply( 1 : nrow( raw.score )  , function(x){
    name = row.names( raw.score )[x]
    op = c( coef = NA , p =  NA )
    if ( name %in% names(models)  ){
      model <- models[[name]]
      ds <- summary(  model )$coefficients  %>% as.data.frame()
      op = c( coef = ds[ time , 'Estimate' ],
              p = ds[ time , 'Pr(>|t|)' ]
      )
    }
    return(   data.table( t(op) )   )
  } , mc.cores = threads ) %>% rbindlist() %>% as.data.frame()

  colnames( res ) <- c(  'coef' , 'p' )
  rownames(res ) <- rownames( raw.score )
  res$p.adj <- p.adjust(  res$p, method = p.adjust.method )
  res <- cbind(  CCC.ID = rownames(res)  , res )
  #
  return( list(
    meta.data = groups,
    result = res ,
    model = models ,
    type = 'time.course',
    parameters = list( data = data , time = time , replicate = replicate , covariance = covariance ,
                       p.adjust.method  =  p.adjust.method , threads = threads
    )
  ) )
  #
}



###
#binary.params,列表对象，list( group = '' , g1 = '' , g2 = ''  )
#anova.column,列名
#glm.column，列名
#time.course.params，list( time  = '' , replicate  = ''  )


#' Detect differential CCC (cell–cell communication) events
#'
#' @param data The object returned by the lr_score function.
#' @param binary.params A list object with the format list(group = '', g1 = '', g2 = ''), where group represents a column name in meta.data, and the fold change (FC) in the results is calculated as g1/g2.
#' @param anova.column A column name in meta.data based on which an ANOVA test is performed.
#' @param glm.column A column name in meta.data used to fit a GLM model.
#' @param time.course.params Provide a list object with the format list(time = '', replicate = ''), where both time and replicate are column names in meta.data of numeric type. time represents the time point, and replicate represents the biological replicate.
#' @param covariance Column names in meta.data representing covariates, used only for GLM modeling and time course analysis.
#' @param p.adjust.method The method for multiple testing correction, as detailed in the p.adjust() function.
#' @param threads The number of threads used for parallel computation.
#'
#' @returns A list object.
#'
#'
#' @examples
#' ccc.binary <- multiCCC( data = LRscore , binary.params = list( group = 'Group' , g1 = 'O' , g2 = 'Y'  ) )
#' ccc.anova <- multiCCC( data = LRscore , anova.column = 'batch' )
#' ccc.glm <- multiCCC( data = LRscore , glm.column = 'weight' )
#' ccc.time <- multiCCC( data = LRscore , time.course.params = list( time  = 'time' , replicate  = 'replicate'  ) )
#'
#'
#' @export
multiCCC <- function( data , binary.params = NULL ,  anova.column = NULL,
                      glm.column = NULL , time.course.params = NULL ,
                      covariance = NULL, p.adjust.method = 'BH', threads = NULL
                    ){

  ###
  suppressMessages({
  library(dplyr)
  library(pbmcapply)
  library(data.table)
  library(stringr)
  library(parallel)
  library(lmerTest)
  })

  ###threads
  if( is.null(threads) ){  threads <- parallel::detectCores()   }

  ###
  oplist <- list( )

  ###
  if( !is.null( binary.params ) ){
    permutation  = 1000
    if( !is.null( binary.params$permutation ) ){  permutation  =  as.integer(binary.params$permutation) }
    oplist[[ 'binary.res' ]] <- get_binary( data = data , group = binary.params$group ,
                                            g1 = binary.params$g1, g2 = binary.params$g2,
                                            permutation  = permutation  , p.adjust.method = p.adjust.method , threads = threads )
  }
  if( !is.null( anova.column ) ){
    oplist[[ 'anova.res' ]] <- get_anova( data = data, group = anova.column,
                                          p.adjust.method = p.adjust.method , threads = threads
    )
  }
  if( !is.null( glm.column ) ){
    oplist[[ 'glm.res' ]] <- get_glm(
      data = data, group = glm.column, covariance = covariance ,
      p.adjust.method = p.adjust.method , threads = threads
    )
  }
  if( !is.null( time.course.params ) ){
    oplist[[ 'time.course.res' ]] <- get_time(
      data = data, time = time.course.params$time , replicate = time.course.params$replicate ,
      covariance = covariance , p.adjust.method = p.adjust.method , threads = threads
    )
  }

  ###
  message( Sys.time() , ' | ','Done.'  )
  return(  oplist )
}




