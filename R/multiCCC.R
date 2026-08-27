
###binary
get_binary <- function( data , group , g1 , g2,  test, permutation , p.adjust.method , threads ){
  #
  message( '[ ', format(Sys.time(), "%Y-%m-%d %H:%M:%S") , ' ] ',
           'Performing ', test ,' test.'  )
  
  #
  new.meta <- data$parameters$meta.data[  , c(  data$parameters$sample  , group   )  ]
  colnames(new.meta) <- c( 'sample'  , 'group'  )
  groups <- dplyr::distinct( new.meta , sample , group ,.keep_all = T   ) %>% arrange( sample, group  )
  groups$sample <- as.character(groups$sample)
  rownames(groups) <- NULL
  
  #
  raw.score <- data$LRscore
  raw.score <- raw.score[  , match(  groups$sample , colnames( raw.score  )    ) ]
  
  #mean,FC
  data1 = raw.score[  , colnames(raw.score) %in% groups$sample[ groups$group == g1  ]  ]
  data2 = raw.score[  , colnames(raw.score) %in% groups$sample[ groups$group == g2  ]  ]
  
  res = suppressWarnings( matrixTests::row_t_equalvar (x = data1 ,  y = data2    ) )
  res <- data.frame( g1.mean = res$mean.x , g2.mean = res$mean.y , row.names = rownames(res) )
  res$log2fc <- log2(  res$g1.mean / res$g2.mean )
  res$p <- 1
  
  colnames( res ) <- c(  paste( 'mean',c(g1, g2), sep='.' ) , 'log2FC' , 'p'  )
  res <- res[ res$log2FC != 'NaN' , ]
  if( nrow(res) < 1 ){ stop( simpleError(  'All input values are zero. No valid data are available for further analysis. Please check the previous analysis step.'  ) )   }
  res$log2FC <- as.numeric( res$log2FC )
  
  temp.raw.score <- raw.score[  rownames(raw.score)  %in% rownames(res)  ,  ]
  temp.raw.score <- temp.raw.score[   match( rownames( res) , rownames(temp.raw.score)   ) , groups$sample  ]
  
  #permutation, U, T
  if( test == 'permutation'  ){
    idx1 <- which( groups$group  == g1 )
    idx2 <- which( groups$group  == g2 )
    n1   <- length(idx1)
    n2   <- length(idx2)
    all_idx <- c(idx1, idx2)
    
    # 
    obs_stat <- rowMeans(temp.raw.score[, idx1, drop = FALSE]) - rowMeans(temp.raw.score[, idx2, drop = FALSE])
    
    #
    perm_stats <- wb.smc( 1:permutation,function(x) {
      set.seed(x)
      shuffled <- sample(all_idx)
      a <- rowMeans(temp.raw.score[, shuffled[seq_len(n1)], drop = FALSE]) - rowMeans(temp.raw.score[, shuffled[(n1 + 1):(n1 + n2)], drop = FALSE])
      return(a)
    }, mc.cores = threads ) %>% data.frame() %>% as.matrix()
    
    pvalues <- rowMeans( abs(perm_stats) >= abs(obs_stat) )
    res$p <- as.numeric( pvalues  )
    #
    
  }else{
    #
    data1 = temp.raw.score[  , colnames(temp.raw.score) %in% groups$sample[ groups$group == g1  ]  ]
    data2 = temp.raw.score[  , colnames(temp.raw.score) %in% groups$sample[ groups$group == g2  ]  ]
    
    #
    if( test == 't' ){
      p.new = suppressWarnings( matrixTests::row_t_equalvar( x = data1 ,  y = data2    ) )
      p.new = p.new[  match( rownames(res) , rownames( p.new  )  )   ,  ]
      res$p <- p.new$pvalue
      
    }else{
      p.new = suppressWarnings( matrixTests::row_wilcoxon_twosample( x = data1 ,  y = data2    ) )
      p.new = p.new[  match( rownames(res) , rownames( p.new  )  )   ,  ]
      res$p <- p.new$pvalue
    }
    #
  }
  
  #padj
  m.info <- data$CCC.info[ data$CCC.info$CCC.ID %in% rownames(res) ,  ]
  m.info <- m.info[  match( rownames(res), m.info$CCC.ID ),  c( "source","target","ligand","receptor","st","lr" )   ]
  
  res$p.adj <- p.adjust(  res$p, method = p.adjust.method )
  res <- cbind(  CCC.ID = rownames(res)    , res  )
  res <- cbind( res , m.info )
  #
  colnames( groups ) <- c(  data$parameters$sample  , group   )
  
  ######LRwFC
  nonzero_means <- c( res[ which(res[,2]>0 ) ,2 ] , res[ which(res[,3]> 0 ) ,3 ] ) %>% unlist() %>% as.numeric() %>% abs()
  pc <- mean(nonzero_means) * ( length( nonzero_means ) / (nrow(res) * 2 ) )
  log2FC.smoothed <- log2( (abs(res[,2]) + pc) / (abs(res[,3]) + pc) )
  max.mean <- pmax(abs(res[,2]), abs(res[,3]) )
  res$LRwFC <- log2FC.smoothed * max.mean / max( nonzero_means )
  
  #
  return( list(
    meta.data = groups,
    result = res ,
    type = 'binary',
    parameters = list( data = data , group = group , g1  = g1 , g2 = g2, test = test, 
                       permutation  = permutation , p.adjust.method = p.adjust.method ,
                       threads = threads  )
  ) )
}



###anova
get_anova <-  function( data , group , p.adjust.method, threads ){
  #
  message( '[ ', format(Sys.time(), "%Y-%m-%d %H:%M:%S") , ' ] ', 'Performing ANOVA test.'  )
  
  #
  new.meta <- data$parameters$meta.data[  , c(  data$parameters$sample  , group   )  ]
  colnames(new.meta) <- c( 'sample'  , 'group'  )
  groups <- dplyr::distinct( new.meta , sample , group ,.keep_all = T   ) %>% arrange( sample, group  )
  groups$sample <- as.character(groups$sample)
  rownames(groups) <- NULL
  gs = unique( groups$group  ) %>% sort()
  
  #
  raw.score <- data$LRscore
  raw.score <- raw.score[  , match(  groups$sample , colnames( raw.score  )    ) ]
  
  #
  res <- wb.smc( 1:nrow( raw.score ) , function(x){
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
    #
    if( is.na(op[['p']]) ){ return(NULL) }else{ return( op ) }
    
  } ,mc.cores = threads )
  #
  null_res <- vapply(res, is.null, logical(1))
  res <- res[  !null_res  ] %>% transpose() %>% as.data.table() %>% setDF()
  
  #
  colnames( res ) <- c(  'F.value' , 'p' , paste( 'mean', gs, sep='.' )  )
  res <- res[ !is.na(res$p) ,  ]
  if( nrow(res) < 1 ){ stop( simpleError(  'All input values are zero. No valid data are available for further analysis. Please check the previous analysis step.'  ) )   }
  rownames(res ) <- rownames( raw.score )[  !null_res  ]
  res$p.adj <- p.adjust(  res$p, method = p.adjust.method )
  res <- dplyr::select( res ,F.value , p , p.adj  , everything() )
  res <- cbind(  CCC.ID = rownames(res)    , res  )
  res <- cbind( res , data$CCC.info[ !null_res , c( "source","target","ligand","receptor","st","lr" )  ]  )
  #
  colnames( groups ) <- c(  data$parameters$sample  , group   )
  
  #
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
  message( '[ ', format(Sys.time(), "%Y-%m-%d %H:%M:%S") , ' ] ', 'Fitting generalized linear model.'  )
  
  #
  new.meta <- data$parameters$meta.data[  , c(  data$parameters$sample  , c(group , covariance )  )  ]
  groups <- dplyr::distinct( new.meta , .keep_all = T   )
  groups[[ group ]] <- as.numeric( as.character( groups[[ group ]] ) )
  rownames(groups) <- NULL
  
  #
  raw.score <- data$LRscore
  raw.score <- raw.score[  , match(  groups[[ data$parameters$sample ]] , colnames( raw.score  )    ) ]
  
  #
  suppressWarnings(
    suppressMessages({
      model_res <- lapply( 1:nrow( raw.score ) , function(x){
        v = raw.score[x,] %>% unlist() %>% as.numeric()
        #
        if( max(v) != 0 ){
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
      })
    })
  )
  model_res <- model_res[ ! unlist( lapply( model_res , is.null   ) )    ]
  
  #
  names <- lapply( model_res, function(x) x[['name']] ) %>% unlist() %>% as.character()
  models <- lapply( model_res, function(x) x[['model']] )
  names(models) <- names
  
  #
  res <- wb.smc( 1 : nrow( raw.score )  , function(x){
    name = row.names( raw.score )[x]
    op = c( coef = NA , p =  NA )
    if ( name %in% names(models)  ){
      model <- models[[name]]
      ds <- summary(  model )$coefficients  %>% as.data.frame()
      op = c( coef = ds[ group , 'Estimate' ],
              p = ds[ group , 'Pr(>|t|)' ]
      )
    }
    #
    if( is.na(op[['p']]) ){ return(NULL) }else{ return( op ) }
    
  } , mc.cores = threads )
  #
  null_res <- vapply(res, is.null, logical(1))
  res <- res[  !null_res  ] %>% transpose() %>% as.data.table() %>% setDF()
  
  colnames( res ) <- c(  'coef' , 'p' )
  res <- res[ !is.na(res$p) ,  ]
  if( nrow(res) < 1 ){ stop( simpleError(  'All input values are zero. No valid data are available for further analysis. Please check the previous analysis step.'  ) ) }
  rownames(res ) <- rownames( raw.score )[  !null_res  ]
  res$p.adj <- p.adjust(  res$p, method = p.adjust.method )
  res <- cbind(  CCC.ID = rownames(res)    , res  )
  res <- cbind( res , data$CCC.info[ !null_res , c( "source","target","ligand","receptor","st","lr" )  ]  )
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
  message( '[ ', format(Sys.time(), "%Y-%m-%d %H:%M:%S") , ' ] ', 'Fitting linear mixed-effects model.'  )
  
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
  suppressWarnings(
    suppressMessages({
      model_res <- lapply( 1:nrow( raw.score ) , function(x){
        v = raw.score[x,] %>% unlist() %>% as.numeric()
        #
        if( max(v) != 0 ){
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
      })
  
    })
  )
  
  model_res <- model_res[ ! unlist( lapply( model_res , is.null   ) )    ]
  
  #
  names <- lapply( model_res, function(x) x[['name']] ) %>% unlist() %>% as.character()
  models <- lapply( model_res, function(x) x[['model']] )
  names(models) <- names
  
  #
  res <- wb.smc( 1 : nrow( raw.score )  , function(x){
    name = row.names( raw.score )[x]
    op = c( coef = NA , p =  NA )
    if ( name %in% names(models)  ){
      model <- models[[name]]
      ds <- summary(  model )$coefficients  %>% as.data.frame()
      op = c( coef = ds[ time , 'Estimate' ],
              p = ds[ time , 'Pr(>|t|)' ]
      )
    }
    #
    if( is.na(op[['p']]) ){ return(NULL) }else{ return( op ) }
    
  } , mc.cores = threads )
  #
  null_res <- vapply(res, is.null, logical(1))
  res <- res[  !null_res  ] %>% transpose() %>% as.data.table() %>% setDF()
  
  colnames( res ) <- c(  'coef' , 'p' )
  res <- res[ !is.na(res$p) ,  ]
  if( nrow(res) < 1 ){ stop( simpleError(  'All input values are zero. No valid data are available for further analysis. Please check the previous analysis step.'  ) ) }
  rownames(res ) <- rownames( raw.score )[  !null_res  ]
  res$p.adj <- p.adjust(  res$p, method = p.adjust.method )
  res <- cbind(  CCC.ID = rownames(res)  , res )
  res <- cbind( res , data$CCC.info[ !null_res , c( "source","target","ligand","receptor","st","lr" )  ]  )
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

#' Detect differential CCC events
#'
#' @description
#' Detect differential CCC (cell–cell communication) events.
#'
#' @param data The object returned by the scoreLR function.
#' @param binary.params 
#' A list defining the comparison. It should contain at least `group`, `g1`, and `g2`, where `group` specifies a column in `meta.data`, and the fold change (FC) is calculated as `g1 / g2`.
#' The optional element `test` specifies the statistical test:
#' 
#' \itemize{
#'   \item Student's t-test (default): list(group = '', g1 = '', g2 = '', test = 't' )
#'   \item Wilcoxon rank-sum test: list(group = '', g1 = '', g2 = '', test = 'wilcoxon' )
#'   \item Permutation test: list(group = '', g1 = '', g2 = '', test = 'permutation' )
#' }
#' 
#' @param anova.column A column name in meta.data based on which an ANOVA test is performed.
#' @param glm.column A column name in meta.data used to fit a GLM model.
#' @param time.course.params Provide a list object with the format list(time = '', replicate = ''), where both time and replicate are column names in meta.data of numeric type. time represents the time point, and replicate represents the biological replicate.
#' @param covariance Column names in meta.data representing covariates, used only for GLM modeling and time course analysis.
#' @param p.adjust.method The method for multiple testing correction, as detailed in the p.adjust() function.
#' @param threads Number of cores used for parallel computation. By default, the maximum computing resources are used.
#'
#' @returns
#' A list object.
#'
#' Note: The multiCCC function supports the simultaneous use of the binary.params, anova.column, glm.column, and time.course.params parameters, and analyses will be performed automatically for each application scenario.
#'
#'
#' @examples
#'
#' data( "CCC.test.data" )
#'
#' #1.scoreLR
#' LRscore <- scoreLR( exp = t( CCC.test.data$exp ) ,
#'                  meta.data = CCC.test.data$meta ,
#'                  LR.species = 'human' ,
#'                  sample = 'orig.ident' ,
#'                  celltype = 'celltype' )
#'
#' #2.multiCCC
#' #Note: The multiCCC function supports the simultaneous use of the binary.params, anova.column, glm.column, and time.course.params parameters, and analyses will be performed automatically for each application scenario.
#' ccc.binary <- multiCCC( data = LRscore , binary.params = list( group = 'Group' , g1 = 'O' , g2 = 'Y'  ) )
#' ccc.anova <- multiCCC( data = LRscore , anova.column = 'batch' )
#' ccc.glm <- multiCCC( data = LRscore , glm.column = 'weight' )
#' ccc.time <- multiCCC( data = LRscore , time.course.params = list( time  = 'time' , replicate  = 'replicate'  ) )
#' 
#' #One-step execution for all scenarios
#' ccc.all <- multiCCC( data = LRscore , 
#' 	   binary.params = list( group = 'Group' , g1 = 'O' , g2 = 'Y'  ),
#' 	   anova.column = 'batch',
#' 	   glm.column = 'weight',
#' 	   time.course.params = list( time  = 'time' , replicate  = 'replicate'  )
#' 	)
#'
#' @export
multiCCC <- function( data , binary.params = NULL ,  anova.column = NULL,
                      glm.column = NULL , time.course.params = NULL ,
                      covariance = NULL, p.adjust.method = 'BH', threads = NULL
){
  
  ###
  if(  is.null( binary.params ) &  is.null( anova.column ) & is.null( glm.column ) & is.null( time.course.params )  ){
    stop( simpleError(  '`binary.params`, `anova.column`, `glm.column`, and `time.course.params` cannot all be `NULL`.'  ) )
  }
  
  #
  run.start = Sys.time()
  suppressMessages({
    library(dplyr)
    library(pbmcapply)
    library(data.table)
    library(stringr)
    library(parallel)
    library(lmerTest)
  })

  
  ###
  oplist <- list( )
  
  ###
  if( !is.null( binary.params ) ){
    test = 't'
    permutation  = 1000
    if( !is.null( binary.params$permutation ) ){  permutation  =  as.integer(binary.params$permutation) }
    if( !is.null( binary.params$test )  ){ test = as.character(  binary.params$test    )   }
    oplist[[ 'binary.res' ]] <- get_binary( data = data , group = binary.params$group ,
                                            g1 = binary.params$g1, g2 = binary.params$g2, test = test , 
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
  run.end = Sys.time()
  message( '[ ', format(Sys.time(), "%Y-%m-%d %H:%M:%S") , ' ] ','Done. Total runtime: ', hms::as_hms( as.numeric( run.end - run.start, units = "secs") ) ,'.' )
  return(  oplist )
}




