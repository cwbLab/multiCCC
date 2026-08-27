# (1)%>%
#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`


# (2):=
#' @importFrom data.table :=
#' @export
data.table::`:=`

rbindlist_n <- function(l, n = 10000, max_chunks = NULL, fill = FALSE, use.names = TRUE, idcol = NULL,  threads = NULL ) {
  
  library(parallel)
  library(data.table)
  library(dplyr)
  
  #detect
  stopifnot(is.list(l))
  if (is.null(threads)) { threads = parallel::detectCores() }
  
  ######
  len <- length(l)
  
  #get n
  if( !is.null(max_chunks) ){  n = max( 1, ceiling(len / max_chunks) )    }
  
  #
  if (len <= n) { return( data.table::rbindlist(l, fill = fill, use.names = use.names, idcol = idcol)) }
  
  #
  starts <- seq(1, len, by = n)
  ends <- pmin(starts + n - 1, len)
  
  #
  intermediate_list <- mclapply(seq_along(starts), function(i){
    data.table::rbindlist( l[starts[i]:ends[i]] , fill = fill, use.names = use.names, idcol = idcol)
  } , mc.cores = threads  )
  
  #
  result <- data.table::rbindlist(intermediate_list, fill = fill, use.names = use.names, idcol = idcol)
  
  ######
  
  return(result)
}

split_vec <- function(x, chunk = 3, min_last = 2) {
  n <- length(x)
  idx <- split(seq_len(n), ceiling(seq_len(n) / chunk))
  
  if (length(idx[[length(idx)]]) < min_last) {
    idx[[length(idx) - 1]] <- c(idx[[length(idx) - 1]], idx[[length(idx)]])
    idx <- idx[-length(idx)]
  }
  
  lapply(idx, function(i) x[i])
}


wb.smc <- function(X, FUN, ..., mc.cores = NULL, mem.ratio.max = 0.8 , mem.max = 16 , pb = T ,time = F ){
  start_time <- Sys.time()
  
  #1
  threads = mc.cores
  is_windows <- .Platform$OS.type == "windows"
  total_cores <- parallel::detectCores()
  
  #windows
  if (is_windows){
    if( time ){ message( wb.log_time_title() , wb.log_time_start_end() , 'Tasks: ' , length(X) ,'.'  )  }
    #
    if(pb){
      res <- pbmcapply::pbmclapply(X, FUN, ... ,mc.cores = 1 )
    }else{
      res <- base::lapply(X, FUN, ...)
    }
    #
    end_time <- Sys.time()
    if( time ){ message( wb.log_time_title() ,
                         wb.log_time_start_end( 'c' ) ,
                         wb.log_time_runtime(  t.minor = start_time ,t.major = end_time  ))
    }
    return(res)
  }
  
  #
  if (  is.null(threads) ){
    target_threads <- total_cores - 1
    target_threads <- max(1, min(target_threads, total_cores))
    
    
    #2
    sample_size <- min(5, length(X))
    set.seed(100)
    sample_idx <- sample(seq_along(X), sample_size)
    
    #
    mean_used <- c()
    for(  temp_idx  in  sample_idx  ){
      temp <- bench::bench_memory(
        suppressMessages( suppressWarnings( capture.output(
          test_results <- base::lapply(X = X[temp_idx], FUN = FUN, ...)
        ) ) )
      )[['mem_alloc']]
      temp <- max( as.numeric( temp ) / 1024^2 , 0.1 )
      
      mean_used <- c( mean_used , temp )
    }
    
    #minimum,0.5 MB
    avg_mem_per_task_mb <- max( median( as.numeric( mean_used ) ), 0.4167 ) * 1.2
    
    #3
    get_total_mem_gb <- function(){
      res <- tryCatch({
        ps::ps_system_memory()[['avail']] / 1024^3
      }, error = function(e) NA )
      return(res)
    }
    
    total_mem_gb <- get_total_mem_gb()
    limit_mem_gb <- if(!is.na(total_mem_gb)) total_mem_gb * mem.ratio.max else mem.max
    
    #4
    max_safe_cores <- floor((limit_mem_gb / (avg_mem_per_task_mb / 1024)) * 0.9)
    max_safe_cores <- max(1, max_safe_cores)
    
    
    if ( max_safe_cores*0.9 > target_threads  ){
      
      target_threads <- min( target_threads  , length(X)  )
      #
      if(pb){
        final_results <- pbmcapply::pbmclapply( X = X, FUN = FUN, ..., mc.cores = target_threads )
      }else{
        final_results <- parallel::mclapply( X = X, FUN = FUN, ..., mc.cores = target_threads )
      }
      #
    }else{
      
      chunk_size <- floor( max(2, min(length(X), max_safe_cores ) ) * 0.9 )
      
      myratio <- chunk_size / target_threads
      if( avg_mem_per_task_mb >= 100 & myratio < 100 ){
        target_threads = floor( max( target_threads / 3,  target_threads / 100   ) )
      }
      
      indices <- split(seq_along(X), ceiling(seq_along(X) / chunk_size))
      
      #5
      final_results <- vector("list", length(X))
      total_chunks <- length(indices)
      
      current_cores <- min( target_threads, max_safe_cores  )
      
      #
      progressr::with_progress({
        a <- 1:total_chunks
        mypb<- progressr::progressor(steps = length(a))
        
        for (i in seq_along(indices)){
          #
          curr_idx <- indices[[i]]
          #
          batch_res <- parallel::mclapply(
            X[curr_idx],
            FUN,
            ...,
            mc.cores = min( current_cores  , length(  X[curr_idx]   )  )
          )
          final_results[curr_idx] <- batch_res
          #
          mypb()
        }
        
      },
      handlers = progressr::handlers(  progressr::handler_progress(
        format = "[:bar] :percent | Elapsed: :elapsed | ETA: :eta",
        clear = FALSE
      )),
      enable = pb
      )
      
    }
  }else{
    threads <- max(1, as.integer( threads ) )
    threads <- min( threads  , length(X) )
    #
    if(pb){
      final_results <- pbmcapply::pbmclapply( X = X, FUN = FUN, ..., mc.cores = threads  )
    }else{
      final_results <- parallel::mclapply( X = X, FUN = FUN, ..., mc.cores = threads  )
    }
    #
  }
  
  #
  end_time <- Sys.time()
  if( time ){ message( wb.log_time_title() ,
                       wb.log_time_start_end( 'c' ) ,
                       wb.log_time_runtime(  t.minor = start_time ,t.major = end_time  ))
  }
  #
  return(final_results)
}

w.pblapply <- function( X , FUN, ... , pb = T , time = F , unlist = F ){
  
  #
  start_time <- Sys.time()
  if( time ){ message( w.log_time_title() , w.log_text_coloured( s.c = 's' ) , 'Tasks: ' , length(X) ,'.'  )  }
  
  #
  result <- NULL
  if( pb ){
    #
    total_length <- length(X) + 1
    mypb <- pbmcapply::progressBar(min = 1 , max = total_length , initial = 1, style = "ETA"  )
    
    result <- base::list()
    for (  i in  seq(  2 , total_length , 1 )  ){
      result <- c( result , base::lapply(  X = X[ i - 1 ] , FUN = FUN , ... ) )
      utils::setTxtProgressBar(mypb, i )
    }
    base::close(mypb)
    #
    names(result) <- names(X)
    
  }else{
    result  <- base::lapply(  X = X , FUN = FUN , ... )
    names(result) <- names(X)
  }
  #
  end_time <- Sys.time()
  if( time ){ message( w.log_time_title() ,
                       w.log_text_coloured( s.c = 'c' ) ,
                       w.log_time_runtime(  t.minor = start_time ,t.major = end_time  ))
  }
  #
  if( unlist ){
    return(  base::unlist( result  )   )
  }else{
    return(result)
  }
}


#
get_database <- function( species  , source ){
  lrdata = NULL
  #
  if ( species == 'human'  ){
    if ( source  == 'CCI' ){
      lrdata = get("mutliCCC.human.lr", envir = asNamespace("multiCCC"))
    }else{
      lrdata = liana::select_resource(  source   ) %>% data.frame()
      lrdata = data.frame( l = lrdata[, grepl("source_genesymbol", colnames(lrdata) ) ],
                           r = lrdata[, grepl("target_genesymbol", colnames(lrdata) ) ]  )
    }
  }else if( species == 'mouse' ){
    if ( source  == 'CCI' ){
      lrdata = get("mutliCCC.mouse.lr", envir = asNamespace("multiCCC"))
    }else{
      lrdata = liana::select_resource(c('MouseConsensus')) %>% data.frame()
      lrdata = data.frame( l = lrdata$MouseConsensus.source_genesymbol , r = lrdata$MouseConsensus.target_genesymbol )
    }
  }
  #
  colnames(lrdata) <- c(  'ligand' , 'receptor'  )
  #
  return( lrdata )
  
}

#
cci_lrscore <- function( exp , meta.data , sample , celltype , LR.ref, lr.database , detect_exp ,threads ){
  
  suppressMessages({library(data.table)})
  #
  samples <-  meta.data[[sample]] %>% unique() %>% as.character()
  celltypes <- meta.data[[celltype]] %>% unique() %>% as.character()
  #
  all_group <- lapply( celltypes , function(x){
    temp <- lapply(celltypes, function(y){
      op <- data.table(   source  = x ,
                          target = y,
                          ligand = lr.database$ligand,
                          receptor = lr.database$receptor
                          
      )
      return(op)
    }) %>% rbindlist()
  }  ) %>% rbindlist() %>% data.frame()
  
  all_group$st <- paste( all_group$source, all_group$target , sep = '_'  )
  all_group$st2 <- paste( all_group$source, all_group$target , sep =  rawToChar(as.raw(c(0xE2, 0x86, 0x92)))  )
  all_group$lr <- paste( all_group$ligand, all_group$receptor , sep = '_'  )
  all_group$CCC.ID <- paste( all_group$st , all_group$lr , sep = '.' )
  #
  all_group <- setDT( all_group  )
  all_group_split <- split(all_group, by = "st2")
  raw_group <- rbindlist(  all_group_split )
  #
  score <- lapply( 1:length( all_group_split ) , function( number  ){
    all_group <- all_group_split[[ number ]]
    score <- wb.smc( 1:nrow(all_group),function(x){
      #
      info <- all_group[x,] %>% unlist() %>% as.character()
      
      res  <- lapply(samples, function(sam){
        #detect
        ld <- subset(detect_exp , sample == sam & celltype == info[1] & gene == info[3]  ) %>% pull(reserved)
        rd <- subset(detect_exp , sample == sam & celltype == info[2] & gene == info[4]  ) %>% pull(reserved)
        if(  length( which( c(ld,rd)  == 'Y'  )) != 2  ){
          LRscore = 0
        }else{
          l.ds <- exp[  meta.data[[sample]] == sam & meta.data[[celltype]] == info[1]   ,  info[3]  ] %>% mean()
          r.ds <- exp[  meta.data[[sample]] == sam & meta.data[[celltype]] == info[2]   ,  info[4]  ] %>% mean()
          LRscore = ( l.ds * r.ds ) / ( l.ds + r.ds )
        }
        #
        return( LRscore )
      })
      #
      return( res )
      
    } , mc.cores = threads , pb = F ) %>% transpose() %>% as.data.table()
    score <- data.table::set( score, j = names(score), value = lapply(score, as.numeric) )
    score <- setDF( score )
    return(score)
  }) %>% rbindlist(use.names = F) %>% setDF()
  #
  rownames(score) <- raw_group$CCC.ID
  colnames(score) <- samples
  score[ is.na(score)  ] <- 0
  #
  return(  list( CCC.info = all_group , LRscore = score , LR.ref = LR.ref  )  )
  
}


#
liana_lrscore <- function( exp,meta.data,sample,celltype, lr.database ,LR.species,LR.method,LR.ref,
                           min.cell  , min.prob , threads , liana_threads = NULL  ){
  #
  if( is.null(threads) ){  liana_threads <- parallel::detectCores()-1  }else{  liana_threads <- threads  }
  #
  suppressWarnings(
    suppressMessages({
      sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = t(as.matrix(exp))),
        metadata = meta.data
      )
      sce$celltype <- meta.data[[ celltype ]] %>% as.character()
      sce@assays@data@listData[["logcounts"]] <- sce@assays@data@listData[["counts"]]
      #
      if( LR.species == 'human' ){
        resource <- lr.database
        lr.database <- liana::select_resource( resource  ) %>% as.data.frame()
        lr.database <- data.frame( ligand = lr.database[, grepl("source_genesymbol", colnames(lr.database) ) ],
                                   receptor = lr.database[, grepl("target_genesymbol", colnames(lr.database) ) ]
        )
      }else if( LR.species == 'mouse' ){
        resource <- 'MouseConsensus'
        lr.database <- liana::select_resource( resource  ) %>% as.data.frame()
        lr.database <- data.frame( ligand = lr.database$MouseConsensus.source_genesymbol,
                                   receptor = lr.database$MouseConsensus.target_genesymbol
        )
      }
    }
    )
  )
  #
  samples <-  meta.data[[sample]] %>% unique() %>% as.character()
  celltypes <- meta.data[[celltype]] %>% unique() %>% as.character()
  #
  all_group <- lapply( celltypes , function(x){
    temp <- lapply(celltypes, function(y){
      op <- data.table(   source  = x ,
                          target = y,
                          ligand = lr.database$ligand,
                          receptor = lr.database$receptor
      )
      return(op)
    }) %>% rbindlist()
  }  ) %>% rbindlist() %>% data.frame()
  
  all_group$st <- paste( all_group$source, all_group$target , sep = '_'  )
  all_group$st2 <- paste( all_group$source, all_group$target , sep =  rawToChar(as.raw(c(0xE2, 0x86, 0x92)))  )
  all_group$lr <- paste( all_group$ligand, all_group$receptor , sep = '_'  )
  all_group$CCC.ID <- paste( all_group$st , all_group$lr , sep = '.' )
  #
  LRscore <- lapply(samples, function(x){
    sce_sub <- sce[, meta.data[[sample]] == x ]
    
    
    ###1.SingleCellSignalR
    if (  LR.method == 'SingleCellSignalR'   ){
      iserror <-tryCatch({
        suppressMessages(
          suppressWarnings(
            liana_res <- liana::liana_wrap( sce_sub,
                                            idents_col = "celltype",
                                            method = c( "sca"),
                                            resource = resource,
                                            min_cells = min.cell,
                                            expr_prop = min.prob,
                                            assay = 'logcounts',
                                            assay.type= 'logcounts',
                                            parallelize =T,  workers = liana_threads
            )
          )
        )
      },error = function(e){
        NULL
      }
      )
      #
      if ( "NULL" %in% class( iserror ) ){  stop( simpleError( paste0( "Failed to evaluate the LR score for the -- ", x , " -- sample using the -- ", LR.method , " -- algorithm. Consider reducing the min.cell, min.exp, and min.prob thresholds." )  ) )  }
      
      liana_res <- iserror
      liana_res$st <- paste( liana_res$source, liana_res$target , sep = '_'  )
      liana_res$lr <- paste( liana_res$ligand, liana_res$receptor , sep = '_'  )
      liana_res$CCC.ID <- paste( liana_res$st , liana_res$lr , sep = '.' )
      #
      myres <- wb.smc( all_group$CCC.ID , function(x){
        op = 0
        if(  x %in% liana_res$CCC.ID ){ op = liana_res$LRscore[ which( liana_res$CCC.ID == x )  ]   }
        return( op )
      }, mc.cores = threads  ,pb = F ) %>% unlist() %>% as.numeric()
      
    }
    
    
    ###2.connectome
    if (  LR.method == 'connectome'   ){
      iserror <-tryCatch({
        suppressMessages(
          suppressWarnings(
            liana_res <- liana::liana_wrap( sce_sub,
                                            idents_col = "celltype",
                                            method = c( "connectome") ,
                                            resource = resource,
                                            min_cells = min.cell,
                                            expr_prop = min.prob,
                                            assay = 'logcounts',
                                            assay.type= 'logcounts',
                                            parallelize =T,  workers = liana_threads
            )
          )
        )
      },error = function(e){
        NULL
      }
      )
      #
      if ( "NULL" %in% class( iserror ) ){  stop( simpleError( paste0( "Failed to evaluate the LR score for the -- ", x , " -- sample using the -- ", LR.method , " -- algorithm. Consider reducing the min.cell, min.exp, and min.prob thresholds." )  ) )  }
      
      liana_res <- iserror
      liana_res$st <- paste( liana_res$source, liana_res$target , sep = '_'  )
      liana_res$lr <- paste( liana_res$ligand, liana_res$receptor , sep = '_'  )
      liana_res$CCC.ID <- paste( liana_res$st , liana_res$lr , sep = '.' )
      #
      myres <- wb.smc( all_group$CCC.ID , function(x){
        op = 0
        if(  x %in% liana_res$CCC.ID ){ op = liana_res$weight_sc[ which( liana_res$CCC.ID == x )  ]   }
        return( op )
      }, mc.cores = threads  ,pb = F  ) %>% unlist() %>% as.numeric()
      
    }
    
    
    ###3.iTALK
    if (  LR.method == 'iTALK'   ){
      iserror <-tryCatch({
        suppressMessages(
          suppressWarnings(
            liana_res <- liana::liana_wrap( sce_sub,
                                            idents_col = "celltype",
                                            method = c( "logfc") ,
                                            resource = resource,
                                            min_cells = min.cell,
                                            expr_prop = min.prob,
                                            assay = 'logcounts',
                                            assay.type= 'logcounts',
                                            parallelize =T,  workers = liana_threads
            )
          )
        )
      },error = function(e){
        NULL
      }
      )
      #
      if ( "NULL" %in% class( iserror ) ){  stop( simpleError( paste0( "Failed to evaluate the LR score for the -- ", x , " -- sample using the -- ", LR.method , " -- algorithm. Consider reducing the min.cell, min.exp, and min.prob thresholds." )  ) )  }
      
      liana_res <- iserror
      liana_res$st <- paste( liana_res$source, liana_res$target , sep = '_'  )
      liana_res$lr <- paste( liana_res$ligand, liana_res$receptor , sep = '_'  )
      liana_res$CCC.ID <- paste( liana_res$st , liana_res$lr , sep = '.' )
      #
      myres <- wb.smc(all_group$CCC.ID , function(x){
        op = 0
        if(  x %in% liana_res$CCC.ID ){ op = liana_res$logfc_comb[ which( liana_res$CCC.ID == x )  ]   }
        return( op )
      }, mc.cores = threads  ,pb = F  ) %>% unlist() %>% as.numeric()
      
    }
    
    
    ###4.NATMI
    if (  LR.method == 'NATMI'   ){
      iserror <-tryCatch({
        suppressMessages(
          suppressWarnings(
            liana_res <- liana::liana_wrap( sce_sub,
                                            idents_col = "celltype",
                                            method = c( "natmi") ,
                                            resource = resource,
                                            min_cells = min.cell,
                                            expr_prop = min.prob,
                                            assay = 'logcounts',
                                            assay.type= 'logcounts',
                                            parallelize =T,  workers = liana_threads
            )
          )
        )
      },error = function(e){
        NULL
      }
      )
      #
      if ( "NULL" %in% class( iserror ) ){  stop( simpleError( paste0( "Failed to evaluate the LR score for the -- ", x , " -- sample using the -- ", LR.method , " -- algorithm. Consider reducing the min.cell, min.exp, and min.prob thresholds." )  ) )  }
      
      liana_res <- iserror
      liana_res$st <- paste( liana_res$source, liana_res$target , sep = '_'  )
      liana_res$lr <- paste( liana_res$ligand, liana_res$receptor , sep = '_'  )
      liana_res$CCC.ID <- paste( liana_res$st , liana_res$lr , sep = '.' )
      #
      myres <- wb.smc( all_group$CCC.ID , function(x){
        op = 0
        if(  x %in% liana_res$CCC.ID ){ op = liana_res$prod_weight[ which( liana_res$CCC.ID == x )  ]   }
        return( op )
      }, mc.cores = threads  ,pb = F  ) %>% unlist() %>% as.numeric()
      
    }
    
    
    ###5.CytoTalk
    if (  LR.method == 'CytoTalk'   ){
      iserror <-tryCatch({
        suppressMessages(
          suppressWarnings(
            liana_res <- liana::liana_wrap( sce_sub,
                                            idents_col = "celltype",
                                            method = c( "cytotalk") ,
                                            resource = resource,
                                            min_cells = min.cell,
                                            expr_prop = min.prob,
                                            assay = 'logcounts',
                                            assay.type= 'logcounts',
                                            parallelize = T,  workers = liana_threads
            )
          )
        )
      },error = function(e){
        NULL
      }
      )
      #
      if ( "NULL" %in% class( iserror ) ){  stop( simpleError( paste0( "Failed to evaluate the LR score for the -- ", x , " -- sample using the -- ", LR.method , " -- algorithm. Consider reducing the min.cell, min.exp, and min.prob thresholds." )  ) )  }
      
      liana_res <- iserror
      liana_res$st <- paste( liana_res$source, liana_res$target , sep = '_'  )
      liana_res$lr <- paste( liana_res$ligand, liana_res$receptor , sep = '_'  )
      liana_res$CCC.ID <- paste( liana_res$st , liana_res$lr , sep = '.' )
      #
      myres <- wb.smc( all_group$CCC.ID , function(x){
        op = 0
        if(  x %in% liana_res$CCC.ID ){ op = liana_res$crosstalk_score[ which( liana_res$CCC.ID == x )  ]   }
        return( op )
      }, mc.cores = threads  ,pb = F  ) %>% unlist() %>% as.numeric()
      
    }
    
    
    ###6.cellchat
    if (  LR.method == 'cellchat'   ){
      
      iserror <-tryCatch({
        suppressMessages(
          suppressWarnings(
            #
            liana_res <- liana::liana_wrap( sce_sub,
                                            idents_col = "celltype",
                                            method = c( "call_cellchat") ,
                                            resource = resource,
                                            min_cells = min.cell,
                                            expr_prop = min.prob,
                                            assay = 'logcounts',
                                            assay.type= 'logcounts',
                                            parallelize = F,  workers = liana_threads,
                                            cellchat.params = list(  organism = LR.species ,
                                                                     de_thresh = 2 ,
                                                                     nboot = 500 ,
                                                                     ligand.pvalues = 2 ,
                                                                     receptor.pvalues = 2
                                            )
            )
            #
          )
        )
      },error = function(e){
        NULL
      }
      )
      #
      if ( "NULL" %in% class( iserror ) ){  stop( simpleError( paste0( "Failed to evaluate the LR score for the -- ", x , " -- sample using the -- ", LR.method , " -- algorithm. Consider reducing the min.cell, min.exp, and min.prob thresholds." )  ) )  }
      
      liana_res <- iserror
      liana_res$st <- paste( liana_res$source, liana_res$target , sep = '_'  )
      liana_res$lr <- paste( liana_res$ligand, liana_res$receptor , sep = '_'  )
      liana_res$CCC.ID <- paste( liana_res$st , liana_res$lr , sep = '.' )
      #
      myres <- wb.smc( all_group$CCC.ID , function(x){
        op = 0
        if(  x %in% liana_res$CCC.ID ){ op = liana_res$prob[ which( liana_res$CCC.ID == x )  ]   }
        return( op )
      }, mc.cores = threads ,pb = F   ) %>% unlist() %>% as.numeric()
      
    }
    
    
    ###7.cellphoneDB
    if (  LR.method == 'cellphoneDB'   ){
      
      iserror <-tryCatch({
        suppressMessages(
          suppressWarnings(
            #
            liana_res <- liana::liana_wrap( sce_sub,
                                            idents_col = "celltype",
                                            method = c( "cellphonedb") ,
                                            resource = resource,
                                            min_cells = min.cell,
                                            expr_prop = min.prob,
                                            assay = 'logcounts',
                                            assay.type= 'logcounts',
                                            parallelize = F,  workers = liana_threads,
                                            cellphonedb.params = list(  score_col = 'LRscore' )
            )
            #
          )
        )
      },error = function(e){
        NULL
      }
      )
      #
      if ( "NULL" %in% class( iserror ) ){  stop( simpleError( paste0( "Failed to evaluate the LR score for the -- ", x , " -- sample using the -- ", LR.method , " -- algorithm. Consider reducing the min.cell, min.exp, and min.prob thresholds." )  ) )  }
      
      liana_res <- iserror
      liana_res$st <- paste( liana_res$source, liana_res$target , sep = '_'  )
      liana_res$lr <- paste( liana_res$ligand, liana_res$receptor , sep = '_'  )
      liana_res$CCC.ID <- paste( liana_res$st , liana_res$lr , sep = '.' )
      #
      myres <- wb.smc( all_group$CCC.ID , function(x){
        op = 0
        if(  x %in% liana_res$CCC.ID ){ op = liana_res$lr.mean[ which( liana_res$CCC.ID == x )  ]   }
        return( op )
      }, mc.cores = threads  ,pb = F  ) %>% unlist() %>% as.numeric()
      
    }
    
    ############
    return( myres )
    
  } ) %>%  transpose() %>% as.data.table() %>%  transpose() %>%  setDF()
  
  #########################################################
  colnames(LRscore) <- samples
  rownames(LRscore) <- all_group$CCC.ID
  #
  return(  list( CCC.info = all_group , LRscore = LRscore ,  LR.ref = LR.ref  )  )
  
}

###metacell
sc.metacell.cluster <- function( object , cluster ,
                                 assay = NULL , layer = NULL,
                                 min.cluster = 50 , min.cluster.reserved = T, min.target = 50 , max.target = NULL ){
  #
  object@meta.data[[ cluster ]] <- as.character( object@meta.data[[ cluster ]]   )
  Idents( object ) <- cluster
  
  #min.cluster
  cell_number <- table( object@meta.data[[ cluster ]] ) %>% sort() %>% as.data.frame()
  cell_number$Var1 <- as.character(  cell_number$Var1  )
  cell_number <- cell_number[ cell_number$Freq > 1 ,  ]
  aba_cells <- cell_number$Var1[ cell_number$Freq < min.cluster  ] %>% as.character()
  cell_number <- cell_number[ ! cell_number$Var1 %in% aba_cells , ]
  if( nrow( cell_number ) < 1  ){  stop( "No clusters are retained for downstream analysis." )  }
  
  #
  if( length( aba_cells ) != 0  ){  aba_obj <- subset( object , idents = aba_cells   )   }
  object <- subset( object  , idents = unique( cell_number$Var1 ) )

  if(  nrow(object) < 1000  ){
    fvf <- rownames(object) %>% as.character()
  }else{
    fvf <- Seurat::VariableFeatures( Seurat::FindVariableFeatures(object, nfeatures = 1000 ,  verbose = F) )
  }

  #
  GE <- Seurat::GetAssayData(  object =  object , assay = assay , layer = layer  )
  sv <- RSpectra::svds(  GE[ which( rownames(GE) %in% fvf  )   ,  ] , k = min(  dim( GE )  ,  max( 30 , nrow( cell_number  ) ) ),
                         opts = list( center = TRUE, scale = TRUE ) )
  pc <- sv$v %*% diag(sv$d)
  rownames( pc ) <- colnames(GE)
  #
  target_ratio <- min(cell_number$Freq) /  min.target
  cell_number$target_n <- round( cell_number$Freq / target_ratio )
  cell_number$target_n[  which( cell_number$target_n < min.target )  ] <- min.target
  if( !is.null( max.target ) ){   cell_number$target_n[  which( cell_number$target_n > max.target )  ] <- max.target  }
  
  mtx_res <- lapply( cell_number$Var1  , function(x){
    #
    s.matrix <- GE[ , object@meta.data[[ cluster ]]  == x  ]
    
    tdata <- pc[ object@meta.data[[ cluster ]]  == x  , ]
    tn <- cell_number$target_n[ cell_number$Var1 == x ]
    
    cluster_res <- mbkmeans::mbkmeans(
      x = t(tdata) ,
      clusters = max( 1 , min( tn , nrow( tdata ) - 2 ) ),
      batch_size = min(tn, nrow( tdata ) ),
      max_iters = 100
    )
    #
    id_trans <- data.frame( raw.id = rownames(tdata)  , cell.id  =  cluster_res[["Clusters"]] )
    groups <- factor( cluster_res[["Clusters"]] , levels = sort(unique(cluster_res[["Clusters"]])  ) )
    
    M <- Matrix::sparseMatrix(i = seq_along(groups), j = as.integer(groups), x = 1)
    sum_expr <- s.matrix %*% M
    n_cell <- Matrix::colSums(M)
    metacell <- base::sweep( sum_expr, 2, n_cell, "/"  )
    
    rownames(metacell) <- rownames(s.matrix)
    colnames(metacell) <- paste0( sprintf( "%s.MC." , x ), levels(groups) )
    #
    id_trans$cell.id <- paste0( sprintf( "%s.MC." , x ), id_trans$cell.id  )
    return(  list( meta = metacell , trans = id_trans  ) )
  } )
  names( mtx_res ) <- as.character(cell_number$Var1)
  #
  id_trans <- lapply(mtx_res, function(x) x$trans  ) %>% bind_rows() 
  mtx_res <- lapply( mtx_res , function(x) x$meta  )
  names( mtx_res ) <- as.character(cell_number$Var1)
  #
  meta.data <- lapply(names(mtx_res), function(x){
    return( data.frame( cell.id = colnames(  mtx_res[[x]] )  , cluster = x )   )
  }) %>% bind_rows()
  exp <- do.call( base::cbind, mtx_res )
  
  ######
  exp <- Matrix::Matrix( exp , sparse = TRUE )
  meta.data <- base::data.frame( meta.data  )
  
  if( length( aba_cells ) != 0 &  min.cluster.reserved  ){
    #
    aba_exp <-  Seurat::GetAssayData(  object = aba_obj , assay = assay , layer = layer  )
    aba_meta <- data.frame( cell.id = colnames(aba_exp)  , cluster = as.character(  Idents( aba_obj ) ) )
    aba_trans <- data.frame(  raw.id =  aba_meta$cell.id   , cell.id  = aba_meta$cell.id  )
    #
    exp <- Matrix::cbind2( exp , aba_exp )
    meta.data <- base::rbind( meta.data , aba_meta )
    id_trans <- base::rbind( id_trans , aba_trans )
    #
  }
  #
  meta.data <- meta.data[ match( colnames(exp) , meta.data$cell.id ) , ]
  id_trans <- id_trans[ match( colnames(exp) , id_trans$cell.id ) , ]
  return( list( mtx = exp , meta.data = meta.data ,
                matched.id = id_trans,
                non.metacell = aba_cells , non.metacell.reserved = min.cluster.reserved  )   )
}

sc.metacell.group <- function( object , cluster , group.by ,
                               assay = NULL , layer = NULL,
                               min.cluster = 50 , min.cluster.reserved = T, min.target = 50 , max.target = NULL ){
  #
  meta.data <- data.frame(object@meta.data)
  
  #
  stat.res <- table(  meta.data[[group.by]] ,  meta.data[[cluster]] ) %>% as.data.frame()
  colnames(stat.res) <- c( 'group' , 'cluster' , 'cell.number'  )
  
  aba_df <- stat.res[ stat.res$cell.number < 2, ]
  stat.res <- stat.res[ stat.res$cell.number > 1, ]
  stat.res$group <- as.character( stat.res$group  )
  stat.res$cluster <- as.character( stat.res$cluster  )
  
  #
  if( nrow( stat.res  ) == 0  ){
    mtx = NULL
    meta.data = NULL
    matched.id = NULL
  }else{
    #
    metacell_res <- lapply( unique( stat.res$group ) , function(xg){
      #
      x  <-  stat.res[  stat.res$group == xg , ]
      #
      scells <- meta.data[ which( meta.data[[group.by]] == xg & meta.data[[cluster]] %in% x$cluster ) , ] %>% rownames()
      
      suppressWarnings(
        #
        res <- sc.metacell.cluster( object = subset( object , cells = scells  )  , cluster = cluster ,
                                    assay = assay , layer =  layer,
                                    min.cluster = min.cluster , min.cluster.reserved = min.cluster.reserved, min.target = min.target ,
                                    max.target = max.target
        )
        #
      )
      
      res$meta.data$group <- xg
      res$meta.data$cell.id <- paste(  xg, res$meta.data$cell.id , sep = '.'  )
      
      res$matched.id$group <- xg
      res$matched.id$cell.id <- paste(  xg, res$matched.id$cell.id , sep = '.'  )
      
      colnames(res$mtx) <- paste(  xg, colnames(res$mtx) , sep = '.'  )
      
      return( res )
      #
    })
    #
    meta.data <- lapply( metacell_res , function(x) x$meta.data ) %>% bind_rows()
    matched.id <- lapply( metacell_res , function(x) x$matched.id ) %>% bind_rows()
    mtx <- lapply( metacell_res, function(x) x$mtx  )
    mtx <- do.call( base::cbind, mtx )
    #
    meta.data <- meta.data[ match( colnames(mtx) , meta.data$cell.id ) , ]
    matched.id <- matched.id[ match( colnames(mtx) , matched.id$cell.id ) , ]
  }
  #
  return(  list(  mtx = mtx , meta.data = meta.data , matched.id = matched.id , non.metacell.reserved = min.cluster.reserved
  ) )
  
}


######comm_score

#' Compute LRscore
#'
#' @description
#' Compute cell–cell communication score (LRscore).
#'
#' @param exp A matrix or dgCMatrix object with cells as row names and genes as column names.
#' 
#' Recommended input: a library-size-normalized data matrix in which the total library size is consistent across individual cells, such as CPM or TPM.
#' @param meta.data A data.frame object whose row names are identical to those of the exp parameter.
#' @param sample The column name in meta.data that represents the sample ID.
#' @param celltype The column name in meta.data that represents the cell type.
#' @param LR.species Species. Options are human or mouse.
#' @param LR.source Two options:
#'
#' (1) Ligand–receptor resource (with CCI or Consensus recommended). Available options include: CCI, Consensus, Baccin2019, CellCall, CellChatDB, Cellinker, CellPhoneDB, CellTalkDB, connectomeDB2020, EMBRACE, Guide2Pharma, HPMR, ICELLNET, iTALK, Kirouac2010, LRdb, Ramilowski2015, OmniPath, and MouseConsensus. 
#' 
#' For mouse, only the CCI and MouseConsensus options are available.
#'
#' (2) A data.frame object can be provided, in which the first column contains ligand gene names and the second column contains receptor gene names.
#' @param LR.method Method for calculating LR score. Available options include: SingleCellSignalR, CCI, connectome, iTALK, NATMI, CytoTalk, cellchat, and cellphoneDB.
#' @param min.cell The minimum number of cells expressing the ligand or receptor.
#' @param min.exp The minimum expression level of ligands and receptors included in the analysis.
#' @param min.prob The minimum proportion of cells expressing the ligand and receptor.
#' @param metacell.min.cells Minimum number of cells required for metacell construction. Clusters containing fewer cells than this threshold will not undergo metacell construction and will retain their original values.
#' @param metacell.max.cells Maximum number of cells allowed in the largest cluster after metacell construction.
#' @param threads Number of cores used for parallel computation. By default, the maximum computing resources are used.
#'
#' @returns A list object.
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
#' @export
#'
scoreLR <- function( exp,meta.data,sample,celltype,
                     LR.species = 'human', LR.source = 'Consensus', LR.method = 'SingleCellSignalR',
                     min.cell = 4, min.exp = 0.1, min.prob = 0.3,
                     metacell.min.cells = 20, metacell.max.cells = 100 ,
                     threads = NULL
){
  ###
  if( !identical( rownames(exp), rownames(meta.data) ) ){
    stop( simpleError( 'Mismatch between row names of exp and meta.data.'  ) )
  }
  run.start = Sys.time()
  suppressWarnings(
    suppressMessages({
      library(dplyr)
      library(pbmcapply)
      library(data.table)
      library(stringr)
      library(parallel)
      library(Seurat)
      library(SingleCellExperiment)
      library(Matrix)
    })
  )
  
  ###LR.source
  if(  is.data.frame( LR.source ) ){
    lr.database <- LR.source
    colnames( lr.database ) <- c(  'ligand' , 'receptor'  )
  }else{
    suppressWarnings(
      suppressMessages({
        lr.database <- get_database( species = LR.species , source = LR.source  )
      })
    )
  }
  lr.genes <- unique(  lr.database$ligand , lr.database$receptor  ) %>% unique()
  
  ###exp,meta.data
  exp <- exp[,  colnames(exp) %in% lr.genes   ] %>% Matrix::Matrix( sparse = T )
  meta.data <- data.frame( meta.data )
  samples <-  meta.data[[sample]] %>% unique() %>% as.character()
  celltypes <-  meta.data[[celltype]] %>% unique() %>% as.character()
  message(  paste0( ">>> Samples ( n = ", length(samples) , ' ): ', paste( sort(samples) ,collapse = ', ' )   )  )
  message(  paste0( ">>> Cell types ( n = ", length( celltypes) , ' ): ', paste( sort(celltypes) ,collapse = ', ' )   )  )
  message( ">>> Running..." )
  message( '[ Step 1/2 | ', format(Sys.time(), "%Y-%m-%d %H:%M:%S") , ' ] ','Checking the expression profiles of ligands and receptors.'  )
  
  ###metacell
  exp <- Matrix::t(exp)
  tsc <- Seurat::CreateSeuratObject( counts = exp , min.cells = 0 , min.features = 0   )
  tsc$group <- meta.data[[sample]]  %>% as.character()
  tsc$cell <- meta.data[[celltype]]  %>% as.character()
  meta_sc <- sc.metacell.group( object = tsc , cluster = 'cell' , group.by = 'group', min.cluster = max( 20 , metacell.min.cells ) ,
                                min.cluster.reserved = T , min.target = max( 20 , metacell.min.cells ) , max.target = metacell.max.cells )
  rm(tsc)
  gc()
  #
  exp <- meta_sc[["mtx"]] %>% as.matrix() %>% t()
  new.meta <- meta_sc[["matched.id"]]
  new.meta <- lapply( rownames(exp)  , function(x){
    #
    ss = new.meta$raw.id[ new.meta$cell.id == x ]
    sm = meta.data[ ss[1] ,  ] 
    #
    return( sm )
  }) %>% bind_rows()
  rownames(new.meta) <- rownames(exp)
  colnames(new.meta) <- colnames(meta.data)
  meta.data <- new.meta
  
  ###detect
  detect_exp <- wb.smc( colnames(exp) ,function(gene){
    level1 <-  lapply(samples, function(m){
      level2 <- lapply(celltypes, function(n){
        ds <- exp[  meta.data[[sample]] == m & meta.data[[celltype]] == n   ,  gene  ] %>% as.numeric()
        #
        my.return='N'
        prob = 0
        if( length(ds) > 0 ){
          prob = length(which(ds >= min.exp)) / length(ds)
          if( length(ds) >= min.cell &  prob >= min.prob   ){
            my.return = 'Y'
          }
        }
        #
        return(   list(sample = m , celltype = n , gene = gene , prob = prob , reserved = my.return )    )
        #
        
      }) %>% rbindlist()
      
    }) %>% rbindlist()
    
  }, mc.cores = threads , pb = T ) %>% rbindlist_n( n = 2000 , use.names = F )
  colnames( detect_exp ) <- c( 'sample' , 'celltype', 'gene', 'prob' ,'reserved'   )
  #
  myfilter <- detect_exp[ which(detect_exp$reserved  == 'Y') , ]
  
  #
  meta.data <- meta.data[  which( meta.data[[sample]] %in% myfilter$sample ) , ]
  exp <- exp[  match( rownames(meta.data) , rownames(exp)  ) ,  
               colnames(exp) %in% myfilter$gene   
  ] %>% as.matrix()
  samples <-  meta.data[[sample]] %>% unique() %>% as.character()
  celltypes <-  meta.data[[celltype]] %>% unique() %>% as.character()
  
  ###LRscore
  message( '[ Step 2/2 | ', format(Sys.time(), "%Y-%m-%d %H:%M:%S") , ' ] ','Calculating communication strength score (LRscore).'  )
  all_lrDB <- data.frame( used.DB = c( "Consensus","Baccin2019","CellCall","CellChatDB",
                                       "Cellinker","CellPhoneDB","CellTalkDB","connectomeDB2020",
                                       "EMBRACE","Guide2Pharma","HPMR","ICELLNET","iTALK","Kirouac2010",
                                       "LRdb","Ramilowski2015","OmniPath","MouseConsensus" )
  )
  all_lrDB$liana_DB <- paste( 'liana' , all_lrDB$used.DB , sep = '_'  )
  all_lrDB <- rbind(all_lrDB , c( 'CCI' , 'CCI' )  )
  
  if( length( samples )  >= 3  ){
    ccc.res.split <- split_vec( 1:length( samples ) , chunk = 3, min_last = 2 )
    #
    ccc.res.split.result <- w.pblapply( 1:length(ccc.res.split), function(chunk){
      #
      s_samples <- samples[ ccc.res.split[[chunk]] ]
      #
      smeta <- meta.data[ meta.data[[sample]] %in% s_samples, ]
      sexp <- exp[ rownames(exp) %in% rownames(smeta) ,   ]
      
      #
      if ( LR.method == 'CCI'  ){
        if(  is.data.frame( LR.source ) ){ LR.ref <- 'custom' }else{ LR.ref <- all_lrDB$liana_DB[ all_lrDB$used.DB == LR.source ]   }
        ccc.res <- cci_lrscore( exp = sexp  , meta.data = smeta, sample =  sample , celltype =  celltype, LR.ref = LR.ref,
                                lr.database =lr.database , detect_exp = detect_exp , threads = threads )
      }else{
        #
        if(  is.data.frame( LR.source ) | LR.source == "CCI"  ){ LR.source <- 'Consensus' }
        LR.ref <- all_lrDB$liana_DB[ all_lrDB$used.DB == LR.source ]
        ccc.res <- liana_lrscore( exp = sexp , meta.data = smeta, sample =  sample , celltype =  celltype,
                                  lr.database = LR.source , LR.species = LR.species,  LR.method = LR.method ,LR.ref = LR.ref,
                                  min.cell = 0 , min.prob = 0 , threads = threads )
      }
      #
      return( ccc.res )
    })
    
    #merge
    all_ccc_info <- lapply(ccc.res.split.result, function(x) data.table( x[['CCC.info']] )) %>% rbindlist( use.names = TRUE )
    final_ccc_info <- unique( all_ccc_info , by = "CCC.ID") %>% data.frame( check.names = F , check.rows = F )
    row_vector <- final_ccc_info$CCC.ID
    
    #
    df_list <- lapply(ccc.res.split.result, function(x) x[['LRscore']] ) 
    all_cols <- unique(unlist(lapply(df_list, colnames)))
    result_matrix <- matrix( 0,  nrow = length(row_vector),  ncol = length(all_cols),
                             dimnames = list(row_vector, all_cols) )
    for (i in seq_along(df_list)) {
      #
      curr_df <- as.matrix(df_list[[i]])
      
      # 
      common_rows <- intersect(rownames(curr_df), row_vector)
      common_cols <- intersect(colnames(curr_df), all_cols)
      
      #
      if (length(common_rows) > 0 && length(common_cols) > 0) {
        block <- curr_df[common_rows, common_cols, drop = FALSE]
        storage.mode(block) <- "numeric"
        #
        result_matrix[common_rows, common_cols] <- as.matrix(block)
      }
    }
    #
    final_ccc_info <- final_ccc_info[ match( rownames(result_matrix) , final_ccc_info$CCC.ID ) ,  ]
    ccc.res <- list(
      CCC.info = final_ccc_info,
      LRscore = data.frame( result_matrix ,check.names = F ,check.rows = F  ) ,
      LR.ref = ccc.res.split.result[[1]][['LR.ref']]
    )
    #
  }else{
    if ( LR.method == 'CCI'  ){
      if(  is.data.frame( LR.source ) ){ LR.ref <- 'custom' }else{ LR.ref <- all_lrDB$liana_DB[ all_lrDB$used.DB == LR.source ]   }
      ccc.res <- cci_lrscore( exp = exp , meta.data = meta.data, sample =  sample , celltype =  celltype, LR.ref = LR.ref,
                              lr.database =lr.database , detect_exp = detect_exp , threads = threads )
      
    }else{
      #
      if(  is.data.frame( LR.source ) | LR.source == "CCI"  ){ LR.source <- 'Consensus' }
      LR.ref <- all_lrDB$liana_DB[ all_lrDB$used.DB == LR.source ]
      ccc.res <- liana_lrscore( exp = exp , meta.data = meta.data, sample =  sample , celltype =  celltype,
                                lr.database = LR.source , LR.species = LR.species,  LR.method = LR.method ,LR.ref = LR.ref,
                                min.cell = 0 , min.prob = 0 , threads = threads )
    }
  }
  
  ###output
  run.end = Sys.time()
  message( '[ Done | ', format(Sys.time(), "%Y-%m-%d %H:%M:%S") , ' ] Total runtime: ', hms::as_hms( as.numeric( run.end - run.start, units = "secs") ) ,'.' )
  #
  return( list(
    CCC.info = ccc.res$CCC.info,
    LRscore = ccc.res$LRscore,
    parameters = list(
      exp = exp ,meta.data = meta.data ,sample = sample,celltype = celltype,
      LR.species = LR.species, LR.source =  ccc.res$LR.ref, LR.method = LR.method,
      min.cell = min.cell, min.exp = min.exp, min.prob = min.prob,
      threads = threads
    )
  ) )
}


