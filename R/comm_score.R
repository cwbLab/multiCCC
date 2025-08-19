

#
get_database <- function( species  , source ){
  lrdata = NULL
  #
  if ( species == 'human'  ){
    if ( source  == 'CCI' ){
      lrdata = get("mutliCCC.human.lr", envir = asNamespace("multiCCC"))
    }else{
      lrdata = liana::select_resource(c('Consensus')) %>% data.frame()
      lrdata = data.frame( l = lrdata$Consensus.source_genesymbol, r = lrdata$Consensus.target_genesymbol   )
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
cci_lrscore <- function( exp , meta.data , sample , celltype , lr.database , detect_exp ,threads ){
  samples <-  meta.data[[sample]] %>% unique() %>% as.character()
  celltypes <- meta.data[[celltype]] %>% unique() %>% as.character()
  #
  all_group <- mclapply( celltypes , function(x){
    temp <- lapply(celltypes, function(y){
      op <- data.table(   source  = x ,
                          target = y,
                          ligand = lr.database$ligand,
                          receptor = lr.database$receptor

      )
      return(op)
    }) %>% rbindlist()
  },mc.cores = threads   ) %>% rbindlist() %>% data.frame()

  all_group$st <- paste( all_group$source, all_group$target , sep = '_'  )
  all_group$st2 <- paste( all_group$source, all_group$target , sep =  rawToChar(as.raw(c(0xE2, 0x86, 0x92)))  )
  all_group$lr <- paste( all_group$ligand, all_group$receptor , sep = '_'  )
  all_group$CCC.ID <- paste( all_group$st , all_group$lr , sep = '.' )
  #
  score <- pbmclapply( 1:nrow(all_group),function(x){
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
    return( data.table( matrix( res , nrow = 1 ) ) )

  } , mc.cores = threads ) %>% rbindlist()
  score <- data.frame( score  )

  score <- apply( score , 2, as.numeric) %>% as.data.frame()
  rownames(score) <- all_group$CCC.ID
  colnames(score) <- samples
  score[ is.na(score)  ] <- 0
  #
  return(  list( CCC.info = all_group , LRscore = score   )  )

}


#
liana_lrscore <- function( exp,meta.data,sample,celltype, lr.database ,LR.species,
                           min.cell  , min.prob , threads ){
  #
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = t(as.matrix(exp))),
    metadata = meta.data
  )
  sce$celltype <- meta.data[[ celltype ]] %>% as.character()
  assay(sce, "logcounts") <- counts(sce)
  #
  if( LR.species == 'human' ){
    resource <- 'Consensus'
    lr.database <- liana::select_resource( resource  ) %>% as.data.frame()
    lr.database <- data.frame( ligand = lr.database$Consensus.source_genesymbol,
                               receptor = lr.database$Consensus.target_genesymbol
    )
  }else if( LR.species == 'mouse' ){
    resource <- 'MouseConsensus'
    lr.database <- liana::select_resource( resource  ) %>% as.data.frame()
    lr.database <- data.frame( ligand = lr.database$MouseConsensus.source_genesymbol,
                               receptor = lr.database$MouseConsensus.target_genesymbol
    )
  }

  #
  samples <-  meta.data[[sample]] %>% unique() %>% as.character()
  celltypes <- meta.data[[celltype]] %>% unique() %>% as.character()
  #
  all_group <- mclapply( celltypes , function(x){
    temp <- lapply(celltypes, function(y){
      op <- data.table(   source  = x ,
                          target = y,
                          ligand = lr.database$ligand,
                          receptor = lr.database$receptor
      )
      return(op)
    }) %>% rbindlist()
  },mc.cores = threads   ) %>% rbindlist() %>% data.frame()

  all_group$st <- paste( all_group$source, all_group$target , sep = '_'  )
  all_group$st2 <- paste( all_group$source, all_group$target , sep =  rawToChar(as.raw(c(0xE2, 0x86, 0x92)))  )
  all_group$lr <- paste( all_group$ligand, all_group$receptor , sep = '_'  )
  all_group$CCC.ID <- paste( all_group$st , all_group$lr , sep = '.' )
  #
  LRscore <- pbmclapply(samples, function(x){
    sce_sub <- sce[, meta.data[[sample]] == x ]
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
                                 parallelize =T,  workers = threads
        )
      )
    )
    liana_res$st <- paste( liana_res$source, liana_res$target , sep = '_'  )
    liana_res$lr <- paste( liana_res$ligand, liana_res$receptor , sep = '_'  )
    liana_res$CCC.ID <- paste( liana_res$st , liana_res$lr , sep = '.' )
    #
    myres <- mclapply( all_group$CCC.ID , function(x){
      op = 0
      if(  x %in% liana_res$CCC.ID ){ op = liana_res$LRscore[ which( liana_res$CCC.ID == x )  ]   }
      return( op )
    }, mc.cores = threads  ) %>% unlist() %>% as.numeric()
    #
    return(  data.table(   matrix( myres , nrow = 1   )  )  )

  } , mc.cores = threads  ) %>% rbindlist() %>% setDF()

  LRscore <- data.frame(t(LRscore))

  colnames(LRscore) <- samples
  rownames(LRscore) <- all_group$CCC.ID
  #
  return(  list( CCC.info = all_group , LRscore = LRscore   )  )

}




######comm_score

#' Compute cell–cell communication score (LRscore)
#'
#' @param exp A matrix object with cells as row names and genes as column names.
#' @param meta.data A data.frame object whose row names are identical to those of the exp parameter.
#' @param sample The column name in meta.data that represents the sample ID.
#' @param celltype The column name in meta.data that represents the cell type.
#' @param LR.species Species. Options are human or mouse.
#' @param LR.source Ligand–receptor resource. Options are CCI or liana.
#' @param LR.method Method for calculating LR score. Options are CCI or SingleCellSignalR.
#' @param min.cell The minimum number of cells expressing the ligand or receptor.
#' @param min.exp The minimum expression level of ligands and receptors included in the analysis.
#' @param min.prob The minimum proportion of cells expressing the ligand and receptor.
#' @param threads The number of threads used for parallel computation.
#'
#' @returns A list object.
#' @examples
#'
#' data( 'test.data' )
#' LRscore <- lr_score( exp = t( test.data$exp ) ,
#'                  meta.data = test.data$meta ,
#'                  LR.species = 'mouse' ,
#'                  sample = 'orig.ident' ,
#'                  celltype = 'celltype' )
#'
#' @export
lr_score <- function( exp,meta.data,sample,celltype,
                        LR.species = 'human', LR.source = 'CCI', LR.method = 'CCI',
                        min.cell = 10, min.exp = 0.1, min.prob = 0.3,
                        threads = NULL
){
  ###
  suppressMessages({
    library(dplyr)
    library(pbmcapply)
    library(data.table)
    library(stringr)
    library(parallel)
    library(SingleCellExperiment)
  })

  ###threads
  if( is.null(threads) ){  threads <- parallel::detectCores()   }

  ###LR.source
  if(  is.data.frame( LR.source ) ){
    lr.database <- LR.source
    colnames( lr.database ) <- c(  'ligand' , 'receptor'  )
  }else{
    lr.database <- get_database( species = LR.species , source = LR.source  )
  }
  lr.genes <- unique(  lr.database$ligand , lr.database$receptor  ) %>% unique()

  ###exp,meta.data
  exp <- exp[,  colnames(exp) %in% lr.genes   ] %>% as.matrix()
  if( !identical( rownames(exp), rownames(meta.data) ) ){
    stop( simpleError( 'Mismatch between row names of exp and meta.data.'  ) )
  }

  ###detect
  meta.data <- data.frame( meta.data )
  samples <-  meta.data[[sample]] %>% unique() %>% as.character()
  celltypes <-  meta.data[[celltype]] %>% unique() %>% as.character()
  message(  paste0( "Samples: ",paste( sort(celltypes) ,collapse = ', ' )   )  )
  message(  paste0( "Cell types: ",paste( sort(celltypes) ,collapse = ', ' )   )  )

  message( Sys.time() , ' | ','Checking the expression profiles of ligands and receptors.'  )


  detect_exp <- pbmclapply( colnames(exp) ,function(gene){
    level1 <-  lapply(samples, function(m){
      level2 <- lapply(celltypes, function(n){
        ds <- exp[  meta.data[[sample]] == m & meta.data[[celltype]] == n   ,  gene  ]
        #
        my.return='N'
        prob = length(which(ds <= min.exp)) / length(ds)
        if( length(ds) >= min.cell &  prob >= min.prob   ){
          my.return = 'Y'
        }
        #
        return(  data.table( matrix( c(sample = m , celltype = n ,
                                       gene = gene , prob = prob , reserved = my.return ) ,
                                     nrow =1 )   )   )
        #

      }) %>% rbindlist()

    }) %>% rbindlist()

  }, mc.cores = threads   ) %>% rbindlist()
  colnames( detect_exp ) <- c( 'sample' , 'celltype', 'gene', 'prob' ,'reserved'   )


  ###LRscore
  message( Sys.time() , ' | ','Calculating communication strength score (LRscore).'  )
  if ( LR.method == 'CCI'  ){
    ccc.res <- cci_lrscore( exp = exp , meta.data = meta.data, sample =  sample , celltype =  celltype,
                            lr.database =lr.database , detect_exp = detect_exp , threads = threads )
  }else{
    ccc.res <- liana_lrscore( exp = exp , meta.data = meta.data, sample =  sample , celltype =  celltype,
                              lr.database =lr.database , LR.species = LR.species,
                              min.cell = min.cell  , min.prob = min.prob , threads = threads )
  }

  ###output
  message( Sys.time() , ' | ','Done.'  )
  return( list(
    CCC.info = ccc.res$CCC.info,
    LRscore = ccc.res$LRscore,
    parameters = list(
      exp = exp ,meta.data = meta.data ,sample = sample,celltype = celltype,
      LR.species = LR.species, LR.source = LR.source, LR.method = LR.method,
      min.cell = min.cell, min.exp = min.exp, min.prob = min.prob,
      threads = threads
    )
  ) )
}
