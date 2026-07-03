filter_ccc <- function( ccc.res , ligand.method = NULL, receptor.method = NULL ){
  #
  suppressMessages(  library( SingleCellExperiment )   )
  #
  methods <- names( ccc.res )
  ccc.res <- lapply(methods, function( method.name ){
    message(  '[ ', format(Sys.time(), "%Y-%m-%d %H:%M:%S") , ' ] ' , method.name    )
    
    x = ccc.res[[ method.name ]]
    #
    sce <- suppressWarnings(
      SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = t(as.matrix( x[["parameters"]][["data"]][["parameters"]][["exp"]] ))),
        metadata = x[["parameters"]][["data"]][["parameters"]][["meta.data"]]
      )
    )
    sce$celltype <- sce@metadata[[ x[["parameters"]][["data"]][["parameters"]][["celltype"]] ]] %>% as.character()
    sce@assays@data@listData[["logcounts"]] <- sce@assays@data@listData[["counts"]]
    
    #markers
    if (  'gene.expression' %in% names(x)  ){
      message(  "'gene.expression' results already exist and will be used directly for filtering."  )
      marker_res <- x[['gene.expression']]$markers
    }else{
      marker_res <- presto::wilcoxauc(sce, group_by = "celltype", assay = "logcounts" )
      x[['gene.expression']] <- list( markers = marker_res , ligand.filter = ligand.method ,
                                      receptor.filter = receptor.method, filtered.ligand = NULL,
                                      filtered.receptor = NULL
      )
    }
    
    #filter
    filtered_result <- x[["result"]]
    filtered_result <- subset(filtered_result , ! is.na( p.adj )  )
    
    CCC.info <-  x[["parameters"]][["data"]][["CCC.info"]]
    CCC.info$ligand.filter <- paste( CCC.info$source , CCC.info$ligand , sep='*'  )
    CCC.info$receptor.filter <- paste( CCC.info$target , CCC.info$receptor , sep='*'  )
    
    if ( ! is.null( ligand.method )  ){
      x$gene.expression$ligand.filter <- ligand.method
      filtered_markers <- subset(  marker_res ,eval(parse(text =  ligand.method )) )
      x$gene.expression$filtered.ligand <- filtered_markers
      filtered_markers$filter <- paste( filtered_markers$group , filtered_markers$feature , sep='*'  )
      CCC.info <- subset( CCC.info ,  ligand.filter %in%  filtered_markers$filter )
    }
    if ( ! is.null( receptor.method )  ){
      x$gene.expression$receptor.filter <- receptor.method
      filtered_markers <- subset(  marker_res  ,eval(parse(text = receptor.method )) )
      x$gene.expression$filtered.receptor <- filtered_markers
      filtered_markers$filter <- paste( filtered_markers$group , filtered_markers$feature , sep='*'  )
      CCC.info <- subset( CCC.info ,  receptor.filter %in%  filtered_markers$filter )
    }
    #
    filtered_result <- subset( filtered_result , CCC.ID %in% CCC.info$CCC.ID  )
    x[[ 'filtered.result' ]] <- filtered_result
    #
    if(  method.name %in% c('glm.res' , 'time.course.res' )  ){
      x <- x[ c(  "meta.data", "result", "gene.expression", "filtered.result", "model","type", "parameters")]
    }else{
      x <- x[ c(  "meta.data", "result", "gene.expression", "filtered.result",  "type", "parameters")]
    }
    return(  x  )
  })
  names( ccc.res ) <- methods
  #
  return( ccc.res )
}

#
#' Filter the CCC results
#'
#' @description
#' Perform the Wilcoxon rank sum test to identify cell-specific markers, and then use the resulting gene list for filtering the CCC results.
#'
#' @param data An object returned by the multiCCC or filterCCC function.
#'
#' Filtering is always applied to the original result table of the input object, preventing cumulative filtering effects.
#' @param ligand.filter Provide a string with the rule defined by the subset function to filter CCC based on the expression of ligands. If NULL, no filtering will be performed.
#'
#' Available filtering metrics include feature, group, avgExpr, logFC, statistic, auc, pval, padj, pct_in, and pct_out.
#' @param receptor.filter Similar to the ligand.filter parameter, provide a string with the rule defined by the subset function to filter CCC based on the expression of receptors. If NULL, no filtering will be performed.
#'
#' @returns
#' The main output content is stored in 'gene.expression' and 'filtered.result'. 'gene.expression' represents the evaluation results of cell-specific markers, and 'filtered.result' represents the filtered CCC results.
#'
#' Note: The filterCCC function simultaneously applies filtering to all analysis results, without requiring separate filtering for different application scenarios.
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
#' ccc.all <- multiCCC( data = LRscore, 
#' 	   binary.params = list( group = 'Group' , g1 = 'O' , g2 = 'Y'  ),
#' 	   anova.column = 'batch',
#' 	   glm.column = 'weight',
#' 	   time.course.params = list( time  = 'time' , replicate  = 'replicate'  )
#' 	)
#'
#' #3.filterCCC
#' #Note: The filterCCC function simultaneously applies filtering to all analysis results, without requiring separate filtering for different application scenarios.
#' filter.binary <- filterCCC(ccc.binary)
#' filter.anova <- filterCCC(ccc.anova)
#' filter.glm <- filterCCC(ccc.glm )
#' filter.time <- filterCCC(ccc.time)
#' 
#' #One-step execution for all scenarios
#' filter.all <- filterCCC(ccc.all)
#'
#' @export
filterCCC <- function( data , ligand.filter = 'logFC > 0.3 & padj < 0.05' , receptor.filter = NULL  ){
  run.start = Sys.time()
  result <- filter_ccc( ccc.res = data , ligand.method = ligand.filter , receptor.method = receptor.filter )
  run.end = Sys.time()
  message( '[ ', format(Sys.time(), "%Y-%m-%d %H:%M:%S") , ' ] ','Done. Total runtime: ', hms::as_hms( as.numeric( run.end - run.start, units = "secs") ) ,'.' )
  #
  return( result )
}
