do_circos <- function(
    res, data.used = 'filtered.result',fill = 'p.adj', threshold =  0.05 ,max.group = NULL,
    ligand = NULL, receptor = NULL, sender = NULL, receiver = NULL ,CCC.ID = NULL,
    title = NULL,plot.title = 14,merge = T,LR=T,colors = NULL,link.alpha = 0.5,labels_rot = F
){
  #
  suppressMessages({
    library( ggplot2 )
    library( dplyr  )
    library( data.table )
    library( stringr )
    library( reshape2  )
    library( sjPlot  )
    library( aplot )
    library( plotthis )
  })
  
  ccc.res <- res
  if(  data.used == 'filtered.result'  ){  ccc.res$result <- ccc.res$filtered.result  }
  ccc.res$result <- ccc.res$result[ !is.na( ccc.res$result$p ) , ]
  
  #filter1
  filter_res1 <- ccc.res$result
  if ( fill == 'p'  ){
    filter_res1 <- filter_res1[ filter_res1$p <  threshold , ]
  }else{
    filter_res1 <- filter_res1[ filter_res1$p.adj <  threshold , ]
  }
  
  #filter2
  filter_res2 <- ccc.res$parameters$data$CCC.info
  if( !is.null(ligand) ){  filter_res2 <- filter_res2[   filter_res2$ligand %in% ligand ,  ]    }
  if( !is.null(receptor) ){  filter_res2 <- filter_res2[   filter_res2$receptor %in% receptor ,  ]    }
  if( !is.null(sender) ){  filter_res2 <- filter_res2[   filter_res2$source %in% sender ,  ]   }
  if( !is.null(receiver) ){  filter_res2 <- filter_res2[   filter_res2$target %in% receiver ,  ]    }
  if( !is.null(CCC.ID) ){  filter_res2 <- filter_res2[   filter_res2$CCC.ID %in% CCC.ID ,  ]    }
  
  #final data
  final_res <- filter_res1[ filter_res1$CCC.ID %in% filter_res2$CCC.ID  ,  ]
  if( nrow( final_res ) == 0  ){ stop(simpleError( 'No data were retained. Please adjust the filtering thresholds.' )  )    }
  
  #plot data
  plot_data <- final_res[  , str_detect(colnames(final_res)  , 'mean'  ) ]
  plot_data <- cbind(plot_data , final_res[ ,c( 'CCC.ID' , fill ) ]  )
  
  plot_data <- reshape2::melt( plot_data , id.vars =c('CCC.ID',fill))
  colnames(plot_data)[3] <- 'Group'
  colnames( plot_data )[4] <- 'LR.score'
  plot_data$Group2 <- str_remove( plot_data$Group , '^mean.'   )
  
  #stat
  ids <- unique( plot_data$CCC.ID  )
  plot_data <- lapply(ids, function(x){
    sd <- subset(plot_data , CCC.ID == x )
    
    op = list(  ccc.id = x ,
                      max.group = sd$Group2[ which.max(  sd$LR.score  )   ],
                      L = filter_res2$ligand[ filter_res2$CCC.ID  == x  ],
                      r = filter_res2$receptor[ filter_res2$CCC.ID  == x  ],
                      LR = filter_res2$lr[ filter_res2$CCC.ID  == x ],
                      s = filter_res2$source[ filter_res2$CCC.ID  == x  ] ,
                      t = filter_res2$target[ filter_res2$CCC.ID  == x  ],
                      st2 = filter_res2$st2[ filter_res2$CCC.ID  == x  ]
    )
    return(  op  )
  }) %>% transpose() %>% as.data.table() %>% setDF()
  
  #
  if( is.null( max.group  ) ){
    max.group = unique(  ccc.res[["parameters"]][["data"]][["parameters"]][["meta.data"]][ , ccc.res[["parameters"]][["group"]]  ]  )[1]
  }
  if( is.null(title ) ){ title =  paste0( max.group , ' vs Others'  ) }
  
  #
  if( LR ){
    final_plot_data <- data.frame(   from = plot_data$L ,  to  = plot_data$r  ,
                                     Group = ifelse( plot_data$max.group == max.group  ,'Up' , 'Down'  )     )
  }else{
    final_plot_data <- data.frame(   from = plot_data$s ,  to  = plot_data$t  ,
                                     Group = ifelse( plot_data$max.group == max.group  ,'Up' , 'Down'  )     )
  }
  
  #
  if (merge){
    final_p <- plotthis::ChordPlot(
      data = final_plot_data,
      from = "from" , to = "to" ,
      links_color = "from",
      palcolor = colors ,
      alpha = link.alpha,
      labels_rot = labels_rot,
      title = title,
      theme_args = list( plot.title = element_text( size = plot.title , hjust = 0.5   ))
    )
  }else{
    final_p_up <- plotthis::ChordPlot(
      data = subset(final_plot_data ,Group == 'Up' ),
      from = "from" , to = "to" ,
      links_color = "from",
      palcolor = colors ,
      alpha = link.alpha,
      labels_rot = labels_rot,
      title = paste0( title,' (Up)' ),
      theme_args = list( plot.title = element_text( size = plot.title , hjust = 0.5   ))
    )
    final_p_down <- plotthis::ChordPlot(
      data = subset(final_plot_data ,Group == 'Down' ),
      from = "from" , to = "to" ,
      links_color = "from",
      palcolor = colors ,
      alpha = link.alpha,
      labels_rot = labels_rot,
      title = paste0( title,' (Down)' ),
      theme_args = list( plot.title = element_text( size = plot.title , hjust = 0.5   ))
    )
    
    #
    final_p <- cowplot::plot_grid( plotlist = list(final_p_up , final_p_down  ) , nrow = 1  )
    
  }
  #
  print( final_p )
  return( final_p )
  #
}


#' Circos plot
#'
#' @description
#' Plot a circos plot for the results of the categorical variable.
#'
#' @param res The object returned by the filterCCC function.
#' @param data.used The data used for plotting. Options are result or filtered.result.
#' @param padj.threshold The threshold for p.adj.
#' @param p.threshold Similar to padj.threshold, this parameter specifies the p-value threshold. When provided, padj.threshold will be ignored.
#' @param max.group Which group is treated as the communication-enhanced group.
#' @param merge Whether upregulated and downregulated CCC events are displayed in a single plot.
#' @param LR Whether to plot ligand–receptor pairs instead of sender–receiver pairs.
#' @param colors Provide a color vector to specify the fill colors of links.
#' @param link.alpha A numeric value specifying the transparency of the links.
#' @param labels.rot Whether to rotate the labels by 90 degrees.
#' @param ligand Character vector for filtering CCC.
#' @param receptor Character vector for filtering CCC.
#' @param sender Character vector for filtering CCC.
#' @param receiver Character vector for filtering CCC.
#' @param CCC.ID Character vector for filtering CCC.
#' @param title Plot title.
#' @param plot.title The font size of the plot title.
#'
#' @returns
#' A ggplot2 object.
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
#' ccc.binary <- multiCCC( data = LRscore , binary.params = list( group = 'Group' , g1 = 'O' , g2 = 'Y'  ) )
#' ccc.anova <- multiCCC( data = LRscore , anova.column = 'batch' )
#'
#' #3.filterCCC
#' filter.binary <- filterCCC(ccc.binary)
#' filter.anova <- filterCCC(ccc.anova)
#'
#' #4.circos
#' plot_circos( res = filter.binary$binary.res , p.threshold = 1   )
#' plot_circos( res = filter.binary$binary.res , p.threshold = 1 , merge = F  )
#' plot_circos( res = filter.binary$binary.res , p.threshold = 1 , LR = T  )
#'
#' plot_circos( res = filter.anova$anova.res , p.threshold = 1 )
#'
#' @export
#'
plot_circos <- function(
    res, data.used = 'filtered.result',padj.threshold =  0.05 , p.threshold =  NULL,max.group = NULL,
    merge = T,LR= F,colors = NULL,link.alpha = 0.7,labels.rot = F,
    ligand = NULL, receptor = NULL, sender = NULL, receiver = NULL ,CCC.ID = NULL,
    title = NULL,plot.title = 14
){
  #
  #
  fill = 'p.adj' ; threshold = padj.threshold
  if( !is.null( p.threshold  )  ){   fill = 'p' ; threshold = p.threshold    }
  
  #
  p <- do_circos(
    res = res, data.used = 'filtered.result',fill = fill, threshold =  threshold ,max.group = max.group,
    ligand = ligand, receptor = receptor, sender = sender, receiver =  receiver,CCC.ID = CCC.ID,
    title = title,plot.title =plot.title,merge = merge ,LR=LR,colors = colors,link.alpha = link.alpha,labels_rot = labels.rot
  )
  #
  return(p)
}




