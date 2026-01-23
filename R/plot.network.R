do_network <- function(res, data.used = 'filtered.result',fill = 'p.adj', threshold =  0.05 ,max.group = NULL,
                       ligand = NULL, receptor = NULL, sender = NULL, receiver = NULL ,CCC.ID = NULL,
                       title = NULL,plot.title = 14,merge = T,
                       link.color = c("#FF0000","#085AFF") , link.alpha = 0.8 , link.weight.range =c(1,5),
                       node.size = NULL , node.border.size = 1 , node.border.color ="grey" , node.fill = "white" ,
                       label.size = 10  , label.color = "black"

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

    op = data.table(  ccc.id = x ,
                      max.group = sd$Group2[ which.max(  sd$LR.score  )   ],
                      L = filter_res2$ligand[ filter_res2$CCC.ID  == x  ],
                      r = filter_res2$receptor[ filter_res2$CCC.ID  == x  ],
                      LR = filter_res2$lr[ filter_res2$CCC.ID  == x ],
                      s = filter_res2$source[ filter_res2$CCC.ID  == x  ] ,
                      t = filter_res2$target[ filter_res2$CCC.ID  == x  ],
                      st2 = filter_res2$st2[ filter_res2$CCC.ID  == x  ]
    )
    return(  op  )
  }) %>% rbindlist() %>% setDF()

  #
  if( is.null( max.group  ) ){
    max.group = unique(  ccc.res[["parameters"]][["data"]][["parameters"]][["meta.data"]][ , ccc.res[["parameters"]][["group"]]  ]  )[1]
  }
  if( is.null(title ) ){ title =  paste0( max.group , ' vs Others'  ) }
  if( is.null( node.size ) ){ node.size = 'number' }
  #
  if( merge ){
    if( nrow( plot_data ) == 0  ){
      stop(simpleError( 'No data were retained. Please adjust the filtering thresholds.' )  )
    }
    #
    plot_links <- plot_data %>%
      group_by(  s, t, st2 ,max.group  ) %>%
      summarise(count = n(), .groups = "drop") %>% as.data.frame()
    plot_links$max.group <- ifelse( plot_links$max.group == max.group , "Up" , "Down"   )
    plot_nodes <- table( c( plot_links$s , plot_links$t  )  ) %>% as.data.frame()
    colnames(plot_nodes) <- c('name' ,'number' )
    plot_nodes$name <- as.character( plot_nodes$name  )
    plot_nodes$number <- as.integer( plot_nodes$number  )
    #
    final_p <- plotthis::Network( links = plot_links ,nodes = plot_nodes ,
                                  from = "s" , to = "t" ,theme = theme_void ,
                                  link_weight_by = "count" , link_weight_name = 'No. of CCCs',link_weight_range = link.weight.range,
                                  link_alpha = link.alpha,link_color_by = "max.group" ,
                                  link_color_name = "Trend",link_palcolor = link.color,
                                  node_by = "name",node_size_by = node.size,node_stroke = node.border.size,
                                  node_color_by = node.border.color,node_fill_by = node.fill ,
                                  label_fg = label.color,label_bg = NA,label_size = label.size,
                                  layout = "circle",title = title,link_curvature = 0.2
    )+theme( plot.title =  element_text( size= plot.title , hjust = 0.5  )  )+
      guides( size = "none"   )
    print(final_p)
    return(final_p)

  }else{
    #####################Up
    splot_data <- plot_data[ plot_data$max.group ==  max.group ,  ]
    if( nrow( splot_data ) == 0  ){
      stop(simpleError( 'No data were retained. Please adjust the filtering thresholds.' )  )
    }
    #
    plot_links <- splot_data %>%
      group_by(  s, t, st2 ,max.group  ) %>%
      summarise(count = n(), .groups = "drop") %>% as.data.frame()
    plot_links$max.group <- ifelse( plot_links$max.group == max.group , "Up" , "Down"   )
    plot_nodes <- table( c( plot_links$s , plot_links$t  )  ) %>% as.data.frame()
    colnames(plot_nodes) <- c('name' ,'number' )
    plot_nodes$name <- as.character( plot_nodes$name  )
    plot_nodes$number <- as.integer( plot_nodes$number  )
    #
    final_p_up <- plotthis::Network( links = plot_links ,nodes = plot_nodes ,
                                     from = "s" , to = "t" ,theme = theme_void ,
                                     link_weight_by = "count" , link_weight_name = 'No. of CCCs',link_weight_range = link.weight.range,
                                     link_alpha = link.alpha,link_color_by = "max.group" ,
                                     link_color_name = "Trend",link_palcolor = link.color[[1]],
                                     node_by = "name",node_size_by = node.size,node_stroke = node.border.size,
                                     node_color_by = node.border.color,node_fill_by = node.fill ,
                                     label_fg = label.color,label_bg = NA,label_size = label.size,
                                     layout = "circle",title = paste0( title,' (Up)' ),
                                     link_curvature = 0.2
    )+theme( plot.title =  element_text( size= plot.title , hjust = 0.5  )  )+
      guides( size = "none"   )

    #down
    splot_data <- plot_data[ plot_data$max.group !=  max.group ,  ]
    if( nrow( splot_data ) == 0  ){
      stop(simpleError( 'No data were retained. Please adjust the filtering thresholds.' )  )
    }
    #
    plot_links <- splot_data %>%
      group_by(  s, t, st2 ,max.group  ) %>%
      summarise(count = n(), .groups = "drop") %>% as.data.frame()
    plot_links$max.group <- ifelse( plot_links$max.group == max.group , "Up" , "Down"   )
    plot_nodes <- table( c( plot_links$s , plot_links$t  )  ) %>% as.data.frame()
    colnames(plot_nodes) <- c('name' ,'number' )
    plot_nodes$name <- as.character( plot_nodes$name  )
    plot_nodes$number <- as.integer( plot_nodes$number  )
    #
    final_p_down <- plotthis::Network( links = plot_links ,nodes = plot_nodes ,
                                       from = "s" , to = "t" ,theme = theme_void ,
                                       link_weight_by = "count" , link_weight_name = 'No. of CCCs',link_weight_range = link.weight.range,
                                       link_alpha = link.alpha,link_color_by = "max.group" ,
                                       link_color_name = "Trend",link_palcolor = link.color[[2]],
                                       node_by = "name",node_size_by = node.size,node_stroke = node.border.size,
                                       node_color_by = node.border.color,node_fill_by = node.fill ,
                                       label_fg = label.color,label_bg = NA,label_size = label.size,
                                       layout = "circle",title = paste0( title,' (Down)' ),
                                       link_curvature = 0.2
    )+theme( plot.title =  element_text( size= plot.title , hjust = 0.5  )  )+
      guides( size = "none"  )

    #
    final_p <- cowplot::plot_grid( plotlist = list( final_p_up , final_p_down ) , nrow = 1  )
    #
    print( final_p )
    return( final_p )
  }
}



#' Network plot of sender-receiver
#'
#' @description
#' A network visualization of cell–cell interactions for categorical variables.
#'
#'
#' @param res The object returned by the filterCCC function.
#' @param data.used The data used for plotting. Options are result or filtered.result.
#' @param padj.threshold The threshold for p.adj.
#' @param p.threshold Similar to padj.threshold, this parameter specifies the p-value threshold. When provided, padj.threshold will be ignored.
#' @param max.group Which group is treated as the communication-enhanced group.
#' @param ligand Character vector for filtering CCC.
#' @param receptor Character vector for filtering CCC.
#' @param sender Character vector for filtering CCC.
#' @param receiver Character vector for filtering CCC.
#' @param CCC.ID Character vector for filtering CCC.
#' @param title Plot title.
#' @param plot.title The font size of the plot title.
#' @param merge Whether upregulated and downregulated CCC events are displayed in a single plot.
#' @param link.color Provide a two-value color vector to specify the fill colors of links.
#' @param link.alpha A numeric value specifying the transparency of the links.
#' @param link.weight.range Provide a two-value numeric vector specifying the range of link widths.
#' @param node.size Node size.
#' @param node.border.size Width of the node outline.
#' @param node.border.color Color of the node outline.
#' @param node.fill Fill color of the node.
#' @param label.size Size of the text inside the node.
#' @param label.color Color of the text inside the node.
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
#' #4.network
#' plot_network( res = filter.binary$binary.res , p.threshold = 1  )
#' plot_network( res = filter.binary$binary.res , p.threshold = 1 , merge = F )
#'
#'
#' @export
#'
plot_network <- function(res, data.used = 'filtered.result',padj.threshold =  0.05 , p.threshold =  NULL,max.group = NULL,
                         ligand = NULL, receptor = NULL, sender = NULL, receiver = NULL ,CCC.ID = NULL,
                         title = NULL,plot.title = 14,merge = T,
                         link.color = c("#FF0000","#085AFF") , link.alpha = 0.8 , link.weight.range =c(1,5),
                         node.size = NULL , node.border.size = 1 , node.border.color ="grey" , node.fill = "white" ,
                         label.size = 5  , label.color = "black"
){
  #
  fill = 'p.adj' ; threshold = padj.threshold
  if( !is.null( p.threshold  )  ){   fill = 'p' ; threshold = p.threshold    }
  #
  p <- do_network(res = res, data.used = data.used,fill = fill, threshold = threshold, max.group = max.group,
                  ligand = ligand, receptor = receptor, sender = sender, receiver = receiver ,CCC.ID = CCC.ID,
                  title = title,plot.title = plot.title,merge = merge,
                  link.color = link.color , link.alpha = link.alpha , link.weight.range = link.weight.range,
                  node.size = node.size , node.border.size = node.border.size , node.border.color = node.border.color , node.fill = node.fill ,
                  label.size = label.size  , label.color = label.color

  )
  return(p)
  #
}
