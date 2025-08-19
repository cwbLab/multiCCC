

######
#' Dot plot
#'
#' @description
#' Plot a dot plot for the results of a categorical variable.
#'
#' @param ccc.res The object returned by the multiCCC function.
#' @param fill The column used for color fill, optionally p or p.adj.
#' @param threshold The threshold for p or p.adj.
#' @param colors A character vector.
#' @param ligand A character vector.
#' @param receptor A character vector.
#' @param source.cell A character vector.
#' @param target.cell A character vector.
#' @param CCC.ID A character vector.
#' @param dot.max The maximum value of the dots.
#' @param dot.min The minimum value of the dots.
#' @param legend.position The position of the legend in the plot.
#' @param strip.text The font size of the facet labels.
#' @param x.text The font size of the X-axis labels.
#' @param y.text The font size of the Y-axis labels.
#' @param facet.row The number of rows in the facets.
#'
#' @returns
#' A ggplot2 object.
#'
#' @export
#'
#' @examples
#' plot_dot( ccc.res = ccc.anova$anova.res ,
#'           ligand = c('Sele' , 'Fgf1' , 'Nrxn1' ,'Nrxn2' ),
#'           threshold =  1 ,strip.text = 8 )
#'
#'
plot_dot <- function( ccc.res , fill = 'p.adj', threshold =  0.05 , colors =  c('#0000FF', "#F7F7F7", '#FF0000'),
                      ligand = NULL, receptor = NULL, source.cell = NULL, target.cell = NULL ,CCC.ID = NULL,
                      dot.max = 3 , dot.min = 1 , legend.position = "right" ,
                      strip.text = 10,x.text = 10 ,y.text = 10 , facet.row = 1
){

  #
  suppressMessages({
  library( ggplot2 )
  library( dplyr  )
  library( data.table )
  library( stringr )
  library( reshape2  )
  library( sjPlot  )
  })

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
  if( !is.null(source.cell) ){  filter_res2 <- filter_res2[   filter_res2$source %in% source.cell ,  ]   }
  if( !is.null(target.cell) ){  filter_res2 <- filter_res2[   filter_res2$target %in% target.cell ,  ]    }
  if( !is.null(CCC.ID) ){  filter_res2 <- filter_res2[   filter_res2$CCC.ID %in% CCC.ID ,  ]    }

  #final data
  final_res <- filter_res1[ filter_res1$CCC.ID %in% filter_res2$CCC.ID  ,  ]


  #plot data
  plot_data <- final_res[  , str_detect(colnames(final_res)  , 'mean'  ) ]
  plot_data <- cbind(plot_data , final_res[ ,c( 'CCC.ID' , fill ) ]  )

  plot_data <- reshape2::melt( plot_data , id.vars =c('CCC.ID',fill))
  colnames(plot_data)[3] <- 'Group'
  colnames( plot_data )[4] <- 'LR.score'

  #shape
  ids <- unique( plot_data$CCC.ID  )
  plot_data <- lapply(ids, function(x){
    sd <- subset(plot_data , CCC.ID == x )
    sd$shape <- 'Middle'
    sd$shape[ which.min( sd$LR.score  ) ] <- 'Min'
    sd$shape[ which.max( sd$LR.score  ) ] <- 'Max'
    #
    sd$LR <- filter_res2$lr[ filter_res2$CCC.ID  == sd$CCC.ID  ]
    sd$st2 <- filter_res2$st2[ filter_res2$CCC.ID  == sd$CCC.ID  ]
    #
    return( sd )
  }) %>% rbindlist()

  #
  plot_data$Group <- str_remove_all( plot_data$Group , 'mean.'  )
  plot_data$LR.score <- as.numeric( plot_data$LR.score )
  plot_data[[fill]] <- as.numeric( plot_data[[fill]]  )
  plot_data$Rank <- plot_data$shape

  #plot
  if( nrow( plot_data ) == 0 ){
    stop( simpleError( 'No CCCs were retained. Please adjust the filtering threshold.'  ) )
  }
  p <- ggplot(plot_data, aes(x = Group, y = LR )) +
    geom_point( aes(  fill = -log10( !!sym( fill ) )  , size = LR.score , shape = Rank , stroke = 0.0001 ) , color = 'white'  ) +
    scale_fill_gradientn( colours = colors ) +
    scale_size( range = c(dot.min  , dot.max )  ) +
    scale_shape_manual(values = c(Max  = 24, Middle = 21 , Min = 25 )) + # 根据 shape 列类别数量调整
    facet_wrap(~st2 ,  scales = "free_x"  , nrow = facet.row  , strip.position = "bottom"  ) +
    guides(
      shape =  guide_legend( override.aes = list( color = 'black' )  ),
      size = guide_legend( override.aes = list( color = 'black' )  )
    )+
    theme_bw(base_size = 14) +
    theme(
      text = element_text( color = 'black'  ),
      panel.grid.major = element_line(color = "grey90" , linetype = "dotted" ),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(  size = x.text ,  angle = 45, vjust = 1 ,hjust = 1),
      axis.text.y = element_text( size = y.text ) ,
      strip.background = element_blank(),
      strip.text = element_text(size = strip.text, face = "bold"),
      strip.placement = "outside",
      legend.position = legend.position
    )
  print(p)
  #
  return( p )
}




######
#' Plot a fitted curve for a continuous variable.
#'
#' @param ccc.res The object returned by the multiCCC function.
#' @param CCC.ID The ID of the CCC (cell–cell communication) event.
#' @param p Whether to display the p-value.
#' @param p.adj Whether to display the p-adjusted.
#' @param colors The color of the line.
#' @param dot.size The size of the dots.
#' @param model.text The font size of the model results.
#' @param axis.text The font size of the axis labels.
#' @param axis.title The font size of the axis titles.
#' @param plot.title The font size of the plot title.
#'
#' @returns
#' A ggplot2 object.
#'
#' @export
#'
#' @examples
#' CCC.ID <- 'Monocyte_Monocyte.C3_Itgb2'
#' plot_line( ccc.res = ccc.glm$glm.res,
#'            CCC.ID = CCC.ID , p.adj = F )
#'
plot_line <- function( ccc.res, CCC.ID, p = T, p.adj =T,
                       colors =  "blue", dot.size = 2, model.text = 3 ,
                       axis.text = 8 , axis.title = 10 , plot.title = 12
){
  #
  suppressMessages({
  library( ggplot2 )
  library( dplyr  )
  library( data.table )
  library( stringr )
  library( reshape2  )
  library( sjPlot  )
  })

  #
  model_res <- ccc.res$model[[CCC.ID]]
  model_summary <- summary(model_res)
  #
  coef_time <- round(model_summary$coefficients[ colnames( ccc.res[["meta.data"]] )[2]  , "Estimate"], 3)
  pval_time <- signif(model_summary$coefficients[colnames( ccc.res[["meta.data"]] )[2], "Pr(>|t|)"], 3)
  padj <- signif( ccc.res$result$p.adj[ ccc.res$result$CCC.ID == CCC.ID ] , 3 )

  # P text
  label_text <- paste0("β = ", coef_time)
  if ( p ){ label_text <- paste0( label_text, "\np = ", pval_time) }
  if ( p.adj ){ label_text <- paste0( label_text, "\np.adj = ", padj) }

  #
  title = ccc.res$parameters$data$CCC.info[ ccc.res$parameters$data$CCC.info$CCC.ID == CCC.ID,  ]
  title = paste0( title$st2 , ' (', title$lr , ')' )

  # line plot
  p <- NULL

  #glm
  if ( ccc.res$type == "glm" ){

    suppressMessages(
      p <- sjPlot::plot_model(model_res,
                              type = "pred",
                              terms = colnames( ccc.res[["meta.data"]] )[2],
                              colors = colors,
                              show.data = T,
                              dot.size = dot.size ,
                              title = title
      )+
        labs( x  = colnames( ccc.res[["meta.data"]] )[2] , y = 'LRscore' )+
        annotate(
          "text",
          x = min(model_res$model[[ colnames( ccc.res[["meta.data"]] )[2] ]], na.rm = TRUE),
          y = max(model_res$model[[ colnames( ccc.res[["meta.data"]] )[2] ]], na.rm = TRUE),
          label = label_text,
          hjust = 0, vjust = 1,
          size = model.text
        )+
        theme_bw()+
        theme(
          plot.title = element_text(  hjust = 0.5 , size = plot.title ),
          axis.text = element_text( size = axis.text ),
          axis.title = element_text( size = axis.title )
        )
    )

  }

  #time course
  if ( ccc.res$type == "time.course" ){

    suppressMessages(
      p <- sjPlot::plot_model(model_res,
                              type = "pred",
                              terms = colnames( ccc.res[["meta.data"]] )[2],
                              colors = colors,
                              show.data = T,
                              dot.size = dot.size ,
                              title = title
      )+
        labs( x  = colnames( ccc.res[["meta.data"]] )[2] , y = 'LRscore' )+
        annotate(
          "text",
          x = min(model_res@frame[[ colnames( ccc.res[["meta.data"]] )[2] ]], na.rm = TRUE),
          y = max(model_res@frame[[ colnames( ccc.res[["meta.data"]] )[2] ]], na.rm = TRUE),
          label = label_text,
          hjust = 0, vjust = 1,
          size = model.text
        )+
        theme_bw()+
        theme(
          plot.title = element_text(  hjust = 0.5 , size = plot.title ),
          axis.text = element_text( size = axis.text ),
          axis.title = element_text( size = axis.title )
        )
    )

  }
  print(p)
  #
  return( p )
}

