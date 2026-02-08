do_heatmap <- function( ccc.res , data.used = 'filtered.result',fill = 'p.adj', threshold =  0.05 ,max.group = NULL,
                        colors =  list( low =  '#002BFF', mid =  "#F7F7F7", high =  '#FF2B00' , midpoint = 0 ),
                        border.color = 'white', border.size = 0.5,
                        ligand = NULL, receptor = NULL, sender = NULL, receiver = NULL ,CCC.ID = NULL,
                        legend.direction = "horizontal" ,x.text = 10 ,y.text = 10 , x.title = 12 , y.title = 12 ,
                        title = NULL, plot.title = 14,
                        aspect.ratio = 1, point.size = c( 3,6 ), point2heatmap = 0.5 ,
                        merge = T,trend = 'Up',LR = F,
                        top.up = NULL, top.down = NULL, top.all = NULL
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
  })
  
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
  plot_data <- final_res[  , str_detect(colnames(final_res)  , 'mean'  ) ] %>% setDT()
  plot_data <- plot_data[, (names(plot_data)) := lapply(.SD, as.numeric)] %>% setDF()
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
  }) %>% rbindlist() %>% setDF()
  
  #
  if( is.null( max.group  ) ){
    max.group = unique(  ccc.res[["parameters"]][["data"]][["parameters"]][["meta.data"]][ , ccc.res[["parameters"]][["group"]]  ]  )[1]
  }
  if( is.null(title ) ){ title =  paste0( max.group , ' vs Others'  )  }
  
  ########################CELL
  if ( !LR  ){
    #
    cells <- unique( as.character( ccc.res[["gene.expression"]][["markers"]][["group"]] ) )
    cells <- sort(cells ,decreasing = F )
    heatmap_data <- expand.grid( cells  , cells ) %>% data.frame()
    colnames(heatmap_data) <- c('source','target')
    heatmap_data <- lapply(1:nrow( heatmap_data  ), function(x){
      #
      d <- heatmap_data[x,] %>% unlist()
      sdata <- subset( plot_data , s == d[[1]] & t == d[[2]]   )
      op <- data.table( target = d[[2]] , source = d[[1]]  ,
                        number = c(
                          nrow( sdata[ sdata$max.group == max.group ,  ]   ),
                          nrow( sdata[ sdata$max.group != max.group ,  ]   )
                        ),
                        trend = c('UP','DOWN')
      )
      return( op )
      
      #
    }) %>% rbindlist()
    colnames(heatmap_data) <- c('X','Y','value','trend')
    
    #X ,sender
    #Y ,receiver
    
    #
    heatmap_data$X <- as.character( heatmap_data$X )
    heatmap_data$Y <- as.character( heatmap_data$Y )
    
    #
    heatmap_data_test <- heatmap_data[ heatmap_data$value != 0 ,  ]
    if( length( unique( heatmap_data_test$X ) ) < 2 | length( unique( heatmap_data_test$Y ) ) < 2 ){
      stop(simpleError( 'No data were retained. Please adjust the filtering thresholds.' )  )
    }
    heatmap_data <- heatmap_data[ heatmap_data$Y %in% heatmap_data_test$Y   ,   ]
    heatmap_data <- heatmap_data[ heatmap_data$X %in% heatmap_data_test$X   ,   ]
    
    #filter
    if( !is.null(top.up) ){
      #
      a <- table( heatmap_data$Y[ heatmap_data$trend == 'UP' & heatmap_data$value != 0   ] )
      if( top.up > length(a) ){ top.up <- length(a) }
      c1 <- names( a[ a >= quantile( a , top.up / length(a) )] )
      #
      a <- table( heatmap_data$X[ heatmap_data$trend == 'UP' & heatmap_data$value != 0   ] )
      if( top.up > length(a) ){ top.up <- length(a) }
      c2 <- names( a[ a >= quantile( a , top.up / length(a) )] )
      
      #
      heatmap_data <- heatmap_data[  heatmap_data$X %in% c2 & heatmap_data$Y %in% c1   , ]
    }
    if( !is.null(top.down) ){
      #
      a <- table( heatmap_data$Y[ heatmap_data$trend == 'DOWN' & heatmap_data$value != 0   ] )
      if( top.down > length(a) ){ top.down <- length(a) }
      c1 <- names( a[ a >= quantile( a , top.down / length(a) )] )
      #
      a <- table( heatmap_data$X[ heatmap_data$trend == 'DOWN' & heatmap_data$value != 0   ] )
      if( top.down > length(a) ){ top.down <- length(a) }
      c2 <- names( a[ a >= quantile( a , top.down / length(a) )] )
      
      #
      heatmap_data <- heatmap_data[  heatmap_data$X %in% c2 & heatmap_data$Y %in% c1   , ]
    }
    if( !is.null(top.all) ){
      #
      a <- table( heatmap_data$Y[  heatmap_data$value != 0   ] )
      if( top.all > length(a) ){ top.all <- length(a) }
      c1 <- names( a[ a >= quantile( a , top.all / length(a) )] )
      #
      a <- table( heatmap_data$X[  heatmap_data$value != 0   ] )
      if( top.all > length(a) ){ top.all <- length(a) }
      c2 <- names( a[ a >= quantile( a , top.all / length(a) )] )
      
      #
      heatmap_data <- heatmap_data[  heatmap_data$X %in% c2 & heatmap_data$Y %in% c1   , ]
    }
    
    
    #
    if( isTRUE( merge )  ){
      #
      get_triangles <- function(df) {
        #
        df$x_num <- as.numeric(as.factor(df$X))
        df$y_num <- as.numeric(as.factor(df$Y))
        
        up_tri <- df %>% filter(trend == "UP") %>% rowwise() %>%
          do(data.frame(
            X = .$X, Y = .$Y, value = .$value, trend = .$trend,
            pos_x = c(.$x_num - 0.5, .$x_num + 0.5, .$x_num + 0.5),
            pos_y = c(.$y_num + 0.5, .$y_num + 0.5, .$y_num - 0.5),
            group = paste0(.$X, .$Y, "UP")
          ))
        
        down_tri <- df %>% filter(trend == "DOWN") %>% rowwise() %>%
          do(data.frame(
            X = .$X, Y = .$Y, value = .$value, trend = .$trend,
            pos_x = c(.$x_num - 0.5, .$x_num - 0.5, .$x_num + 0.5),
            pos_y = c(.$y_num + 0.5, .$y_num - 0.5, .$y_num - 0.5),
            group = paste0(.$X, .$Y, "DOWN")
          ))
        
        return(rbind(up_tri, down_tri))
      }
      #
      triangles_data <- get_triangles(heatmap_data)
      triangles_data$value2 <- triangles_data$value
      triangles_data$value2[ which( triangles_data$trend == 'DOWN'  )  ] <- -triangles_data$value2[ which( triangles_data$trend == 'DOWN'  )  ]
      
      #
      heatmap_main <- ggplot(triangles_data, aes(x = pos_x, y = pos_y, fill = value2, group = group)) +
        geom_polygon(color = border.color, size = border.size ) +
        scale_fill_gradient2(low = colors$low, mid = colors$mid, high = colors$high, midpoint = colors$midpoint ) +
        scale_x_discrete( expand = c( 0 , 0 ) ,  limits = cells[ cells %in% heatmap_data$X  ] ) +
        scale_y_discrete( expand = c( 0 , 0 ) , limits = cells[ cells %in% heatmap_data$Y  ]  ) +
        labs( x = "Sender", y = "Receiver") +
        coord_fixed( ratio = aspect.ratio   ) +
        theme_void(   )+
        guides( fill=guide_colorbar( title = 'No. of CCCs'  )  )+
        theme(
          text = element_text( color = 'black'  ),
          plot.title = element_text( size = plot.title , hjust = 0.5  ),
          axis.title.x = element_text(  size = x.title , hjust = 0.5 ),
          axis.title.y = element_text( size = y.title , angle = 90 ) ,
          axis.text.x = element_text(  size = x.text , vjust = 0.5 ,angle = 90 ,hjust = 1  ),
          axis.text.y = element_text( size = y.text ,hjust = 1 ) ,
          plot.margin = margin( 1, 1, 1, 1),
          legend.position = "right",
          legend.direction = legend.direction
        )
      
      #
      margin_x <- heatmap_data %>% group_by(X, trend) %>% summarise(total = sum(value))
      margin_x$alpha <- ifelse( margin_x$trend == 'UP', 1 , 0.9 )
      margin_y <- heatmap_data %>% group_by(Y, trend) %>% summarise(total = sum(value))
      margin_y$alpha <- ifelse( margin_y$trend == 'UP', 1 , 0.9  )
      
      #receiver
      p_top <- ggplot(margin_x, aes(x = X, y = total, color = trend)) +
        geom_line(aes(group = X), color = "black" ) +
        geom_point( aes(size = total , alpha = alpha  ) ) + scale_size( range = point.size )+
        scale_alpha( range = c( 1,1 ) )+
        scale_x_discrete( limits = cells[ cells %in% heatmap_data$X  ]   ) +
        scale_color_manual(values = c("UP" = colors$high, "DOWN" = colors$low)) +
        labs(  y = 'Total CCCs (Receiver)' ,
               title = title  )+
        theme_classic()+
        guides( color=guide_legend( title = 'Trend' , override.aes = list( size  = max( point.size ) ) ) ,
                size = "none" ,
                alpha = "none"
        )+
        theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.y = element_line( linetype = "dashed" ,color ="grey80"   ),
          axis.title.y = element_text( size = y.title , angle = 90 ) ,
          axis.text.y = element_text( size = y.text ,hjust = 1 ) ,
          aspect.ratio = point2heatmap * aspect.ratio * ( length(unique( margin_y$Y )) / length(unique( margin_x$X ))  ) ,
          plot.margin = margin(0, 0, 0, 0 ),
          legend.position = "right",
          legend.direction = legend.direction,
          plot.title = element_text( hjust = 0.5 , size = plot.title   )
        )
      
      
      #sender
      right_ratio <- 1 / point2heatmap
      colnames(margin_y)[1] <- 'X'
      p_right <- ggplot(margin_y, aes(x = total, y = X, color = trend)) +
        geom_line(aes(group = X), color = "black" ) +
        geom_point( aes(size = total , alpha = alpha  ) ) + scale_size( range = point.size )+
        scale_alpha( range = c( 1,1 ) )+
        scale_y_discrete( limits = cells[ cells %in% heatmap_data$Y  ]  ) +
        scale_color_manual(values = c("UP" = colors$high, "DOWN" = colors$low)) +
        labs(  x = 'Total CCCs (Sender)' )+
        theme(legend.position = "none") +
        theme_classic()+
        guides( color=guide_legend( title = 'Trend' , override.aes = list( size  = max( point.size ) ) ) ,
                size = "none" ,
                alpha = "none"
        )+
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.x = element_line( linetype = "dashed" ,color ="grey80"   ),
          axis.title.x = element_text( size = x.title , hjust = 0.5 ) ,
          axis.text.x = element_text( size = x.text ,hjust = 0.5 ) ,
          aspect.ratio = right_ratio  ,
          plot.margin = margin( 0, 0, 0, 0),
          legend.position = "right",
          legend.direction = legend.direction
        )
      
      #
      final_p <- heatmap_main %>%  insert_top(p_top,height= point2heatmap ) %>% insert_right( p_right )
      
      #
      print( final_p  )
      return( final_p )
      
    }else{
      #split,UP,DOWN
      
      if( trend == 'Up' ){
        heatmap_data2 <- heatmap_data[ heatmap_data$trend == 'UP' ,  ]
        heatmap_data_test <- heatmap_data2[ heatmap_data2$value != 0 ,  ]
        if( length( unique( heatmap_data_test$X ) ) < 2 | length( unique( heatmap_data_test$Y ) ) < 2 ){
          stop(simpleError( 'No data were retained. Please adjust the filtering thresholds.' )  )
        }
        heatmap_data2 <- heatmap_data2[ heatmap_data2$Y %in% heatmap_data_test$Y   ,   ]
        heatmap_data2 <- heatmap_data2[ heatmap_data2$X %in% heatmap_data_test$X   ,   ]
        
        
        heatmap_main2 <- ggplot( heatmap_data2, ) +
          geom_tile( aes( x = X , y = Y , fill =  value  ), color = border.color, size = border.size ) +
          scale_fill_gradient2(low = colors$mid, high = colors$high ) +
          scale_x_discrete( expand = c( 0 , 0 ) ,  limits = cells[ cells %in% heatmap_data2$X  ] ) +
          scale_y_discrete( expand = c( 0 , 0 ) , limits = cells[ cells %in% heatmap_data2$Y  ]  ) +
          labs( x = "Receiver", y = "Sender") +
          coord_fixed( ratio = aspect.ratio   ) +
          theme_void(   )+
          guides( fill=guide_colorbar( title = 'No. of CCCs'  )  )+
          theme(
            text = element_text( color = 'black'  ),
            plot.title = element_text( size = plot.title , hjust = 0.5  ),
            axis.title.x = element_text(  size = x.title , hjust = 0.5 ),
            axis.title.y = element_text( size = y.title , angle = 90 ) ,
            axis.text.x = element_text(  size = x.text , vjust = 0.5 ,angle = 90 ,hjust = 1  ),
            axis.text.y = element_text( size = y.text ,hjust = 1 ) ,
            plot.margin = margin( 1, 1, 1, 1),
            legend.position = "left",
            legend.direction = legend.direction
          )
        
        #
        margin_x <- heatmap_data2 %>% group_by(X, trend) %>% summarise(total = sum(value))
        margin_y <- heatmap_data2 %>% group_by(Y, trend) %>% summarise(total = sum(value))
        
        #receiver
        p_top <- ggplot(margin_x, aes(x = X, y = total)) +
          geom_point( aes(size = total , alpha = total  ) , color = colors$high  ) + scale_size( range = point.size )+
          scale_x_discrete( limits = cells[ cells %in% heatmap_data2$X  ]   ) +
          labs(  y = 'Total CCCs (Receiver)' ,
                 title = paste0( title,' (Up)' )  )+
          theme_classic()+
          guides( color="none" ,
                  size = "none"  ,
                  alpha = "none"
          )+
          theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major.y = element_line( linetype = "dashed" ,color ="grey80"   ),
            axis.title.y = element_text( size = y.title , angle = 90 ) ,
            axis.text.y = element_text( size = y.text ,hjust = 1 ) ,
            aspect.ratio = point2heatmap * aspect.ratio * ( length(unique( margin_y$Y )) / length(unique( margin_x$X ))  )   ,
            plot.margin = margin( 0, 0, 0, 0 ),
            legend.position = "right",
            legend.direction = legend.direction,
            plot.title = element_text( hjust = 0.5 , size = plot.title   )
          )
        
        
        #sender
        right_ratio <- 1 / point2heatmap
        colnames(margin_y)[1] <- 'X'
        p_right <- ggplot(margin_y, aes(x = total, y = X)) +
          geom_point( aes(size = total , alpha = total  ) , color = colors$high  ) + scale_size( range = point.size )+
          scale_y_discrete( limits = cells[ cells %in% heatmap_data2$Y  ]   ) +
          labs(  x = 'Total CCCs (Sender)' )+
          theme_classic()+
          guides( color="none" ,
                  size = "none" ,
                  alpha = "none"
          )+
          theme(
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major.x = element_line( linetype = "dashed" ,color ="grey80"   ),
            axis.title.x = element_text( size = x.title , hjust = 0.5 ) ,
            axis.text.x = element_text( size = x.text ,hjust = 0.5 ) ,
            aspect.ratio = right_ratio  ,
            plot.margin = margin( 0, 0, 0, 0),
            legend.position = "right",
            legend.direction = legend.direction
          )
        
        #
        up_p <- heatmap_main2 %>%  insert_top(p_top, height= point2heatmap ) %>% insert_right( p_right )
        print( up_p )
        return( up_p )
      }
      
      
      #down
      if( trend == 'Down' ){
        heatmap_data2 <- heatmap_data[ heatmap_data$trend == 'DOWN' ,  ]
        heatmap_data_test <- heatmap_data2[ heatmap_data2$value != 0 ,  ]
        if( length( unique( heatmap_data_test$X ) ) < 2 | length( unique( heatmap_data_test$Y ) ) < 2 ){
          stop(simpleError( 'No data were retained. Please adjust the filtering thresholds.' )  )
        }
        heatmap_data2 <- heatmap_data2[ heatmap_data2$Y %in% heatmap_data_test$Y   ,   ]
        heatmap_data2 <- heatmap_data2[ heatmap_data2$X %in% heatmap_data_test$X   ,   ]
        
        
        heatmap_main2 <- ggplot( heatmap_data2, ) +
          geom_tile( aes( x = X , y = Y , fill =  value  ), color = border.color, size = border.size ) +
          scale_fill_gradient2(low = colors$mid, high = colors$low ) +
          scale_x_discrete( expand = c( 0 , 0 ) ,  limits = cells[ cells %in% heatmap_data2$X  ] ) +
          scale_y_discrete( expand = c( 0 , 0 ) , limits = cells[ cells %in% heatmap_data2$Y  ]  ) +
          labs( x = "Receiver", y = "Sender") +
          coord_fixed( ratio = aspect.ratio   ) +
          theme_void(   )+
          guides( fill=guide_colorbar( title = 'No. of CCCs'  )  )+
          theme(
            text = element_text( color = 'black'  ),
            plot.title = element_text( size = plot.title , hjust = 0.5  ),
            axis.title.x = element_text(  size = x.title , hjust = 0.5 ),
            axis.title.y = element_text( size = y.title , angle = 90 ) ,
            axis.text.x = element_text(  size = x.text , vjust = 0.5 ,angle = 90 ,hjust = 1  ),
            axis.text.y = element_text( size = y.text ,hjust = 1 ) ,
            plot.margin = margin( 1, 1, 1, 1),
            legend.position = "left",
            legend.direction = legend.direction
          )
        
        #
        margin_x <- heatmap_data2 %>% group_by(X, trend) %>% summarise(total = sum(value))
        margin_y <- heatmap_data2 %>% group_by(Y, trend) %>% summarise(total = sum(value))
        
        #receiver
        p_top <- ggplot(margin_x, aes(x = X, y = total)) +
          geom_point( aes(size = total , alpha = total  ) , color = colors$low  ) + scale_size( range = point.size )+
          scale_x_discrete( limits = cells[ cells %in% heatmap_data2$X  ]   ) +
          labs(  y = 'Total CCCs (Receiver)' ,
                 title = paste0( title,' (Down)' )  )+
          theme_classic()+
          guides( color="none" ,
                  size = "none"  ,
                  alpha = "none"
          )+
          theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major.y = element_line( linetype = "dashed" ,color ="grey80"   ),
            axis.title.y = element_text( size = y.title , angle = 90 ) ,
            axis.text.y = element_text( size = y.text ,hjust = 1 ) ,
            aspect.ratio = point2heatmap * aspect.ratio * ( length(unique( margin_y$Y )) / length(unique( margin_x$X ))  ) ,
            plot.margin = margin( 0, 0, 0, 0 ),
            legend.position = "right",
            legend.direction = legend.direction,
            plot.title = element_text( hjust = 0.5 , size = plot.title   )
          )
        
        
        #sender
        right_ratio <- 1 / point2heatmap
        colnames(margin_y)[1] <- 'X'
        p_right <- ggplot(margin_y, aes(x = total, y = X)) +
          geom_point( aes(size = total , alpha = total  ) , color = colors$low  ) + scale_size( range = point.size )+
          scale_y_discrete( limits = cells[ cells %in% heatmap_data2$Y  ]   ) +
          labs(  x = 'Total CCCs (Sender)' )+
          theme_classic()+
          guides( color="none" ,
                  size = "none" ,
                  alpha = "none"
          )+
          theme(
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major.x = element_line( linetype = "dashed" ,color ="grey80"   ),
            axis.title.x = element_text( size = x.title , hjust = 0.5 ) ,
            axis.text.x = element_text( size = x.text ,hjust = 0.5 ) ,
            aspect.ratio = right_ratio   ,
            plot.margin = margin( 0, 0, 0, 0),
            legend.position = "right",
            legend.direction = legend.direction
          )
        
        #
        down_p <- heatmap_main2 %>%  insert_top(p_top, height= point2heatmap ) %>% insert_right( p_right )
        print( down_p  )
        return( down_p )
        
      }
      
    }
    
    
  }else{
    ######################################################LR
    
    #
    cells <- unique( c( as.character( plot_data$L ) ,  as.character( plot_data$r   ) ) )
    cells <- sort(cells ,decreasing = F )
    heatmap_data <- expand.grid( unique(plot_data$L)  , unique( plot_data$r  )) %>% data.frame()
    colnames(heatmap_data) <- c('source','target')
    heatmap_data <- lapply(1:nrow( heatmap_data  ), function(x){
      #
      d <- heatmap_data[x,] %>% unlist()
      sdata <- subset( plot_data , L == d[[1]] & r == d[[2]]   )
      op <- data.table( target = d[[2]] , source = d[[1]]  ,
                        number = c(
                          nrow( sdata[ sdata$max.group == max.group ,  ]   ),
                          nrow( sdata[ sdata$max.group != max.group ,  ]   )
                        ),
                        trend = c('UP','DOWN')
      )
      return( op )
      
      #
    }) %>% rbindlist()
    colnames(heatmap_data) <- c('X','Y','value','trend')
    
    
    #X ,receptor
    #Y ,ligand
    
    #
    heatmap_data$X <- as.character( heatmap_data$X )
    heatmap_data$Y <- as.character( heatmap_data$Y )
    
    #
    heatmap_data_test <- heatmap_data[ heatmap_data$value != 0 ,  ]
    if( length( unique( heatmap_data_test$X ) ) < 2 | length( unique( heatmap_data_test$Y ) ) < 2 ){
      stop(simpleError( 'No data were retained. Please adjust the filtering thresholds.' )  )
    }
    heatmap_data <- heatmap_data[ heatmap_data$Y %in% heatmap_data_test$Y   ,   ]
    heatmap_data <- heatmap_data[ heatmap_data$X %in% heatmap_data_test$X   ,   ]
    
    
    #filter
    if( !is.null(top.up) ){
      #
      a <- table( heatmap_data$Y[ heatmap_data$trend == 'UP' & heatmap_data$value != 0   ] )
      if( top.up > length(a) ){ top.up <- length(a) }
      c1 <- names( a[ a >= quantile( a , top.up / length(a) )] )
      #
      a <- table( heatmap_data$X[ heatmap_data$trend == 'UP' & heatmap_data$value != 0   ] )
      if( top.up > length(a) ){ top.up <- length(a) }
      c2 <- names( a[ a >= quantile( a , top.up / length(a) )] )
      
      #
      heatmap_data <- heatmap_data[  heatmap_data$X %in% c2 & heatmap_data$Y %in% c1   , ]
    }
    if( !is.null(top.down) ){
      #
      a <- table( heatmap_data$Y[ heatmap_data$trend == 'DOWN' & heatmap_data$value != 0   ] )
      if( top.down > length(a) ){ top.down <- length(a) }
      c1 <- names( a[ a >= quantile( a , top.down / length(a) )] )
      #
      a <- table( heatmap_data$X[ heatmap_data$trend == 'DOWN' & heatmap_data$value != 0   ] )
      if( top.down > length(a) ){ top.down <- length(a) }
      c2 <- names( a[ a >= quantile( a , top.down / length(a) )] )
      
      #
      heatmap_data <- heatmap_data[  heatmap_data$X %in% c2 & heatmap_data$Y %in% c1   , ]
    }
    if( !is.null(top.all) ){
      #
      a <- table( heatmap_data$Y[  heatmap_data$value != 0   ] )
      if( top.all > length(a) ){ top.all <- length(a) }
      c1 <- names( a[ a >= quantile( a , top.all / length(a) )] )
      #
      a <- table( heatmap_data$X[  heatmap_data$value != 0   ] )
      if( top.all > length(a) ){ top.all <- length(a) }
      c2 <- names( a[ a >= quantile( a , top.all / length(a) )] )
      
      #
      heatmap_data <- heatmap_data[  heatmap_data$X %in% c2 & heatmap_data$Y %in% c1   , ]
    }
    
    
    #
    if( isTRUE( merge )  ){
      #
      get_triangles <- function(df) {
        #
        df$x_num <- as.numeric(as.factor(df$X))
        df$y_num <- as.numeric(as.factor(df$Y))
        
        up_tri <- df %>% filter(trend == "UP") %>% rowwise() %>%
          do(data.frame(
            X = .$X, Y = .$Y, value = .$value, trend = .$trend,
            pos_x = c(.$x_num - 0.5, .$x_num + 0.5, .$x_num + 0.5),
            pos_y = c(.$y_num + 0.5, .$y_num + 0.5, .$y_num - 0.5),
            group = paste0(.$X, .$Y, "UP")
          ))
        
        down_tri <- df %>% filter(trend == "DOWN") %>% rowwise() %>%
          do(data.frame(
            X = .$X, Y = .$Y, value = .$value, trend = .$trend,
            pos_x = c(.$x_num - 0.5, .$x_num - 0.5, .$x_num + 0.5),
            pos_y = c(.$y_num + 0.5, .$y_num - 0.5, .$y_num - 0.5),
            group = paste0(.$X, .$Y, "DOWN")
          ))
        
        return(rbind(up_tri, down_tri))
      }
      #
      triangles_data <- get_triangles(heatmap_data)
      triangles_data$value2 <- triangles_data$value
      triangles_data$value2[ which( triangles_data$trend == 'DOWN'  )  ] <- -triangles_data$value2[ which( triangles_data$trend == 'DOWN'  )  ]
      
      #
      heatmap_main <- ggplot(triangles_data, aes(x = pos_x, y = pos_y, fill = value2, group = group)) +
        geom_polygon(color = border.color, size = border.size ) +
        scale_fill_gradient2(low = colors$low, mid = colors$mid, high = colors$high, midpoint = colors$midpoint ) +
        scale_x_discrete( expand = c( 0 , 0 ) ,  limits = cells[ cells %in% heatmap_data$X  ] ) +
        scale_y_discrete( expand = c( 0 , 0 ) , limits = cells[ cells %in% heatmap_data$Y  ]  ) +
        labs( x = "Receptor", y = "Ligand") +
        coord_fixed( ratio = aspect.ratio   ) +
        theme_void(   )+
        guides( fill=guide_colorbar( title = 'No. of CCCs'  )  )+
        theme(
          text = element_text( color = 'black'  ),
          plot.title = element_text( size = plot.title , hjust = 0.5  ),
          axis.title.x = element_text(  size = x.title , hjust = 0.5 ),
          axis.title.y = element_text( size = y.title , angle = 90 ) ,
          axis.text.x = element_text(  size = x.text , vjust = 0.5 ,angle = 90 ,hjust = 1  ),
          axis.text.y = element_text( size = y.text ,hjust = 1 ) ,
          plot.margin = margin( 1, 1, 1, 1),
          legend.position = "right",
          legend.direction = legend.direction
        )
      
      #
      margin_x <- heatmap_data %>% group_by(X, trend) %>% summarise(total = sum(value))
      margin_x$alpha <- ifelse( margin_x$trend == 'UP', 1 , 0.9 )
      margin_y <- heatmap_data %>% group_by(Y, trend) %>% summarise(total = sum(value))
      margin_y$alpha <- ifelse( margin_y$trend == 'UP', 1 , 0.9  )
      
      #receptor
      p_top <- ggplot(margin_x, aes(x = X, y = total, color = trend)) +
        geom_line(aes(group = X), color = "black" ) +
        geom_point( aes(size = total , alpha = alpha  ) ) + scale_size( range = point.size )+
        scale_alpha( range = c( 1,1 ) )+
        scale_x_discrete( limits = cells[ cells %in% heatmap_data$X  ]   ) +
        scale_color_manual(values = c("UP" = colors$high, "DOWN" = colors$low)) +
        labs(  y = 'Total CCCs (Receptor)' ,
               title = title  )+
        theme_classic()+
        guides( color=guide_legend( title = 'Trend' , override.aes = list( size  = max( point.size ) ) ) ,
                size = "none" ,
                alpha = "none"
        )+
        theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.y = element_line( linetype = "dashed" ,color ="grey80"   ),
          axis.title.y = element_text( size = y.title , angle = 90 ) ,
          axis.text.y = element_text( size = y.text ,hjust = 1 ) ,
          aspect.ratio = point2heatmap * aspect.ratio * ( length(unique( margin_y$Y )) / length(unique( margin_x$X ))  )   ,
          plot.margin = margin(0, 0, 0, 0 ),
          legend.position = "right",
          legend.direction = legend.direction,
          plot.title = element_text( hjust = 0.5 , size = plot.title   )
        )
      
      
      #ligand
      right_ratio <- 1 / point2heatmap
      colnames(margin_y)[1] <- 'X'
      p_right <- ggplot(margin_y, aes(x = total, y = X, color = trend)) +
        geom_line(aes(group = X), color = "black" ) +
        geom_point( aes(size = total , alpha = alpha  ) ) + scale_size( range = point.size )+
        scale_alpha( range = c( 1,1 ) )+
        scale_y_discrete( limits = cells[ cells %in% heatmap_data$Y  ]  ) +
        scale_color_manual(values = c("UP" = colors$high, "DOWN" = colors$low)) +
        labs(  x = 'Total CCCs (Ligand)' )+
        theme(legend.position = "none") +
        theme_classic()+
        guides( color=guide_legend( title = 'Trend' , override.aes = list( size  = max( point.size ) ) ) ,
                size = "none" ,
                alpha = "none"
        )+
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.x = element_line( linetype = "dashed" ,color ="grey80"   ),
          axis.title.x = element_text( size = x.title , hjust = 0.5 ) ,
          axis.text.x = element_text( size = x.text ,hjust = 0.5 ) ,
          aspect.ratio = right_ratio  ,
          plot.margin = margin( 0, 0, 0, 0),
          legend.position = "right",
          legend.direction = legend.direction
        )
      
      #
      final_p <- heatmap_main %>%  insert_top(p_top,height = point2heatmap ) %>% insert_right( p_right )
      
      #
      print( final_p  )
      return( final_p )
      
    }else{
      #split,UP,DOWN
      
      if( trend == 'Up' ){
        heatmap_data2 <- heatmap_data[ heatmap_data$trend == 'UP' ,  ]
        heatmap_data_test <- heatmap_data2[ heatmap_data2$value != 0 ,  ]
        if( length( unique( heatmap_data_test$X ) ) < 2 | length( unique( heatmap_data_test$Y ) ) < 2 ){
          stop(simpleError( 'No data were retained. Please adjust the filtering thresholds.' )  )
        }
        heatmap_data2 <- heatmap_data2[ heatmap_data2$Y %in% heatmap_data_test$Y   ,   ]
        heatmap_data2 <- heatmap_data2[ heatmap_data2$X %in% heatmap_data_test$X   ,   ]
        
        
        heatmap_main2 <- ggplot( heatmap_data2, ) +
          geom_tile( aes( x = X , y = Y , fill =  value  ), color = border.color, size = border.size ) +
          scale_fill_gradient2(low = colors$mid, high = colors$high ) +
          scale_x_discrete( expand = c( 0 , 0 ) ,  limits = cells[ cells %in% heatmap_data2$X  ] ) +
          scale_y_discrete( expand = c( 0 , 0 ) , limits = cells[ cells %in% heatmap_data2$Y  ]  ) +
          labs( x = "Receptor", y = "Ligand") +
          coord_fixed( ratio = aspect.ratio   ) +
          theme_void(   )+
          guides( fill=guide_colorbar( title = 'No. of CCCs'  )  )+
          theme(
            text = element_text( color = 'black'  ),
            plot.title = element_text( size = plot.title , hjust = 0.5  ),
            axis.title.x = element_text(  size = x.title , hjust = 0.5 ),
            axis.title.y = element_text( size = y.title , angle = 90 ) ,
            axis.text.x = element_text(  size = x.text , vjust = 0.5 ,angle = 90 ,hjust = 1  ),
            axis.text.y = element_text( size = y.text ,hjust = 1 ) ,
            plot.margin = margin( 1, 1, 1, 1),
            legend.position = "left",
            legend.direction = legend.direction
          )
        
        #
        margin_x <- heatmap_data2 %>% group_by(X, trend) %>% summarise(total = sum(value))
        margin_y <- heatmap_data2 %>% group_by(Y, trend) %>% summarise(total = sum(value))
        
        #receptor
        p_top <- ggplot(margin_x, aes(x = X, y = total)) +
          geom_point( aes(size = total , alpha = total  ) , color = colors$high  ) + scale_size( range = point.size )+
          scale_x_discrete( limits = cells[ cells %in% heatmap_data2$X  ]   ) +
          labs(  y = 'Total CCCs (Receptor)' ,
                 title = paste0( title,' (Up)' )  )+
          theme_classic()+
          guides( color="none" ,
                  size = "none"  ,
                  alpha = "none"
          )+
          theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major.y = element_line( linetype = "dashed" ,color ="grey80"   ),
            axis.title.y = element_text( size = y.title , angle = 90 ) ,
            axis.text.y = element_text( size = y.text ,hjust = 1 ) ,
            aspect.ratio = point2heatmap * aspect.ratio * ( length(unique( margin_y$Y )) / length(unique( margin_x$X ))  ) ,
            plot.margin = margin( 0, 0, 0, 0 ),
            legend.position = "right",
            legend.direction = legend.direction,
            plot.title = element_text( hjust = 0.5 , size = plot.title   )
          )
        
        
        #ligand
        right_ratio <- 1 / point2heatmap
        colnames(margin_y)[1] <- 'X'
        p_right <- ggplot(margin_y, aes(x = total, y = X)) +
          geom_point( aes(size = total , alpha = total  ) , color = colors$high  ) + scale_size( range = point.size )+
          scale_y_discrete( limits = cells[ cells %in% heatmap_data2$Y  ]   ) +
          labs(  x = 'Total CCCs (Ligand)' )+
          theme_classic()+
          guides( color="none" ,
                  size = "none" ,
                  alpha = "none"
          )+
          theme(
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major.x = element_line( linetype = "dashed" ,color ="grey80"   ),
            axis.title.x = element_text( size = x.title , hjust = 0.5 ) ,
            axis.text.x = element_text( size = x.text ,hjust = 0.5 ) ,
            aspect.ratio = right_ratio  ,
            plot.margin = margin( 0, 0, 0, 0),
            legend.position = "right",
            legend.direction = legend.direction
          )
        
        #
        up_p <- heatmap_main2 %>%  insert_top(p_top, height =  point2heatmap ) %>% insert_right( p_right )
        print( up_p )
        return( up_p )
      }
      
      
      #down
      if( trend == 'Down' ){
        heatmap_data2 <- heatmap_data[ heatmap_data$trend == 'DOWN' ,  ]
        heatmap_data_test <- heatmap_data2[ heatmap_data2$value != 0 ,  ]
        if( length( unique( heatmap_data_test$X ) ) < 2 | length( unique( heatmap_data_test$Y ) ) < 2 ){
          stop(simpleError( 'No data were retained. Please adjust the filtering thresholds.' )  )
        }
        heatmap_data2 <- heatmap_data2[ heatmap_data2$Y %in% heatmap_data_test$Y   ,   ]
        heatmap_data2 <- heatmap_data2[ heatmap_data2$X %in% heatmap_data_test$X   ,   ]
        
        heatmap_main2 <- ggplot( heatmap_data2 ) +
          geom_tile( aes( x = X , y = Y , fill =  value  ), color = border.color, size = border.size ) +
          scale_fill_gradient2(low = colors$mid, high = colors$low ) +
          scale_x_discrete( expand = c( 0 , 0 ) ,  limits = cells[ cells %in% heatmap_data2$X  ] ) +
          scale_y_discrete( expand = c( 0 , 0 ) , limits = cells[ cells %in% heatmap_data2$Y  ]  ) +
          labs( x = "Receptor", y = "Ligand") +
          coord_fixed( ratio = aspect.ratio   ) +
          theme_void(   )+
          guides( fill=guide_colorbar( title = 'No. of CCCs'  )  )+
          theme(
            text = element_text( color = 'black'  ),
            plot.title = element_text( size = plot.title , hjust = 0.5  ),
            axis.title.x = element_text(  size = x.title , hjust = 0.5 ),
            axis.title.y = element_text( size = y.title , angle = 90 ) ,
            axis.text.x = element_text(  size = x.text , vjust = 0.5 ,angle = 90 ,hjust = 1  ),
            axis.text.y = element_text( size = y.text ,hjust = 1 ) ,
            plot.margin = margin( 1, 1, 1, 1),
            legend.position = "left",
            legend.direction = legend.direction
          )
        
        #
        margin_x <- heatmap_data2 %>% group_by(X, trend) %>% summarise(total = sum(value))
        margin_y <- heatmap_data2 %>% group_by(Y, trend) %>% summarise(total = sum(value))
        
        #receptor
        p_top <- ggplot(margin_x, aes(x = X, y = total)) +
          geom_point( aes(size = total , alpha = total  ) , color = colors$low  ) + scale_size( range = point.size )+
          scale_x_discrete( limits = cells[ cells %in% heatmap_data2$X  ]   ) +
          labs(  y = 'Total CCCs (Receptor)' ,
                 title = paste0( title,' (Down)' )  )+
          theme_classic()+
          guides( color="none" ,
                  size = "none"  ,
                  alpha = "none"
          )+
          theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major.y = element_line( linetype = "dashed" ,color ="grey80"   ),
            axis.title.y = element_text( size = y.title , angle = 90 ) ,
            axis.text.y = element_text( size = y.text ,hjust = 1 ) ,
            aspect.ratio = point2heatmap * aspect.ratio * ( length(unique( margin_y$Y )) / length(unique( margin_x$X ))  ) ,
            plot.margin = margin( 0, 0, 0, 0 ),
            legend.position = "right",
            legend.direction = legend.direction,
            plot.title = element_text( hjust = 0.5 , size = plot.title   )
          )
        
        
        #ligand
        right_ratio <- 1 / point2heatmap
        colnames(margin_y)[1] <- 'X'
        p_right <- ggplot(margin_y, aes(x = total, y = X)) +
          geom_point( aes(size = total , alpha = total  ) , color = colors$low  ) + scale_size( range = point.size )+
          scale_y_discrete( limits = cells[ cells %in% heatmap_data2$Y  ]   ) +
          labs(  x = 'Total CCCs (Ligand)' )+
          theme_classic()+
          guides( color="none" ,
                  size = "none" ,
                  alpha = "none"
          )+
          theme(
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major.x = element_line( linetype = "dashed" ,color ="grey80"   ),
            axis.title.x = element_text( size = x.title , hjust = 0.5 ) ,
            axis.text.x = element_text( size = x.text ,hjust = 0.5 ) ,
            aspect.ratio = right_ratio  ,
            plot.margin = margin( 0, 0, 0, 0),
            legend.position = "right",
            legend.direction = legend.direction
          )
        
        #
        down_p <- heatmap_main2 %>%  insert_top(p_top, height = point2heatmap ) %>% insert_right( p_right )
        print( down_p  )
        return( down_p )
        
      }
    }
  }
}

#
#' Heatmap of sender-receiver
#'
#' @description
#' Plot a heatmap for the results of the categorical variable.
#'
#' @param res The object returned by the filterCCC function.
#' @param max.group Which group is treated as the communication-enhanced group.
#' @param merge Whether upregulated and downregulated CCC events are displayed in a single plot.
#' @param trend When merge = FALSE, selects upregulated or downregulated CCC events in max.group. Options are Up and Down.
#' @param data.used The data used for plotting. Options are result or filtered.result.
#' @param padj.threshold The threshold for p.adj.
#' @param p.threshold Similar to padj.threshold, this parameter specifies the p-value threshold. When provided, padj.threshold will be ignored.
#' @param ligand Character vector for filtering CCC.
#' @param receptor Character vector for filtering CCC.
#' @param sender Character vector for filtering CCC.
#' @param receiver Character vector for filtering CCC.
#' @param CCC.ID Character vector for filtering CCC.
#' @param top.up Selects the top n cell types with the highest number of upregulated CCC events for plotting.
#' @param top.down Selects the top n cell types with the highest number of downregulated CCC events for plotting.
#' @param top.all Selects the top n cell types with the most CCC events for plotting.
#' @param colors A list specifying the colors.
#' @param border.color Color of grid borders in the plot area.
#' @param border.size Size of grid borders in the plot area.
#' @param x.text The font size of the X-axis labels.
#' @param y.text The font size of the Y-axis labels.
#' @param x.title The font size of the X-axis title.
#' @param y.title The font size of the Y-axis title.
#' @param title Plot title.
#' @param plot.title The font size of the plot title.
#' @param aspect.ratio aspect ratio, expressed as y / x.
#' @param point.size Specifies the range of point sizes.
#' @param point2heatmap Ratio of the point plot to the heatmap.
#' @param legend.direction Legend arrangement, options: 'horizontal' or 'vertical'.
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
#' #4.heatmap
#' plot_heatmap_cell( res = filter.binary$binary.res , max.group = "O" , merge = T , p.threshold = 1 )
#' plot_heatmap_cell( res = filter.binary$binary.res , max.group = "O" , merge = F , trend = 'Up' , p.threshold = 1 )
#' plot_heatmap_cell( res = filter.binary$binary.res , max.group = "O" , merge = F , trend = 'Down' , p.threshold = 1 )
#'
#' plot_heatmap_cell( res = filter.anova$anova.res , max.group = "O" , merge = T , p.threshold = 1 )
#'
#' @export
#'

plot_heatmap_cell <- function(res , max.group = NULL, merge = T,trend = 'Up',
                              data.used = 'filtered.result',padj.threshold =  0.05 , p.threshold =  NULL,
                              ligand = NULL, receptor = NULL, sender = NULL, receiver = NULL ,CCC.ID = NULL,
                              top.up = NULL, top.down = NULL, top.all = NULL,
                              colors =  list( low =  '#002BFF', mid =  "#F7F7F7", high =  '#FF2B00' , midpoint = 0 ),
                              border.color = 'white', border.size = 0.5,
                              x.text = 10 ,y.text = 10 , x.title = 12 , y.title = 12 ,
                              title = NULL, plot.title = 14,
                              aspect.ratio = 1, point.size = c( 3,6 ), point2heatmap = 0.5 ,
                              legend.direction = "horizontal"
){
  #
  fill = 'p.adj' ; threshold = padj.threshold
  if( !is.null( p.threshold  )  ){   fill = 'p' ; threshold = p.threshold    }
  #
  suppressMessages(
    p <- do_heatmap( ccc.res = res , data.used = data.used, fill = fill, threshold =  threshold ,max.group = max.group,
                     colors =  colors,
                     border.color = border.color , border.size = border.size,
                     ligand = ligand, receptor = receptor, sender = sender, receiver = receiver ,CCC.ID = CCC.ID ,
                     legend.direction = legend.direction ,x.text = x.text ,y.text = y.text , x.title = x.title, y.title = y.title ,
                     title = title, plot.title = plot.title,
                     aspect.ratio = aspect.ratio, point.size = point.size, point2heatmap = point2heatmap ,
                     merge = merge,trend = trend,
                     LR = F,
                     top.up = top.up, top.down = top.down, top.all = top.all
    )
  )
  
  return( p )
}



#
#' Heatmap of ligand-receptor
#'
#' @description
#' Plot a heatmap for the results of the categorical variable.
#'
#' @param res The object returned by the filterCCC function.
#' @param max.group Which group is treated as the communication-enhanced group.
#' @param merge Whether upregulated and downregulated CCC events are displayed in a single plot.
#' @param trend When merge = FALSE, selects upregulated or downregulated CCC events in max.group. Options are Up and Down.
#' @param data.used The data used for plotting. Options are result or filtered.result.
#' @param padj.threshold The threshold for p.adj.
#' @param p.threshold Similar to padj.threshold, this parameter specifies the p-value threshold. When provided, padj.threshold will be ignored.
#' @param ligand Character vector for filtering CCC.
#' @param receptor Character vector for filtering CCC.
#' @param sender Character vector for filtering CCC.
#' @param receiver Character vector for filtering CCC.
#' @param CCC.ID Character vector for filtering CCC.
#' @param top.up Selects the top n cell types with the highest number of upregulated CCC events for plotting.
#' @param top.down Selects the top n cell types with the highest number of downregulated CCC events for plotting.
#' @param top.all Selects the top n cell types with the most CCC events for plotting.
#' @param colors A list specifying the colors.
#' @param border.color Color of grid borders in the plot area.
#' @param border.size Size of grid borders in the plot area.
#' @param x.text The font size of the X-axis labels.
#' @param y.text The font size of the Y-axis labels.
#' @param x.title The font size of the X-axis title.
#' @param y.title The font size of the Y-axis title.
#' @param title Plot title
#' @param plot.title The font size of the plot title.
#' @param aspect.ratio aspect ratio, expressed as y / x.
#' @param point.size Specifies the range of point sizes.
#' @param point2heatmap Ratio of the point plot to the heatmap.
#' @param legend.direction Legend arrangement, options: 'horizontal' or 'vertical'.
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
#' #4.heatmap
#' plot_heatmap_lr( res = filter.binary$binary.res , max.group = "O" , merge = T , p.threshold = 1 )
#' plot_heatmap_lr( res = filter.binary$binary.res , max.group = "O" , merge = F , trend = 'Up' , p.threshold = 1 )
#' plot_heatmap_lr( res = filter.binary$binary.res , max.group = "O" , merge = F , trend = 'Down' , p.threshold = 1 )
#'
#' plot_heatmap_lr( res = filter.anova$anova.res , max.group = "O" , merge = T , p.threshold = 1 )
#'
#' @export
#'
plot_heatmap_lr <- function(res , max.group = NULL, merge = T,trend = 'Up',
                            data.used = 'filtered.result',padj.threshold =  0.05 , p.threshold =  NULL,
                            ligand = NULL, receptor = NULL, sender = NULL, receiver = NULL ,CCC.ID = NULL,
                            top.up = NULL, top.down = NULL, top.all = NULL,
                            colors =  list( low =  '#002BFF', mid =  "#F7F7F7", high =  '#FF2B00' , midpoint = 0 ),
                            border.color = 'white', border.size = 0.5,
                            x.text = 10 ,y.text = 10 , x.title = 12 , y.title = 12 ,
                            title = NULL, plot.title = 14,
                            aspect.ratio = 1, point.size = c( 3,6 ), point2heatmap = 0.5 ,
                            legend.direction = "horizontal"
){
  #
  fill = 'p.adj' ; threshold = padj.threshold
  if( !is.null( p.threshold  )  ){   fill = 'p' ; threshold = p.threshold    }
  #
  suppressMessages(
    p <- do_heatmap( ccc.res = res , data.used = data.used, fill = fill, threshold =  threshold ,max.group = max.group,
                     colors =  colors,
                     border.color = border.color , border.size = border.size,
                     ligand = ligand, receptor = receptor, sender = sender, receiver = receiver ,CCC.ID = CCC.ID ,
                     legend.direction = legend.direction ,x.text = x.text ,y.text = y.text , x.title = x.title, y.title = y.title ,
                     title = title, plot.title = plot.title,
                     aspect.ratio = aspect.ratio, point.size = point.size, point2heatmap = point2heatmap ,
                     merge = merge,trend = trend,
                     LR = T,
                     top.up = top.up, top.down = top.down, top.all = top.all
    )
  )
  
  return( p )
}


