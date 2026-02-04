# multiCCC
`multiCCC` is an R package designed to test the association between cell–cell communication (CCC) and phenotypes across multiple contexts.

![Screenshot](https://github.com/cwbLab/multiCCC/blob/main/data/pipeline.png)

## Install
    devtools::install_github("cwbLab/multiCCC")

## Quick Start
`multiCCC` currently supports scRNA-seq data from both human and mouse.

### 1. Evaluate cell–cell communication scores (LRscore)
	data( "CCC.test.data" )
	LRscore <- scoreLR( exp = t( CCC.test.data$exp ) ,
						meta.data = CCC.test.data$meta ,
						LR.species = 'human' ,
						sample = 'orig.ident' ,
						celltype = 'celltype' )

### 2. Assess the relationship between cell–cell communication events and phenotypes
Note: The multiCCC function supports the simultaneous use of the binary.params, anova.column, glm.column, and time.course.params parameters, and analyses will be performed automatically for each application scenario.

	ccc.binary <- multiCCC( data = LRscore , binary.params = list( group = 'Group' , g1 = 'O' , g2 = 'Y'  ) )
	ccc.anova <- multiCCC( data = LRscore , anova.column = 'batch' )
	ccc.glm <- multiCCC( data = LRscore , glm.column = 'weight' )
	ccc.time <- multiCCC( data = LRscore , time.course.params = list( time  = 'time' , replicate  = 'replicate'  ) )
	
	#One-step execution for all scenarios
	ccc.all <- multiCCC( data = LRscore , 
                     binary.params = list( group = 'Group' , g1 = 'O' , g2 = 'Y'  ),
                     anova.column = 'batch',
                     glm.column = 'weight',
                     time.course.params = list( time  = 'time' , replicate  = 'replicate'  )
                     )
	

### 3. Filter CCC events based on marker genes
Note: The filterCCC function simultaneously applies filtering to all analysis results, without requiring separate filtering for different application scenarios.

	filter.binary <- filterCCC(ccc.binary)
	filter.anova <- filterCCC(ccc.anova)
	filter.glm <- filterCCC(ccc.glm )
	filter.time <- filterCCC(ccc.time)
	
	#One-step execution for all scenarios
	filter.all <- filterCCC(ccc.all)
	

### 4. Visualize

#### 4.1. Circos plot
	plot_circos( res = filter.binary$binary.res , p.threshold = 1 ,LR =T  )
![Screenshot](https://github.com/cwbLab/multiCCC/blob/main/data/test_circos.png)

#### 4.2. Heatmap plot
	plot_heatmap_cell( res = filter.binary$binary.res , 
					   max.group = "O" , merge = F ,trend = "Up",
					   p.threshold = 1 ,x.title = 8, y.title = 8 )
	plot_heatmap_cell( res = filter.binary$binary.res , 
					   max.group = "O" , merge = F ,trend = "Down",
					   p.threshold = 1 ,x.title = 8, y.title = 8)
![Screenshot](https://github.com/cwbLab/multiCCC/blob/main/data/test_heatmap_cell.png)

	plot_heatmap_lr( res = filter.binary$binary.res , 
					 max.group = "O" , merge = T , 
					 p.threshold = 1 ,aspect.ratio = 0.6 )
![Screenshot](https://github.com/cwbLab/multiCCC/blob/main/data/test_heatmap_lr.png)

#### 4.3. Network plot
	plot_network( res = filter.binary$binary.res , p.threshold = 1 , merge = F ,
				  node.fill = "#FFFF3E",link.color = c("#E61A33", "#1AB2FF"),
				  link.alpha = 0.5 )
![Screenshot](https://github.com/cwbLab/multiCCC/blob/main/data/test_network.png)

#### 4.4. Dot plot
	plot_dot( res = filter.anova$anova.res ,
          threshold =  1 ,strip.text = 12 ,dot.max = 7 , legend.position = "right",
          x.text = 12,y.text = 12 )
![Screenshot](https://github.com/cwbLab/multiCCC/blob/main/data/test_dot.png)

#### 4.5. Line plot
	ccc.glm <- multiCCC( data = LRscore , glm.column = 'time' )
	filter.glm <- filterCCC(ccc.glm ,ligand.filter = NULL,receptor.filter = NULL )
	plot_line( res = filter.glm$glm.res,
			   CCC.ID = "cell1_cell3.ADAM15_ITGAV" , p.adj = F )
![Screenshot](https://github.com/cwbLab/multiCCC/blob/main/data/test_line.png)

## Contact
For any technical questions, please contact Wenbo Chen (wenbo_chen@bjmu.edu.cn).
