# multiCCC
`multiCCC` is an R package designed to test the association between cell–cell communication (CCC) and phenotypes across multiple contexts
![Screenshot](https://github.com/cwbLab/multiCCC/blob/main/data/pipeline.jpg)

## Install
    devtools::install_github("cwbLab/multiCCC")

## Quick Start
`multiCCC` currently supports scRNA-seq data from both human and mouse.

### 1. Evaluate cell–cell communication scores (LRscore)
	data( 'test.data' )
	LRscore <- lr_score( exp = t( test.data$exp ) ,
						 meta.data = test.data$meta ,
						 LR.species = 'mouse' ,
						 sample = 'orig.ident' ,
						 celltype = 'celltype',threads = 10  
						)

### 2. Assess the relationship between cell–cell communication events and phenotypes
	ccc.binary <- multiCCC( data = LRscore , binary.params = list( group = 'Group' , g1 = 'O' , g2 = 'Y'  ) ,threads = 10 )
	ccc.anova <- multiCCC( data = LRscore , anova.column = 'batch' ,threads = 10 )
	ccc.glm <- multiCCC( data = LRscore , glm.column = 'weight' , threads = 10 )
	ccc.time <- multiCCC( data = LRscore , time.course.params = list( time  = 'time' , replicate  = 'replicate'  ) ,threads = 10   )

### 3. Visualize 
#### Dot plot
	plot_dot( ccc.res = ccc.anova$anova.res , ligand = c('Sele' , 'Fgf1' , 'Nrxn1' ,'Nrxn2' ), threshold =  1 ,strip.text = 8 )

#### Line plot
	CCC.ID <- 'Monocyte_Monocyte.C3_Itgb2'
	plot_line( ccc.res = ccc.glm$glm.res, CCC.ID = CCC.ID , p.adj = F )

## Contact
Any technical question please contact Wenbo Chen (cwb528@whu.edu.cn).


