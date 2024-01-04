
###PACKAGES###

library(geomorph)
library(rgl)
library(ape)
library(phytools)
library(geiger)
library(ggplot2)
library(graphics)
library(NbClust)
library(svgViewR)
library(formattable)
library(htmlwidgets)
library(scatterplot3d)
library(RColorBrewer)
library(WGCNA)
library(geometry)
library(dispRity)
library(RRPP)
library(landvR)
library(caper)
library(StereoMorph)
library(Morpho)
library(matlib)
library(geiger)
library(ouch)
library(surface)
library(convevol)


#### COMO ABRIR .csv CON NOMBRES, CLASIFICADORES Y VARIABLES, Y .txt CON CLASIFICADORES ####

Rawcoords_todos<-read.csv("Rawcoords150_todos381.csv",header=T,sep=";",row.names=1)
coords.todos<-arrayspecs(Rawcoords_todos[,1:150],50,3)

Rawcoords_aves<-read.csv("Rawcoords150_aves199.csv",header=T,sep=";",row.names=1)
coords.aves<-arrayspecs(Rawcoords_aves[,1:150],50,3)
Rawcoords_mamos<-read.csv("Rawcoords150_mamos150.csv",header=T,sep=";",row.names=1)
coords.mamos<-arrayspecs(Rawcoords_mamos[,1:150],50,3)
Rawcoords_reptiles<-read.csv("Rawcoords150_reptiles32.csv",header=T,sep=";",row.names=1)
coords.reptiles<-arrayspecs(Rawcoords_reptiles[,1:150],50,3)
Rawcoords_cocos<-read.csv("coords.cocos.csv",header=T,sep=";",row.names=1)  
coords.cocos<-arrayspecs(Rawcoords_cocos[,1:93],31,3) 

Rawcoords_fossils<-read.csv("Rawcoords150_fossils394.csv",header=T,sep=";",row.names=1)
coords.fossils<-arrayspecs(Rawcoords_fossils[,1:150],50,3)


#SIN PREMAXILA
Rawcoords_todos<-read.csv("Rawcoords150_todos381.csv",header=T,sep=";",row.names=1)
coords.todos.sinpremax<-arrayspecs(Rawcoords_todos[,4:150],49,3)

Rawcoords_aves<-read.csv("Rawcoords150_aves199.csv",header=T,sep=";",row.names=1)
coords.aves.sinpremax<-arrayspecs(Rawcoords_aves[,4:150],49,3)
Rawcoords_mamos<-read.csv("Rawcoords150_mamos150.csv",header=T,sep=";",row.names=1)
coords.mamos.sinpremax<-arrayspecs(Rawcoords_mamos[,4:150],49,3)
Rawcoords_reptiles<-read.csv("Rawcoords150_reptiles27.csv",header=T,sep=";",row.names=1)
coords.reptiles.sinpremax<-arrayspecs(Rawcoords_reptiles[,4:150],49,3)    
Rawcoords_fossils<-read.csv("Rawcoords150_fossils394.csv",header=T,sep=";",row.names=1)
coords.fossils.sinpremax<-arrayspecs(Rawcoords_fossils[,4:150],49,3)


#VARIABLES Y CLASIFICADORES
variables.aves<-read.csv("clasif.variables.aves.csv",header=T,sep=";",row.names=1)
variables.mamos<-read.csv("clasif.variables.mamos.csv",header=T,sep=";",row.names=1)
variables.reptiles<-read.csv("clasif.variables.reptiles.csv",header=T,sep=";",row.names=1)  

variables.amniotas <- read.csv("variables.amniotas.csv",header=T,sep=";",row.names=1)


#### COMO ABRIR ARCHIVO DE ARBOLES FILOGENETICOS .nex ####

tree.aves<-read.nexus("BirdsTree_consensus.nex")   #phylo pasada por TreeAnnotator de BEAST
tree.mamos<-read.nexus("MammalsTree_consensus.nex")   #phylo pasada por TreeAnnotator de BEAST
tree.reptiles<-read.nexus("SquamatesTree_consensus.nex")   #phylo pasada por TreeAnnotator de BEAST
tree.crocs<-read.nexus("CrocodileTree.nex")
  intervalTimes <- read.table('IntervalTimes.txt',header=TRUE,row.names=1)
  taxonTimes <- read.table('TaxonInterval.txt',header=TRUE,row.names=1)
    timeList_cocos <- list(intervalTimes,taxonTimes)
      bintrees <- bin_timePaleoPhy(tree.crocs, timeList=timeList_cocos, type = "mbl", vartime = 1, ntrees = 10, nonstoch.bin = FALSE, randres = TRUE, timeres = F,sites = NULL, point.occur = FALSE, add.term = T, inc.term.adj = T, dateTreatment = "firstLast", node.mins = NULL, noisyDrop = TRUE, plot = TRUE)
        writeNexus(bintrees, 'cocos_10_mlb.nex')
          cocos.consensus<-read.nexus('cocos_consensus.nex')
            cocos.consensus[["edge.length"]]
          y <- as.numeric(c(29.401962,62.006324,4.7708891,4.7708891,41.000000,25.777215,25.777215,96.179177)) ##cambiar en función de las longitudes que quieras meterle
          cocos.consensus$edge.length = y
          writeNexus(cocos.consensus, 'CrocsTree_consensus.nex')

tree.amniota <- read.nexus("AmnioteTree.full.nex")


##Para combinar filogenias seleccionando la rama y la posición temporal de los nodos, creando una filo compuesta calibrada

plot(tree_mesozoicbirds, show.node.label = TRUE, font = 1, root.edge = TRUE)
  max(vcv.phylo(tree.aves)) #esto es para conocer la edad máxima del árbol
  max(vcv.phylo(tree.mamos))
  max(vcv.phylo(tree.reptiles))

tree.amniota$edge.length = as.numeric(c(65.5,259.5,3.85,255.65,255.65,325)) #con esto le otorgamos longitudes a las ramas, basadas en registro fosil de los grades clados (amniota, sauria, archosauria)
  tree.amniota.filled <- bind.tree(tree.amniota,tree.aves, interactive = TRUE) #con esto localizamos la rama donde tenemos que incrustar el arbol
  tree.amniota.filled <- bind.tree(tree.amniota,tree.aves, where = 2, position = 108.5918) #ahora metemos el arbol diciendo directamente la rama (2) y la edad del nodo (la hemos sacado antes con el max(vcv.phylo(tree.aves)))
    tree.amniota.filled2 <- bind.tree(tree.amniota.filled,cocos.consensus, where = 3, position = 96.17918)
      tree.amniota.filled3 <- bind.tree(tree.amniota.filled2,tree.reptiles, where = 1, position = 243.2984)
        tree.amniota.filled4 <- bind.tree(tree.amniota.filled3,tree.mamos, where = 4, position = 188.3652)
          plot(tree.amniota.full, show.node.label = TRUE, font = 0.1, root.edge = TRUE)

          tip<-c("Aves", "Mammalia", "Crocodylia", "Squamates") #ahora cortamos las ramas que sobran, las de los grandes clados que hemos utilizado para anclar las filos calibradas
            tree.amniota.full <- drop.tip(tree.amniota.filled4,c("Aves", "Mammalia", "Crocodylia", "Squamates"))  
              writeNexus(tree.amniota.full, 'AmnioteTree.full.nex')

# this two lines check if tree names and taxa names are the same

data.names <- rownames(variables.reptiles)
  drop.taxa <- tree.amniota$tip.label[ ! tree.amniota$tip.label %in% data.names ] 
    drop.taxa #if this line returns: 'character(0)', then you are good to go





#### VARIABLES CONTINUAS EN FILOGENIAS ####
    
tree<-tree.amniota
tree.amniota<-ladderize(tree.amniota)
  x<-variables.amniotas$SPH1.res
    names(x)<-row.names(variables.amniotas)
    
    phylo<-contMap(tree.amniota, 
                   x, 
                   res=50, 
                   lwd=0.5, 
                   legend = 40, 
                   outline=FALSE, 
                   sig=3, 
                   type="fan",
                   plot=FALSE)
    
    obj<-setMap(phylo,colors=c("#D83E39","#D83E39","#E8E37F","#B9D1A1","#B9D1A1","#64B3E0","#64B3E0")) #B9D1A1 otro color intermedio #gradiente de siempre: 
    #obj<-setMap(phylo,colors=c("#1a2a6c","#1a2a6c","#1a2a6c", "#f5f5f5", "#f5f5f5", "#b21f1f","#fdbb2d","#fdbb2d")) colores PHEvsPHY
    #gradiente de fucsia para cerebro
    obj <- setMap(phylo.headbrain.res,colors=c("#fff1f7","#f3d7e5","#e0b0c9","#392D69"))
    #gradiente de azul para cerebro
    obj <- setMap(phylo.PC2,colors=c("#e8e8e8","#c9d4e9","#9eb4e0","#160070"))
    plot(obj, type = "fan", fsize=0.2, legend = 100, outline = FALSE, lwd=3,xlim=c(-1.8,1.8),leg.txt= "Relative brain size over the cranial base (SPH)")
    
    
    ## PHENOGRAMS
    
    x<-variables.amniotas$clivus.foramen.angle
      names(x)<-row.names(variables.amniotas)
        phenogram.clivus.ang<-phenogram(tree.amniota, x, fsize = 0.1)
        phenogram.clivus.ang.uncert<-fancyTree(tree.amniota,type="phenogram95",x=x, fsize = 0.1)
    
    #ARCOSAURIOS
    variables.angulos.arcos<-read.csv("variables.angulos.arcos.csv",header=T,sep=";",row.names=1)
      data.names <- dimnames(variables.angulos.arcos) [[1]]
        drop.taxa <- tree.amniota$tip.label[ ! tree.amniota$tip.label %in% data.names ]
          pruntree.arcos <- drop.tip( tree.amniota , drop.taxa )
    
          x<-variables.angulos.arcos$clivus.foramen.angle
            names(x)<-row.names(variables.angulos.arcos)
              phenogram.clivus.ang<-phenogram(pruntree.arcos, x, fsize = 0.1)
              phenogram.clivus.ang.uncert<-fancyTree(pruntree.arcos,type="phenogram95",x=x, fsize = 0.1)
    
    #MAMOS
    variables.angulos.mamos<-read.csv("variables.angulos.mamos.csv",header=T,sep=";",row.names=1)
      data.names <- dimnames(variables.angulos.mamos) [[1]]
        drop.taxa <- tree.amniota$tip.label[ ! tree.amniota$tip.label %in% data.names ]
          pruntree.mamos <- drop.tip( tree.amniota , drop.taxa )
    
          x<-variables.angulos.mamos$clivus.foramen.angle
            names(x)<-row.names(variables.angulos.mamos)
              phenogram.clivus.ang<-phenogram(pruntree.mamos, x, fsize = 0.1)
              phenogram.clivus.ang.uncert<-fancyTree(pruntree.mamos,type="phenogram95",x=x, fsize = 0.1)
    
#### EXTRACCIÓN DE SUBCLADOS EN UNA FILOGENIA Y CALCULO DE VALORES POR NODO ####    
    
    get.descendants <- function( phy , node ) { extract.clade( phy , node )$tip.label } #extract.clade te genera un arbol a partir de todos los tips de un feterminado nodo, podando todos los demás tips
    subset.rows <- function( data , rownames ) { data[ rownames , ] }
    
    phy <- tree.amniota
    
    ContVar.by.node <- function( data , phy ) {
      nodes.temp <- unique( phy$edge[,1] ) #para generar un vector que contenga los valores de los nodos internos de la filogenia
      descendants.temp <- lapply( nodes.temp , FUN = get.descendants , phy = phy ) #esta línea devuelve una lista de la longitud de x que se llena de los resultados obtenidos a partir de aplicar la función FUN a cada uno de los elementos de x
      names( descendants.temp ) <- nodes.temp			
      node.data.list <- lapply( descendants.temp , FUN = subset.rows , data = data )
      headbrain.temp<-list() 
      
      for (i in 1: length( node.data.list ) ) {
        headbrain.indiv.values <- node.data.list[[i]]
        headbrain.temp[i]<- mean(headbrain.indiv.values)                         
      }
      
      names(headbrain.temp)<-nodes.temp
      headbrain.temp<-unlist(headbrain.temp)
      headbrain.temp
      #modularity.temp.scaled<-modularity.temp / modularity.temp[1] ## activar si queremos que el valor para cada nodo esté escalado y se dé de manera relativa (sobre 1)
      #modularity.temp.scaled
      
    }
    
    headbrain <- variables.amniotas$headbrain.res
      names(headbrain) <- row.names(variables.amniotas)
    headbrain <- as.data.frame(headbrain)
      headbrain.by.node <- ContVar.by.node( data = headbrain , phy = phy )
        headbrain.by.node.mod <- headbrain.by.node + 3
          log.headbrain.by.node.mod <- log(headbrain.by.node.mod)
    
    
    clivus.angle <- variables.amniotas$clivus.foramen.angle
      names(clivus.angle) <- row.names(clivus.angle)
    clivus.angle <- as.data.frame(clivus.angle)
      angle.clivus.by.node <- ContVar.by.node( data = clivus.angle , phy = phy )
    
    
    
    

#### LOOP BUENO PARA REFLEJAR/TRANSFORMAR/EXPORTAR LANDMARKS ####
      
      
    LMs.aves.R<-read.csv('LMs_aves_rawcoords_R.csv', header = T, sep = ";", row.names = 1)
      coords<-arrayspecs(LMs.aves.R[,1:93],31,3)
      
    LMs.aves.L<-read.csv('LMs_aves_rawcoords_L.csv', header = T, sep = ";", row.names = 1)
      coords<-arrayspecs(LMs.aves.L[,1:93],31,3)
      
    LMs.mamos.R<-read.csv('RawCoords_mamos_todos_R.csv', header = T, sep = ";", row.names = 1)
      coords<-arrayspecs(LMs.mamos.R[,1:93],31,3)
      
    LMs.reptiles.R<-read.csv('RawCoords150_reptiles27.csv', header = T, sep = ";", row.names = 1)
      coords<-arrayspecs(LMs.reptiles.R[,1:93],31,3)
      
    LMs.cocos.R<-read.csv("gavialis.csv",header=T,sep=";",row.names=1)  
      coords<-arrayspecs(LMs.cocos.R[,1:93],31,3)
      
    LMs.fossils.R<-read.csv("RawCoords_fossils12.csv",header=T,sep=";",row.names=1)  
      coords<-arrayspecs(LMs.fossils.R[,1:93],31,3)
      
    LMs.fossils.L<-read.csv("RawCoords_fossils1.csv",header=T,sep=";",row.names=1)  
      coords<-arrayspecs(LMs.fossils.L[,1:93],31,3)
      
      plot3d(coords[,,1]); aspect3d("iso")
      
      
      mirrored.coords <- list() #haces una lista vacía donde vas a meter los cráneos mirroreados, puede ser una lista un array una matriz lo que sea.
      
      for (i in 1: length(coords)) { #para probar cada línea, empezar con i = 1 e ir corriendo cada línea. Si funcionan todas, quitar la i = 1 y correr la función
        craneo.individual <- coords[,,i] #coge el elemento del número que vaya tocando en cada iteración, 1 para la primera, 2 para la segunda hasta la 250 vez.
          craneo.mirrored <- mirror2plane(craneo.individual, v1 = craneo.individual[12,], v2 = craneo.individual[15,], v3 = craneo.individual[30,])
            mirrored.coords[[i]]<-craneo.mirrored #metes el craneo mirroreado en una lista nueva
      }
      
      plot3d(mirrored.coords[[1]]); aspect3d("iso") #visualización 3D sin deformación de ejes
      
      
      ## CAMINO ALTERNATIVO PARA COMBINAR EL ARRAY DEL LADO LANDMARKEADO Y EL ARRAY DEL LADO REFLEJADO ##
      
      array.mirrored.coords <- simplify2array(mirrored.coords) #transformamos la lista de cráneos reflejados creada en el loop en un 3D array
        pruned.array.mirrored.coords <- array.mirrored.coords[-c(1,2,12:17,25:26,30:31),,] #borramos las filas de los landmarks repetidos (los mediales)
          concatenated.arrays.with.mirrored.coords <- bindArr(coords, pruned.array.mirrored.coords, along = 1) #concatenamos los dos 3D arrays
            LMs.fossils2 <- two.d.array(concatenated.arrays.with.mirrored.coords) #transformamos el 3D array combinado en una matriz 2D como las que tenemos en excel con las coordenadas de los landmarks para guardarlo así
              write.csv(LMs.fossils2, 'Rawcoords150_fossils1.csv')
                plot3d(concatenated.arrays.with.mirrored.coords[,,1]); aspect3d("iso") #ploteamos para comprobar que han salido bien
  

#### LOOP PARA CALCULAR VOLUMENES A PARTIR DE CONFIGURACIONES GEOMETRICAS ####
  
  #BRAIN
  Rawcoords_brain.reptiles<-read.csv("brain.coords.reptiles.csv",header=T,sep=";",row.names=1)
    coords.brain.reptiles<-arrayspecs(Rawcoords_brain.reptiles[,1:63],21,3)
  
  gpa.reptiles <- gpagen(coords.brain.reptiles)
    coords.brain.reptiles <- gpa.reptiles$coords
  
  
  brain.convhulls <- list() #haces una lista vacía donde vas a meter los convex hull de cada cerebro, puede ser una lista un array una matriz lo que sea.
  
  for (i in 1: length(coords.brain.aves)) { #para probar cada línea, empezar con i = 1 e ir corriendo cada línea. Si funcionan todas, quitar la i = 1 y correr la función
    coords.individual <- coords.brain.aves[,,i] #coge el elemento del número que vaya tocando en cada iteración, 1 para la primera, 2 para la segunda hasta la 250 vez.
      convhull.individual <- convhulln(coords.individual, output.options = TRUE) #función simple
        brain.convhulls[[i]]<-convhull.individual$vol #metes el convex hull de cada cerebro en una lista nueva
  }
  
  brain.convhulls<-as.data.frame(brain.convhulls, row.names = row.names(Rawcoords_brain.aves))
    write.csv(brain.convhulls$V1, "brain.convhulls.csv")
  
    
  ##FUNCION BUENA ALTERNATIVA A LA ANTERIOR 
    
  brain.alphashape <- list() #haces una lista vacía donde vas a meter las alpha shape de cada cerebro, puede ser una lista un array una matriz lo que sea.
  
  for (i in 1: length(coords.brain.reptiles)) { #para probar cada línea, empezar con i = 1 e ir corriendo cada línea. Si funcionan todas, quitar la i = 1 y correr la función
    coords.individual <- coords.brain.reptiles[,,i] #coge el elemento del número que vaya tocando en cada iteración, 1 para la primera, 2 para la segunda hasta la 250 vez.
      alphashape.individual <- ashape3d(coords.individual, alpha = 100, pert = FALSE, eps = 1e-09) #función simple
        vol.ashape.individual <- volume_ashape3d(alphashape.individual, byComponents = FALSE, indexAlpha = 1)
          brain.alphashape[[i]]<-vol.ashape.individual #metes la alpha-shape de cada cerebro en una lista nueva
  }
  
  plot(alphashape.individual, indexAlpha = "all")
    write.csv(brain.alphashape, 'brain.convhull.reptiles.csv')
  
  
      
#### LOOP PARA CALCULAR ANGULOS A PARTIR DE CONFIGURACIONES DE LMs ####

angulos.orientation<-list()

for (i in 1: length(gpa.todos$coords)) { #para probar cada línea, empezar con i = 1 e ir corriendo cada línea. Si funcionan todas, quitar la i = 1 y correr la función
  craneo.individual <- gpa.todos$coords[,,i] #coge el elemento del número que vaya tocando en cada iteración, 1 para la primera, 2 para la segunda hasta la 250 vez.
    vect1 <- matrix(c((gpa.todos$coords[15,,i])-(gpa.todos$coords[13,,i])), nrow = 3)
    vect2 <- matrix(c((gpa.todos$coords[12,,i])-(gpa.todos$coords[13,,i])), nrow = 3)
      angulo.individual <- acos(sum(vect1*vect2)/(sqrt(sum(vect1*vect1))*sqrt(sum(vect2*vect2))))
        angulo.individual <- rad2deg(angulo.individual)
  angulos.orientation[[i]]<-angulo.individual #metes el craneo mirroreado en una lista nueva
}

angulos.orientation <- as.matrix(angulos.orientation)
  write.csv2(angulos.orientation, 'facial.flexion.csv')


#############  
#### GPA ####
#############

gpa.aves<-gpagen(coords.aves)
gpa.mamos<-gpagen(coords.mamos)
gpa.reptiles<-gpagen(coords.reptiles)
gpa.todos<-gpagen(coords.todos)
gpa.fossils<-gpagen(coords.fossils)

#SIN PREMAXILA
gpa.aves.sinpremax<-gpagen(coords.aves.sinpremax)
gpa.mamos.sinpremax<-gpagen(coords.mamos.sinpremax)
gpa.reptiles.sinpremax<-gpagen(coords.reptiles.sinpremax)
gpa.todos.sinpremax<-gpagen(coords.todos.sinpremax)
gpa.fossils.sinpremax<-gpagen(coords.fossils.sinpremax)


#para plotear todo viendo la media
plotAllSpecimens(gpa.todos.sinpremax$coords, mean = TRUE, links = NULL, label = FALSE, 
                 plot.param = list(pt.cex = 0.1, mean.cex = 5)) 

#############
#### PCA ####
#############

#correr antes los geomorph.data.frame de arriba

phyloPCA.todos<-gm.prcomp(gpa.todos$coords, phy = tree.amniota)
  summary(phyloPCA.todos)
    colours.clades.b <- setNames(c( "#FFEC50", "#FF5050", "#FFBA50",  "#517B1B", "#94E161",  "#75F597",
                                "#95E9E3", "#84C0F7", "#8558B1"),levels(variables.amniotas$clades.b))
      plot(phyloPCA.todos, phylo = TRUE, main = "phyloPCA", phylo.par = list(node.labels = FALSE, anc.states = FALSE, node.cex = 0.2, edge.color = "grey", edge.width = 0.5, tip.txt.cex = 0.5))
  phyloPCs<-phyloPCA.todos$x
  phylomorph.amniotas <- phylomorphospace(tree.amniota, phyloPCs[,1:2], colors = colours.clades.b, xlab = "PC1", ylab = "PC2")

plotRefToTarget(phyloPCA.todos$shapes$shapes.comp1$min, phyloPCA.todos$shapes$shapes.comp1$max,
                method = "vector", mag = 0.75, pt.size = 0.5, label = TRUE)
plotRefToTarget(phyloPCA.todos$shapes$shapes.comp2$min, phyloPCA.todos$shapes$shapes.comp2$max,
                method = "vector", mag = 0.75, pt.size = 1, label = TRUE)
plotRefToTarget(phyloPCA.todos$shapes$shapes.comp5$min, phyloPCA.todos$shapes$shapes.comp5$max,
                method = "vector", mag = 0.75, pt.size = 1, label = TRUE)


PCA.fossils <- gm.prcomp(gpa.fossils$coords)
  plot(PCA.fossils, main = "PCA")


### PCA SIN PREMAXILA ###

phyloPCA.todos.sinpmx<-gm.prcomp(gpa.todos.sinpremax$coords, phy = tree.amniota)
  summary(phyloPCA.todos.sinpmx)
    plot(phyloPCA.todos.sinpmx, phylo = TRUE, main = "phyloPCA", phylo.par = list(node.labels = FALSE, anc.states = FALSE, node.cex = 0.2, edge.color = "grey", edge.width = 0.5, tip.txt.cex = 0.5))
  phyloPCs.sinpmx<-phyloPCA.todos.sinpmx$x

phylomorph.PC23<-phylomorphospace(tree.amniota, phyloPCs.sinpmx[,2:3], colors = "grey", label = "horizontal", lwd = 1, xlab = "PC2", ylab = "PC3")


phylomorphospace3d(tree = tree.amniota, phyloPCs.sinpmx[,1:3], label = FALSE, 
                   control = list(box = TRUE), method = "dynamic", angle = -30)

plotRefToTarget(phyloPCA.todos.sinpmx$shapes$shapes.comp1$min, phyloPCA.todos.sinpmx$shapes$shapes.comp1$max, 
                method = "vector", mag = 0.75, label = TRUE)
plotRefToTarget(phyloPCA.todos.sinpmx$shapes$shapes.comp2$min, phyloPCA.todos.sinpmx$shapes$shapes.comp2$max, 
                method = "vector", mag = 0.75, label = TRUE)
plotRefToTarget(phyloPCA.todos.sinpmx$shapes$shapes.comp5$min, phyloPCA.todos.sinpmx$shapes$shapes.comp5$max, 
                method = "vector", mag = 0.75, label = TRUE)


#PCA 3D#

pca<-prcomp(PCs,scale. = TRUE)
  gr<-classifier2.1
    summary(gr)
      pca3d(pca,group=gr,legend="topleft")

phylomorphospace3d(tree = newzombie.tree, scores(PCs)[,1:3], label = TRUE, 
                   control = list(box = TRUE), method = "static", angle = -30)

#para mapear una variable continua sobre la filogenia del morfoespacio y ver su distribución en filo y morfo
x<-variables.zombie$IRE_vol
y<-variables.zombie$brain.res
  names(x)<-row.names(variables.zombie)
    tree<-newzombie.tree
phylo.IREvol<-contMap(tree,x,res=100,lwd=3,legend = 40,outline=FALSE,sig=3,
                      type="phylogram",direction="rightwards",plot=FALSE)
obj<-setMap(phylo.IREvol,colors=c("#64B3E0","#E8E37F","#D83E39"))
  plot(obj, outline = FALSE)

phylomorph_cont_IREvol<-phylomorphospace(obj$tree, PCs[,1:2], colors = obj$cols, lwd = 4, xlab = "PC1", ylab = "PC2")


#FILOMORFOESPACIO 3D (cronofilomorfoespacio) del PCA PLOTEANDO TIEMPO#

cronophylomorphospace<-plotGMPhyloMorphoSpace(newzombie.tree, newzombie.gdf$coords, tip.labels = FALSE,
                                              ancStates=TRUE,plot.param = list(t.bg=as.integer(newzombie.gdf$superorders),
                                                                               l.col="grey",lwd=1),
                                              node.labels = TRUE,zaxis="time", shadow = TRUE)

#Shape predictor#

PC12<-PCA$x[,1:2]
  M <- mshape(newzombie.gpa$coords)
    predshape.node88<-shape.predictor(newzombie.gpa$coords, x = PC12, Intercept = FALSE, pred1 = c(-0.213822,0.5785432))
    predshape.node86<-shape.predictor(newzombie.gpa$coords, x = PC12, Intercept = FALSE, pred1 = c(-0.2243689,0.3807889))
    predshape.node90<-shape.predictor(newzombie.gpa$coords, x = PC12, Intercept = FALSE, pred1 = c(-0.3562053,0.1408467))




df<-data.frame(clades.a = variables.amniotas$clades.a, clades.b = variables.amniotas$clades.b,
               clades.c = variables.amniotas$clades.c, clades.d = variables.amniotas$clades.comp, 
               log.skullCS = log(variables.amniotas$skull.CS.sinpmx), log.cranialbase.CS = log(variables.amniotas$cranialbase.CS),
               log.brainconvhull = log(variables.amniotas$brain.convhull), headbrain.res = variables.amniotas$headbrain.res, 
               headbrain.res.conpmx = variables.amniotas$headbrain.res.conpmx, encephalization = variables.amniotas$encephalization,
               skullbase.res = variables.amniotas$skullbase.res, foramen = variables.amniotas$foramen.angule, base.ang = variables.amniotas$base.angule, 
               flexion = variables.amniotas$flexion.angule, sph.res = variables.amniotas$SPH1.res, 
               clivus.foramen.ang = variables.amniotas$clivus.foramen.angle, facial.ang = variables.amniotas$facial.angle,
               PC1 = phyloPCs.sinpmx[,1], PC2 = phyloPCs.sinpmx[,2], PC3 =phyloPCs.sinpmx[,3], PC4 =phyloPCs.sinpmx[,4],
               RegScores.skullCS.sinpremax = RegScores.skullCS.sinpremax, 
               RegScores.headbrain.sinpremax = RegScores.headbrain.sinpremax,
               RegScores.skullbase.sinpmx = RegScores.skullbase.sinpmx, 
               row.names = row.names(variables.amniotas))


# VIOLIN PLOTS

plot1 <- ggplot(df, aes(x = clades.b,y = headbrain.res.conpmx))
  plot2 <- print(plot1 + geom_violin(trim = FALSE, adjust = .7) + geom_boxplot(width = .1, fill = "black", outlier.colour = NA) + stat_summary(fun = median, geom = "point", fill = "white", shape = 21, size = 2.5))

plot1 <- ggplot(df, aes(x = clades.a, y = headbrain.res))
  plot2 <- print(plot1 + geom_violin(trim = FALSE, adjust = .5) + geom_boxplot(width = .1, fill = "black", outlier.colour = NA) + stat_summary(fun = median, geom = "point", fill = "white", shape = 21, size = 2.5))

  
############################
#### SEÑAL FILOGENÉTICA ####
############################
  
  physignal(gpa.aves$coords, phy = tree.aves, iter = 9999) #K multivariante
  physignal(gpa.aves.sinpremax$coords, phy = tree.aves, iter = 9999)
  physignal(gpa.mamos$coords, phy = tree.mamos, iter = 9999)
  physignal(gpa.mamos.sinpremax$coords, phy = tree.mamos, iter = 9999)
  physignal(gpa.reptiles$coords, phy = tree.reptiles, iter = 9999)
  physignal(gpa.reptiles.sinpremax$coords, phy = tree.reptiles, iter = 9999)
  
  phylosig(phyloPCA.reptiles$x[,1], tree = tree.reptiles, method ="K", test = T, nsim = 9999)
  phylosig(phyloPCA.reptiles$x[,2], tree = tree.reptiles, method ="K", test = T, nsim = 9999)
  phylosig(phyloPCA.reptiles$x[,3], tree = tree.reptiles, method ="K", test = T, nsim = 9999)
  phylosig(phyloPCA.reptiles.sinpremax$x[,1], tree = tree.reptiles, method ="K", test = T, nsim = 9999)
  phylosig(phyloPCA.reptiles.sinpremax$x[,2], tree = tree.reptiles, method ="K", test = T, nsim = 9999)
  phylosig(phyloPCA.reptiles.sinpremax$x[,3], tree = tree.reptiles, method ="K", test = T, nsim = 9999)
  
  headbrain.res.aves<-variables.aves$headbrain.res
  names(headbrain.res.aves)<-row.names(variables.aves)
  phylosig(headbrain.res.aves, tree = tree.aves, method = "K", test = T, nsim = 9999)
  
  headbrain.res.mamos<-variables.mamos$headbrain.res
  names(headbrain.res.mamos)<-row.names(variables.mamos)
  phylosig(headbrain.res.mamos, tree = tree.mamos, method = "K", test = T, nsim = 9999) 
  
  
  physignal(gpa.todos$coords, phy = tree.amniota, iter = 999) #K multivariante
  physignal(gpa.todos.sinpremax$coords, phy = tree.amniota, iter = 999)
  
  phylosig(phyloPCA.todos$x[,1], tree = tree.amniota, method ="K", test = T, nsim = 9999)
  phylosig(phyloPCA.todos$x[,2], tree = tree.amniota, method ="K", test = T, nsim = 9999)
  phylosig(phyloPCA.todos$x[,3], tree = tree.amniota, method ="K", test = T, nsim = 9999)
  phylosig(phyloPCA.todos$x[,4], tree = tree.amniota, method ="K", test = T, nsim = 9999)
  phylosig(phyloPCA.todos.sinpmx$x[,1], tree = tree.amniota, method ="K", test = T, nsim = 9999)
  phylosig(phyloPCA.todos.sinpmx$x[,2], tree = tree.amniota, method ="K", test = T, nsim = 9999)
  phylosig(phyloPCA.todos.sinpmx$x[,3], tree = tree.amniota, method ="K", test = T, nsim = 9999)
  phylosig(phyloPCA.todos.sinpmx$x[,4], tree = tree.amniota, method ="K", test = T, nsim = 9999)
  
  foramen<-variables.amniotas$foramen.angule
  names(foramen)<-row.names(variables.amniotas)
  phylosig(foramen, tree = tree.amniota, method = "K", test = T, nsim = 9999)
  
  base.ang<-variables.amniotas$base.angule
  names(base.ang)<-row.names(variables.amniotas)
  phylosig(base.ang, tree = tree.amniota, method = "K", test = T, nsim = 9999)
  
  flexion<-variables.amniotas$flexion.angule
  names(flexion)<-row.names(variables.amniotas)
  phylosig(flexion, tree = tree.amniota, method = "K", test = T, nsim = 9999)
  
  skull.CS.sinpmx<-variables.amniotas$skull.CS.sinpmx
  names(skull.CS.sinpmx)<-row.names(variables.amniotas)
  phylosig(skull.CS.sinpmx, tree = tree.amniota, method = "K", test = T, nsim = 9999)
  
  brainconvhull<-variables.amniotas$brain.convhull
  names(brainconvhull)<-row.names(variables.amniotas)
  phylosig(log(brainconvhull), tree = tree.amniota, method = "K", test = T, nsim = 9999)
  
  headbrain<-variables.amniotas$headbrain.res
  names(headbrain)<-row.names(variables.amniotas)
  phylosig(headbrain, tree = tree.amniota, method = "K", test = T, nsim = 9999)
  
  SPH<-variables.amniotas$SPH1.res
  names(SPH)<-row.names(variables.amniotas)
  phylosig(SPH, tree = tree.amniota, method = "K", test = T, nsim = 9999)
  
  
  LD3.LD2<-variables.newzombie$LD3.LD2
  names(LD3.LD2)<-row.names(variables.newzombie)
  LD3.LD2 <- LD3.LD2[ , ,  unlist( pruntree$tip.label ) ] #subsetting by another element re-orders the data with the element's order
  phylosig(variables.newzombie$LD3.LD2, tree = pruntree, method ="K", test = T, nsim = 9999) #k univariante
  
  
#### RIDGELINE PLOTS ####

library(ggridges)
library(ggplot2)
library(viridis)
library(hrbrthemes)

plot1 <- ggplot(df, aes(x = headbrain.res,y = clades.b, fill = stat(x))) + geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.01) +
  scale_fill_gradientn(#colours = c("#1a2a6c","#1a2a6c", "#838bae", "#f5f5f5","#b21f1f","#b21f1f","#fdbb2d"),
    colours = c("#64B3E0","#64B3E0","#64B3E0","#E8E37F","#E8E37F","#E8E37F","#D83E39","#D83E39"),
    values = NULL, space = "Lab", guide = "colourbar", aesthetics = "fill") +
  labs(title = 'encephalization') +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )
plot2 <- print(plot1)



##############################################
#### DISTANCIAS PROCRUSTES Y DISPARIDADES ####
##############################################

#con la función gpagen#

pairwise.ProcDist.matrix<-as.matrix(dist(gpa.todos$coords))
  write.csv(proc.coords.2darray.sinpmx, 'proc.coords.2darray.sinpmx.csv')
    proc.coords.2darray.sinpmx.ordenada <- read.csv('proc.coords.2darray.sinpmx.csv',header=T,sep=";",row.names=1)
    proc.coords.2darray.sinpmx.ordenada <- as.matrix(proc.coords.2darray.sinpmx.ordenada)
      proc.dist.sinpmx <- as.matrix(dist(proc.coords.2darray.sinpmx.ordenada))


proc.dist.sinpmx.list <- t(combn(colnames(proc.dist.sinpmx), 2))
proc.dist.sinpmx.list <- data.frame(proc.dist.sinpmx.list, dist=proc.dist.sinpmx[proc.dist.sinpmx.list])
  write.csv(proc.dist.sinpmx.list, 'proc.dist.sinpmx.list.csv')


#con una función hecha a mano#

mean_shape <- landvR::select.procrustes(newzombie.gpa, selector = mean, factors = groups) #Para seleccionar valores específicos de una superposición procrustes (forma media, mediana, mínima, máxima, etc.)
mean_shape <- geomorph::mshape(newzombie.gpa$coords)

p.dist <- function (a) {
  #### a is an array of ALIGNED Procrusted coordinates
  #### (e.g. after GPA). 
  #### An array is of the form (p x k x n) where
  #### p is number of landmarks,
  #### k is dimensionality of the data (2D or 3D)
  #### n is number of individuals.
  #### Function p.dist calculates and returns procrustes distance(s)
  #### between each individual and mean shape of the data.
  
  p <- dim (a)[3]
  ms <- apply (a, c(1, 2), mean)
  dists <- vector ("numeric", p)
  for(i in 1:p) {
    dists[i] <- sqrt (sum ((a[,,i] - ms)^2)) # this is the formula for PD
  }
  dists
}

dim(newzombie.gpa$coords)
p.dist(newzombie.gpa$coords)


#Otra función para calcular distancias procrustes#

proc.distance <- function(x, mean) {
  return(dist(rbind(as.vector(x), as.vector(mean)), "euclidean"))
}
proc_distances <- unlist(lapply(newzombie.gpa$coords, proc.distance, mean = mean_shape[[1]]))



## DISPARIDADES

gdf.disp <- geomorph.data.frame(gpa.todos, clades.a = variables.amniotas$clades.a,
                                clades.b = variables.amniotas$clades.b, clades.c = variables.amniotas$clades.c)
gdf.disp.sinpmx <- geomorph.data.frame(gpa.todos.sinpremax, clades.a = variables.amniotas$clades.a,
                                       clades.b = variables.amniotas$clades.b, clades.c = variables.amniotas$clades.c)

disparity <- morphol.disparity(coords ~ clades.c, groups = ~ clades.c, data = gdf.disp, iter=999)
  summary(disparity)

disparity.sinpmx <- morphol.disparity(coords ~ clades.c, groups = ~ clades.c, data = gdf.disp.sinpmx, iter=999)
  summary(disparity.sinpmx)


##DISPARITY THROUGH TIME

matrix.coords <- two.d.array(gpa.todos.sinpremax$coords)
  PCs.con <- phyloPCs[,1:25]
  PCs.sin <- phyloPCs.sinpmx[,1:30]
    DTT.all <- dtt(phy = tree.amniota, data = matrix.coords, index = c("avg.sq"), mdi.range = c(325,0), nsim = 1000, CI = 0.95, plot = TRUE)      
    DTT.all.PCs <- dtt(phy = tree.amniota, data = PCs.sin, index = c("avg.sq"), mdi.range = c(325,0), nsim = 1000, CI = 0.95, plot = TRUE)

#por grupos
data.names <- row.names(grouping1.conpmx[[4]])
  drop.taxa <- tree.amniota$tip.label[ ! tree.amniota$tip.label %in% data.names ]
    pruntree.mamos <- drop.tip( tree.amniota , drop.taxa )

coords.birds <- as.matrix(grouping1.conpmx[[1]][,5:154])
  DTT.birds <- dtt(phy = pruntree.birds, data = coords.birds, index = c("avg.sq"), mdi.range = c(0,1), nsim = 1000, CI = 0.95, plot = TRUE)      

coords.crocs <- as.matrix(grouping1.conpmx[[2]][,5:154])
  DTT.crocs <- dtt(phy = pruntree.crocs, data = coords.crocs, index = c("avg.sq"), mdi.range = c(0,1), nsim = 1000, CI = 0.95, plot = TRUE)      

coords.lepids <- as.matrix(grouping1.conpmx[[3]][,5:154])
  DTT.lepids <- dtt(phy = pruntree.lepids, data = coords.lepids, index = c("avg.sq"), mdi.range = c(0,1), nsim = 1000, CI = 0.95, plot = TRUE)      

coords.mamos <- as.matrix(grouping1.conpmx[[4]][,5:154])
  DTT.mamos <- dtt(phy = pruntree.mamos, data = coords.mamos, index = c("avg.sq"), mdi.range = c(0,1), nsim = 1000, CI = 0.95, plot = TRUE)      



##DISPARITY THROUGH PHYLOGENETIC EVOLUTION (NODE BY NODE)  

disparity.by.node.conpmx <- disparity(phy = phy, data = twodarray.proccoords.conpmx[,5:154], index = "avg.sq")
disparity.by.node.sinpmx <- disparity(phy = phy, data = twodarray.proccoords.sinpmx[,5:151], index = "avg.sq")

gpa.coords.temp <- gpa.todos$coords
gpa.coords.temp <- gpa.todos.sinpremax$coords


disparity.phylo <- function( data , phy ) {
  nodes.temp <- unique( phy$edge[,1] )
    descendants.temp <- lapply( nodes.temp , FUN = get.descendants , phy = phy )
      names( descendants.temp ) <- nodes.temp	
      #data<- two.d.array(data)
      node.data.list <- lapply( descendants.temp , FUN = subset.rows , data = data )
        disparity.temp<-list() 
  for (i in 1: length( node.data.list ) ) {
    distances<-dist(node.data.list[[i]], method = "euclidean")
      disparity.temp[i]<- mean(distances)      # loop does the same process as Procrustes variance:
    
  }
  
  names(disparity.temp)<-nodes.temp
    disparity.temp<-unlist(disparity.temp)
      disparity.temp
    disparity.temp.scaled<-disparity.temp / disparity.temp[1]
      disparity.temp.scaled
  
}

disparity.nodes.temp <- disparity.phylo( data = gpa.coords.temp , phy = tree.temp ) #si no funciona, cambiar el gpa.coords.temp a una matriz 2D en lugar de 3D


#simulated disparity data using mvMORPH fully multivariate functions

library(mvMORPH)

fit1 <- mvgls( two.d.array(gpa.coords.temp) ~ 1 , tree = tree.temp , model = "BM" , method = "H&L" ) #primero ajustamos los datos al modelo evolutivo elegido, en este caso BM no direccional
  sim.temp.new <- mvSIM(tree.temp, nsim = 200, model = "BM1", param = list(sigma = fit1$sigma$Pinv, theta = fit1$coefficients )) #con mvSIM simulamos nuestros datos ajustados a un modelo evolutivo determinado, en este caso BM no direccional
    x <- array(unlist(sim.temp.new), c(381,147,200), dimnames = list(dimnames(gpa.coords.temp)[[3]],as.character(c(1:147)),as.character(c(1:200)))) #con esto transformamos la lista de las 500 matrices de 381x150 que hemos simulado, en una sola matriz de 3 dimensiones (especies x coordenadas x simulaciones)

sim.disparities <- apply( x , 3 , disparity.phylo , phy = tree.temp )

### Generates sizes and symbols for deviations of empirical disparities from BM simulations

disparities.deviation <- log10( disparity.nodes.temp ) - log10( apply (sim.disparities, 1, median))
  max(disparities.deviation)
  mean(disparities.deviation)
  min(disparities.deviation)

edge.values <-  disparities.deviation - min( disparities.deviation )
  colours.rgb <- c("#fdf02d", "#d5cfce",  "#22c1c3" ) ##because we want people to look at outliers from BMw, grey is a good colour for unwanted attention
    colour.scale <- colorRamp( rev( colours.rgb ) , space = "rgb", bias = 0.45) #bias ensures grey corresponds always in good BM fitted values (adjust conveniently and check, more sophisticated in the final plots)
      BG.node.table <- colour.scale( edge.values / max( edge.values ) )
BG.disp <- c()
for( i in 1:nrow( BG.node.table ) ) {
  BG.disp[i] <- rgb( BG.node.table[i,1], BG.node.table[i,2], BG.node.table[i,3] , alpha = 255 , max = 255 )
}	

point.size.temp <- disparity.nodes.temp[-c(125)] * 2

#gets nodes only for the daughter-clade comparisons, i.e., only picks daughter clades with 2 or more descendent tips

tree.temp <- tree
  variable = edge.values
    tree <- tree.temp						

nodes.ages.temp <- tree.age(tree.temp) 
  node.ages.temp.mod <- nodes.ages.temp[,1]
    names(node.ages.temp.mod)<- as.character(nodes.ages.temp[,2])
      node.ages.temp.mod<-node.ages.temp.mod[ names(variable)]						

#plot to check nodes
plot(edge.values, node.ages.temp.mod, xlim = c(0, 2), ylim = c(325, 0) , type = "p", cex = point.size.temp, pch = 21, bg = BG.disp, frame.plot = T)
  text( y = node.ages.temp.mod , x =  edge.values, labels = names( edge.values ), pos = 4, cex =0.7, col = "black")	      

#tree with pls.angles
plot( tree.temp , type = "fan", show.tip.label = T, no.margin = T, y.lim =  c(-5, 72.907 ))
  nodelabels( node = as.numeric(names(edge.values)), frame = "n" , cex = point.size.temp , pch = 21 , bg = BG.disp )  


n.descendants.nodes <- c()
for (i in 1:length(descendants.temp)) {
  length.node.indiv <- length(descendants.temp[[i]])  
    n.descendants.nodes[i] <- length.node.indiv 
}
n.descendants.nodes <- as.numeric(n.descendants.nodes)
  names(n.descendants.nodes) <- names(descendants.temp)
    point.size.temp <- log10(n.descendants.nodes)


nodes.ages.temp <- tree.age(tree.temp) 
  node.ages.temp.mod <- nodes.ages.temp[,1]
    names(node.ages.temp.mod)<- as.character(nodes.ages.temp[,2])
      node.ages.temp.mod<-node.ages.temp.mod[ names(variable)]						

# generates colours for node ages 
edge.values <- node.ages.temp.mod
  colours.rgb <- c("#3664c7","#22b357","#FAFF7E") #azul-verde(life) ##because we want people to look at outliers from BMw, grey is a good colour for unwanted attention
  colour.scale <- colorRamp( rev( colours.rgb ) , space = "rgb", bias = 0.5) #bias ensures grey corresponds always in good BM fitted values (adjust conveniently and check, more sophisticated in the final plots)
    BG.node.table <- colour.scale( edge.values / max( edge.values ) )
BG.node.age<- c()
for( i in 1:nrow( BG.node.table ) ) {
  BG.node.age[i]<- rgb( BG.node.table[i,1], BG.node.table[i,2], BG.node.table[i,3] , alpha = 255 , max = 255 )
}	
disp.modul.nodes <- read.csv("DisparityVSModularity.perNode.csv",header=T,sep=";",row.names=1)
  plot(disparities.deviation, modularities.deviation, xlim = c(-0.75, 0.5), ylim = c(-3.5, 4) , type = "p", cex = point.size.temp, pch = 21, bg = BG.node.age, frame.plot = T)
    text( y = modularities.deviation , x =  disparities.deviation, labels = names( disparities.deviation ), pos = 4, cex =0.5, col = "black")	  

## REVISAR Y PROBAR
Plot1 <- ggplot(data, aes(log(Procrustes.variance), log(Convex.hulls), 
                          colour = Procrustes.variance,label = row.names(data))) + theme_bw() + geom_point(alpha = 1, size = 10)
Plot2 <-print(Plot1 + scale_colour_gradient(low = "gold2", high = "#063F10") 
              + geom_text(color = "black", size = 4, check_overlap = T))


# subset node.ages to data, do this before plotting scores 
tree.temp <- tree
  variable = disparity.by.node.conpmx
    tree <- tree.temp						

nodes.ages.temp <- tree.age(tree.temp) 
  node.ages.temp.mod <- nodes.ages.temp[,1]
    names(node.ages.temp.mod)<- as.character(nodes.ages.temp[,2])
      node.ages.temp.mod<-node.ages.temp.mod[ names(variable)]						

# generates colours for node ages 
edge.values <- node.ages.temp.mod
  colours.rgb <- c("#003949","#078148","#FFF835") #rojo-azul ##because we want people to look at outliers from BMw, grey is a good colour for unwanted attention
    colour.scale <- colorRamp( rev( colours.rgb ) , space = "rgb", bias = 0.5) #bias ensures grey corresponds always in good BM fitted values (adjust conveniently and check, more sophisticated in the final plots)
      BG.node.table <- colour.scale( edge.values / max( edge.values ) )
BG.node.age<- c()
for( i in 1:nrow( BG.node.table ) ) {
  BG.node.age[i]<- rgb( BG.node.table[i,1], BG.node.table[i,2], BG.node.table[i,3] , alpha = 255 , max = 255 )
}	

#plot to check nodes
plot(node.ages.temp.mod, modularity.per.node.temp.absval.sinpmx, xlim = c(max(node.ages.temp.mod), 0), ylim = c(min(modularity.per.node.temp.absval.sinpmx), max(modularity.per.node.temp.absval.sinpmx)) , type = "p", cex = 2, pch = 21, col = "transparent", bg = BG.node.age, frame.plot = F)
  text( y = modularity.per.node.temp.absval.conpmx , x =  node.ages.temp.mod, labels = names( node.ages.temp.mod ), pos = 4, cex =0.3, col = "black")	      




n.descendants.nodes <- list()

for (i in 1:length(descendants.daughters.temp)) {
  length.indiv <- length(descendants.daughters.temp[[i]])  
    n.descendants.nodes[i] <- length.indiv
}

n.descendants.nodes <- as.data.frame(unlist(n.descendants.nodes))
  write.csv(n.descendants.nodes$V1, 'temp.csv')

disp.modul.nodes.red <- read.csv("DisparityVSModularity.perNode.red.csv",header=T,sep=";",row.names=1)        
nodes.ages.temp <- tree.age(tree.temp) 
  node.ages.temp.mod <- nodes.ages.temp[,1]
    names(node.ages.temp.mod)<- as.character(nodes.ages.temp[,2])
      node.ages.temp.mod<-node.ages.temp.mod[ row.names(disp.modul.nodes.red)]
df.disp.modul.nodes.red <- data.frame(row.names = row.names(disp.modul.nodes.red),
                                      groups = disp.modul.nodes.red$node.clasif,
                                      disparity.conpmx = disp.modul.nodes.red$disparity.conpmx,
                                      disparity.sinpmx = disp.modul.nodes.red$disparity.sinpmx,
                                      modularity.conpmx = disp.modul.nodes.red$modularity.conpmx,
                                      modularity.sinpmx = disp.modul.nodes.red$modularity.sinpmx,
                                      deviation.modularity = disp.modul.nodes.red$modularity.deviation.sinpmx,
                                      deviation.disparity.conpmx = disp.modul.nodes.red$disparity.deviation.conpmx,
                                      deviation.disparity.sinpmx = disp.modul.nodes.red$disparity.deviation.sinpmx,
                                      node.ages = node.ages.temp.mod)

Plot1 <- ggplot(df.disp.modul.nodes.red, aes(deviation.disparity.conpmx, node.ages, label = rownames(df.disp.modul.nodes.red), colour = groups)) + 
  geom_point(aes(alpha = 1, size = disparity.conpmx)) 
  #geom_smooth(method = lm) 
  #geom_line(aes( x = reg.lines.x, y = reg.lines.y, color =  clades.comp), size = 2, alpha = 1)
    Plot2 <- print(Plot1 + geom_text(color = "black", size = 1, check_overlap = T))

# VIOLIN PLOTS // NO SE VE NADA CLARO

plot1 <- ggplot(df.disp.modul.nodes.red, aes(x = groups,y = modularity.sinpmx, color = groups))
  plot2 <- print(plot1 + geom_violin(trim = FALSE, adjust = .2) + geom_boxplot(width = .1, fill = "transparent", outlier.colour = "black") + stat_summary(fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) + theme_minimal() + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5))



## DISPARIDAD POR MÓDULOS
# PRECORDAL
Rawcoords_block1<-read.csv("Rawcoords78_todos381_block1.csv",header=T,sep=";",row.names=1)
  coords.block1<-arrayspecs(Rawcoords_block1[,1:78],26,3)
    gpa.block1<-gpagen(coords.block1)
      gdf.disp.block1 <- geomorph.data.frame(gpa.block1, clades.a = variables.amniotas$clades.a, clades.b = variables.amniotas$clades.b, 
                                       clades.c = variables.amniotas$clades.c, clades.comp = variables.amniotas$clades.comp)
disparity.block1 <- morphol.disparity(coords ~ 1, data = gdf.disp.block1, iter=999) # Morphological disparity for entire data set
  summary(disparity.block1)
disparity.block1 <- morphol.disparity(coords ~ Csize, data = gdf.disp.block1, iter=999) # Morphological disparity for entire data set, accounting for allometry
  summary(disparity.block1)
disparity.block1 <- morphol.disparity(coords ~ clades.comp, # Morphological disparity without covariates, using group means
                                      groups = ~ clades.comp, 
                                      data = gdf.disp.block1, iter=999)
  summary(disparity.block1)

matrix.coords.block1 <- two.d.array(gpa.block1$coords)           
  DTT.block1 <- dtt(phy = tree.amniota, data = matrix.coords.block1, index = c("avg.sq"), mdi.range = c(325,0), nsim = 1000, CI = 0.95, plot = TRUE)      

coords.block1.sinpmx<-arrayspecs(Rawcoords_block1[,4:78],25,3)
  gpa.block1.sinpmx<-gpagen(coords.block1.sinpmx)
    gdf.disp.block1.sinpmx <- geomorph.data.frame(gpa.block1.sinpmx, clades.a = variables.amniotas$clades.a, clades.b = variables.amniotas$clades.b, 
                                              clades.c = variables.amniotas$clades.c, clades.comp = variables.amniotas$clades.comp)
disparity.block1.sinpmx <- morphol.disparity(coords ~ 1, data = gdf.disp.block1.sinpmx, iter=999) # Morphological disparity for entire data set
  summary(disparity.block1.sinpmx)
disparity.block1.sinpmx <- morphol.disparity(coords ~ Csize, data = gdf.disp.block1.sinpmx, iter=999) # Morphological disparity for entire data set, accounting for allometry
  summary(disparity.block1.sinpmx)
disparity.block1.sinpmx <- morphol.disparity(coords ~ clades.comp, # Morphological disparity without covariates, using group means
                                             groups = ~ clades.comp, 
                                             data = gdf.disp.block1.sinpmx, iter=999)
  summary(disparity.block1.sinpmx)

matrix.coords.block1.sinpmx <- two.d.array(gpa.block1.sinpmx$coords)           
  DTT.block1.sinpmx <- dtt(phy = tree.amniota, data = matrix.coords.block1.sinpmx, index = c("avg.sq"), mdi.range = c(325,0), nsim = 1000, CI = 0.95, plot = TRUE)


# CORDAL
Rawcoords_block2<-read.csv("Rawcoords72_todos381_block2.csv",header=T,sep=";",row.names=1)
  coords.block2<-arrayspecs(Rawcoords_block2[,1:72],24,3)
    gpa.block2<-gpagen(coords.block2)
      gdf.disp.block2 <- geomorph.data.frame(gpa.block2, clades.a = variables.amniotas$clades.a, clades.b = variables.amniotas$clades.b, 
                                       clades.c = variables.amniotas$clades.c, clades.comp = variables.amniotas$clades.comp)
disparity.block2 <- morphol.disparity(coords ~ 1, data = gdf.disp.block2, iter=999) # Morphological disparity for entire data set
  summary(disparity.block2)
disparity.block2 <- morphol.disparity(coords ~ Csize, data = gdf.disp.block2, iter=999) # Morphological disparity for entire data set, accounting for allometry
  summary(disparity.block2)
disparity.block2 <- morphol.disparity(coords ~ clades.comp, # Morphological disparity without covariates, using group means
                                      groups = ~ clades.comp, 
                                      data = gdf.disp.block2, iter=999)
  summary(disparity.block2)

matrix.coords.block2 <- two.d.array(gpa.block2$coords)           
  DTT.block2 <- dtt(phy = tree.amniota, data = matrix.coords.block2, index = c("avg.sq"), mdi.range = c(325,0), nsim = 1000, CI = 0.95, plot = TRUE)



## DIFERENCIAS POR LANDMARK Y ENTRE ESPECIMENES MAX Y MIN

differences_from_mean <- coordinates.difference(coordinates = GPA$coords,reference = GPA$consensus, type = "spherical")  
procrustes_2d_array <- geomorph::two.d.array(GPA$coords)
ordination <- stats::prcomp(procrustes_2d_array)
pca<-gm.prcomp(GPA$coords)
PC1_variation <- variation.range(procrustes, axis = 1, ordination = ordination, type = "spherical")
hypothetical_1<-pca$shapes$shapes.comp1$min
hypothetical_2<-pca$shapes$shapes.comp1$max
mean<-newzombie.gpa$consensus
procrustes.var.plot(hypothetical_1, hypothetical_2, col = heat.colors,col.val = PC1_variation[, "radius"], labels = TRUE)
procrustes.var.plot(mean, hypothetical_2, col = heat.colors,col.val = PC1_variation[, "radius"], labels = TRUE)

  

#################################     
#### CONTRASTE DE DISTANCIAS ####
#################################     

##DATOS OBSERVADOS
#calcular la distancia filogenetica

phy.dist<-cophenetic.phylo(tree.amniota)

phy.dist.sinpmx.list <- t(combn(colnames(phy.dist), 2))
  phy.dist.sinpmx.list <- data.frame(phy.dist.sinpmx.list, dist=phy.dist[phy.dist.sinpmx.list]) #estas dos para transformar en lista uno de los triangulos de la pairwise matrix
    write.csv(phy.dist.sinpmx.list, 'phy.dist.sinpmx.list.csv')

#calcular la distancia fenotipica (distancia procrustes) // cambiar la base de datos fuente, la matriz de coordenadas con o sin premaxila en 2D

gpa.coords.temp <- gpa.brain$coords 
  gpa.coords.temp <- gpa.coords.temp[ , ,  unlist( tree.amniota$tip.label ) ] #reordenar la matriz para darle el orden en el que salen los tips de la filogenia
    proc.coords.2darray.ordenada <- as.matrix(two.d.array(gpa.coords.temp))
      write.csv(proc.coords.2darray, 'proc.coords.2darray.csv') #la exportamos para ponerla en orden alfabetico (mismo orden que la lista de distancias filogeneticas)
        proc.coords.2darray.ordenada <- read.csv('proc.coords.2darray.csv',header=T,sep=";",row.names=1)
          proc.coords.2darray.ordenada <- as.matrix(proc.coords.2darray.ordenada) #la volvemos a cargar y transformamos el dataframe en matriz
            proc.dist <- as.matrix(dist(proc.coords.2darray.ordenada)) #calculamos la distancia y la guardamos como matriz

              proc.dist.list <- t(combn(colnames(proc.dist), 2))
                proc.dist.list <- data.frame(proc.dist.list, dist=proc.dist[proc.dist.list])
                  write.csv(proc.dist.list, 'proc.dist.list.brain.csv')


##DATOS SIMULADOS (NO MULTIVARIATE) // NO SIRVE PARA DATOS DE ALTA DIMENSIONALIDAD, LINEAS CORRECTAS MAS ABAJO

vcv.procdist.sim <- ratematrix(tree.amniota, proc.coords.2darray.sinpmx.ordenada)
  procdist.sim <- sim.char(tree.amniota, vcv.procdist.sim, nsim = 1000, model = "BM")
    procdist.sim.mean <- apply(procdist.sim, c(1,2), mean)
      proc.dist.sinpmx.sim.mean <- as.matrix(dist(procdist.sim.mean))

        proc.dist.sinpmx.sim.mean.list <- t(combn(colnames(proc.dist.sinpmx.sim.mean), 2))
          proc.dist.sinpmx.sim.mean.list <- data.frame(proc.dist.sinpmx.sim.mean.list, dist=proc.dist.sinpmx.sim.mean[proc.dist.sinpmx.sim.mean.list])
            write.csv(proc.dist.sinpmx.sim.mean.list, 'proc.dist.sinpmx.sim.mean.list.csv')



## PHENOTYPIC DISTANCE VS PHYLOGENETIC DISTANCE (EMPIRICAL VS SIMULATED (mvMORPH)) // MORE ACCURATE FOR HIGH-DIMENSIONAL DATA

gpa.coords.temp <- gpa.todos$coords
gpa.coords.temp <- gpa.todos.sinpremax$coords 

block1.conpmx <- integration.conpmx$A1
block2.conpmx <- integration.conpmx$A2
block1.sinpmx <- integration.sinpmx$A1
block2.sinpmx <- integration.sinpmx$A2

gpa.coords.temp <- block2.sinpmx

gpa.coords.temp <- gpa.coords.temp[ , ,  unlist( tree.temp$tip.label ) ] #reordenar la matriz para darle el orden en el que salen los tips de la filogenia

proc.coords.2darray.ordenada <- as.matrix(two.d.array(gpa.coords.temp))
  proc.dist <- as.matrix(dist(proc.coords.2darray.ordenada)) #calculamos la distancia y la guardamos como matriz
    proc.dist.list <- t(combn(colnames(proc.dist), 2))
      proc.dist.list <- data.frame(proc.dist.list, dist=proc.dist[proc.dist.list])

#BM
library(mvMORPH)
fit1 <- mvgls( two.d.array(gpa.coords.temp) ~ 1 , tree = tree.temp , model = "BM" , method = "H&L" ) #primero ajustamos los datos al modelo evolutivo elegido, en este caso BM no direccional
  sim.temp.new <- mvSIM(tree.temp, nsim = 500, model = "BM1", param = list(sigma = fit1$sigma$Pinv, theta = fit1$coefficients )) #con mvSIM simulamos nuestros datos ajustados a un modelo evolutivo determinado, en este caso BM no direccional

#x <- array(unlist(sim.temp.new), c(381,147,500), dimnames = list(dimnames(gpa.coords.temp)[[3]],as.character(c(1:147)),as.character(c(1:500)))) #con esto transformamos la lista de las 500 matrices de 381x150 que hemos simulado, en una sola matriz de 3 dimensiones (especies x coordenadas x simulaciones)
#cambiar el valor de coordenadas de 150 a 147 

# vamos a probarlo al reves, haciendo distancias procrustes de todas las simulaciones
proc.dist.sim.list <- list() #hacemos una lista donde se irán almacenando los resultados hasta llegar a 500 elementos de 72390 comparaciones cada uno

for (i in 1: length(sim.temp.new)) { 
  
  simulated.data.individual <- sim.temp.new[[i]] 
  
  proc.dist.individual <- as.matrix(dist(simulated.data.individual)) #calculamos la distancia para la simulacion individual y transformamos el resultado en matriz
    proc.dist.individual.ordered <- t(combn(colnames(proc.dist.individual), 2)) #ordenamos los resultados en una sola columna
      proc.dist.individual.ordered <- data.frame(proc.dist.individual.ordered, dist=proc.dist.individual[proc.dist.individual.ordered]) #los transformamos en dataframe
        proc.dist.individual.num <- proc.dist.individual.ordered$dist #seleccionamos la columna que nos interesa
          names(proc.dist.individual.num) <- row.names(phydist.procdist.df) #le damos el nombre de las comparaciones a pares (taxon1:taxon2)
  
  proc.dist.sim.list[[i]]<-proc.dist.individual.num #metes el vector numerico con sus nombres en los 500 elementos de la lista
}

#x <- array(unlist(proc.dist.sim.list), c(72390,1,500), dimnames = list(row.names(phydist.procdist.df),as.character(c(1)),as.character(c(1:500)))) #con esto transformamos la lista de las 500 matrices de 381x150 que hemos simulado, en una sola matriz de 3 dimensiones (especies x coordenadas x simulaciones)

x <- array(unlist(proc.dist.sim.list), c(72390,500), dimnames = list(row.names(phydist.procdist.df),as.character(c(1:500)))) #hacemos una matriz de nºcomparacions x nºsimulaciones

#con las siguientes lineas hacemos estadistica de la matrix de 72390 x 500,
#pero aun no sale bien porque combina valores de una y otra simulacion, hay que seguir bicheando  

proc.dist.mean <- x[,134]                  
  write.csv(proc.dist.mean, 'proc.dist.sim.mean.brain.csv')
proc.dist.quan95 <- apply(x, 1, quantile, probs = 0.95)                  
  write.csv(proc.dist.quan95, 'proc.dist.sim.quan95.brain.csv')
proc.dist.quan05 <- apply(x, 1, quantile, probs = 0.05)                  
  write.csv(proc.dist.quan05, 'proc.dist.sim.quan05.brain.csv') 
proc.dist.quan <- apply(x, 1, quantile, probs = c(0.05, 0.95))                  
  write.csv(proc.dist.quan95, 'proc.dist.sim.quan05.brain.csv')

proc.dist.max <- apply(x , 1 , max)                  
proc.dist.min <- apply(x , 1 , min) 


#con la siguiente linea sacamos un vector con la media de cada simulacion, 
#estoy podria valer para seleccionar la simulacion que tenga los valores 
#mas cercanos a la media de todo, y los dos extremos max y min    
z <- apply(x , 2 , mean) 
mean(z)
z
#simulacion mas cercana a la media en los datos sin premaxila -> #57  
#simulacion mas cercana a la media en los datos con premaxila -> #471 
#simulacion mas cercana a la media en los datos del bloque 1 con premaxila -> #123
#simulacion mas cercana a la media en los datos del bloque 2 con premaxila -> #365
#simulacion mas cercana a la media en los datos del bloque 1 sin premaxila -> #62
#simulacion mas cercana a la media en los datos del bloque 2 sin premaxila -> #470
#simulacion mas cercana a la media en los datos endocraneales -> #134


{zplot(phydist.procdist.df$phy.dist, x[,57], type = "p", pch = 19, cex = 1, col = alpha("grey", 0.3))
  points(phydist.procdist.df$phy.dist, phydist.procdist.df$proc.dist.obs.rel, type = "p", pch = 19, cex = 2, col = alpha("red", 0.1)) #ploteamos encima los empíricos
    points(phydist.procdist.df$phy.dist, proc.dist.min, type = "p", pch = 19, cex = 0.5, col = alpha("blue", 0.3))
      lines(phydist.procdist.df$phy.dist, proc.dist.min, type = "l", col = "blue")
  
  
  phydist.procdist.df <- read.csv('PHYDIST.PROCDIST.csv',header=T,sep=";",row.names=1) #cambiar datos en excel fuente
    plot(phydist.procdist.df$phy.dist, phydist.procdist.df$proc.dist.sim.mean..rel, type = "p", pch = 19, cex = 2, col = alpha("grey", 0.3))
      points(phydist.procdist.df$phy.dist, phydist.procdist.df$proc.dist.sim.quan95..rel, type = "p", pch = 19, col = alpha("grey", 0.3))   
        points(phydist.procdist.df$phy.dist, phydist.procdist.df$proc.dist.sim.quan05..rel, type = "p", pch = 19, col = alpha("grey", 0.3))
          points(phydist.procdist.df$phy.dist, phydist.procdist.df$proc.dist.obs.rel, type = "p", pch = 19, cex = 0.5, col = alpha("red", 0.2))} #INTENTOS FALLIDOS DE PLOTEAR, PERO EL CÓDIGO ESTÁ BIEN


## PLOT EN BINS DE EMPIRICOS VS SIMULADOS DE PHE VS PHY  [FUNCIONA]
df1 <- data.frame(x = phydist.procdist.df$phy.dist, y = proc.dist.list$dist) #usar siempre valores absolutos de distancia procrustes
df2 <- data.frame(x = phydist.procdist.df$phy.dist, y = proc.dist.mean)

library(dplyr)
bothDF <- dplyr::bind_rows(A = df1, B = df2, .id = "df")
bothHex <- hexbin::hexbin(x = bothDF$x, y = bothDF$y, IDs = TRUE, xbins = 50)
counts <- hexbin::hexTapply(bothHex, factor(bothDF$df), table) %>% 
  simplify2array %>% 
  t %>% 
  data.frame() %>% 
  mutate(id = as.numeric(row.names(.)), diff = log(A+1) - log(B+1)) %>%   #logaritmizamos cada elemento de la diferencia entre empiricos y simulados, para suavizar el big data y poder mantener los valores negativos
  dplyr::left_join(data.frame(id = bothHex@cell, hexbin::hcell2xy(bothHex)))
head(counts)


counts %>%
  ggplot(aes(x = x, y = y, fill = diff)) + 
  geom_hex(stat = "identity") +
  scale_fill_gradientn(values = c(1, 0.91, 0.64, 0.37, 0), 
                       colours = c("#fdbb2d", "#b21f1f", "#1a2a6c", "#f5f5f5", "#f5f5f5")) +  #utilizamos rangos de entre 0 y 1 para que los cambios de color sean a partir de ellos y las tonalidades contribuyan a observar las diferencias mas reseñables
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

##################################################
####  COMMON AXES OF MORPHOLOGICAL VARIATION  ####
##################################################

### NEW FUNCTION THAT COMPARES ALL NODES between them ###

# Funcion creada desde cero para sacar los valores de angulos entre los vectores de variacion de todos los clados (13695 comparaciones a pares)
# Tambien saca el nombre de cada comparacion (nodo1 x nodo2), la edad de su MRCA, y la distancia filogenética (patristica) que los separa

phy$node.label <- nodes.temp #le damos nombre a los nodos para que sepa reconocerlos al hacer los analisis

eigendirections.transverse <- function( data , phy ) {
  
  coords.matrix<-data
  nodes.temp <- unique( phy$edge[,1] )
  descendants.temp <- lapply( nodes.temp , FUN = get.descendants , phy = phy )
  names( descendants.temp ) <- nodes.temp
  discarded.descendants <- as.vector(names(which(lapply(descendants.temp, length)<5)))
  descendants.temp.mod <- descendants.temp[names(descendants.temp) %in% discarded.descendants == FALSE]
  node.data.list <- lapply( descendants.temp.mod , FUN = subset.rows , data = coords.matrix )
  
  
  vectors.matrix1 <- matrix( ncol = length(node.data.list) , nrow = 147)
  
  for (i in 1:length(node.data.list)) {
    node1.data <- node.data.list[[i]]
    GM.node1<-arrayspecs(as.matrix(node1.data), dim(as.matrix(node1.data))[2]/3, 3)
    pc.node1.vector<-gm.prcomp(GM.node1)$rotation[,3]
    name.node1 <- names(node.data.list)[i]
    
    vectors.matrix1[,i] <- pc.node1.vector
  }
  
  colnames(vectors.matrix1) <- names(node.data.list)
  
  
  vectors.matrix2 <- matrix( ncol = length(node.data.list) , nrow = 147)
  
  for (i in 1:length(node.data.list)) {
    node2.data <- node.data.list[[i]]
    GM.node2<-arrayspecs(as.matrix(node2.data), dim(as.matrix(node2.data))[2]/3, 3)
    pc.node2.vector<-gm.prcomp(GM.node2)$rotation[,3]
    name.node2 <- names(node.data.list)[i]
    
    vectors.matrix2[,i] <- pc.node2.vector      
  }
  
  colnames(vectors.matrix2) <- names(node.data.list)
  
  angles.vectors <- matrix(ncol = length(node.data.list) , nrow = length(node.data.list))  
  
  for (x in 1:length(node.data.list)) {
    vector1 <- vectors.matrix1[,x]
    angles.vectors[x,] <- angle( vector1, vectors.matrix2[,1:length(node.data.list)], degree = T)
  }  
  
  colnames(angles.vectors) <- names(node.data.list)
  rownames(angles.vectors) <- names(node.data.list)
  
  angles.vectors.list <- t(combn(colnames(angles.vectors), 2))
  angles.vectors.list <- data.frame(angles.vectors.list, angle=angles.vectors[angles.vectors.list]) #estas dos para transformar en lista uno de los triangulos de la pairwise matrix
  
  
  pc.angles.temp <- as.matrix(angles.vectors.list)
  
  library(castor)
  age.mrca.comparisons.list <- matrix(ncol = 1 , nrow = length(pc.angles.temp[,1]))
  for (i in 1:length(pc.angles.temp[,1])) {
    name.node1 <- as.numeric(pc.angles.temp[i,1])
    name.node2 <- as.numeric(pc.angles.temp[i,2])
    node.mrca <- get_pairwise_mrcas( phy , name.node1, name.node2)
    age.mrca.comparisons.list[i,] <- 325 - (nodeheight(phy, node.mrca))
  }
  
  phy.dist.comparisons.list <- matrix(ncol = 1 , nrow = length(pc.angles.temp[,1]))
  for (i in 1:length(pc.angles.temp[,1])) {
    name.node1 <- as.numeric(pc.angles.temp[i,1])
    name.node2 <- as.numeric(pc.angles.temp[i,2])
    node.mrca <- get_pairwise_distances( phy , name.node1, name.node2)
    phy.dist.comparisons.list[i,] <- node.mrca
  }
  
  names.compared.nodes <- matrix(ncol = 1 , nrow = length(pc.angles.temp[,1]))
  for (i in 1:length(pc.angles.temp[,1])) {
    name.node1 <- as.numeric(pc.angles.temp[i,1])
    name.node2 <- as.numeric(pc.angles.temp[i,2])
    names.compared.nodes[i] <- paste( name.node1, "x", name.node2)
  }
  
  
  pc.angles <- data.frame(node1 = pc.angles.temp[,1], 
                          node2 = pc.angles.temp[,2], 
                          angle = pc.angles.temp[,3], 
                          age.mrca.comparison = age.mrca.comparisons.list,
                          patris.dist.comparison = phy.dist.comparisons.list,
                          row.names = names.compared.nodes)
  
  #this block makes angles all angles below 90 degrees: because PLS vectors are arbitrary in direction the maximum angle is orthogonal	
  pc.angles.mod <- as.numeric(pc.angles$angle)
  names(pc.angles.mod) <- row.names(pc.angles)
  pc.angles.mod.wrong <-  pc.angles.mod [which(  pc.angles.mod > 90 )]
  pc.angles.mod.corrected <- (pc.angles.mod.wrong -180) * -1
  pc.angles.mod[names(pc.angles.mod.corrected)] <- pc.angles.mod.corrected
  pc.angles$angle <- pc.angles.mod
  
  pc.angles 
  
}  #### FUNCIONAAAAAAAAA VAAAAMOOOOOS

PC1.angles.df <- eigendirections.transverse( data = gpa.todos.sinpremax$coords , phy = phy)
write.csv(PC1.angles.df, 'pc1.angles.sinpmx.allnodescombinations.csv')

PC2.angles.df <- eigendirections.transverse( data = gpa.todos.sinpremax$coords , phy = phy)
write.csv(PC2.angles.df, 'pc2.angles.sinpmx.allnodescombinations.csv')

PC3.angles.df <- eigendirections.transverse( data = gpa.todos.sinpremax$coords , phy = phy)
write.csv(PC3.angles.df, 'pc3.angles.sinpmx.allnodescombinations.csv')


### SIMULADOS ###

library(mvMORPH)

gpa.coords.temp <- data
fit1 <- mvgls( two.d.array(gpa.coords.temp) ~ 1 , tree = phy , model = "BM" , method = "H&L" ) #primero ajustamos los datos al modelo evolutivo elegido, en este caso BM no direccional
sim.temp.new <- mvSIM(phy, nsim = 50, model = "BM1", param = list(sigma = fit1$sigma$Pinv, theta = fit1$coefficients )) #con mvSIM simulamos nuestros datos ajustados a un modelo evolutivo determinado, en este caso BM no direccional
x <- array(unlist(sim.temp.new), c(381,147,50), dimnames = list(dimnames(gpa.coords.temp)[[3]],as.character(c(1:147)),as.character(c(1:50)))) #con esto transformamos la lista de las 500 matrices de 381x150 que hemos simulado, en una sola matriz de 3 dimensiones (especies x coordenadas x simulaciones)

eigendirections.transverse.4simul <- function( data , phy ) {
  
  coords.matrix<-data
  nodes.temp <- unique( phy$edge[,1] )
  descendants.temp <- lapply( nodes.temp , FUN = get.descendants , phy = phy )
  names( descendants.temp ) <- nodes.temp
  discarded.descendants <- as.vector(names(which(lapply(descendants.temp, length)<5)))
  descendants.temp.mod <- descendants.temp[names(descendants.temp) %in% discarded.descendants == FALSE]
  node.data.list <- lapply( descendants.temp.mod , FUN = subset.rows , data = coords.matrix )
  
  
  vectors.matrix1 <- matrix( ncol = length(node.data.list) , nrow = 147)
  
  for (i in 1:length(node.data.list)) {
    node1.data <- node.data.list[[i]]
    GM.node1<-arrayspecs(as.matrix(node1.data), dim(as.matrix(node1.data))[2]/3, 3)
    pc.node1.vector<-gm.prcomp(GM.node1)$rotation[,2]
    name.node1 <- names(node.data.list)[i]
    
    vectors.matrix1[,i] <- pc.node1.vector
  }
  colnames(vectors.matrix1) <- names(node.data.list)
  
  
  vectors.matrix2 <- matrix( ncol = length(node.data.list) , nrow = 147)
  
  for (i in 1:length(node.data.list)) {
    node2.data <- node.data.list[[i]]
    GM.node2<-arrayspecs(as.matrix(node2.data), dim(as.matrix(node2.data))[2]/3, 3)
    pc.node2.vector<-gm.prcomp(GM.node2)$rotation[,2]
    name.node2 <- names(node.data.list)[i]
    
    vectors.matrix2[,i] <- pc.node2.vector      
  }
  colnames(vectors.matrix2) <- names(node.data.list)
  
  
  angles.vectors <- matrix(ncol = length(node.data.list) , nrow = length(node.data.list))  
  
  for (x in 1:length(node.data.list)) {
    vector1 <- vectors.matrix1[,x]
    angles.vectors[x,] <- angle( vector1, vectors.matrix2[,1:length(node.data.list)], degree = T)
  } 
  colnames(angles.vectors) <- names(node.data.list)
  rownames(angles.vectors) <- names(node.data.list)
  
  angles.vectors.list <- t(combn(colnames(angles.vectors), 2))
  angles.vectors.list <- data.frame(angles.vectors.list, angle=angles.vectors[angles.vectors.list]) #estas dos para transformar en lista uno de los triangulos de la pairwise matrix
  
  
  pc.angles.temp <- as.matrix(angles.vectors.list)
  
  names.compared.nodes <- matrix(ncol = 1 , nrow = length(pc.angles.temp[,1]))
  for (i in 1:length(pc.angles.temp[,1])) {
    name.node1 <- as.numeric(pc.angles.temp[i,1])
    name.node2 <- as.numeric(pc.angles.temp[i,2])
    names.compared.nodes[i] <- paste( name.node1, "x", name.node2)
  }
  
  pc.angles <- data.frame(node1 = pc.angles.temp[,1], 
                          node2 = pc.angles.temp[,2], 
                          angle = pc.angles.temp[,3], 
                          row.names = names.compared.nodes)
  
  #this block makes angles all angles below 90 degrees: because PLS vectors are arbitrary in direction the maximum angle is orthogonal	
  pc.angles.mod <- as.numeric(pc.angles$angle)
  names(pc.angles.mod) <- row.names(pc.angles)
  pc.angles.mod.wrong <-  pc.angles.mod [which(  pc.angles.mod > 90 )]
  pc.angles.mod.corrected <- (pc.angles.mod.wrong -180) * -1
  pc.angles.mod[names(pc.angles.mod.corrected)] <- pc.angles.mod.corrected
  pc.angles$angle <- pc.angles.mod
  pc.angles.vec <- as.numeric(pc.angles$angle)
  pc.angles.vec
}

sim.PC1.angles.matrix <- apply( x , 3 , eigendirections.transverse.4simul , phy = phy ) #cambiar a $rotation[,1]
sim.PC1.angles.mean <- apply( sim.PC1.angles.matrix , 1 , mean )

sim.PC2.angles.matrix <- apply( x , 3 , eigendirections.transverse.4simul , phy = phy ) #cambiar a $rotation[,2]
sim.PC2.angles.mean <- apply( sim.PC2.angles.matrix , 1 , mean )
sim.PC2.angles.indiv <- sim.PC2.angles.matrix[,49]

sim.PC3.angles.matrix <- apply( x , 3 , eigendirections.transverse.4simul , phy = phy ) #cambiar a $rotation[,3]
sim.PC3.angles.mean <- apply( sim.PC3.angles.matrix , 1 , mean )



sim.PC1.angles.df <- eigendirections.transverse.BMsim( data = gpa.todos.sinpremax$coords , phy = phy)
write.csv(sim.PC1.angles.df, 'pc1.angles.sinpmx.allnodescombinations.csv')

sim.PC2.angles.df <- eigendirections.transverse.BMsim( data = gpa.todos.sinpremax$coords , phy = phy)
write.csv(sim.PC2.angles.df, 'pc2.angles.sinpmx.allnodescombinations.csv')

sim.PC3.angles.df <- eigendirections.transverse.BMsim( data = gpa.todos.sinpremax$coords , phy = phy)
write.csv(sim.PC3.angles.df, 'pc3.angles.sinpmx.allnodescombinations.csv')


### Generates sizes and symbols for deviations of empirical disparities from BM simulations

angles.deviation <- log10( PC2.angles.df$angle ) - log10( sim.PC2.angles.mean )
max(angles.deviation)
mean(angles.deviation)
min(angles.deviation)

edge.values <-  angles.deviation - min( angles.deviation )
names(edge.values) <- rownames( PC1.angles.df )
colours.rgb <- c("#c765c5", "#d5cfce",  "#6cac5c" ) ##because we want people to look at outliers from BMw, grey is a good colour for unwanted attention
colour.scale <- colorRamp( rev( colours.rgb ) , space = "rgb", bias = 0.60) #bias ensures grey corresponds always in good BM fitted values (adjust conveniently and check, more sophisticated in the final plots)
BG.node.table <- colour.scale( edge.values / max( edge.values ) )
BG.ang.dev <- c()
for( i in 1:nrow( BG.node.table ) ) {
  BG.ang.dev[i] <- rgb( BG.node.table[i,1], BG.node.table[i,2], BG.node.table[i,3] , alpha = 255 , max = 255 )
}	

#plot to check nodes

plot(edge.values, PC1.angles.df$age.mrca.comparison , xlim = c(0, 3.3), ylim = c(330, 0) , type = "p", cex = 2, pch = 21, bg = BG.ang.dev, frame.plot = T)
text( y = PC2.angles.df$age.mrca.comparison , x =  edge.values, labels = names( edge.values ), pos = 4, cex = 0.3, col = "black")	      

# generates colours for angles 

x <- sim.PC1.angles.mean
y1 <- PC3.angles.df$age.mrca.comparison
y2 <- PC3.angles.df$patris.dist.comparison
edge.values <- x
colours.rgb <- c("#96271C","#FBF399","#2f74a2") #naranja rojizo-amarillo-gris ##because we want people to look at outliers from BMw, grey is a good colour for unwanted attention
colour.scale <- colorRamp( rev( colours.rgb ) , space = "rgb") #bias ensures grey corresponds always in good BM fitted values (adjust conveniently and check, more sophisticated in the final plots)
BG.node.table <- colour.scale( edge.values / max( edge.values ) )
BG.angles<- c()
for( i in 1:nrow( BG.node.table ) ) {
  BG.angles[i]<- rgb( BG.node.table[i,1], BG.node.table[i,2], BG.node.table[i,3] , alpha = 255 , max = 255 )
}	

phydist.point.size <- y2  ### AUN NO REFLEJA CORRECTAMENTE LA IDEA DE LOS TAMAÑOS EN FUNCION DE LA DIST PHY
log.phydist.point.size <- log(phydist.point.size)
log.phydist.point.size.mod <- (log.phydist.point.size + (min(log.phydist.point.size)*-1))/3
log.phydist.point.size.mod

# plot to check nodes

plot(x, y1, xlim = c(-10, 100), ylim = c(330, 0) , type = "p", cex = 1, pch = 21, bg = BG.angles, frame.plot = T)
text( y = y1 , x =  x, labels = rownames( PC1.angles.df ), pos = 4, cex =0.2, col = "black")	      



## NEW FUNCTION THAT COMPARES DAUGHTER NODES between them ##

eigendirections.horizontal.phylo <- function( data , phy ) {
  
  coords.matrix<-two.d.array(data)
  internal.node.pairs <- phy$edge[!phy$edge[, 2] <= length(phy$tip.label),] #subsets a matrix with all the node comparisons we want: first column inmmediately parent nodes, second column immediately daughter nodes
  
  descendants.daughters.temp <- lapply( internal.node.pairs[,2] , FUN = get.descendants , phy = phy ) #this block makes the ata list for daugther clades
  names( descendants.daughters.temp ) <- internal.node.pairs[, 2]
  node.data.daughters.list <- lapply( descendants.daughters.temp , FUN = subset.rows , data = coords.matrix )
  
  pc.angles<-list()
  
  for (i in 1:length(internal.node.pairs[,1])) { # loop reconverts data to 3D array so morphol.disparity can work
    
    node <- as.character(internal.node.pairs [i, 1])
    daughters <- internal.node.pairs[internal.node.pairs[,1] == node,]
    
    if (length(daughters) == 4 ) { 
      
      
      daughter.1.node <- as.character(daughters[1, 2])
      daughter.2.node <- as.character(daughters[2, 2])
      
      GM.node.data.daughters.1.list<-arrayspecs(as.matrix(node.data.daughters.list[[daughter.1.node]]), dim(as.matrix(node.data.daughters.list[[daughter.1.node]]))[2]/3, 3)
      pc.daughter.1.vector<-gm.prcomp(GM.node.data.daughters.1.list)$rotation[,1]
      
      GM.node.data.daughters.2.list<-arrayspecs(as.matrix(node.data.daughters.list[[daughter.2.node]]), dim(as.matrix(node.data.daughters.list[[daughter.2.node]]))[2]/3, 3)
      pc.daughter.2.vector<-gm.prcomp(GM.node.data.daughters.2.list)$rotation[,1]
      
      pc.angles[i]<-angle( pc.daughter.1.vector, pc.daughter.2.vector, degree = T)
      
    }
    
    else {
      daughter.1.node <- as.character(daughters[2])
      daughter.2.node <- as.character(daughters[2])
      
      GM.node.data.daughters.1.list<-arrayspecs(as.matrix(node.data.daughters.list[[daughter.1.node]]), dim(as.matrix(node.data.daughters.list[[daughter.1.node]]))[2]/3, 3)
      pc.daughter.1.vector<-gm.prcomp(GM.node.data.daughters.1.list)$rotation[,1]
      
      GM.node.data.daughters.2.list<-arrayspecs(as.matrix(node.data.daughters.list[[daughter.2.node]]), dim(as.matrix(node.data.daughters.list[[daughter.2.node]]))[2]/3, 3)
      pc.daughter.2.vector<-gm.prcomp(GM.node.data.daughters.2.list)$rotation[,1]
      
      pc.angles[i]<-angle( pc.daughter.1.vector, pc.daughter.2.vector, degree = T)
      
    }
    
  }
  
  
  pc.angles<-unlist(pc.angles)
  names(pc.angles)<- ###as.character(internal.node.pairs[,1])
    pc.angles
  
  
}


#ANGLES BETWEEN PC1 VECTORS OF DESCENDANT CLADES
PC1.angles.daugthers.per.node.conpmx <- eigendirections.horizontal.phylo( data = gpa.todos$coords , phy = phy )
PC1.angles.daugthers.per.node.sinnpmx <- eigendirections.horizontal.phylo( data = gpa.todos.sinpremax$coords , phy = phy )

#ANGLES BETWEEN PC2 VECTORS OF DESCENDANT CLADES (modify $rotation[,2])
PC2.angles.daugthers.per.node.conpmx <- eigendirections.horizontal.phylo( data = gpa.todos$coords , phy = phy )
PC2.angles.daugthers.per.node.sinnpmx <- eigendirections.horizontal.phylo( data = gpa.todos.sinpremax$coords , phy = phy )

#ANGLES BETWEEN PC3 VECTORS OF DESCENDANT CLADES (modify $rotation[,3])
PC3.angles.daugthers.per.node.conpmx <- eigendirections.horizontal.phylo( data = gpa.todos$coords , phy = phy )
PC3.angles.daugthers.per.node.sinnpmx <- eigendirections.horizontal.phylo( data = gpa.todos.sinpremax$coords , phy = phy )


## PLS vector angles  between clades 
# PLS angles between daughters clades, plots value in parent clade

tree.temp <- phy
pc.angles <- PC1.angles.daugthers.per.node.conpmx #cambiar al elemento que corresponda

# prepares angles to be plotted
pc.angles.mod <- pc.angles[!is.na(pc.angles)]# removes NaNs product of clade-tip daughter comparisons
pc.angles.mod <- pc.angles.mod[!duplicated(pc.angles.mod)] # two daughters per clade are obtained but we want only to keep one to plot in the immediate parent node

#this block makes all angles below 90 degrees: because PC vectors are arbitrary in direction the maximum angle is orthogonal	
pc.angles.mod
pc.angles.mod.wrong <-  pc.angles.mod [which(  pc.angles.mod > 90 )]
pc.angles.mod.corrected <-  (pc.angles.mod.wrong -180) * -1
pc.angles.mod[names(pc.angles.mod.corrected)] <- pc.angles.mod.corrected
pc.angles.mod<-pc.angles.mod[- which(pc.angles.mod >= 0 & pc.angles.mod <= 2)]

#prunes modularity values to match only the represented nodes
#variance.percentage.point.size <- modularity.per.node.temp.absval.sinpmx [names (pc.angles.mod)]
#variance.percentage.point.size.mod <- variance.percentage.point.size/10
#variance.percentage.point.size.mod

# generates colours for angles 

edge.values <- pc.angles.mod
colours.rgb <- c("#ff6200","#fff367","#e8e8e8") #naranja rojizo-amarillo-gris ##because we want people to look at outliers from BMw, grey is a good colour for unwanted attention
colour.scale <- colorRamp( rev( colours.rgb ) , space = "rgb", bias = 0.6) #bias ensures grey corresponds always in good BM fitted values (adjust conveniently and check, more sophisticated in the final plots)
BG.node.table <- colour.scale( edge.values / max( edge.values ) )
BG.angles<- c()
for( i in 1:nrow( BG.node.table ) ) {
  BG.angles[i]<- rgb( BG.node.table[i,1], BG.node.table[i,2], BG.node.table[i,3] , alpha = 255 , max = 255 )
}	

#gets nodes only for the daughter-clade comparisons, i.e., only picks daughter clades with 2 or more descendent tips
tree.temp <- tree
variable = pc.angles.mod
tree <- tree.temp						

nodes.ages.temp <- tree.age(tree.temp) 
node.ages.temp.mod <- nodes.ages.temp[,1]
names(node.ages.temp.mod)<- as.character(nodes.ages.temp[,2])
node.ages.temp.mod<-node.ages.temp.mod[ names(variable)]						

#plot to check nodes
plot(pc.angles.mod, node.ages.temp.mod, xlim = c(0, 100), ylim = c(120, 0) , type = "p", cex = variance.percentage.point.size.mod, pch = 21, bg = BG.angles, frame.plot = T)
text( y = node.ages.temp.mod , x =  pc.angles.mod, labels = names( pc.angles.mod ), pos = 4, cex =0.7, col = "black")	      

#tree with pls.angles
plot( phy , type = "fan", show.tip.label = T, no.margin = T, y.lim =  c(-5, 72.907 ))
nodelabels( node = as.numeric(names(pc.angles.mod)), frame = "n" , cex = variance.percentage.point.size.mod , pch = 21 , bg = BG.angles )





###########################
#### INTEGRATION TESTS ####
###########################

#hipótesis a partir del origen de desarrollo (precordal (cresta neural) vs. cordal (mesodérmico))
modules.LM.list.conpmx <- c("A","A","A","A","A","A","A","A","A","A","A","A","B","B","B","B","B","B","B","B","B","B","B","B","B","B","A","B","A","A","B","A","A","A","A","A","A","A","A","A","B","B","B","B","B","B","B","A","B","A") 
integration.conpmx <- phylo.integration(gpa.todos$coords, partition.gp = modules.LM.list.conpmx, phy = tree.amniota, iter = 999)
  summary(integration.conpmx)
    max.shapes.conpmx <- picknplot.shape(plot(integration.conpmx))

preds.min.pls.conpmx <- read.csv("preds.min.pls.conpmx.csv",header=T,sep=";",row.names=1)
  preds.min.pls.conpmx <- as.matrix(preds.min.pls.conpmx)
preds.max.pls.conpmx <- read.csv("preds.max.pls.conpmx.csv",header=T,sep=";",row.names=1)
  preds.max.pls.conpmx <- as.matrix(preds.max.pls.conpmx)

library(tibble)
per_lm_distance <- function (shape.data.min , shape.data.max) #para modificar los colores de los LMs, meter los codigos que gusten en cols1
{
  df <- as.data.frame(cbind(shape.data.min, shape.data.max))
  distances <- apply(df, 1, function(x) dist(matrix(x, nrow = 2, byrow = TRUE)))
  cols1 <- colorRampPalette(c("#1a2a6c", "#7B90E3", "#fdbb2d", "#b21f1f"))
  cols <- cols1(100)
  x = (log10(distances))
  xlims <- NULL
  tol <- 1e-06
  xlims <- range(x) + c(-tol, tol)
  nbin = 100
  breaks <- 0:nbin/nbin * (xlims[2] - xlims[1]) + xlims[1]
  whichColor <- function(p, cols, breaks) {
    i <- 1
    while (p >= breaks[i] && p > breaks[i + 1]) i <- i +
        1
    cols[i]
  }
  distancecolors <- sapply(x, whichColor, cols = cols, breaks = breaks)
  distance_table <- tibble(Per_Lm_Distance = distances, Log_Distance = x,Distance_Colors = distancecolors)
  return(distance_table)
}

my.distances.preds.pls.conpmx <- per_lm_distance(shape.data.min = preds.min.pls.conpmx , shape.data.max = preds.max.pls.conpmx)
coords.original <- coords.todos[,,"Ptilonorhynchus_violaceus"] #Ptilonorhynchus_violaceus , Mephitis_mephitis
open3d()
par3d(windowRect = c(0,0,500,500)) #para dimensiones del cuadro RGL
Sys.sleep(1)
mfrow3d(nr = 1, nc = 2, byrow = TRUE, sharedMouse = TRUE)
mesh.min <- tps3d(mesh.original, refmat = coords.original, tarmat = preds.min.pls.conpmx) #maps landmarks or a triangular mesh via thin plate spline based on a reference and a target configuration in 2D and 3D
shade3d(mesh.min, col= 8, alpha = 0.3) #dibuja la mesh con la forma predicha en la línea anterior #alpha cambia la transparencia de la surface
spheres3d(preds.min.pls.conpmx, col = my.distances.preds.pls.conpmx$Distance_Colors, radius = 0.01) #radius cambia el tamaño del LM
next3d()
mesh.max <- tps3d(mesh.original, refmat = coords.original, tarmat = preds.max.pls.conpmx, lambda = 0.5)
shade3d(mesh.max, col= 8, alpha = 0.3) 
spheres3d(preds.max.pls.conpmx, col = my.distances.preds.pls.conpmx$Distance_Colors, radius = 0.01)


df.integration.conpmx <- data.frame(XScores.conpmx = integration.conpmx$XScores[,1], 
                                    YScores.conpmx = integration.conpmx$YScores[,1],
                                    log.skullCS = log(variables.amniotas$skull.CS), 
                                    log.cranialbase.CS = log(variables.amniotas$cranialbase.CS),
                                    log.brainconvhull = log(variables.amniotas$brain.convhull),
                                    headbrain.res = variables.amniotas$headbrain.res,
                                    skullbase.res = variables.amniotas$skullbase.res,
                                    foramen = variables.amniotas$foramen.angule,
                                    base.ang = variables.amniotas$base.angule,
                                    flexion = variables.amniotas$flexion.angule,
                                    sph.res = variables.amniotas$SPH1.res,
                                    clivus.foramen.ang = variables.amniotas$clivus.foramen.angle,
                                    facial.ang = variables.amniotas$facial.angle,
                                    clades.a = variables.amniotas$clades.a,
                                    clades.b = variables.amniotas$clades.b,
                                    clades.c = variables.amniotas$clades.c,
                                    clades.comp = variables.amniotas$clades.comp,
                                    row.names = row.names(variables.amniotas))

Plot1 <- ggplot(df.integration.conpmx, aes(XScores.conpmx, YScores.conpmx, label = rownames(df.integration.conpmx), colour = clades.a))+ geom_point(alpha = 1, size = 5) + geom_smooth(method = glm) #+ geom_line(aes( x = reg.lines.x, y = reg.lines.y, color =  clades.b), size = 2, alpha = 1)
  Plot2 <- print(Plot1 + geom_text(color = "black", size = 1, check_overlap = T) + theme(aspect.ratio = 1))

Plot1 <- ggplot(df.integration.conpmx, aes(XScores.conpmx, YScores.conpmx, label = row.names(df.integration.conpmx),colour = log.skullCS)) + geom_point(alpha = 1, size=5) + geom_smooth(method = glm)
  Plot2 <- print(Plot1 + scale_colour_gradientn(colours = c("#638B95", "#AC7E5A", "#F4711F"))
               + geom_text(color = "black", size = 0.5, check_overlap = T)
               + theme(aspect.ratio = 1))



modules.LM.list.sinpmx <- c("A","A","A","A","A","A","A","A","A","A","A","B","B","B","B","B","B","B","B","B","B","B","B","B","B","A","B","A","A","B","A","A","A","A","A","A","A","A","A","B","B","B","B","B","B","B","A","B","A")
integration.sinpmx <- phylo.integration(gpa.todos.sinpremax$coords, partition.gp = modules.LM.list.sinpmx, phy = tree.amniota, iter = 999)
  summary(integration.sinpmx)
    plot1 <- plot(integration.sinpmx) 
    plot2 <- print(plot1 + theme(aspect.ratio = 1))

preds.pls.skull.block1 <- shape.predictor(integration.sinpmx$A1, integration.sinpmx$XScores[,1], 
                                          method = "PLS", predmin = min(integration.sinpmx$XScores[,1]),
                                          predmax = max(integration.sinpmx$XScores[,1]),  pref  = 0)
  write.csv(preds.pls.skull.block1$predmin, 'preds.min.block1.sinpmx.csv')
  write.csv(preds.pls.skull.block1$predmax, 'preds.max.block1.sinpmx.csv')
preds.pls.skull.block2 <- shape.predictor(integration.sinpmx$A2, integration.sinpmx$YScores[,1], 
                                          method = "PLS", predmin = min(integration.sinpmx$YScores[,1]),
                                          predmax = max(integration.sinpmx$YScores[,1]),  pref  = 0)
  write.csv(preds.pls.skull.block2$predmin, 'preds.min.block2.sinpmx.csv')
  write.csv(preds.pls.skull.block2$predmax, 'preds.max.block2.sinpmx.csv')

preds.min.pls.sinpmx <- read.csv("preds.min.pls.sinpmx.csv",header=T,sep=";",row.names=1)
preds.min.pls.sinpmx <- as.matrix(preds.min.pls.sinpmx)
preds.max.pls.sinpmx <- read.csv("preds.max.pls.sinpmx.csv",header=T,sep=";",row.names=1)
preds.max.pls.sinpmx <- as.matrix(preds.max.pls.sinpmx)

library(tibble)
per_lm_distance <- function (shape.data.min , shape.data.max) #para modificar los colores de los LMs, meter los codigos que gusten en cols1
{
  df <- as.data.frame(cbind(shape.data.min, shape.data.max))
  distances <- apply(df, 1, function(x) dist(matrix(x, nrow = 2, byrow = TRUE)))
  cols1 <- colorRampPalette(c("#1a2a6c", "#7B90E3", "#fdbb2d", "#b21f1f"))
  cols <- cols1(100)
  x = (log10(distances))
  xlims <- NULL
  tol <- 1e-06
  xlims <- range(x) + c(-tol, tol)
  nbin = 100
  breaks <- 0:nbin/nbin * (xlims[2] - xlims[1]) + xlims[1]
  whichColor <- function(p, cols, breaks) {
    i <- 1
    while (p >= breaks[i] && p > breaks[i + 1]) i <- i +
        1
    cols[i]
  }
  distancecolors <- sapply(x, whichColor, cols = cols, breaks = breaks)
  distance_table <- tibble(Per_Lm_Distance = distances, Log_Distance = x,Distance_Colors = distancecolors)
  return(distance_table)
}

my.distances.preds.pls.sinpmx <- per_lm_distance(shape.data.min = preds.min.pls.sinpmx , shape.data.max = preds.max.pls.sinpmx)
coords.original <- coords.todos.sinpremax[,,"Mephitis_mephitis"] #Ptilonorhynchus_violaceus , Mephitis_mephitis
open3d()
par3d(windowRect = c(0,0,500,500)) #para dimensiones del cuadro RGL
Sys.sleep(1)
mfrow3d(nr = 1, nc = 2, byrow = TRUE, sharedMouse = TRUE)
mesh.min <- tps3d(mesh.original2, refmat = coords.original, tarmat = preds.min.pls.sinpmx) #maps landmarks or a triangular mesh via thin plate spline based on a reference and a target configuration in 2D and 3D
shade3d(mesh.min, col= 8, alpha = 0.3) #dibuja la mesh con la forma predicha en la línea anterior #alpha cambia la transparencia de la surface
spheres3d(preds.min.pls.sinpmx, col = my.distances.preds.pls.sinpmx$Distance_Colors, radius = 0.01) #radius cambia el tamaño del LM
next3d()
mesh.max <- tps3d(mesh.original2, refmat = coords.original, tarmat = preds.max.pls.sinpmx, lambda = 0.5)
shade3d(mesh.max, col= 8, alpha = 0.3) 
spheres3d(preds.max.pls.sinpmx, col = my.distances.preds.pls.sinpmx$Distance_Colors, radius = 0.01)


df.integration.sinpmx <- data.frame(XScores.sinpmx = integration.sinpmx$XScores[,1], 
                                    YScores.sinpmx = integration.sinpmx$YScores[,1],
                                    log.skullCS = log(variables.amniotas$skull.CS), 
                                    log.cranialbase.CS = log(variables.amniotas$cranialbase.CS),
                                    log.brainconvhull = log(variables.amniotas$brain.convhull),
                                    headbrain.res = variables.amniotas$headbrain.res,
                                    skullbase.res = variables.amniotas$skullbase.res,
                                    foramen = variables.amniotas$foramen.angule,
                                    base.ang = variables.amniotas$base.angule,
                                    flexion = variables.amniotas$flexion.angule,
                                    sph.res = variables.amniotas$SPH1.res,
                                    clivus.foramen.ang = variables.amniotas$clivus.foramen.angle,
                                    facial.ang = variables.amniotas$facial.angle,
                                    clades.a = variables.amniotas$clades.a,
                                    clades.b = variables.amniotas$clades.b,
                                    clades.c = variables.amniotas$clades.c,
                                    clades.comp = variables.amniotas$clades.comp,
                                    row.names = row.names(variables.amniotas))

Plot1 <- ggplot(df.integration.sinpmx, aes(XScores.sinpmx, YScores.sinpmx, label = rownames(df.integration.sinpmx), colour = clades.comp))+ geom_point(alpha = 1, size = 5) + geom_smooth(method = glm) #+ geom_line(aes( x = reg.lines.x, y = reg.lines.y, color =  clades.b), size = 2, alpha = 1)
Plot2 <- print(Plot1 + geom_text(color = "black", size = 0.4, check_overlap = T) + theme(aspect.ratio = 1))

Plot1 <- ggplot(df.integration.sinpmx, aes(XScores.sinpmx, YScores.sinpmx, label = row.names(df.integration.sinpmx),colour = log.skullCS)) + geom_point(alpha = 1, size=5) #+ geom_smooth(method = glm)
Plot2 <- print(Plot1 + scale_colour_gradientn(colours = c("#638B95", "#AC7E5A", "#F4711F"))
               + geom_text(color = "black", size = 0.5, check_overlap = T)
               + theme(aspect.ratio = 1))



## FUNCION PARA CALCULAR VALORES DE INTEGRACIÓN POR NODOS

get.descendants <- function( phy , node ) { extract.clade( phy , node )$tip.label } #extract.clade te genera un arbol a partir de todos los tips de un feterminado nodo, podando todos los demás tips
subset.rows <- function( data , rownames ) { data[ rownames , ] } 

phy <- tree.amniota

modularity.per.node <- function( data , phy ) {
  nodes.temp <- unique( phy$edge[,1] ) #para generar un vector que contenga los valores de los nodos internos de la filogenia
  descendants.temp <- lapply( nodes.temp , FUN = get.descendants , phy = phy ) #esta línea devuelve una lista de la longitud de x que se llena de los resultados obtenidos a partir de aplicar la función FUN a cada uno de los elementos de x
  names( descendants.temp ) <- nodes.temp			
  node.data.list <- lapply( descendants.temp , FUN = subset.rows , data = data )
  modularity.temp<-list() 
  for (i in 1: length( node.data.list ) ) {
    node.2darray.data.list <- node.data.list[[i]]
    node.3darray.data.list <- arrayspecs(node.2darray.data.list[,1:150],50,3) #cambiar [,1:150] si la matriz de origen es la que contiene los clasificadores [,5:154]
    data.names <- dimnames(node.3darray.data.list) [[3]]
    drop.taxa <- phy$tip.label[ ! phy$tip.label %in% data.names ]
    prunned.phy <- drop.tip( phy , drop.taxa )
    pls<-phylo.integration(node.3darray.data.list, partition.gp = modules.LM.list.conpmx, phy = prunned.phy, iter = 999)
    modularity.temp[i]<- pls$Z                         
    
  }
  
  names(modularity.temp)<-nodes.temp
  modularity.temp<-unlist(modularity.temp)
  modularity.temp
  #modularity.temp.scaled<-modularity.temp / modularity.temp[1] ## activar si queremos que el valor para cada nodo esté escalado y se dé de manera relativa (sobre 1)
  #modularity.temp.scaled
  
}

modularity.per.node.temp.conpmx <- modularity.per.node( data = twodarray.proccoords.conpmx , phy = phy)
modularity.per.node.temp.absval.conpmx <- modularity.per.node( data = twodarray.proccoords.conpmx , phy = phy)

#FUNCIONA: hay que trabajar con una 2darray de inicio para que lea bien los rownames, 
#y luego hay que podar el arbol específicamente para cada nodo, 
#y que asi la funcion phylo.integration() pueda correr 
#sin dar el aviso de que hay más taxones en el arbol que el los datos

## COMO PLOTEAR LOS VALORES DE MODULARIDAD POR NODO EN CADA NODO DE LA FILOGENIA:
# PLS z.scores

tree.temp <- phy
z.scores <- modularity.per.node.temp.absval.conpmx #cambiar al elemento que corresponda


# prepares z.scores to be plotted
z.scores.mod <- z.scores[!is.na(z.scores)]# removes NaNs product of clade-tip daughter comparisons
z.scores.mod <- z.scores.mod[!duplicated(z.scores.mod)] # two daughters per clade are obtained but we want only to keep one to plot in the immediate parent node


#prunes modularity values to match only the represented nodes
zscores.point.size <- modularity.per.node.temp.absval.sinpmx [names (z.scores.mod)]
log.zscores.point.size <- log(zscores.point.size)
log.zscores.point.size

disparities.deviation.plot <- disparities.deviation.plot [names (z.scores.mod)]

# generates colours for angles 

edge.values <- z.scores
resta <- min(edge.values)
resta.positivo <- -1*resta
edge.values <- edge.values + resta.positivo
colours.rgb <- c("#ff0000","#8f0441","#1f0781") #rojo-azul ##because we want people to look at outliers from BMw, grey is a good colour for unwanted attention
colour.scale <- colorRamp( rev( colours.rgb ) , space = "rgb", bias = 0.6) #bias ensures grey corresponds always in good BM fitted values (adjust conveniently and check, more sophisticated in the final plots)
BG.node.table <- colour.scale( edge.values / max( edge.values ) )
BG.zscores<- c()
for( i in 1:nrow( BG.node.table ) ) {
  BG.zscores[i]<- rgb( BG.node.table[i,1], BG.node.table[i,2], BG.node.table[i,3] , alpha = 255 , max = 255 )
}	


# subset node.ages to data, do this before plotting z.scores 

tree.temp <- tree
variable = z.scores
tree <- tree.temp						

nodes.ages.temp <- tree.age(tree.temp) 
node.ages.temp.mod <- nodes.ages.temp[,1]
names(node.ages.temp.mod)<- as.character(nodes.ages.temp[,2])
node.ages.temp.mod<-node.ages.temp.mod[ names(variable)]						

#plot to check nodes

plot(z.scores, node.ages.temp.mod, xlim = c(-2, max(z.scores)), ylim = c(max(node.ages.temp.mod), 0) , type = "p", cex = 2, pch = 21, col = "transparent", bg = BG.zscores, frame.plot = F)
text( y = node.ages.temp.mod , x =  z.scores, labels = names( pls.angles.mod ), pos = 4, cex =0.7, col = "black")	      


#tree with z.scores

plot( tree.temp , type = "fan", show.tip.label = T, no.margin = T, y.lim =  c(-5, 72.907 ))
nodelabels( node = as.numeric(names(z.scores)), frame = "n" , cex = log.zscores.point.size , pch = 21 , bg = BG.zscores )





## NEW FUNCTION THAT COMPARES DAUGHTER NODES WITH IMMEDIATELY PARENTAL NODES ##

eigendirections.vertical.phylo <- function( data , phy ) {
  
  coords.matrix<-two.d.array(data)
  internal.node.pairs <- phy$edge[!phy$edge[, 2] <= length(phy$tip.label),] #subsets a matrix with all the node comparisons we want: first column inmmediately parent nodes, second column immediately daughter nodes
  
  descendants.parents.temp <- lapply( internal.node.pairs[,1] , FUN = get.descendants , phy = phy ) # this block makes the data list for parent clades
  names( descendants.parents.temp ) <- internal.node.pairs[,1 ]
  node.data.parents.list <- lapply( descendants.parents.temp , FUN = subset.rows , data = coords.matrix )
  
  descendants.daughters.temp <- lapply( internal.node.pairs[,2] , FUN = get.descendants , phy = phy ) #this block makes the ata list for daugther clades
  names( descendants.daughters.temp ) <- internal.node.pairs[, 2]
  node.data.daughters.list <- lapply( descendants.daughters.temp , FUN = subset.rows , data = coords.matrix )
  
  
  pls1.angles<-list()
  
  for (i in 1:dim(internal.node.pairs)[1] ) { # loop reconverts data to 3D array so phylo.integration can work
    
    parent.node <- as.character(internal.node.pairs[i, 1])
    daughter.node <- as.character(internal.node.pairs[i, 2])
    
    GM.node.data.parents.list<-arrayspecs(as.matrix(node.data.parents.list[[parent.node]]), dim(as.matrix(node.data.parents.list[[parent.node]]))[2]/3, 3)
    data.names <- dimnames(GM.node.data.parents.list) [[3]]
    drop.taxa <- phy$tip.label[ ! phy$tip.label %in% data.names ]
    prunned.phy <- drop.tip( phy , drop.taxa )
    pls1.parent.vector<-phylo.integration(GM.node.data.parents.list, partition.gp = modules.LM.list.conpmx, phy = prunned.phy, iter = 999)$right.pls.vectors[,1]
    
    GM.node.data.daughters.list<-arrayspecs(as.matrix(node.data.daughters.list[[daughter.node]]), dim(as.matrix(node.data.daughters.list[[daughter.node]]))[2]/3, 3)
    data.names <- dimnames(GM.node.data.daughters.list) [[3]]
    drop.taxa <- phy$tip.label[ ! phy$tip.label %in% data.names ]
    prunned.phy <- drop.tip( phy , drop.taxa )  
    pls1.daughter.vector<-phylo.integration(GM.node.data.daughters.list, partition.gp = modules.LM.list.conpmx, phy = prunned.phy, iter = 999)$right.pls.vectors[,1]
    
    pls1.angles[i]<-angle(pls1.parent.vector, pls1.daughter.vector, degree = T)
  }
  
  pls1.angles<-unlist(pls1.angles)
  names(pls1.angles)<- internal.node.pairs[,2]
  pls1.angles
  
  
}

#cambiar dentro de la funcion $left.pls.vectors[,1]
diff.pls1.left.angles.conpmx <- eigendirections.vertical.phylo( data = gpa.todos$coords , phy = phy )
#cambiar dentro de la funcion $right.pls.vectors[,1]
diff.pls1.right.angles.conpmx <- eigendirections.vertical.phylo( data = gpa.todos$coords , phy = phy )



## NEW FUNCTION THAT COMPARES DAUGHTER (SISTER) NODES BETWEEN THEM ##

eigendirections.horizontal.phylo <- function( data , phy ) {
  
  coords.matrix<-two.d.array(data)
  internal.node.pairs <- phy$edge[!phy$edge[, 2] <= length(phy$tip.label),] #subsets a matrix with all the node comparisons we want: first column inmmediately parent nodes, second column immediately daughter nodes
  
  descendants.daughters.temp <- lapply( internal.node.pairs[,2] , FUN = get.descendants , phy = phy ) #this block makes the ata list for daugther clades
  names( descendants.daughters.temp ) <- internal.node.pairs[, 2]
  node.data.daughters.list <- lapply( descendants.daughters.temp , FUN = subset.rows , data = coords.matrix )
  
  
  pls1.angles<-list()
  
  for (i in 1:length(internal.node.pairs[,1])) { # loop reconverts data to 3D array so phylo.integration can work
    
    node <- as.character(internal.node.pairs[i, 1])
    daughters <- internal.node.pairs[internal.node.pairs[,1] == node,]
    
    if (length(daughters) == 4 ) {
      
      
      daughter.1.node <- as.character(daughters[1, 2])
      daughter.2.node <- as.character(daughters[2, 2])
      
      GM.node.data.daughters.1.list<-arrayspecs(as.matrix(node.data.daughters.list[[daughter.1.node]]), dim(as.matrix(node.data.daughters.list[[daughter.1.node]]))[2]/3, 3)
      data.names <- dimnames(GM.node.data.daughters.1.list) [[3]]
      drop.taxa <- phy$tip.label[ ! phy$tip.label %in% data.names ]
      prunned.phy <- drop.tip( phy , drop.taxa )  
      pls1.daughter.1.vector<-phylo.integration(GM.node.data.daughters.1.list, partition.gp = modules.LM.list.sinpmx, phy = prunned.phy, iter = 999)$right.pls.vectors[,1]
      
      GM.node.data.daughters.2.list<-arrayspecs(as.matrix(node.data.daughters.list[[daughter.2.node]]), dim(as.matrix(node.data.daughters.list[[daughter.2.node]]))[2]/3, 3)
      data.names <- dimnames(GM.node.data.daughters.2.list) [[3]]
      drop.taxa <- phy$tip.label[ ! phy$tip.label %in% data.names ]
      prunned.phy <- drop.tip( phy , drop.taxa )  
      pls1.daughter.2.vector<-phylo.integration(GM.node.data.daughters.2.list, partition.gp = modules.LM.list.sinpmx, phy = prunned.phy, iter = 999)$right.pls.vectors[,1]
      
      
      
      pls1.angles[i]<-angle(pls1.daughter.1.vector, pls1.daughter.2.vector, degree = T)
      
    }
    
    else {
      
      daughter.1.node <- as.character(daughters[2])
      daughter.2.node <- as.character(daughters[2])
      
      GM.node.data.daughters.1.list<-arrayspecs(as.matrix(node.data.daughters.list[[daughter.1.node]]), dim(as.matrix(node.data.daughters.list[[daughter.1.node]]))[2]/3, 3)
      data.names <- dimnames(GM.node.data.daughters.1.list) [[3]]
      drop.taxa <- phy$tip.label[ ! phy$tip.label %in% data.names ]
      prunned.phy <- drop.tip( phy , drop.taxa )  
      pls1.daughter.1.vector<-phylo.integration(GM.node.data.daughters.1.list, partition.gp = modules.LM.list.sinpmx, phy = prunned.phy, iter = 999)$right.pls.vectors[,1]
      
      GM.node.data.daughters.2.list<-arrayspecs(as.matrix(node.data.daughters.list[[daughter.2.node]]), dim(as.matrix(node.data.daughters.list[[daughter.2.node]]))[2]/3, 3)
      data.names <- dimnames(GM.node.data.daughters.2.list) [[3]]
      drop.taxa <- phy$tip.label[ ! phy$tip.label %in% data.names ]
      prunned.phy <- drop.tip( phy , drop.taxa )  
      pls1.daughter.2.vector<-phylo.integration(GM.node.data.daughters.2.list, partition.gp = modules.LM.list.sinpmx, phy = prunned.phy, iter = 999)$right.pls.vectors[,1]
      
      pls1.angles[i]<-angle(pls1.daughter.1.vector, pls1.daughter.2.vector, degree = T)
      
    }
    
    
  }
  
  pls1.angles<-unlist(pls1.angles)
  names(pls1.angles)<- internal.node.pairs[,1]
  pls1.angles
  
  
}

#cambiar dentro de la funcion $left.pls.vectors[,1]
diff.pls1.left.sister.angles.conpmx <- eigendirections.horizontal.phylo( data = gpa.todos$coords , phy = phy )
#cambiar dentro de la funcion $right.pls.vectors[,1]
diff.pls1.right.sister.angles.conpmx <- eigendirections.horizontal.phylo( data = gpa.todos$coords , phy = phy )


#cambiamos ahora a sin premaxila (modificar -> partition.gp = modules.LM.list.sinpmx)
#cambiar dentro de la funcion $left.pls.vectors[,1]
diff.pls1.left.sister.angles.sinpmx <- eigendirections.horizontal.phylo( data = gpa.todos.sinpremax$coords , phy = phy )
#cambiar dentro de la funcion $right.pls.vectors[,1]
diff.pls1.right.sister.angles.sinpmx <- eigendirections.horizontal.phylo( data = gpa.todos.sinpremax$coords , phy = phy )



# PLS vector angles  between clades 
# PLS angles between daughters clades, plots value in parent clade
tree.temp <- phy
pls.angles <- diff.pls1.right.sister.angles.sinpmx #cambiar al elemento que corresponda

# prepares angles to be plotted
pls.angles.mod <- pls.angles[!is.na(pls.angles)]# removes NaNs product of clade-tip daughter comparisons
pls.angles.mod <- pls.angles.mod[!duplicated(pls.angles.mod)] # two daughters per clade are obtained but we want only to keep one to plot in the immediate parent node

#this block makes angles all angles below 90 degrees: because PLS vectors are arbitrary in direction the maximum angle is orthogonal	
pls.angles.mod
pls.angles.mod.wrong <-  pls.angles.mod [which(  pls.angles.mod > 90 )]
pls.angles.mod.corrected <-  (pls.angles.mod.wrong -180) * -1
pls.angles.mod[names(pls.angles.mod.corrected)] <- pls.angles.mod.corrected
pls.angles.mod<-pls.angles.mod[- which(pls.angles.mod >= 0 & pls.angles.mod <= 2)]

#prunes modularity values to match only the represented nodes
zscores.point.size <- modularity.per.node.temp.absval.sinpmx [names (pls.angles.mod)]
log.zscores.point.size <- log(zscores.point.size)
log.zscores.point.size

disparities.deviation.plot <- disparities.deviation.plot [names (pls.angles.mod)]

# generates colours for angles 
edge.values <- pls.angles.mod
colours.rgb <- c("#ff6200","#fff367","#e8e8e8") #naranja rojizo-amarillo-gris ##because we want people to look at outliers from BMw, grey is a good colour for unwanted attention
colour.scale <- colorRamp( rev( colours.rgb ) , space = "rgb", bias = 0.6) #bias ensures grey corresponds always in good BM fitted values (adjust conveniently and check, more sophisticated in the final plots)
BG.node.table <- colour.scale( edge.values / max( edge.values ) )
BG.angles<- c()
for( i in 1:nrow( BG.node.table ) ) {
  BG.angles[i]<- rgb( BG.node.table[i,1], BG.node.table[i,2], BG.node.table[i,3] , alpha = 255 , max = 255 )
}	


#gets nodes only for the daughter-clade comparisons, i.e., only picks daughter clades with 2 or more descendent tips
tree.temp <- tree
variable = pls.angles.mod
tree <- tree.temp						

nodes.ages.temp <- tree.age(tree.temp) 
node.ages.temp.mod <- nodes.ages.temp[,1]
names(node.ages.temp.mod)<- as.character(nodes.ages.temp[,2])
node.ages.temp.mod<-node.ages.temp.mod[ names(variable)]						

#plot to check nodes
plot(pls.angles.mod, node.ages.temp.mod, xlim = c(0, 100), ylim = c(120, 0) , type = "p", cex = log.zscores.point.size, pch = 21, bg = BG.angles, frame.plot = T)
text( y = node.ages.temp.mod , x =  pls.angles.mod, labels = names( pls.angles.mod ), pos = 4, cex =0.7, col = "black")	      

#tree with pls.angles
plot( tree.temp , type = "fan", show.tip.label = T, no.margin = T, y.lim =  c(-5, 72.907 ))
nodelabels( node = as.numeric(names(pls.angles.mod)), frame = "n" , cex = log.zscores.point.size , pch = 21 , bg = BG.angles )



## NEW FUNCTION THAT PERFORMS COMPARE.PLS BETWEEN DAUGHTER (SISTER) NODES ##

compare.pls.horizontal.phylo <- function( data , phy ) {
  
  coords.matrix<-two.d.array(data)
  internal.node.pairs <- phy$edge[!phy$edge[, 2] <= length(phy$tip.label),] #subsets a matrix with all the node comparisons we want: first column inmmediately parent nodes, second column immediately daughter nodes
  
  descendants.daughters.temp <- lapply( internal.node.pairs[,2] , FUN = get.descendants , phy = phy ) #this block makes the Data list for daugther clades
  names( descendants.daughters.temp ) <- internal.node.pairs[, 2]
  node.data.daughters.list <- lapply( descendants.daughters.temp , FUN = subset.rows , data = coords.matrix )
  
  
  pairwise.pls<-list()
  
  for (i in 1:length(internal.node.pairs[,1])) { # loop reconverts data to 3D array so phylo.integration can work
    
    node <- as.character(internal.node.pairs[i, 1])
    daughters <- internal.node.pairs[internal.node.pairs[,1] == node,]
    
    if (length(daughters) == 4 ) {
      
      daughter.1.node <- as.character(daughters[1, 2])
      daughter.2.node <- as.character(daughters[2, 2])
      
      GM.node.data.daughters.1.list<-arrayspecs(as.matrix(node.data.daughters.list[[daughter.1.node]]), dim(as.matrix(node.data.daughters.list[[daughter.1.node]]))[2]/3, 3)
      data.names <- dimnames(GM.node.data.daughters.1.list) [[3]]
      drop.taxa <- phy$tip.label[ ! phy$tip.label %in% data.names ]
      prunned.phy <- drop.tip( phy , drop.taxa )  
      pls.daughter.1<-phylo.integration(GM.node.data.daughters.1.list, partition.gp = modules.LM.list.sinpmx, phy = prunned.phy, iter = 999)
      
      GM.node.data.daughters.2.list<-arrayspecs(as.matrix(node.data.daughters.list[[daughter.2.node]]), dim(as.matrix(node.data.daughters.list[[daughter.2.node]]))[2]/3, 3)
      data.names <- dimnames(GM.node.data.daughters.2.list) [[3]]
      drop.taxa <- phy$tip.label[ ! phy$tip.label %in% data.names ]
      prunned.phy <- drop.tip( phy , drop.taxa )  
      pls.daughter.2<-phylo.integration(GM.node.data.daughters.2.list, partition.gp = modules.LM.list.sinpmx, phy = prunned.phy, iter = 999)
      
      
      pairwise.pls[i] <- compare.pls(pls.daughter.1, pls.daughter.2)$pairwise.z[1,2]
      
    }
    
    else {
      
      daughter.1.node <- as.character(daughters[2])
      daughter.2.node <- as.character(daughters[2])
      
      GM.node.data.daughters.1.list<-arrayspecs(as.matrix(node.data.daughters.list[[daughter.1.node]]), dim(as.matrix(node.data.daughters.list[[daughter.1.node]]))[2]/3, 3)
      data.names <- dimnames(GM.node.data.daughters.1.list) [[3]]
      drop.taxa <- phy$tip.label[ ! phy$tip.label %in% data.names ]
      prunned.phy <- drop.tip( phy , drop.taxa )  
      pls.daughter.1<-phylo.integration(GM.node.data.daughters.1.list, partition.gp = modules.LM.list.sinpmx, phy = prunned.phy, iter = 999)
      
      GM.node.data.daughters.2.list<-arrayspecs(as.matrix(node.data.daughters.list[[daughter.2.node]]), dim(as.matrix(node.data.daughters.list[[daughter.2.node]]))[2]/3, 3)
      data.names <- dimnames(GM.node.data.daughters.2.list) [[3]]
      drop.taxa <- phy$tip.label[ ! phy$tip.label %in% data.names ]
      prunned.phy <- drop.tip( phy , drop.taxa )  
      pls.daughter.2<-phylo.integration(GM.node.data.daughters.2.list, partition.gp = modules.LM.list.sinpmx, phy = prunned.phy, iter = 999)
      
      pairwise.pls[i] <- compare.pls(pls.daughter.1, pls.daughter.2)$pairwise.z[1,2]
      
    }
    
    
  }
  
  pairwise.pls<-unlist(pairwise.pls)
  names(pairwise.pls)<- internal.node.pairs[,1]
  pairwise.pls
  
}


pls.z.comparison.sisters.conpmx <- compare.pls.horizontal.phylo( data = gpa.todos$coords , phy = phy )
pls.P.comparison.sisters.conpmx <- compare.pls.horizontal.phylo( data = gpa.todos$coords , phy = phy ) #cambiar a $pairwise.P[1,2] en compare.pls

#cambiamos ahora a sin premaxila (modificar -> partition.gp = modules.LM.list.sinpmx)
pls.z.comparison.sisters.sinpmx <- compare.pls.horizontal.phylo( data = gpa.todos.sinpremax$coords , phy = phy )
pls.P.comparison.sisters.sinpmx <- compare.pls.horizontal.phylo( data = gpa.todos.sinpremax$coords , phy = phy ) #cambiar a $pairwise.P[1,2] en compare.pls                        



# PLS pairwise.Z between clades 
# PLS pairwise.Z between daughters clades, plots value in parent clade
tree.temp <- phy
pairwise.z <- pls.z.comparison.sisters.sinpmx #cambiar al elemento que corresponda
pairwise.P <- pls.P.comparison.sisters.sinpmx #cambiar al elemento que corresponda  


# prepares pairwise.z to be plotted discarding those whose p-value is higher than 0.05
pairwise.z.mod <- pairwise.z[!duplicated(pairwise.z)] # two daughters per clade are obtained but we want only to keep one to plot in the immediate parent node
pairwise.P.mod <- pairwise.P[!duplicated(pairwise.P)]

pairwise.pls.sig <- pairwise.P.mod [which(  pairwise.P.mod <= 0.05 )]
pairwise.pls.sig


#prunes pairwise.z values to match only the significative nodes
pairwise.z.sig <- pairwise.z.mod [names (pairwise.pls.sig)]
pairwise.z.sig

zscores.point.size <- modularity.per.node.temp.absval.sinpmx [names (pairwise.z.sig)]
log.zscores.point.size <- log(zscores.point.size)
log.zscores.point.size


##disparities.deviation.plot <- disparities.deviation.plot [names (pls.angles.mod)]

# generates colours for angles 
edge.values <- pairwise.z.sig
colours.rgb <- c("#f4d03f","#16a085") #naranja rojizo-amarillo-gris ##because we want people to look at outliers from BMw, grey is a good colour for unwanted attention
colour.scale <- colorRamp( rev( colours.rgb ) , space = "rgb", bias = 0.8) #bias ensures grey corresponds always in good BM fitted values (adjust conveniently and check, more sophisticated in the final plots)
BG.node.table <- colour.scale( edge.values / max( edge.values ) )
BG.pairwise.z <- c()
for( i in 1:nrow( BG.node.table ) ) {
  BG.pairwise.z[i]<- rgb( BG.node.table[i,1], BG.node.table[i,2], BG.node.table[i,3] , alpha = 255 , max = 255 )
}	


#gets nodes only for the daughter-clade comparisons, i.e., only picks daughter clades with 2 or more descendent tips
tree.temp <- tree
variable = pairwise.z.sig
tree <- tree.temp						

nodes.ages.temp <- tree.age(tree.temp) 
node.ages.temp.mod <- nodes.ages.temp[,1]
names(node.ages.temp.mod)<- as.character(nodes.ages.temp[,2])
node.ages.temp.mod<-node.ages.temp.mod[ names(variable)]						

#plot to check nodes
plot(pairwise.z.sig, node.ages.temp.mod, xlim = c(0, 8), ylim = c(120, 0) , type = "p", cex = log.zscores.point.size, pch = 21, bg = BG.pairwise.z, frame.plot = T)
text( y = node.ages.temp.mod , x =  pairwise.z.sig, labels = names( pairwise.z.sig ), pos = 4, cex =0.7, col = "black")	      


#tree with pls.angles
plot( tree.temp , type = "fan", show.tip.label = T, no.margin = T, y.lim =  c(-5, 72.907 ))
nodelabels( node = as.numeric(names(pairwise.z.sig)), frame = "n" , cex = log.zscores.point.size , pch = 21 , bg = BG.pairwise.z )



## MODULARIDAD POR NODO A PARTIR DE LOS DATOS SIMULADOS ##        
#BM
library(mvMORPH)

gpa.coords.temp <- gpa.todos$coords
gpa.coords.temp <- gpa.todos.sinpremax$coords 

fit1 <- mvgls( two.d.array(gpa.coords.temp) ~ 1 , tree = tree.temp , model = "BM" , method = "H&L" ) #primero ajustamos los datos al modelo evolutivo elegido, en este caso BM no direccional
sim.temp.new <- mvSIM(tree.temp, nsim = 20, model = "BM1", param = list(sigma = fit1$sigma$Pinv, theta = fit1$coefficients )) #con mvSIM simulamos nuestros datos ajustados a un modelo evolutivo determinado, en este caso BM no direccional
x <- array(unlist(sim.temp.new), c(381,147,20), dimnames = list(dimnames(gpa.coords.temp)[[3]],as.character(c(1:147)),as.character(c(1:20)))) #con esto transformamos la lista de las 500 matrices de 381x150 que hemos simulado, en una sola matriz de 3 dimensiones (especies x coordenadas x simulaciones)

modularity.per.node <- function( data , phy ) {
  nodes.temp <- unique( phy$edge[,1] ) #para generar un vector que contenga los valores de los nodos internos de la filogenia
  descendants.temp <- lapply( nodes.temp , FUN = get.descendants , phy = phy ) #esta línea devuelve una lista de la longitud de x que se llena de los resultados obtenidos a partir de aplicar la función FUN a cada uno de los elementos de x
  names( descendants.temp ) <- nodes.temp			
  node.data.list <- lapply( descendants.temp , FUN = subset.rows , data = data )
  modularity.temp<-list() 
  for (i in 1: length( node.data.list ) ) {
    node.2darray.data.list <- node.data.list[[i]]
    node.3darray.data.list <- arrayspecs(node.2darray.data.list[,1:147],49,3) #cambiar [,1:150] si la matriz de origen es la que contiene los clasificadores [,5:154]
    data.names <- dimnames(node.3darray.data.list) [[3]]
    drop.taxa <- phy$tip.label[ ! phy$tip.label %in% data.names ]
    prunned.phy <- drop.tip( phy , drop.taxa )
    pls<-phylo.integration(node.3darray.data.list, partition.gp = modules.LM.list.sinpmx, phy = prunned.phy, iter = 499)
    modularity.temp[i]<- pls$Z                         
    
  }
  
  names(modularity.temp)<-nodes.temp
  modularity.temp<-unlist(modularity.temp)
  modularity.temp
  #modularity.temp.scaled<-modularity.temp / modularity.temp[1] ## activar si queremos que el valor para cada nodo esté escalado y se dé de manera relativa (sobre 1)
  #modularity.temp.scaled
  
}

sim.modularities.nodes <- apply( x , 3 , modularity.per.node , phy = tree.temp ) 


### Generates sizes and symbols for deviations of empirical disparities from BM simulations

modularities.deviation <- modularity.per.node.temp.absval.sinpmx - apply (sim.modularities.nodes, 1, median)
max(modularities.deviation)
mean(modularities.deviation)
min(modularities.deviation)  


edge.values <-  modularities.deviation - min( modularities.deviation )
log.zscores.point.size <- log(disp.modul.nodes.red$modularity.sinpmx+1)
names(log.zscores.point.size) <- row.names(disp.modul.nodes.red)
log.zscores.point.size <- na.omit(log.zscores.point.size)
edge.values<-edge.values[ names(log.zscores.point.size)]
colours.rgb <- c("#C8963E", "#d5cfce",  "#55185D" ) ##because we want people to look at outliers from BMw, grey is a good colour for unwanted attention
colour.scale <- colorRamp( rev( colours.rgb ) , space = "rgb", bias = 1.1) #bias ensures grey corresponds always in good BM fitted values (adjust conveniently and check, more sophisticated in the final plots)
BG.node.table <- colour.scale( edge.values / max( edge.values ) )
BG.modul <- c()
for( i in 1:nrow( BG.node.table ) ) {
  BG.modul[i] <- rgb( BG.node.table[i,1], BG.node.table[i,2], BG.node.table[i,3] , alpha = 255 , max = 255 )
}	

#gets nodes only for the daughter-clade comparisons, i.e., only picks daughter clades with 2 or more descendent tips
tree.temp <- tree
variable = edge.values
tree <- tree.temp						

nodes.ages.temp <- tree.age(tree.temp) 
node.ages.temp.mod <- nodes.ages.temp[,1]
names(node.ages.temp.mod)<- as.character(nodes.ages.temp[,2])
node.ages.temp.mod<-node.ages.temp.mod[ names(log.zscores.point.size)]						

#plot to check nodes
plot(edge.values, node.ages.temp.mod, xlim = c(0, 7), ylim = c(325, 0) , type = "p", cex = log.zscores.point.size, pch = 21, bg = BG.modul, frame.plot = T)
text( y = node.ages.temp.mod , x =  edge.values, labels = names( edge.values ), pos = 4, cex =0.7, col = "black")	      


#tree with pls.angles
plot( tree.temp , type = "fan", show.tip.label = T, no.margin = T, y.lim =  c(-5, 72.907 ))
nodelabels( node = as.numeric(names(edge.values)), frame = "n" , cex = log.zscores.point.size , pch = 21 , bg = BG.modul ) 



## FUNCION PARA CALCULAR VALORES DE INTEGRACIÓN POR NODOS

get.descendants <- function( phy , node ) { extract.clade( phy , node )$tip.label } #extract.clade te genera un arbol a partir de todos los tips de un feterminado nodo, podando todos los demás tips
subset.rows <- function( data , rownames ) { data[ rownames , ] } 

phy <- tree.amniota

modularity.per.node <- function( data , phy ) {
  nodes.temp <- unique( phy$edge[,1] ) #para generar un vector que contenga los valores de los nodos internos de la filogenia
  descendants.temp <- lapply( nodes.temp , FUN = get.descendants , phy = phy ) #esta línea devuelve una lista de la longitud de x que se llena de los resultados obtenidos a partir de aplicar la función FUN a cada uno de los elementos de x
  names( descendants.temp ) <- nodes.temp			
  node.data.list <- lapply( descendants.temp , FUN = subset.rows , data = data )
  modularity.temp<-list() 
  for (i in 1: length( node.data.list ) ) {
    node.2darray.data.list <- node.data.list[[i]]
    node.3darray.data.list <- arrayspecs(node.2darray.data.list[,1:147],49,3)
    data.names <- dimnames(node.3darray.data.list) [[3]]
    drop.taxa <- phy$tip.label[ ! phy$tip.label %in% data.names ]
    prunned.phy <- drop.tip( phy , drop.taxa )
    pls<-phylo.integration(node.3darray.data.list, partition.gp = modules.LM.list.sinpmx, phy = prunned.phy, iter = 999)
    modularity.temp[i]<- pls$Z                         
    
  }
  
  names(modularity.temp)<-nodes.temp
  modularity.temp<-unlist(modularity.temp)
  modularity.temp
  #modularity.temp.scaled<-modularity.temp / modularity.temp[1] ## activar si queremos que el valor para cada nodo esté escalado y se dé de manera relativa (sobre 1)
  #modularity.temp.scaled
  
}

modularity.per.node.temp.sinpmx <- modularity.per.node( data = twodarray.proccoords.sinpmx , phy = phy)
modularity.per.node.temp.absval.sinpmx <- modularity.per.node( data = twodarray.proccoords.sinpmx[,5:151] , phy = phy)



## NEW FUNCTION THAT COMPARES DAUGHTER NODES WITH IMMEDIATELY PARENTAL NODES ##

eigendirections.vertical.phylo <- function( data , phy ) {
  
  coords.matrix<-two.d.array(data)
  internal.node.pairs <- phy$edge[!phy$edge[, 2] <= length(phy$tip.label),] #subsets a matrix with all the node comparisons we want: first column inmmediately parent nodes, second column immediately daughter nodes
  
  descendants.parents.temp <- lapply( internal.node.pairs[,1] , FUN = get.descendants , phy = phy ) # this block makes the data list for parent clades
  names( descendants.parents.temp ) <- internal.node.pairs[,1 ]
  node.data.parents.list <- lapply( descendants.parents.temp , FUN = subset.rows , data = coords.matrix )
  
  descendants.daughters.temp <- lapply( internal.node.pairs[,2] , FUN = get.descendants , phy = phy ) #this block makes the ata list for daugther clades
  names( descendants.daughters.temp ) <- internal.node.pairs[, 2]
  node.data.daughters.list <- lapply( descendants.daughters.temp , FUN = subset.rows , data = coords.matrix )
  
  
  pls1.angles<-list()
  
  for (i in 1:dim(internal.node.pairs)[1] ) { # loop reconverts data to 3D array so phylo.integration can work
    
    parent.node <- as.character(internal.node.pairs[i, 1])
    daughter.node <- as.character(internal.node.pairs[i, 2])
    
    GM.node.data.parents.list<-arrayspecs(as.matrix(node.data.parents.list[[parent.node]]), dim(as.matrix(node.data.parents.list[[parent.node]]))[2]/3, 3)
    data.names <- dimnames(GM.node.data.parents.list) [[3]]
    drop.taxa <- phy$tip.label[ ! phy$tip.label %in% data.names ]
    prunned.phy <- drop.tip( phy , drop.taxa )
    pls1.parent.vector<-phylo.integration(GM.node.data.parents.list, partition.gp = modules.LM.list.sinpmx, phy = prunned.phy, iter = 999)$left.pls.vectors[,1]
    
    GM.node.data.daughters.list<-arrayspecs(as.matrix(node.data.daughters.list[[daughter.node]]), dim(as.matrix(node.data.daughters.list[[daughter.node]]))[2]/3, 3)
    data.names <- dimnames(GM.node.data.daughters.list) [[3]]
    drop.taxa <- phy$tip.label[ ! phy$tip.label %in% data.names ]
    prunned.phy <- drop.tip( phy , drop.taxa )  
    pls1.daughter.vector<-phylo.integration(GM.node.data.daughters.list, partition.gp = modules.LM.list.sinpmx, phy = prunned.phy, iter = 999)$left.pls.vectors[,1]
    
    pls1.angles[i]<-angle(pls1.parent.vector, pls1.daughter.vector, degree = T)
  }
  
  pls1.angles<-unlist(pls1.angles)
  names(pls1.angles)<- internal.node.pairs[,2]
  pls1.angles
  
  
}

#cambiar dentro de la funcion $left.pls.vectors[,1]
diff.pls1.left.angles.sinpmx <- eigendirections.vertical.phylo( data = gpa.todos.sinpremax$coords , phy = phy )
#cambiar dentro de la funcion $right.pls.vectors[,1]
diff.pls1.right.angles.sinpmx <- eigendirections.vertical.phylo( data = gpa.todos.sinpremax$coords , phy = phy )




#######################
##### REGRESIONES #####
#######################

##PHYLOGENETIC REGRESSION (PGLS)## OKKKKKKK
#AVES#   
gdf.pgls<-geomorph.data.frame(phy = tree.amniota, 
                              coords = gpa.todos$coords,
                              log.skullCS = log(variables.amniotas$skull.CS), 
                              log.cranialbase.CS = log(variables.amniotas$cranialbase.CS),
                              log.brainconvhull = log(variables.amniotas$brain.convhull),
                              headbrain.res = variables.amniotas$headbrain.res,
                              skullbase.res = variables.amniotas$skullbase.res,
                              foramen = variables.amniotas$foramen.angule,
                              base.ang = variables.amniotas$base.angule,
                              flexion = variables.amniotas$flexion.angule,
                              sph.res = variables.amniotas$SPH1.res,
                              clivus.foramen.ang = variables.amniotas$clivus.foramen.angle,
                              facial.ang = variables.amniotas$facial.angle,
                              clades.a = variables.amniotas$clades.a,
                              clades.b = variables.amniotas$clades.b,
                              clades.c = variables.amniotas$clades.c,
                              clades.comp = variables.amniotas$clades.comp,
                              row.names = row.names(variables.amniotas))

gdf.pgls.sinpmx<-geomorph.data.frame(phy = tree.amniota, 
                                     coords = gpa.todos.sinpremax$coords,
                                     log.skullCS = log(variables.amniotas$skull.CS.sinpmx), 
                                     log.cranialbase.CS = log(variables.amniotas$cranialbase.CS),
                                     log.brainconvhull = log(variables.amniotas$brain.convhull),
                                     headbrain.res = variables.amniotas$headbrain.res.conpmx,
                                     skullbase.res = variables.amniotas$skullbase.res,
                                     foramen = variables.amniotas$foramen.angule,
                                     base.ang = variables.amniotas$base.angule,
                                     flexion = variables.amniotas$flexion.angule,
                                     sph.res = variables.amniotas$SPH1.res,
                                     clivus.foramen.ang = variables.amniotas$clivus.foramen.angle,
                                     facial.ang = variables.amniotas$facial.angle,
                                     clades.a = variables.amniotas$clades.a,
                                     clades.b = variables.amniotas$clades.b,
                                     clades.c = variables.amniotas$clades.c,
                                     clades.comp = variables.amniotas$clades.comp,
                                     row.names = row.names(variables.amniotas))

df.plots<-data.frame(row.names = row.names(variables.zombie),
                     superorders=superorders,
                     ence.coeff.BM=variables.zombie$brain.res,
                     IRE_vol=variables.zombie$IRE_vol,
                     IRE_mid_large=variables.zombie$IRE_classic,
                     IRE_mid_short=variables.zombie$IRE_midplane,
                     Csize=variables.zombie$cs_skull,
                     RegScores.pgls.ence.IREok=RegScores.pgls.ence.IREok,
                     RegScores.pgls.allom=RegScores.pgls.allom)


##CEREBRALIZACION Y CEFALIZACION##  

x <- log(variables.amniotas$body.mass)
names(x) <- row.names(variables.amniotas)
y <- log(variables.amniotas$brain.convhull)
names(y) <- row.names(variables.amniotas)
phylores.encephalization <- phyl.resid(tree.amniota, x, y, method = "BM")
summary(phylores.encephalization)

x <- log(variables.aves$body.mass)
names(x) <- row.names(variables.aves)
y <- log(variables.aves$skull.CS.sinpremax)
names(y) <- row.names(variables.aves)
pgls.aves.headres <- phyl.resid(tree.aves, x, y, method = "BM")
summary(pgls.aves.headres)    

x <- log(variables.amniotas$skull.CS)
names(x) <- row.names(variables.amniotas)
y <- log(variables.amniotas$brain.convhull)
names(y) <- row.names(variables.amniotas)
phylresid.headbrainres <- phyl.resid(tree.amniota, x, y, method = "BM") #### ESTA PARECE LA BUENA, LA ENCEFALIZACION ES REAL

x <- log(variables.amniotas$skull.CS.sinpmx)
names(x) <- row.names(variables.amniotas)
y <- log(variables.amniotas$brain.convhull)
names(y) <- row.names(variables.amniotas)
phylresid.headbrainres.sinpmx <- phyl.resid(tree.amniota, x, y, method = "BM")

x <- log(variables.amniotas$skull.CS.sinpmx)
names(x) <- row.names(variables.amniotas)
y <- log(variables.amniotas$cranialbase.CS)
names(y) <- row.names(variables.amniotas)
phylresid.skullbaseres <- phyl.resid(tree.amniota, x, y, method = "BM") #### TAMAÑO RELATIVO DE LA BASE CRANEAL DENTRO DEL CRANEO

x <- variables.amniotas$headbrain.res
names(x) <- row.names(variables.amniotas)
y <- variables.amniotas$skullbase.res
names(y) <- row.names(variables.amniotas)
phylresid.SPH <- phyl.resid(tree.amniota, x, y, method = "BM")

x <- log(variables.amniotas$foramen.CS)
names(x) <- row.names(variables.amniotas)
y1 <- log(variables.amniotas$skull.CS.sinpmx)
y2 <- log(variables.amniotas$brain.convhull)
names(y2) <- row.names(variables.amniotas)
phylresid.foram <- phyl.resid(tree.amniota, x, y1, method = "BM")
phylresid.foram.brain <- phyl.resid(tree.amniota, x, y2, method = "BM")


##ANOVAS ALOMETRIAS CRANEALES##
pgls.skullCS <- procD.pgls(coords ~ log.skullCS, phy = phy, data = gdf.pgls, iter = 9999, SS.type = "II")
plot(pgls.skullCS, type = "regression", reg.type = "RegScore", predictor = gdf.pgls$log.skullCS)
pgls.skullCS.sinpmx <- procD.pgls(coords ~ log.skullCS, phy = phy, data = gdf.pgls.sinpmx, iter = 9999, SS.type = "II")  
plot.pgls.skullCS.sinpmx<-plot(pgls.skullCS.sinpmx, type = "regression", reg.type = "RegScore", predictor = gdf.pgls.sinpmx$log.skullCS) ## ploteas regscores para ver la distribucion general
RegScores.skullCS.sinpmx<-plot.pgls.skullCS.sinpmx$RegScore

##ANOVAS ALOMETRIAS CRANEO~CEREBRO##
pgls.headbrain <- procD.pgls(coords ~ headbrain.res, phy = phy, data = gdf.pgls, iter = 9999, SS.type = "II")
plot(pgls.headbrain, type = "regression", reg.type = "RegScore", predictor = gdf.pgls$headbrain.res)
pgls.headbrain.sinpmx <- procD.pgls(coords ~ headbrain.res, phy = phy, data = gdf.pgls.sinpmx, iter = 9999, SS.type = "II")  
plot.pgls.headbrain.sinpmx<-plot(pgls.headbrain.sinpmx, type = "regression", reg.type = "RegScore", predictor = gdf.pgls.sinpmx$headbrain.res)
RegScores.headbrain.sinpmx<-plot.pgls.headbrain.sinpmx$RegScore

pgls.skullbase <- procD.pgls(coords ~ skullbase.res, phy = phy, data = gdf.pgls, iter = 9999, SS.type = "II")
plot(pgls.skullbase, type = "regression", reg.type = "RegScore", predictor = gdf.pgls$skullbase.res)
pgls.skullbase.sinpmx <- procD.pgls(coords ~ skullbase.res, phy = phy, data = gdf.pgls.sinpmx, iter = 9999, SS.type = "II")  
plot.pgls.skullbase.sinpmx<-plot(pgls.skullbase.sinpmx, type = "regression", reg.type = "RegScore", predictor = gdf.pgls.sinpmx$skullbase.res)
RegScores.skullbase.sinpmx<-plot.pgls.skullbase.sinpmx$RegScore


##ANGULOS Y ORIENTACIONES DE FORAMEN Y BASE    
pgls.foramen.ang <- procD.pgls(coords ~ foramen, phy = phy, data = gdf.pgls, iter = 9999, SS.type = "II")
plot(pgls.foramen.ang, type = "regression", reg.type = "RegScore", predictor = gdf.pgls$foramen)
pgls.foramen.ang.sinpmx <- procD.pgls(coords ~ foramen, phy = phy, data = gdf.pgls.sinpmx, iter = 9999, SS.type = "II")  
plot.pgls.foramen.sinpmx<-plot(pgls.foramen.ang.sinpmx, type = "regression", reg.type = "RegScore", predictor = gdf.pgls.sinpmx$foramen)

pgls.base.ang <- procD.pgls(coords ~ base.ang, phy = phy, data = gdf.pgls, iter = 9999, SS.type = "II")
plot(pgls.base.ang, type = "regression", reg.type = "RegScore", predictor = gdf.pgls$base.ang)
pgls.base.ang.sinpmx <- procD.pgls(coords ~ base.ang, phy = phy, data = gdf.pgls.sinpmx, iter = 9999, SS.type = "II")  
plot(pgls.base.ang.sinpmx, type = "regression", reg.type = "RegScore", predictor = gdf.pgls.sinpmx$base.ang)

pgls.flexion <- procD.pgls(coords ~ flexion, phy = phy, data = gdf.pgls, iter = 9999, SS.type = "II")
plot(pgls.flexion, type = "regression", reg.type = "RegScore", predictor = gdf.pgls$flexion)
pgls.flexion.sinpmx <- procD.pgls(coords ~ flexion, phy = phy, data = gdf.pgls.sinpmx, iter = 9999, SS.type = "II")  
plot(pgls.flexion.sinpmx, type = "regression", reg.type = "RegScore", predictor = gdf.pgls.sinpmx$flexion)    

pgls.clivus.foramen.ang.sinpmx <- procD.pgls(coords ~ clivus.foramen.ang, phy = phy, data = gdf.pgls.sinpmx, iter = 9999, SS.type = "II")  
plot.pgls.clivus.foramen.ang.sinpmx<-plot(pgls.clivus.foramen.ang.sinpmx, type = "regression", reg.type = "RegScore", predictor = gdf.pgls.sinpmx$clivus.foramen.ang)    
RegScores.clivus.foramen.ang.sinpmx<-plot.pgls.clivus.foramen.ang.sinpmx$RegScore


##ANOVAS CON SPH -> INFLUENCIA SOBRE LA FORMA DEL CRANEO QUE TIENE EL CEREBRO AL CRECER DENTRO DEL CRANEO (SPATIAL PACKING INDEX)
pgls.SPH <- procD.pgls(coords ~ sph.res, phy = phy, data = gdf.pgls, iter = 9999, SS.type = "II")
plot(pgls.SPH, type = "regression", reg.type = "RegScore", predictor = gdf.pgls$sph.res)
pgls.SPH.sinpmx <- procD.pgls(coords ~ sph.res, phy = phy, data = gdf.pgls.sinpmx, iter = 9999, SS.type = "II")  
plot.pgls.SPH.sinpmx<-plot(pgls.SPH.sinpmx, type = "regression", reg.type = "RegScore", predictor = gdf.pgls.sinpmx$sph.res)


##MODELOS Y ANOVAS DE MODELOS

pgls.modelos.indep <- procD.pgls(coords ~ headbrain.res+skullbase.res+foramen, phy = phy, data = gdf.pgls.sinpmx, iter = 9999, SS.type = "II")
pgls.modelos.interac <- procD.pgls(coords ~ headbrain.res*skullbase.res*foramen, phy = phy, data = gdf.pgls.sinpmx, iter = 9999, SS.type = "II")

anova(pgls.modelos.indep,pgls.modelos.interac)


pgls.aves.SPH.indep <- procD.pgls(coords ~ log.skullCS+headbrain.res+skullbase.res, phy = phy, data = gdf.aves.sinpremax.pgls, iter = 9999, SS.type = "II")

pgls.aves.SPH.indep2 <- procD.pgls(coords ~ log.skullCS+headbrain.res+skullbase.res+foramen.angule, phy = phy, data = gdf.aves.sinpremax.pgls, iter = 9999, SS.type = "II")


##ANOVAS ALOMETRIAS POR GRUPOS##
#skullCS
#modelo complejo, el que incluye interaccion de las variables
pgls.skullCS.sinpremax.unique <- procD.pgls(coords ~ log.skullCS*clades.b, phy = phy, data = gdf.pgls.sinpmx, iter = 999, SS.type = "II")
pgls.skullCS.sinpremax.unique$aov.table
plot.pgls.skullCS.sinpremax<-plot(pgls.skullCS.sinpremax.unique, type="regression", reg.type="RegScore", predictor = gdf.pgls.sinpmx$log.skullCS) ## ploteas regscores para ver la distribucion general
RegScores.skullCS.sinpremax<-plot.pgls.skullCS.sinpremax$RegScore
comparisons.allometry<-pairwise(pgls.skullCS.sinpremax.unique, groups = gdf.pgls.sinpmx$clades.b, covariate = gdf.pgls.sinpmx$log.brainconvhull)
summ.comparisons.allom<-summary(comparisons.allometry, test.type = "dist") ##cambiar entre "DL"(para diferencia absoluta entre slopes), "VC"(para correlación entre slopes) y "var"(para la dispersion de cada slope)
summ.comparisons.allom
pred.pgls.skullCS.sinpremax<-plot(pgls.skullCS.sinpremax.unique, type="regression", reg.type="PredLine", predictor = gdf.pgls.sinpmx$log.skullCS)
PredLine.skullCS <- pred.pgls.skullCS.sinpremax$PredLine  
reg.lines.x <- as.vector(pred.pgls.skullCS.sinpremax$plot.args$x)
names(reg.lines.x) <- gdf.pgls.sinpmx$row.names
reg.lines.y <- pred.pgls.skullCS.sinpremax$plot.args$y
names(reg.lines.y) <- gdf.pgls.sinpmx$row.names
df.skullCS.plots<-data.frame(RegScores.skullCS.sinpremax = RegScores.skullCS.sinpremax, 
                             log.skullCS = gdf.pgls.sinpmx$log.skullCS,
                             clades.b = gdf.pgls.sinpmx$clades.b,
                             reg.lines.x = reg.lines.x, reg.lines.y = reg.lines.y)
rownames(df.skullCS.plots) <- row.names(variables.amniotas)
Plot1 <- ggplot(df.skullCS.plots, aes(log.skullCS, RegScores.skullCS.sinpremax, label = rownames(df.skullCS.plots), colour = clades.b))+ geom_point(alpha = 1, size = 5) + geom_smooth(method = glm)#+ geom_line(aes( x = reg.lines.x, y = reg.lines.y, color =  clades.b), size = 2, alpha = 1)
Plot2 <- print(Plot1+ geom_text(color = "black", size = 1, check_overlap = T))

#modelo nulo      
pgls.skullCS.sinpremax.common <- procD.pgls(coords ~ log.skullCS+clades.comp, phy = phy, data = gdf.pgls.sinpmx, iter = 999, SS.type = "II")
pgls.skullCS.sinpremax$aov.table
plot.pgls.skullCS.sinpremax<-plot(pgls.skullCS.sinpremax, type="regression", reg.type="RegScore", predictor = gdf.pgls.sinpmx$log.skullCS) ## ploteas regscores para ver la distribucion general
RegScores.skullCS.sinpremax<-plot.pgls.skullCS.sinpremax$RegScore
comparisons.allometry<-pairwise(pgls.skullCS.sinpremax, groups = gdf.pgls.sinpmx$clades.comp, covariate = gdf.pgls.sinpmx$log.skullCS)
summ.comparisons.allom<-summary(comparisons.allometry, test.type = "DL") ##cambiar entre "DL"(para diferencia absoluta entre slopes), "VC"(para correlación entre slopes) y "var"(para la dispersion de cada slope)
summ.comparisons.allom
pred.pgls.skullCS.sinpremax<-plot(pgls.skullCS.sinpremax, type="regression", reg.type="PredLine", predictor = gdf.pgls.sinpmx$log.skullCS)
PredLine.skullCS <- pred.pgls.skullCS.sinpremax$PredLine  
reg.lines.x <- as.vector(pred.pgls.skullCS.sinpremax$plot.args$x)
names(reg.lines.x) <- gdf.pgls.sinpmx$row.names
reg.lines.y <- pred.pgls.skullCS.sinpremax$plot.args$y
names(reg.lines.y) <- gdf.pgls.sinpmx$row.names
df.skullCS.plots<-data.frame(RegScores.skullCS.sinpremax = RegScores.skullCS.sinpremax, 
                             log.skullCS = gdf.pgls.sinpmx$log.skullCS,
                             clades.comp = gdf.pgls.sinpmx$clades.comp,
                             reg.lines.x = reg.lines.x, reg.lines.y = reg.lines.y)
rownames(df.skullCS.plots) <- row.names(variables.amniotas)
Plot1 <- ggplot(df.skullCS.plots, aes(log.skullCS, RegScores.skullCS.sinpremax, label = rownames(df.skullCS.plots), colour = clades.comp))+ geom_point(alpha = 1, size = 5) + geom_smooth(method = lm)#+ geom_line(aes( x = reg.lines.x, y = reg.lines.y, color =  clades.comp), size = 2, alpha = 1)
Plot2 <- print(Plot1+ geom_text(color = "black", size = 1, check_overlap = T))

#comparativa de modelos para ver cual es mas explicativo
anova(pgls.skullCS.sinpremax.common,pgls.skullCS.sinpremax.unique,print.progress = FALSE)



#brainconvhull
pgls.aves.brainconvhull.sinpremax <- procD.pgls(coords ~ log.brainconvhull*clades.comp, phy = phy, data = gdf.aves.sinpremax.pgls, iter = 9999, SS.type = "II")
pgls.aves.brainconvhull.sinpremax$aov.table
plot.pgls.aves.brainconvhull.sinpremax<-plot(pgls.aves.brainconvhull.sinpremax, type="regression", reg.type="RegScore", predictor = gdf.aves.sinpremax.pgls$log.brainconvhull) ## ploteas regscores para ver la distribucion general
RegScores.aves.brainconvhull.sinpremax<-plot.pgls.aves.brainconvhull.sinpremax$RegScore
comparisons.brain.aves<-pairwise(pgls.aves.brainconvhull.sinpremax, groups = gdf.aves.sinpremax.pgls$clades.comp, covariate = gdf.aves.sinpremax.pgls$log.brainconvhull)
summ.comparisons.brain.aves<-summary(comparisons.brain.aves)
summ.comparisons.brain.aves
pred.pgls.aves.brainconvhull.sinpremax<-plot(pgls.aves.brainconvhull.sinpremax, type="regression", reg.type="PredLine", predictor = gdf.aves.sinpremax.pgls$log.brainconvhull)
PredLine.aves.brainconvhull <- pred.pgls.aves.brainconvhull.sinpremax$PredLine  
reg.lines.x <- as.vector(pred.pgls.aves.brainconvhull.sinpremax$plot.args$x)
names(reg.lines.x) <- gdf.aves.sinpremax.pgls$row.names
reg.lines.y <- pred.pgls.aves.brainconvhull.sinpremax$plot.args$y
names(reg.lines.y) <- gdf.aves.sinpremax.pgls$row.names
df.aves.brainconvhull.plots<-data.frame(RegScores.aves.brainconvhull.sinpremax = RegScores.aves.brainconvhull.sinpremax, 
                                        log.brainconvhull = gdf.aves.sinpremax.pgls$log.brainconvhull,
                                        clades.e = gdf.aves.sinpremax.pgls$clades.e,
                                        reg.lines.x = reg.lines.x, reg.lines.y = reg.lines.y)
rownames(df.aves.brainconvhull.plots) <- row.names(variables.aves)
Plot1 <- ggplot(df.aves.brainconvhull.plots, aes(log.brainconvhull, RegScores.aves.brainconvhull.sinpremax, label = rownames(df.aves.brainconvhull.plots), colour = clades.e)) + geom_line(aes( x = reg.lines.x, y = reg.lines.y, color =  clades.e), size = 2, alpha = 1)
Plot2 <- print(Plot1)


#headbrain.res -> INFLUENCIA SOBRE LA FORMA DEL CRANEO QUE TIENE EL CEREBRO AL CRECER DENTRO DEL CRANEO (SPATIAL PACKING INDEX)
pgls.brain.sinpremax.unique <- procD.pgls(coords ~ headbrain.res*clades.b, phy = phy, data = gdf.pgls.sinpmx, iter = 999, SS.type = "II")
pgls.headbrain.sinpremax.unique$aov.table
plot.pgls.headbrain.sinpremax.unique<-plot(pgls.headbrain.sinpremax.unique, type="regression", reg.type="RegScore", predictor = gdf.pgls.sinpmx$headbrain.res) ## ploteas regscores para ver la distribucion general
RegScores.headbrain.sinpremax<-plot.pgls.headbrain.sinpremax.unique$RegScore
comparisons.brain<-pairwise(pgls.brain.sinpremax.unique, groups = gdf.pgls.sinpmx$clades.b, covariate = gdf.pgls.sinpmx$headbrain.res)
summ.comparisons.brain<-summary(comparisons.brain, test.type = "VC") ##cambiar entre "DL"(para diferencia absoluta entre slopes), "VC"(para correlación entre slopes) y "var"(para la dispersion de cada slope)
summ.comparisons.brain
pred.pgls.headbrain.sinpremax<-plot(pgls.headbrain.sinpremax.unique, type="regression", reg.type="PredLine", predictor = gdf.pgls.sinpmx$headbrain.res)
PredLine.headbrain <- pred.pgls.headbrain.sinpremax$PredLine  
reg.lines.x <- as.vector(pred.pgls.headbrain.sinpremax$plot.args$x)
names(reg.lines.x) <- gdf.pgls.sinpmx$row.names
reg.lines.y <- pred.pgls.headbrain.sinpremax$plot.args$y
names(reg.lines.y) <- gdf.pgls.sinpmx$row.names
df.headbrain.plots<-data.frame(RegScores.headbrain.sinpremax = RegScores.headbrain.sinpremax, 
                               headbrain.res = gdf.pgls.sinpmx$headbrain.res,
                               clades.b = gdf.pgls.sinpmx$clades.b,
                               reg.lines.x = reg.lines.x, reg.lines.y = reg.lines.y)
rownames(df.headbrain.plots) <- row.names(variables.amniotas)
Plot1 <- ggplot(df.headbrain.plots, aes(headbrain.res, RegScores.headbrain.sinpremax, label = rownames(df.headbrain.plots), colour = clades.b))+ geom_point(alpha = 1, size = 5) + geom_smooth(method = glm)#+ geom_line(aes( x = reg.lines.x, y = reg.lines.y, color =  clades.b), size = 2, alpha = 1)
Plot2 <- print(Plot1+ geom_text(color = "black", size = 1, check_overlap = T))

#comparativa de modelos para ver cual es mas explicativo    
pgls.headbrain.sinpremax.common <- procD.pgls(coords ~ headbrain.res*clades.b, phy = phy, data = gdf.pgls.sinpmx, iter = 999, SS.type = "II")    
pgls.headbrain.sinpremax.common$aov.table    
plot.pgls.headbrain.sinpremax.common<-plot(pgls.headbrain.sinpremax.common, type="regression", reg.type="RegScore", predictor = gdf.pgls.sinpmx$headbrain.res) ## ploteas regscores para ver la distribucion general

anova(pgls.headbrain.sinpremax.common,pgls.headbrain.sinpremax.unique,print.progress = FALSE)  



#clivus.foramen.ang -> INFLUENCIA SOBRE LA FORMA DEL CRANEO QUE TIENE EL FORAMEN AL FLEXIONARSE
pgls.clivus.foramen.ang.unique <- procD.pgls(coords ~ clivus.foramen.ang*clades.b, phy = phy, data = gdf.pgls.sinpmx, iter = 999, SS.type = "II")
pgls.clivus.foramen.ang.unique$aov.table
plot.pgls.clivus.foramen.ang.unique<-plot(pgls.clivus.foramen.ang.unique, type="regression", reg.type="RegScore", predictor = gdf.pgls.sinpmx$clivus.foramen.ang) ## ploteas regscores para ver la distribucion general
RegScores.clivus.foramen.ang<-plot.pgls.clivus.foramen.ang.unique$RegScore
comparisons.clivus.foramen.ang<-pairwise(pgls.clivus.foramen.ang.unique, groups = gdf.pgls.sinpmx$clades.comp, covariate = gdf.pgls.sinpmx$clivus.foramen.ang)
summ.comparisons.clivus.foramen.ang<-summary(comparisons.clivus.foramen.ang, test.type = "VC") ##cambiar entre "DL"(para diferencia absoluta entre slopes), "VC"(para correlación entre slopes) y "var"(para la dispersion de cada slope)
summ.comparisons.clivus.foramen.ang
pred.pgls.clivus.foramen.ang<-plot(pgls.clivus.foramen.ang.unique, type="regression", reg.type="PredLine", predictor = gdf.pgls.sinpmx$clivus.foramen.ang)
PredLine.clivus.foramen.ang <- pred.pgls.clivus.foramen.ang$PredLine  
reg.lines.x <- as.vector(pred.pgls.clivus.foramen.ang$plot.args$x)
names(reg.lines.x) <- gdf.pgls.sinpmx$row.names
reg.lines.y <- pred.pgls.clivus.foramen.ang$plot.args$y
names(reg.lines.y) <- gdf.pgls.sinpmx$row.names
df.clivus.foramen.ang.plots<-data.frame(RegScores.clivus.foramen.ang = RegScores.clivus.foramen.ang, 
                                        clivus.foramen.ang = gdf.pgls.sinpmx$clivus.foramen.ang,
                                        clades.comp = gdf.pgls.sinpmx$clades.comp,
                                        reg.lines.x = reg.lines.x, reg.lines.y = reg.lines.y)
rownames(df.clivus.foramen.ang.plots) <- row.names(variables.amniotas)
Plot1 <- ggplot(df.clivus.foramen.ang.plots, aes(clivus.foramen.ang, RegScores.clivus.foramen.ang, label = rownames(df.clivus.foramen.ang.plots), colour = clades.comp))+ geom_point(alpha = 1, size = 5) + geom_smooth(method = glm)#+ geom_line(aes( x = reg.lines.x, y = reg.lines.y, color =  clades.a), size = 2, alpha = 1)
Plot2 <- print(Plot1+ geom_text(color = "black", size = 1, check_overlap = T))

#comparativa de modelos para ver cual es mas explicativo    
pgls.clivus.foramen.ang.common <- procD.pgls(coords ~ clivus.foramen.ang+clades.comp, phy = phy, data = gdf.pgls.sinpmx, iter = 999, SS.type = "II")    
anova(pgls.clivus.foramen.ang.common,pgls.clivus.foramen.ang.unique,print.progress = FALSE)  



################################
#### PLOTTING SHAPE CHANGES ####
################################

#checking extreme shapes for main axes

GPA.fit <- gpa.todos
GPA.fit.sinpmx <- gpa.todos.sinpremax

PCA.fit <- phyloPCA.todos
PCA.fit.sinpmx <- phyloPCA.todos.sinpmx

plot(PCA.fit, axis1 = 1, axis2 = 2)

PC <- 3 #meter el pc sobre el que queramos hacer los morphing
PC.temp<-phyloPCs.sinpmx[,PC]
preds<-shape.predictor(GPA.fit.sinpmx$coords, x = as.numeric( phyloPCs.sinpmx[,PC] ), Intercept = F, 
                       pred1 = quantile(PC.temp, probs = 0.05), pred2 = quantile(PC.temp, probs = 0.95)) #modificar los valores del probs para exagerar o suavizar los cambios de forma

preds.allometry<-shape.predictor(pgls.skullCS.sinpmx$GM$pgls.fitted, x = plot.pgls.skullCS.sinpmx$RegScore, Intercept = F, 
                                 predmin = min(plot.pgls.skullCS.sinpmx$RegScore), predmax = max(plot.pgls.skullCS.sinpmx$RegScore))

preds.headbrain.res<-shape.predictor(pgls.headbrain.sinpmx$GM$pgls.fitted, x = plot.pgls.headbrain.sinpmx$RegScore, Intercept = F, 
                                     predmin = min(plot.pgls.headbrain.sinpmx$RegScore), predmax = max(plot.pgls.headbrain.sinpmx$RegScore))     

preds.clivus.foramen.ang<-shape.predictor(pgls.clivus.foramen.ang.sinpmx$GM$pgls.fitted, x = plot.pgls.clivus.foramen.ang.sinpmx$RegScore, Intercept = F, 
                                          predmin = min(plot.pgls.clivus.foramen.ang.sinpmx$RegScore), predmax = max(plot.pgls.clivus.foramen.ang.sinpmx$RegScore))

preds.foramen<-shape.predictor(pgls.foramen.ang.sinpmx$GM$pgls.fitted, x = plot.pgls.foramen.sinpmx$RegScore, Intercept = F, 
                               predmin = min(plot.pgls.foramen.sinpmx$RegScore), predmax = max(plot.pgls.foramen.sinpmx$RegScore))



#cargar meshes y configuraciones de referencia
mesh.original <- read.ply("Ptilonorhynchus_violaceus.ply", ShowSpecimen = F, addNormals = F ) #asegurarse de que no tiene el binary encoding marcado en MeshLab
coords.original <- coords.todos.sinpremax[,,"Ptilonorhynchus_violaceus"]
mesh.original2 <- read.ply("Mephitis_mephitis-fmnh_mammals_129327_M81787-157915_red.ply", ShowSpecimen = F, addNormals = F )
coords.original2 <- coords.todos.sinpremax[,,"Mephitis_mephitis"]
mesh.original3 <- read.ply("Physignathus_cocincinus.ply", ShowSpecimen = F, addNormals = F )
coords.original3 <- coords.todos.sinpremax[,,"Physignathus_cocincinus"]

M <- mshape(GPA.fit.sinpmx$coords)
M.allometry <- pgls.skullCS.sinpmx$GM$pgls.mean


#calculates per landmarks variance 
library(hot.dots)
library(tibble)

per_lm_variance <- function (shape.data) #para modificar los colores de los LMs, meter los codigos que gusten en cols1
{
  variances <- rowSums(apply(shape.data, c(1, 2), var))
  cols1 <- colorRampPalette(c("gray88", "yellow", "red"))
  cols <- cols1(100)
  x = (log10(variances))
  xlims <- NULL
  tol <- 1e-06
  xlims <- range(x) + c(-tol, tol)
  nbin = 100
  breaks <- 0:nbin/nbin * (xlims[2] - xlims[1]) + xlims[1]
  whichColor <- function(p, cols, breaks) {
    i <- 1
    while (p >= breaks[i] && p > breaks[i + 1]) i <- i +
        1
    cols[i]
  }
  variancecolors <- sapply(x, whichColor, cols = cols, breaks = breaks)
  variance_table <- tibble(Per_Lm_Variance = variances, Log_Variance = x,Variance_Colors = variancecolors)
  return(variance_table)
}

my.variances <- per_lm_variance(shape.data = GPA.fit$coords)
my.variances.sinpmx <- per_lm_variance(shape.data = GPA.fit.sinpmx$coords)

my.variances.allometry <- per_lm_variance(shape.data = pgls.skullCS.sinpmx$GM$pgls.fitted)
my.variances.headbrain.res <- per_lm_variance(shape.data = pgls.headbrain.sinpmx$GM$pgls.fitted)
my.variances.clivus.foramen.ang <- per_lm_variance(shape.data = pgls.clivus.foramen.ang.sinpmx$GM$pgls.fitted)
my.variances.foramen <- per_lm_variance(shape.data = pgls.foramen.ang.sinpmx$GM$pgls.fitted)

per_lm_distance <- function (shape.data) #para modificar los colores de los LMs, meter los codigos que gusten en cols1
{
  shape.data <- preds
  df <- as.data.frame(cbind(shape.data$pred1, shape.data$pred2))
  distances <- apply(df, 1, function(x) dist(matrix(x, nrow = 2, byrow = TRUE)))
  cols1 <- colorRampPalette(c("gray88", "yellow", "red"))
  cols <- cols1(100)
  x = (log10(distances))
  xlims <- NULL
  tol <- 1e-06
  xlims <- range(x) + c(-tol, tol)
  nbin = 100
  breaks <- 0:nbin/nbin * (xlims[2] - xlims[1]) + xlims[1]
  whichColor <- function(p, cols, breaks) {
    i <- 1
    while (p >= breaks[i] && p > breaks[i + 1]) i <- i +
        1
    cols[i]
  }
  distancecolors <- sapply(x, whichColor, cols = cols, breaks = breaks)
  distance_table <- tibble(Per_Lm_Distance = distances, Log_Distance = x,Distance_Colors = distancecolors)
  return(distance_table)
}
my.distances <- per_lm_distance(shape.data = preds)

# now we need to put all shapes in the same space to plot warped mesh and hot.dots on top of it
#PCA
open3d()
par3d(windowRect = c(0,0,500,500)) #para dimensiones del cuadro RGL
Sys.sleep(1)
mfrow3d(nr = 1, nc = 2, byrow = TRUE, sharedMouse = TRUE)

mesh.min <- tps3d(mesh.original, refmat = coords.original, tarmat = preds$pred1) #maps landmarks or a triangular mesh via thin plate spline based on a reference and a target configuration in 2D and 3D
shade3d(mesh.min, col= 8, alpha = 0.3) #dibuja la mesh con la forma predicha en la línea anterior #alpha cambia la transparencia de la surface
spheres3d(preds$pred1, col = my.distances$Distance_Colors, radius = 0.01) #radius cambia el tamaño del LM

next3d()
mesh.max <- tps3d(mesh.original, refmat = coords.original, tarmat = preds$pred2, lambda = 0.5)
shade3d(mesh.max, col= 8, alpha = 0.3) 
spheres3d(preds$pred2, col = my.distances$Distance_Colors, radius = 0.01)


mesh2ply(mesh.max, filename = "mesh.max", col = NULL, writeNormals = FALSE)



#ALLOMETRY
mesh.mean <- tps3d(mesh.original, refmat = coords.original, tarmat = M.allometry)
shade3d(mesh.mean, col= 8, alpha = 0.3)

open3d()
par3d(windowRect = c(0,0,500,500)) #para dimensiones del cuadro RGL
Sys.sleep(1)
mfrow3d(nr = 1, nc = 2, byrow = TRUE, sharedMouse = TRUE)

mesh.min <- tps3d(mesh.original2, refmat = coords.original2, tarmat = preds.allometry$predmin) #maps landmarks or a triangular mesh via thin plate spline based on a reference and a target configuration in 2D and 3D
shade3d(mesh.min, col= 8, alpha = 0.3) #dibuja la mesh con la forma predicha en la línea anterior #alpha cambia la transparencia de la surface
spheres3d(preds.allometry$predmin, col = my.variances.allometry$Variance_Colors, radius = 0.01) #radius cambia el tamaño del LM

next3d()
mesh.max <- tps3d(mesh.original2, refmat = coords.original2, tarmat = preds.allometry$predmax, lambda = 0.5)
shade3d(mesh.max, col= 8, alpha = 0.3) 
spheres3d(preds.allometry$predmax, col = my.variances.allometry$Variance_Colors, radius = 0.01)


#HEADBRAIN.RES
open3d()
par3d(windowRect = c(0,0,500,500)) #para dimensiones del cuadro RGL
Sys.sleep(1)
mfrow3d(nr = 1, nc = 2, byrow = TRUE, sharedMouse = TRUE)

mesh.min <- tps3d(mesh.original, refmat = coords.original, tarmat = preds.headbrain.res$predmin) #maps landmarks or a triangular mesh via thin plate spline based on a reference and a target configuration in 2D and 3D
shade3d(mesh.min, col= 8, alpha = 0.3) #dibuja la mesh con la forma predicha en la línea anterior #alpha cambia la transparencia de la surface
spheres3d(preds.headbrain.res$predmin, col = my.variances.headbrain.res$Variance_Colors, radius = 0.01) #radius cambia el tamaño del LM

next3d()
mesh.max <- tps3d(mesh.original, refmat = coords.original, tarmat = preds.headbrain.res$predmax, lambda = 0.5)
shade3d(mesh.max, col= 8, alpha = 0.3) 
spheres3d(preds.headbrain.res$predmax, col = my.variances.headbrain.res$Variance_Colors, radius = 0.01)


#CLIVUS.FORAMEN.ANG ###ESTE NO TERMINA DE EXPRESAR BIEN LOS CAMBIOS DE FORMA
open3d()
par3d(windowRect = c(0,0,500,500)) #para dimensiones del cuadro RGL
Sys.sleep(1)
mfrow3d(nr = 1, nc = 2, byrow = TRUE, sharedMouse = TRUE)

mesh.min <- tps3d(mesh.original, refmat = coords.original, tarmat = preds.clivus.foramen.ang$predmin) #maps landmarks or a triangular mesh via thin plate spline based on a reference and a target configuration in 2D and 3D
shade3d(mesh.min, col= 8, alpha = 0.3) #dibuja la mesh con la forma predicha en la línea anterior #alpha cambia la transparencia de la surface
spheres3d(preds.clivus.foramen.ang$predmin, col = my.variances.clivus.foramen.ang$Variance_Colors, radius = 0.01) #radius cambia el tamaño del LM

next3d()
mesh.max <- tps3d(mesh.original, refmat = coords.original, tarmat = preds.clivus.foramen.ang$predmax, lambda = 0.5)
shade3d(mesh.max, col= 8, alpha = 0.3) 
spheres3d(preds.clivus.foramen.ang$predmax, col = my.variances.clivus.foramen.ang$Variance_Colors, radius = 0.01)


#FORAMEN
open3d()
par3d(windowRect = c(0,0,500,500)) #para dimensiones del cuadro RGL
Sys.sleep(1)
mfrow3d(nr = 1, nc = 2, byrow = TRUE, sharedMouse = TRUE)

mesh.min <- tps3d(mesh.original2, refmat = coords.original2, tarmat = preds.foramen$predmin) #maps landmarks or a triangular mesh via thin plate spline based on a reference and a target configuration in 2D and 3D
shade3d(mesh.min, col= 8, alpha = 0.3) #dibuja la mesh con la forma predicha en la línea anterior #alpha cambia la transparencia de la surface
spheres3d(preds.foramen$predmin, col = my.variances.foramen$Variance_Colors, radius = 0.01) #radius cambia el tamaño del LM

next3d()
mesh.max <- tps3d(mesh.original2, refmat = coords.original2, tarmat = preds.foramen$predmax, lambda = 0.5)
shade3d(mesh.max, col= 8, alpha = 0.3) 
spheres3d(preds.foramen$predmax, col = my.variances.foramen$Variance_Colors, radius = 0.01)



##########################################
#### ANALISIS DE FLEXION CRANEOFACIAL ####
##########################################


division.LM.list.CB.sinpmx <- c("A","A","A","A","A","A","A","A","A","A","A","B","B","B","B","B","B","B","B","B","B","B","A","B","B","A","A","A","A","B","A","A","A","A","A","A","A","A","A","B","B","B","B","B","B","A","A","A","A")
integration.CB.sinpmx <- phylo.integration(gpa.todos.sinpremax$coords, partition.gp = division.LM.list.CB.sinpmx, phy = tree.amniota, iter = 999)
summary(integration.CB.sinpmx)
plot1 <- plot(integration.CB.sinpmx)


rawcoords.cranialbase.todos <- read.csv("rawcoords.cranialbase.todos.csv",header=T,sep=";",row.names=1)
coords.cranialbase.todos<-arrayspecs(rawcoords.cranialbase.todos[,1:60],20,3)
gpa.cranialbase<-gpagen(coords.cranialbase.todos)
proc.coords.cranialbase.2darray <- two.d.array(gpa.cranialbase$coords)

rawcoords.brain.todos <- read.csv("rawcoords.brain.todos.csv",header=T,sep=";",row.names=1)
coords.brain.todos<-arrayspecs(rawcoords.brain.todos[,1:60],20,3)
gpa.brain<-gpagen(coords.brain.todos)
proc.coords.brain.2darray <- two.d.array(gpa.brain$coords)

rawcoords.foramen.todos <- read.csv("rawcoords.foramen.todos.csv",header=T,sep=";",row.names=1)
coords.foramen.todos<-arrayspecs(rawcoords.foramen.todos[,1:12],4,3)
gpa.foramen<-gpagen(coords.foramen.todos)
pca.foramen<-gm.prcomp(gpa.cranialbase$coords, phy = tree.amniota)
summary(pca.foramen)
plot(pca.foramen, phylo = TRUE, main = "phyloPCA", phylo.par = list(node.labels = FALSE, anc.states = FALSE, node.cex = 0.2, edge.color = "grey", edge.width = 0.5, tip.txt.cex = 0.5))
proc.coords.foramen.2darray <- two.d.array(gpa.foramen$coords)     
variables.amniotas$foramen.CS <- gpa.foramen$Csize
variables.amniotas.temp[,21] 

foram.res <- as.data.frame(read.csv("foram.res.csv",header=F,sep=";"))
variables.amniotas$foram.res <- foram.res$V1
foram.brain.res <- as.data.frame(read.csv("foram.brain.res.csv",header=F,sep=";"))
variables.amniotas$foram.brain.res <- foram.brain.res$V1

df.variables.amniotas <- as.data.frame(variables.amniotas)  
angles <- as.matrix(df.variables.amniotas[,18:22])
names(angles) <- row.names(df.variables.amniotas)

angles.paleo <- angles.extantxfossils[,5:8]
names(angles.paleo) <- row.names(angles.extantxfossils)
coords.angles.paleo <- arrayspecs(angles.paleo[,1:4],1,4)


#### PCAs ###

phyloPCA.cranialbase<-gm.prcomp(gpa.cranialbase$coords, phy = tree.amniota)
summary(phyloPCA.cranialbase)
plot(phyloPCA.cranialbase, axis1 = 1, axis2 = 2)
plot(phyloPCA.cranialbase, phylo = TRUE, main = "phyloPCA", phylo.par = list(node.labels = FALSE, anc.states = FALSE, node.cex = 0.2, edge.color = "grey", edge.width = 0.5, tip.txt.cex = 0.5))

phyloPCA.brain<-gm.prcomp(gpa.brain$coords, phy = tree.amniota)
summary(phyloPCA.brain)
plot(phyloPCA.brain, axis1 = 1, axis2 = 2)
plot(phyloPCA.brain, phylo = TRUE, main = "phyloPCA", phylo.par = list(node.labels = FALSE, anc.states = FALSE, node.cex = 0.2, edge.color = "grey", edge.width = 0.5, tip.txt.cex = 0.5))



phyloPCA.angles <- gm.prcomp(angles, phy = tree.amniota)   
summary(phyloPCA.angles)
plot(phyloPCA.angles, axis1 = 1, axis2 = 2)
plot(phyloPCA.angles, phylo = TRUE, main = "phyloPCA", loadings = TRUE , phylo.par = list(node.labels = FALSE, anc.states = FALSE, node.cex = 0.2, edge.color = "grey", edge.width = 0.5, tip.txt.cex = 0.5))
arrows(0, 0, phyloPCA.angles$rotation[,1]*100, phyloPCA.angles$rotation[,2]*100, length = 0.2)

phyloPCA.angles.paleo <- gm.prcomp(angles.paleo, phy = tree.amniotafossils.full)


# Extract PC axes for plotting
PCAvalues <- data.frame(phyloPCA.angles$x)

# Extract loadings of the variables
PCAloadings <- data.frame(Variables = rownames(phyloPCA.angles$rotation), phyloPCA.angles$rotation)

df<-data.frame(clades.a = variables.amniotas$clades.a, clades.b = variables.amniotas$clades.b,
               clades.c = variables.amniotas$clades.c, clades.d = variables.amniotas$clades.comp, 
               log.skullCS = log(variables.amniotas$skull.CS.sinpmx), log.cranialbase.CS = log(variables.amniotas$cranialbase.CS),
               log.brainconvhull = log(variables.amniotas$brain.convhull), headbrain.res = variables.amniotas$headbrain.res, 
               headbrain.res.conpmx = variables.amniotas$headbrain.res.conpmx, encephalization = variables.amniotas$encephalization,
               skullbase.res = variables.amniotas$skullbase.res, foramen = variables.amniotas$foramen.angule, base.ang = variables.amniotas$base.angule, 
               flexion = variables.amniotas$flexion.angule, sph.res = variables.amniotas$SPH1.res, 
               clivus.foramen.ang = variables.amniotas$clivus.foramen.angle, facial.ang = variables.amniotas$facial.angle,
               PC1 = phyloPCA.angles$x[,1], PC2 = phyloPCA.angles$x[,2], PC3 = phyloPCA.angles$x[,3], PC4 = phyloPCA.angles$x[,4],
               RegScores.skullCS.sinpremax = RegScores.skullCS.sinpremax, 
               RegScores.headbrain.sinpremax = RegScores.headbrain.sinpremax,
               RegScores.skullbase.sinpmx = RegScores.skullbase.sinpmx, 
               row.names = row.names(variables.amniotas))

# Plot
Plot1 <- ggplot(df, aes(PC1, PC2, label = row.names(df),  colour = clades.b)) + geom_point(aes(size = log.cranialbase.CS))
Plot2 <- print(Plot1 
               + geom_text(color = "black", size = 0.5, check_overlap = F)) 
+ geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (Comp1*100),
                                       yend = (Comp2*100)), arrow = arrow(length = unit(1/2, "picas")), color = "black") 
+ annotate("text", x = (PCAloadings$Comp1*100), y = (PCAloadings$Comp2*100), label = PCAloadings$Variables)



# Extract loadings of the variables
PCAloadings <- data.frame(Variables = rownames(PCA$rotation), PCA$rotation)

# Plot
ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = Species)) +
  geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (PC1*5),
                                       yend = (PC2*5)), arrow = arrow(length = unit(1/2, "picas")),
               color = "black") +
  geom_point(size = 3) +
  annotate("text", x = (PCAloadings$PC1*5), y = (PCAloadings$PC2*5),
           label = PCAloadings$Variables)


arrows(0, 0, phyloPCA.angles$rotation[,1]*100, phyloPCA.angles$rotation[,2]*100, length = 0.2)


autoplot(phyloPCA.angles, loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)
library(ggfortify)
pca_angles <- prcomp(angles, scale. = TRUE)
autoplot(pca_angles, loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)  



#### MORPHOLOGICAL INTEGRATION OF THE CRANIAL BASE ft. ANGLES ###

integration.cranialbase.angles <- phylo.integration(A = proc.coords.cranialbase.2darray , A2 = angles, phy = tree.amniota, iter = 999)
summary(integration.cranialbase.angles)
plot(integration.cranialbase.angles)
ord <- ordinate(integration.cranialbase.angles$YScores, A = integration.cranialbase.angles$XScores)
round(ord$RV, 4)

preds.pls.base.angles <- shape.predictor(arrayspecs(integration.cranialbase.angles$A1[,1:60],20,3), integration.cranialbase.angles$XScores[,1], 
                                         method = "PLS", predmin = min(integration.cranialbase.angles$XScores[,1]),
                                         predmax = max(integration.cranialbase.angles$XScores[,1]),  pref  = 0) 
#modificar los valores del probs para exagerar o suavizar los cambios de forma
plotRefToTarget(preds.pls.base.angles$pref, preds.pls.base.angles$predmin)
plotRefToTarget(preds.pls.base.angles$pref, preds.pls.base.angles$predmax)
library(hot.dots)
library(tibble)
per_lm_distance <- function (shape.data) #para modificar los colores de los LMs, meter los codigos que gusten en cols1
{
  df <- as.data.frame(cbind(shape.data$predmin, shape.data$predmax))
  distances <- apply(df, 1, function(x) dist(matrix(x, nrow = 2, byrow = TRUE)))
  cols1 <- colorRampPalette(c("#1a2a6c", "#7B90E3", "#fdbb2d", "#b21f1f"))
  cols <- cols1(100)
  x = (log10(distances))
  xlims <- NULL
  tol <- 1e-06
  xlims <- range(x) + c(-tol, tol)
  nbin = 100
  breaks <- 0:nbin/nbin * (xlims[2] - xlims[1]) + xlims[1]
  whichColor <- function(p, cols, breaks) {
    i <- 1
    while (p >= breaks[i] && p > breaks[i + 1]) i <- i +
        1
    cols[i]
  }
  distancecolors <- sapply(x, whichColor, cols = cols, breaks = breaks)
  distance_table <- tibble(Per_Lm_Distance = distances, Log_Distance = x,Distance_Colors = distancecolors)
  return(distance_table)
}
#sale el pico metido para dentro en la pred minima
my.distances.preds.pls.base.ang <- per_lm_distance(shape.data = preds.pls.base.angles)
open3d()
par3d(windowRect = c(0,0,500,500)) #para dimensiones del cuadro RGL
Sys.sleep(1)
mfrow3d(nr = 1, nc = 2, byrow = TRUE, sharedMouse = TRUE)
spheres3d(preds.pls.base.angles$predmin, col = my.distances.preds.pls.base.ang$Distance_Colors, radius = 0.01) #radius cambia el tamaño del LM
next3d()
spheres3d(preds.pls.base.angles$predmax, col = my.distances.preds.pls.base.ang$Distance_Colors, radius = 0.01)

integration.skull.angles <- phylo.integration(A = proc.coords.2darray.sinpmx , A2 = angles, phy = tree.amniota, iter = 999)
summary(integration.skull.angles)
plot(integration.skull.angles)
ord <- ordinate(integration.skull.angles$YScores, A = integration.skull.angles$XScores)
round(ord$RV, 4)

preds.pls.skull.angles <- shape.predictor(arrayspecs(integration.skull.angles$A1[,1:147],49,3), integration.skull.angles$XScores[,1], 
                                          method = "PLS", predmin = min(integration.skull.angles$XScores[,1]),
                                          predmax = max(integration.skull.angles$XScores[,1]),  pref  = 0) 
#modificar los valores del probs para exagerar o suavizar los cambios de forma
plotRefToTarget(preds.pls.skull.angles$pref, preds.pls.skull.angles$predmin)
plotRefToTarget(preds.pls.skull.angles$pref, preds.pls.skull.angles$predmax)
library(hot.dots)
library(tibble)
per_lm_distance <- function (shape.data) #para modificar los colores de los LMs, meter los codigos que gusten en cols1
{
  df <- as.data.frame(cbind(shape.data$predmin, shape.data$predmax))
  distances <- apply(df, 1, function(x) dist(matrix(x, nrow = 2, byrow = TRUE)))
  cols1 <- colorRampPalette(c("#1a2a6c", "#7B90E3", "#fdbb2d", "#b21f1f"))
  cols <- cols1(100)
  x = (log10(distances))
  xlims <- NULL
  tol <- 1e-06
  xlims <- range(x) + c(-tol, tol)
  nbin = 100
  breaks <- 0:nbin/nbin * (xlims[2] - xlims[1]) + xlims[1]
  whichColor <- function(p, cols, breaks) {
    i <- 1
    while (p >= breaks[i] && p > breaks[i + 1]) i <- i +
        1
    cols[i]
  }
  distancecolors <- sapply(x, whichColor, cols = cols, breaks = breaks)
  distance_table <- tibble(Per_Lm_Distance = distances, Log_Distance = x,Distance_Colors = distancecolors)
  return(distance_table)
}

my.distances.preds.pls.skull.ang <- per_lm_distance(shape.data = preds.pls.skull.angles)
open3d()
par3d(windowRect = c(0,0,500,500)) #para dimensiones del cuadro RGL
Sys.sleep(1)
mfrow3d(nr = 1, nc = 2, byrow = TRUE, sharedMouse = TRUE)
mesh.min <- tps3d(mesh.original2, refmat = coords.original2, tarmat = preds.pls.skull.angles$predmin) #maps landmarks or a triangular mesh via thin plate spline based on a reference and a target configuration in 2D and 3D
shade3d(mesh.min, col= 8, alpha = 0.3) #dibuja la mesh con la forma predicha en la línea anterior #alpha cambia la transparencia de la surface
spheres3d(preds.pls.skull.angles$predmin, col = my.distances.preds.pls.skull.ang$Distance_Colors, radius = 0.01) #radius cambia el tamaño del LM
next3d()
mesh.max <- tps3d(mesh.original2, refmat = coords.original2, tarmat = preds.pls.skull.angles$predmax, lambda = 0.5)
shade3d(mesh.max, col= 8, alpha = 0.3) 
spheres3d(preds.pls.skull.angles$predmax, col = my.distances.preds.pls.skull.ang$Distance_Colors, radius = 0.01)


integration.skull.con.angles <- phylo.integration(A = proc.coords.2darray , A2 = angles, phy = tree.amniota, iter = 999)
summary(integration.skull.con.angles)
plot(integration.skull.con.angles)
ord <- ordinate(integration.skull.con.angles$YScores, A = integration.skull.con.angles$XScores)
round(ord$RV, 4)

preds.pls.skull.con.angles <- shape.predictor(arrayspecs(integration.skull.con.angles$A1[,1:150],50,3), integration.skull.con.angles$XScores[,1], 
                                              method = "PLS", predmin = min(integration.skull.con.angles$XScores[,1]),
                                              predmax = max(integration.skull.con.angles$XScores[,1]),  pref  = 0) 
#modificar los valores del probs para exagerar o suavizar los cambios de forma
plotRefToTarget(preds.pls.skull.con.angles$pref, preds.pls.skull.con.angles$predmin)
plotRefToTarget(preds.pls.skull.con.angles$pref, preds.pls.skull.con.angles$predmax)
library(hot.dots)
library(tibble)
per_lm_distance <- function (shape.data) #para modificar los colores de los LMs, meter los codigos que gusten en cols1
{
  df <- as.data.frame(cbind(shape.data$predmin, shape.data$predmax))
  distances <- apply(df, 1, function(x) dist(matrix(x, nrow = 2, byrow = TRUE)))
  cols1 <- colorRampPalette(c("#1a2a6c", "#7B90E3", "#fdbb2d", "#b21f1f"))
  cols <- cols1(100)
  x = (log10(distances))
  xlims <- NULL
  tol <- 1e-06
  xlims <- range(x) + c(-tol, tol)
  nbin = 100
  breaks <- 0:nbin/nbin * (xlims[2] - xlims[1]) + xlims[1]
  whichColor <- function(p, cols, breaks) {
    i <- 1
    while (p >= breaks[i] && p > breaks[i + 1]) i <- i +
        1
    cols[i]
  }
  distancecolors <- sapply(x, whichColor, cols = cols, breaks = breaks)
  distance_table <- tibble(Per_Lm_Distance = distances, Log_Distance = x,Distance_Colors = distancecolors)
  return(distance_table)
}
#sale el pico metido para dentro en la pred minima
my.distances.preds.pls.skull.con.ang <- per_lm_distance(shape.data = preds.pls.skull.con.angles)
open3d()
par3d(windowRect = c(0,0,500,500)) #para dimensiones del cuadro RGL
Sys.sleep(1)
mfrow3d(nr = 1, nc = 2, byrow = TRUE, sharedMouse = TRUE)
mesh.min <- tps3d(mesh.original, refmat = coords.original, tarmat = preds.pls.skull.con.angles$predmin) #maps landmarks or a triangular mesh via thin plate spline based on a reference and a target configuration in 2D and 3D
shade3d(mesh.min, col= 8, alpha = 0.3) #dibuja la mesh con la forma predicha en la línea anterior #alpha cambia la transparencia de la surface
spheres3d(preds.pls.skull.con.angles$predmin, col = my.distances.preds.pls.skull.con.ang$Distance_Colors, radius = 0.01) #radius cambia el tamaño del LM
next3d()
mesh.max <- tps3d(mesh.original, refmat = coords.original, tarmat = preds.pls.skull.con.angles$predmax, lambda = 0.5)
shade3d(mesh.max, col= 8, alpha = 0.3) 
spheres3d(preds.pls.skull.con.angles$predmax, col = my.distances.preds.pls.skull.con.ang$Distance_Colors, radius = 0.01)

integration.brain.angles <- phylo.integration(A = proc.coords.brain.2darray , A2 = angles, phy = tree.amniota, iter = 999)
summary(integration.brain.angles)
plot(integration.brain.angles)
ord <- ordinate(integration.brain.angles$YScores, A = integration.brain.angles$XScores)
round(ord$RV, 4)

preds.pls.brain.angles <- shape.predictor(arrayspecs(integration.brain.angles$A1[,1:60],20,3), integration.brain.angles$XScores[,1], 
                                          method = "PLS", predmin = min(integration.brain.angles$XScores[,1]),
                                          predmax = max(integration.brain.angles$XScores[,1]),  pref  = 0) 
#modificar los valores del probs para exagerar o suavizar los cambios de forma
plotRefToTarget(preds.pls.brain.angles$pref, preds.pls.brain.angles$predmin)
plotRefToTarget(preds.pls.brain.angles$pref, preds.pls.brain.angles$predmax)
library(hot.dots)
library(tibble)
per_lm_distance <- function (shape.data) #para modificar los colores de los LMs, meter los codigos que gusten en cols1
{
  df <- as.data.frame(cbind(shape.data$predmin, shape.data$predmax))
  distances <- apply(df, 1, function(x) dist(matrix(x, nrow = 2, byrow = TRUE)))
  cols1 <- colorRampPalette(c("#1a2a6c", "#7B90E3", "#fdbb2d", "#b21f1f"))
  cols <- cols1(100)
  x = (log10(distances))
  xlims <- NULL
  tol <- 1e-06
  xlims <- range(x) + c(-tol, tol)
  nbin = 100
  breaks <- 0:nbin/nbin * (xlims[2] - xlims[1]) + xlims[1]
  whichColor <- function(p, cols, breaks) {
    i <- 1
    while (p >= breaks[i] && p > breaks[i + 1]) i <- i +
        1
    cols[i]
  }
  distancecolors <- sapply(x, whichColor, cols = cols, breaks = breaks)
  distance_table <- tibble(Per_Lm_Distance = distances, Log_Distance = x,Distance_Colors = distancecolors)
  return(distance_table)
}
#sale el pico metido para dentro en la pred minima
my.distances.preds.pls.brain.ang <- per_lm_distance(shape.data = preds.pls.brain.angles)
open3d()
par3d(windowRect = c(0,0,500,500)) #para dimensiones del cuadro RGL
Sys.sleep(1)
mfrow3d(nr = 1, nc = 2, byrow = TRUE, sharedMouse = TRUE)
spheres3d(preds.pls.brain.angles$predmin, col = my.distances.preds.pls.brain.ang$Distance_Colors, radius = 0.01) #radius cambia el tamaño del LM
next3d()
spheres3d(preds.pls.brain.angles$predmax, col = my.distances.preds.pls.brain.ang$Distance_Colors, radius = 0.01)

df.integration.shape.angles <- data.frame(XScores.base = integration.cranialbase.angles$XScores[,1],
                                          YScores.base = integration.cranialbase.angles$YScores[,1],
                                          XScores.skull = integration.skull.con.angles$XScores[,1],
                                          YScores.skull = integration.skull.con.angles$YScores[,1],
                                          XScores.skullsin = integration.skull.angles$XScores[,1],
                                          YScores.skullsin = integration.skull.angles$YScores[,1],
                                          XScores.brain = integration.brain.angles$XScores[,1],
                                          YScores.brain = integration.brain.angles$YScores[,1],
                                          clades.a = variables.amniotas$clades.a, clades.b = variables.amniotas$clades.b,
                                          clades.c = variables.amniotas$clades.c, clades.comp = variables.amniotas$clades.comp,
                                          foramen.angle = variables.amniotas$foramen.angule, base.angle = variables.amniotas$base.angule,
                                          flexion.angle = variables.amniotas$flexion.angule, clivus.foramen.angle = variables.amniotas$clivus.foramen.angle,
                                          facial.angle = variables.amniotas$facial.angle,
                                          bipedalism = as.factor(variables.amniotas$biped), flight = as.factor(variables.amniotas$flight),
                                          log.skullCS = log(variables.amniotas$skull.CS.sinpmx),
                                          log.cranialbase.CS = log(variables.amniotas$cranialbase.CS),
                                          log.brainconvhull = log(variables.amniotas$brain.convhull),
                                          headbrain.res = variables.amniotas$headbrain.res, 
                                          headbrain.res.conpmx = variables.amniotas$headbrain.res.conpmx,
                                          encephalization = variables.amniotas$encephalization,
                                          skullbase.res = variables.amniotas$skullbase.res,
                                          row.names = row.names(variables.amniotas))

Plot1 <- ggplot(df.integration.shape.angles, aes(XScores.brain, YScores.brain, label = rownames(df.integration.shape.angles), colour = clades.b))+ geom_point(aes(size = log.cranialbase.CS)) + geom_smooth(method = glm) #+ geom_line(aes( x = reg.lines.x, y = reg.lines.y, color =  clades.b), size = 2, alpha = 1)
Plot2 <- print(Plot1 + geom_text(color = "black", size = 2, check_overlap = T) + theme(aspect.ratio = 1))

#PLS con el X.res marcado como gradiente de rosa (cerebros pequeños) a morado (cerebros grandes) y el tamaño de los puntos como logCS
Plot1 <- ggplot(df.integration.CB.angles, aes(XScores, YScores, label = row.names(df.integration.CB.angles),  colour = log.cranialbase.CS)) + geom_point(aes(size = log.cranialbase.CS))
Plot2 <- print(Plot1 
               + scale_colour_gradient2(low = "#F9D9E7", mid = "#E9B7CE", high = "#392D69") 
               + geom_text(color = "black", size = 2, check_overlap = T))

#PLS con ángulos marcados como gradiente de azul grisaceo (foramen caudal) a naranja (foramen ventral) y el tamaño de los puntos como logCS
Plot1 <- ggplot(df.integration.CB.angles, aes(XScores, YScores, label = row.names(df.integration.CB.angles),  colour = facial.angle)) + geom_point(aes(size = log.cranialbase.CS))
Plot2 <- print(Plot1 
               + scale_colour_gradientn(colours = c("#638B95", "#AC7E5A", "#F4711F")) #utilizar scale_colour_gradientn para variables con amplio rango
               + geom_text(color = "black", size = 2, check_overlap = T))



## FUNCION PARA CALCULAR VALORES DE INTEGRACIÓN POR NODOS

get.descendants <- function( phy , node ) { extract.clade( phy , node )$tip.label } #extract.clade te genera un arbol a partir de todos los tips de un feterminado nodo, podando todos los demás tips
subset.rows <- function( data , rownames ) { data[ rownames , ] } 

phy <- tree.amniota

modularity.per.node <- function( block1 , block2 , phy ) {
  nodes.temp <- unique( phy$edge[,1] ) #para generar un vector que contenga los valores de los nodos internos de la filogenia
  descendants.temp <- lapply( nodes.temp , FUN = get.descendants , phy = phy ) #esta línea devuelve una lista de la longitud de x que se llena de los resultados obtenidos a partir de aplicar la función FUN a cada uno de los elementos de x
  names( descendants.temp ) <- nodes.temp			
  node.data.block1.list <- lapply( descendants.temp , FUN = subset.rows , data = block1 )
  node.data.block2.list <- lapply( descendants.temp , FUN = subset.rows , data = block2 )
  modularity.temp<-list() 
  for (i in 1: length( node.data.block1.list ) ) {
    node.2darray.data.list <- node.data.block1.list[[i]]
    node.3darray.data.list <- arrayspecs(node.2darray.data.list[,1:60],20,3) #cambiar [,1:150] si la matriz de origen es la que contiene los clasificadores [,5:154]
    data.names <- dimnames(node.3darray.data.list) [[3]]
    node.angles.data.list <- node.data.block2.list[[i]]
    drop.taxa <- phy$tip.label[ ! phy$tip.label %in% data.names ]
    prunned.phy <- drop.tip( phy , drop.taxa )
    pls<-phylo.integration(A = node.3darray.data.list, A2 = node.angles.data.list, phy = prunned.phy, iter = 999)
    modularity.temp[i]<- pls$Z                         
    
  }
  
  names(modularity.temp)<-nodes.temp
  modularity.temp<-unlist(modularity.temp)
  modularity.temp
  #modularity.temp.scaled<-modularity.temp / modularity.temp[1] ## activar si queremos que el valor para cada nodo esté escalado y se dé de manera relativa (sobre 1)
  #modularity.temp.scaled
  
}

modularity.per.node.base.angles <- modularity.per.node( block1 = proc.coords.cranialbase.2darray , block2 = angles , phy = phy)
#modularity.per.node.temp.absval.conpmx <- modularity.per.node( data = twodarray.proccoords.conpmx , phy = phy)


## COMO PLOTEAR LOS VALORES DE MODULARIDAD POR NODO EN CADA NODO DE LA FILOGENIA:
# PLS z.scores

tree.temp <- phy
z.scores <- modularity.per.node.base.angles #cambiar al elemento que corresponda


# prepares z.scores to be plotted
z.scores.mod <- z.scores[!is.na(z.scores)]# removes NaNs product of clade-tip daughter comparisons
z.scores.mod <- z.scores.mod[!duplicated(z.scores.mod)] # two daughters per clade are obtained but we want only to keep one to plot in the immediate parent node


#prunes modularity values to match only the represented nodes
zscores.point.size <- modularity.per.node.base.angles [names (z.scores)]
log.zscores.point.size <- log(zscores.point.size)+2
log.zscores.point.size

disparities.deviation.plot <- disparities.deviation.plot [names (z.scores)]

# generates colours for angles 

edge.values <- z.scores
resta <- min(edge.values)
resta.positivo <- -1*resta
edge.values <- edge.values + resta.positivo
colours.rgb <- c("#ff0000","#8f0441","#1f0781") #rojo-azul ##because we want people to look at outliers from BMw, grey is a good colour for unwanted attention
colour.scale <- colorRamp( rev( colours.rgb ) , space = "rgb", bias = 0.6) #bias ensures grey corresponds always in good BM fitted values (adjust conveniently and check, more sophisticated in the final plots)
BG.node.table <- colour.scale( edge.values / max( edge.values ) )
BG.zscores<- c()
for( i in 1:nrow( BG.node.table ) ) {
  BG.zscores[i]<- rgb( BG.node.table[i,1], BG.node.table[i,2], BG.node.table[i,3] , alpha = 255 , max = 255 )
}	


# subset node.ages to data, do this before plotting z.scores 

tree.temp <- tree
variable = edge.values
tree <- tree.temp						

nodes.ages.temp <- tree.age(tree.temp) 
node.ages.temp.mod <- nodes.ages.temp[,1]
names(node.ages.temp.mod)<- as.character(nodes.ages.temp[,2])
node.ages.temp.mod<-node.ages.temp.mod[ names(variable)]						

#plot to check nodes

plot(z.scores, node.ages.temp.mod, xlim = c(-2, max(z.scores)), ylim = c(max(node.ages.temp.mod), 0) , type = "p", cex = log.zscores.point.size, pch = 21, col = "transparent", bg = BG.zscores, frame.plot = F)
text( y = node.ages.temp.mod , x =  z.scores, labels = names( node.ages.temp.mod ), pos = 4, cex =0.7, col = "black")	      


#tree with z.scores

plot( tree.temp , type = "fan", show.tip.label = F, no.margin = T, y.lim =  c(-5, 72.907 ))
nodelabels( node = as.numeric(names(z.scores)), frame = "n" , cex = log.zscores.point.size, pch = 21 , bg = BG.zscores )



#### DEVIATION FROM BM OF THE MORPHOLOGICAL INTEGRATION ft. ANGLES ###

## MODULARIDAD POR NODO A PARTIR DE LOS DATOS SIMULADOS ##        
#BM
library(mvMORPH)

gpa.coords.temp <- gpa.cranialbase$coords

fit1 <- mvgls( two.d.array(gpa.coords.temp) ~ 1 , tree = tree.temp , model = "BM" , method = "H&L" ) #primero ajustamos los datos al modelo evolutivo elegido, en este caso BM no direccional
sim.temp.new <- mvSIM(tree.temp, nsim = 20, model = "BM1", param = list(sigma = fit1$sigma$Pinv, theta = fit1$coefficients )) #con mvSIM simulamos nuestros datos ajustados a un modelo evolutivo determinado, en este caso BM no direccional
block1.sim <- array(unlist(sim.temp.new), c(381,60,20), dimnames = list(dimnames(gpa.coords.temp)[[3]],as.character(c(1:60)),as.character(c(1:20)))) #con esto transformamos la lista de las 500 matrices de 381x150 que hemos simulado, en una sola matriz de 3 dimensiones (especies x coordenadas x simulaciones)

fit2 <- mvgls( angles ~ 1 , tree = tree.temp , model = "BM" , method = "H&L" ) #primero ajustamos los datos al modelo evolutivo elegido, en este caso BM no direccional
sim2.temp.new <- mvSIM(tree.temp, nsim = 20, model = "BM1", param = list(sigma = fit1$sigma$Pinv, theta = fit1$coefficients )) #con mvSIM simulamos nuestros datos ajustados a un modelo evolutivo determinado, en este caso BM no direccional
block2.sim <- array(unlist(sim2.temp.new), c(381,5,20), dimnames = list(dimnames(gpa.coords.temp)[[3]],as.character(c(1:5)),as.character(c(1:20)))) #con esto transformamos la lista de las 500 matrices de 381x150 que hemos simulado, en una sola matriz de 3 dimensiones (especies x coordenadas x simulaciones)

blocks.sim.merged <- abind::abind( block1.sim , block2.sim , along = 2)
dim(blocks.sim.merged) #[1] 381  65  20
blocks.id.base.angles <- c("A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","B","B","B","B","B")    

modularity.per.node.for.simulations <- function( data , phy ) {
  nodes.temp <- unique( phy$edge[,1] ) #para generar un vector que contenga los valores de los nodos internos de la filogenia
  descendants.temp <- lapply( nodes.temp , FUN = get.descendants , phy = phy ) #esta línea devuelve una lista de la longitud de x que se llena de los resultados obtenidos a partir de aplicar la función FUN a cada uno de los elementos de x
  names( descendants.temp ) <- nodes.temp			
  node.data.list <- lapply( descendants.temp , FUN = subset.rows , data = data )
  modularity.temp<-list() 
  for (i in 1: length( node.data.list ) ) {
    node.2darray.data.list <- node.data.list[[i]]
    data.names <- dimnames(node.2darray.data.list) [[1]]
    drop.taxa <- phy$tip.label[ ! phy$tip.label %in% data.names ]
    prunned.phy <- drop.tip( phy , drop.taxa )
    pls<-phylo.integration(node.2darray.data.list, partition.gp = blocks.id.base.angles, phy = prunned.phy, iter = 499)
    modularity.temp[i]<- pls$Z                         
    
  }
  
  names(modularity.temp)<-nodes.temp
  modularity.temp<-unlist(modularity.temp)
  modularity.temp
  #modularity.temp.scaled<-modularity.temp / modularity.temp[1] ## activar si queremos que el valor para cada nodo esté escalado y se dé de manera relativa (sobre 1)
  #modularity.temp.scaled
  
}

sim.modularities.nodes.BASE.ANGLES <- apply( blocks.sim.merged , 3 , modularity.per.node.for.simulations , phy = tree.temp )

### Generates sizes and symbols for deviations of empirical modularity from BM simulations

modularities.BASE.ANGLES.deviation <- modularity.per.node.base.angles - apply (sim.modularities.nodes.BASE.ANGLES, 1, mean)
max(modularities.BASE.ANGLES.deviation)
mean(modularities.BASE.ANGLES.deviation)
min(modularities.BASE.ANGLES.deviation)  


edge.values <-  modularities.BASE.ANGLES.deviation - min( modularities.BASE.ANGLES.deviation )
zscores.point.size <- modularity.per.node.base.angles [names (z.scores)]
log.zscores.point.size <- log(zscores.point.size)+2
names(log.zscores.point.size) <- row.names(disp.modul.nodes)
log.zscores.point.size <- na.omit(log.zscores.point.size)
edge.values<-edge.values[ names(log.zscores.point.size)]
colours.rgb <- c("#a4306b", "#d5cfce",  "#29956e" ) ##because we want people to look at outliers from BMw, grey is a good colour for unwanted attention
colour.scale <- colorRamp( rev( colours.rgb ) , space = "rgb", bias = 0.7) #bias ensures grey corresponds always in good BM fitted values (adjust conveniently and check, more sophisticated in the final plots)
BG.node.table <- colour.scale( edge.values / max( edge.values ) )
BG.modul <- c()
for( i in 1:nrow( BG.node.table ) ) {
  BG.modul[i] <- rgb( BG.node.table[i,1], BG.node.table[i,2], BG.node.table[i,3] , alpha = 255 , max = 255 )
}	

#gets nodes only for the daughter-clade comparisons, i.e., only picks daughter clades with 2 or more descendent tips

tree.temp <- tree
variable = edge.values
tree <- tree.temp						

nodes.ages.temp <- tree.age(tree.temp) 
node.ages.temp.mod <- nodes.ages.temp[,1]
names(node.ages.temp.mod)<- as.character(nodes.ages.temp[,2])
node.ages.temp.mod<-node.ages.temp.mod[ names(log.zscores.point.size)]						

#plot to check nodes

plot(edge.values, node.ages.temp.mod, xlim = c(0, 6), ylim = c(325, 0) , type = "p", cex = log.zscores.point.size, pch = 21, bg = BG.modul, frame.plot = T)
text( y = node.ages.temp.mod , x =  edge.values, labels = names( edge.values ), pos = 4, cex =0.7, col = "black")	      


#tree with pls.angles

plot( tree.temp , type = "fan", show.tip.label = F, no.margin = T)
nodelabels( node = as.numeric(names(edge.values)), frame = "n" , cex = log.zscores.point.size , pch = 21 , bg = BG.modul ) 



#### PGLSs BASICRANIAL SHAPE / BRAIN SHAPE / ANGLES ~ ANGLES ###

variables.amniotas.temp <- as.matrix(df.variables.amniotas[,5:27])

gdf.pgls.angles<-geomorph.data.frame(phy = tree.amniota, 
                                     coords.skull = gpa.todos$coords,
                                     coords.skull.sinpmx = gpa.todos.sinpremax$coords,
                                     coords.base = gpa.cranialbase$coords,
                                     coords.brain = gpa.brain$coords,
                                     variables.amniotas = as.matrix(variables.amniotas[,5:25]),
                                     log.BM = log(variables.amniotas$body.mass),
                                     log.skullCS = log(variables.amniotas$skull.CS),
                                     log.skullCS.sinpmx = log(variables.amniotas$skull.CS.sinpmx),
                                     log.brainCS = log(variables.amniotas$brain.CS),
                                     log.cranialbase.CS = log(variables.amniotas$cranialbase.CS),
                                     log.brainconvhull = log(variables.amniotas$brain.convhull),
                                     log.foramenCS = log(variables.amniotas$foramen.CS),
                                     headbrain.res = variables.amniotas$headbrain.res,
                                     headbrain.res.conpmx = variables.amniotas$headbrain.res.conpmx,
                                     encephalization = variables.amniotas$encephalization,
                                     skullbase.res = variables.amniotas$skullbase.res,
                                     foram.res = variables.amniotas$foram.res,
                                     foram.brain.res = variables.amniotas$foram.brain.res,
                                     SPH1 = variables.amniotas$SPH1.res,
                                     SPH2 = variables.amniotas$SPH2.res,
                                     SPH3 = variables.amniotas$SPH3.res,
                                     foramen = variables.amniotas$foramen.angule,
                                     base.ang = variables.amniotas$base.angule,
                                     flexion = variables.amniotas$flexion.angule,
                                     clivus.foramen.ang = variables.amniotas$clivus.foramen.angle,
                                     facial.ang = variables.amniotas$facial.angle,
                                     clades.a = variables.amniotas$clades.a,
                                     clades.b = variables.amniotas$clades.b,
                                     clades.c = variables.amniotas$clades.c,
                                     clades.comp = variables.amniotas$clades.comp,
                                     bipedal = variables.amniotas$biped,
                                     flight = variables.amniotas$flight,
                                     row.names = row.names(variables.amniotas))

## UNIVARIATES

angles.univar.models <- procD.pgls(log(variables.amniotas[,17]) ~ foram.brain.res, phy = phy, data = gdf.pgls.angles, iter = 9999, SS.type = "II")
summary(angles.univar.models)
plot.angles.univar.models <- plot(angles.univar.models, type = "regression", reg.type = "RegScore", predictor = gdf.pgls.angles$log.foramenCS)
angles.univar.models.cladeindep <- procD.pgls(log(variables.amniotas[,6]) ~ log.foramenCS + clades.b , phy = phy, data = gdf.pgls.angles, iter = 9999, SS.type = "II")
summary(angles.univar.models.cladeindep)
plot.pgls.angles.univar.clades <- plot(angles.univar.models.cladeindep, type = "regression", reg.type = "RegScore", predictor = gdf.pgls.angles$log.foramenCS)
RegScores.angles.univar.clades<-plot.pgls.angles.univar.clades$RegScore
pred.pgls.angles.univar.clades<-plot(angles.univar.models.cladeindep, type="regression", reg.type="PredLine", predictor = gdf.pgls.angles$log.foramenCS)
PredLine.angles.univar <- pred.pgls.angles.univar.clades$PredLine  
reg.lines.x <- as.vector(pred.pgls.angles.univar.clades$plot.args$x)
names(reg.lines.x) <- gdf.pgls.angles$row.names
reg.lines.y <- pred.pgls.angles.univar.clades$plot.args$y
names(reg.lines.y) <- gdf.pgls.angles$row.names
df.angles.univar.plots<-data.frame(RegScores.angles.univar.clades = RegScores.angles.univar.clades, 
                                   log.brainconvhull = gdf.pgls.angles$log.brainconvhull,
                                   headbrain.res.conpmx = gdf.pgls.angles$headbrain.res.conpmx,
                                   log.foramenCS = gdf.pgls.angles$log.foramenCS,
                                   clades.b = gdf.pgls.angles$clades.b,
                                   reg.lines.x = reg.lines.x, reg.lines.y = reg.lines.y)
rownames(df.angles.univar.plots) <- row.names(variables.amniotas)
Plot1 <- ggplot(df.angles.univar.plots, aes(log.foramenCS, RegScores.angles.univar.clades, label = rownames(df.angles.univar.plots), colour = clades.b))+ geom_point(aes(size = log.foramenCS)) + geom_smooth(method = glm)#+ geom_line(aes( x = reg.lines.x, y = reg.lines.y, color =  clades.b), size = 2, alpha = 1)
Plot2 <- print(Plot1+ geom_text(color = "black", size = 1, check_overlap = T))


## SHAPES ~ VARIABLE 1

angles.shape.models <- procD.pgls(coords.base ~ SPH1, phy = phy, data = gdf.pgls.angles, iter = 9999, SS.type = "II")
summary(angles.shape.models)
plot.angles.shape.models <- plot(angles.shape.models, type = "regression", reg.type = "RegScore", predictor = gdf.pgls.angles$SPH1)
angles.shape.models.cladeindep <- procD.pgls(coords.base ~ SPH1 + clades.b , phy = phy, data = gdf.pgls.angles, iter = 9999, SS.type = "II")
summary(angles.shape.models.cladeindep)
plot.pgls.angles.shape.clades <- plot(angles.shape.models.cladeindep, type = "regression", reg.type = "RegScore", predictor = gdf.pgls.angles$SPH1)
RegScores.angles.shape.clades<-plot.pgls.angles.shape.clades$RegScore
pred.pgls.angles.shape.clades<-plot(angles.shape.models.cladeindep, type="regression", reg.type="PredLine", predictor = gdf.pgls.angles$SPH1)
PredLine.angles.shape <- pred.pgls.angles.shape.clades$PredLine  
reg.lines.x <- as.vector(pred.pgls.angles.shape.clades$plot.args$x)
names(reg.lines.x) <- gdf.pgls.angles$row.names
reg.lines.y <- pred.pgls.angles.shape.clades$plot.args$y
names(reg.lines.y) <- gdf.pgls.angles$row.names
df.angles.shape.plots<-data.frame(RegScores.angles.shape.clades = RegScores.angles.shape.clades, 
                                  SPH1 = gdf.pgls.angles$SPH1,
                                  cranialbase.CS = gdf.pgls.angles$log.cranialbase.CS,
                                  headbrain.res.conpmx = gdf.pgls.angles$headbrain.res.conpmx,
                                  clades.b = gdf.pgls.angles$clades.b,
                                  reg.lines.x = reg.lines.x, reg.lines.y = reg.lines.y)
rownames(df.angles.shape.plots) <- row.names(variables.amniotas)
Plot1 <- ggplot(df.angles.shape.plots, aes(SPH1, RegScores.angles.shape.clades, label = rownames(df.angles.shape.plots), colour = clades.b))+ geom_point(aes(size = headbrain.res.conpmx+3)) + geom_smooth(method = glm)#+ geom_line(aes( x = reg.lines.x, y = reg.lines.y, color =  clades.b), size = 2, alpha = 1)
Plot2 <- print(Plot1+ geom_text(color = "black", size = 1, check_overlap = T))


## SHAPES ~ VARIABLE 1 + VARIABLE 2

angles.shape.models <- procD.pgls(coords.brain ~ log.brainconvhull : log.cranialbase.CS : foramen, phy = phy, data = gdf.pgls.angles, iter = 9999, SS.type = "II")
summary(angles.shape.models)
plot.angles.shape.models <- plot(angles.shape.models, type = "regression", reg.type = "RegScore", predictor = gdf.pgls.angles$foramen)
angles.shape.models.cladeindep <- procD.pgls(coords.brain ~ log.brainconvhull : log.cranialbase.CS + clades.b , phy = phy, data = gdf.pgls.angles, iter = 9999, SS.type = "II")
summary(angles.shape.models.cladeindep)
plot.pgls.angles.shape.clades <- plot(angles.shape.models.cladeindep, type = "regression", reg.type = "RegScore", predictor = gdf.pgls.angles$foramen)
RegScores.angles.shape.clades<-plot.pgls.angles.shape.clades$RegScore
pred.pgls.angles.shape.clades<-plot(angles.shape.models.cladeindep, type="regression", reg.type="PredLine", predictor = gdf.pgls.angles$foramen)
PredLine.angles.shape <- pred.pgls.angles.shape.clades$PredLine  
reg.lines.x <- as.vector(pred.pgls.angles.shape.clades$plot.args$x)
names(reg.lines.x) <- gdf.pgls.angles$row.names
reg.lines.y <- pred.pgls.angles.shape.clades$plot.args$y
names(reg.lines.y) <- gdf.pgls.angles$row.names
df.angles.shape.plots<-data.frame(RegScores.angles.shape.clades = RegScores.angles.shape.clades, 
                                  log.brainconvhull = gdf.pgls.angles$log.brainconvhull,
                                  cranialbase.CS = gdf.pgls.angles$log.cranialbase.CS,
                                  foramen = gdf.pgls.angles$foramen,
                                  headbrain.res.conpmx = gdf.pgls.angles$headbrain.res.conpmx,
                                  clades.b = gdf.pgls.angles$clades.b,
                                  reg.lines.x = reg.lines.x, reg.lines.y = reg.lines.y)
rownames(df.angles.shape.plots) <- row.names(variables.amniotas)
Plot1 <- ggplot(df.angles.shape.plots, aes(foramen, RegScores.angles.shape.clades, label = rownames(df.angles.shape.plots), colour = clades.b))+ geom_point(aes(size = headbrain.res.conpmx-min(headbrain.res.conpmx))) + geom_smooth(method = glm)#+ geom_line(aes( x = reg.lines.x, y = reg.lines.y, color =  clades.b), size = 2, alpha = 1)
Plot2 <- print(Plot1+ geom_text(color = "black", size = 1, check_overlap = T))


#### DTT FOR ANGLES ###

DTT.angles.all <- dtt(phy = tree.amniota, data = angles, index = c("avg.sq"), mdi.range = c(325,0), nsim = 1000, CI = 0.95, plot = TRUE)      
#raro

#por grupos
data.names <- row.names(grouping1.conpmx[[4]])
drop.taxa <- tree.amniota$tip.label[ ! tree.amniota$tip.label %in% data.names ]
pruntree.mamos <- drop.tip( tree.amniota , drop.taxa )

coords.birds <- as.matrix(grouping1.conpmx[[1]][,5:154])
DTT.birds <- dtt(phy = pruntree.birds, data = coords.birds, index = c("avg.sq"), mdi.range = c(0,1), nsim = 1000, CI = 0.95, plot = TRUE)      

coords.crocs <- as.matrix(grouping1.conpmx[[2]][,5:154])
DTT.crocs <- dtt(phy = pruntree.crocs, data = coords.crocs, index = c("avg.sq"), mdi.range = c(0,1), nsim = 1000, CI = 0.95, plot = TRUE)      

coords.lepids <- as.matrix(grouping1.conpmx[[3]][,5:154])
DTT.lepids <- dtt(phy = pruntree.lepids, data = coords.lepids, index = c("avg.sq"), mdi.range = c(0,1), nsim = 1000, CI = 0.95, plot = TRUE)      

coords.mamos <- as.matrix(grouping1.conpmx[[4]][,5:154])
DTT.mamos <- dtt(phy = pruntree.mamos, data = coords.mamos, index = c("avg.sq"), mdi.range = c(0,1), nsim = 1000, CI = 0.95, plot = TRUE)      


#### ANCESTRAL STATE ESTIMATION FOR ANGLES AND COMPARISON TO MODULARITY DEVIATIONS ###

x<-variables.amniotas$clivus.foramen.angle
x<-variables.amniotas$headbrain.res.conpmx
x<-variables.amniotas$encephalization
names(x)<-row.names(variables.amniotas)
ace.clivus <- ace(x , tree.amniota , type = "continuous" , method = "pic" , CI = TRUE )
ace.headbraincon <- ace(x , tree.amniota , type = "continuous" , method = "pic" , CI = TRUE )
ace.encephalization <- ace(x , tree.amniota , type = "continuous" , method = "pic" , CI = TRUE )

edge.values <-  modularities.deviation - min( modularities.deviation )
nodes.ages.temp <- tree.age(tree.temp) 
node.ages.temp.mod <- nodes.ages.temp[,1]
names(node.ages.temp.mod)<- as.character(nodes.ages.temp[,2])
node.ages.temp.mod<-node.ages.temp.mod[ names(log.zscores.point.size)]
log.nodeage.point.size <- log(node.ages.temp.mod)
names(log.zscores.point.size) <- row.names(disp.modul.nodes.red)
log.zscores.point.size <- na.omit(log.zscores.point.size)
edge.values<-edge.values[ names(log.zscores.point.size)]


colours.rgb <- c("#C8963E", "#d5cfce",  "#55185D" ) ##because we want people to look at outliers from BMw, grey is a good colour for unwanted attention
colour.scale <- colorRamp( rev( colours.rgb ) , space = "rgb", bias = 1.1) #bias ensures grey corresponds always in good BM fitted values (adjust conveniently and check, more sophisticated in the final plots)
BG.node.table <- colour.scale( edge.values / max( edge.values ) )
BG.modul <- c()
for( i in 1:nrow( BG.node.table ) ) {
  BG.modul[i] <- rgb( BG.node.table[i,1], BG.node.table[i,2], BG.node.table[i,3] , alpha = 255 , max = 255 )
}	

tree.temp <- tree
variable = edge.values
tree <- tree.temp						

nodes.ages.temp <- tree.age(tree.temp) 
node.ages.temp.mod <- nodes.ages.temp[,1]
names(node.ages.temp.mod)<- as.character(nodes.ages.temp[,2])
node.ages.temp.mod<-node.ages.temp.mod[ names(log.zscores.point.size)]						

ace.clivus.nodesvalues <- ace.clivus$ace[ names(log.zscores.point.size)]
ace.headbraincon.nodesvalues <- ace.headbraincon$ace[ names(log.zscores.point.size)]
ace.encephalization.nodesvalues <- ace.encephalization$ace[ names(log.zscores.point.size)]

#plot to check nodes
plot(edge.values, ace.clivus.nodesvalues, xlim = c(0, 7), ylim = c(60, 140) , cex = log.zscores.point.size+1, pch = 21, bg = BG.modul, frame.plot = T)
text( y = ace.clivus.nodesvalues , x =  edge.values, labels = names( edge.values ), pos = 4, cex =0.7, col = "black")	      

library(plot3D)
library(plotly)
library(rgl)

scatter3D(edge.values , ace.clivus.nodesvalues , node.ages.temp.mod ,  pch = 20 , cex = log.zscores.point.size+1 , ticktype = "detailed" , phi = 0 , theta = 0)
x = edge.values/max(edge.values)
y = (ace.headbraincon.nodesvalues-min(ace.headbraincon.nodesvalues))/max(ace.headbraincon.nodesvalues-min(ace.headbraincon.nodesvalues))
z = (node.ages.temp.mod-min(node.ages.temp.mod))/max(node.ages.temp.mod-min(node.ages.temp.mod))
rgl.open()
rgl.bg(color = c("white"))
rgl.spheres(-x,-z,y,color=BG.modul,r=y/35)
rgl.bbox(color = "#333377")
#

#### ANÁLISIS CON FÓSILES ###

# tree construction binding molecular phylogenies and fossil-based calibrated trees of extinct species  
library(paleotree)
tree.fossils.base<-read.nexus("TreeFossils_base.nex")
plot(tree.fossils.base)
is.ultrametric(tree.fossils.base)

#tree.fossils.base<-read.tree("TreeFossils_base.tre")   
#row.names(intervalTimes.fossils)<-intervalTimes.fossils[,1]
#setdiff(row.names(intervalTimes.fossils),tree.fossils.base$tip.label)
#intervalTimes.fossils[,1]<-as.character(intervalTimes.fossils[,1])
#timeData<-intervalTimes.fossils=intervalTimes.fossils[,2:3]  

intervalTimes.fossils <- read.table('IntervalTimes.txt',header=TRUE,row.names=1)
taxonTimes.fossils <- read.table('TaxonInterval_TreeFossils.txt',header=TRUE,row.names=1)

timeList_fossils <- list(intervalTimes.fossils , taxonTimes.fossils)

#tree.fossils.calib <- timePaleoPhy(tree.fossils.base, intervalTimes.fossils, type = "basic", vartime = NULL, ntrees = 1, randres=FALSE, timeres=FALSE, randres = TRUE, add.term=FALSE, inc.term.adj=FALSE, dateTreatment="firstLast", node.mins=NULL, noisyDrop=TRUE, plot=FALSE)

bintrees.fossils <- bin_timePaleoPhy(tree.fossils.base, timeList = timeList_fossils, type = "mbl", vartime = 1, ntrees = 100, nonstoch.bin = FALSE, randres = TRUE, timeres = F,sites = NULL, point.occur = FALSE, add.term = T, inc.term.adj = T, dateTreatment = "firstLast", node.mins = NULL, noisyDrop = TRUE, plot = TRUE)
writeNexus(bintrees.fossils, 'amniota_fossils_100_mlb.nex')

tree.fossils.base.calib <- read.nexus("amniota_fossils_consensus.nex")
plot( tree.fossils.base.calib , show.tip.label = T, no.margin = T, y.lim =  c(-5, 72.907 ))

##Para combinar filogenias seleccionando la rama y la posición temporal de los nodos, creando una filo compuesta calibrada

plot(tree_mesozoicbirds, show.node.label = TRUE, font = 1, root.edge = TRUE)
max(vcv.phylo(tree.aves)) #esto es para conocer la edad máxima del árbol
max(vcv.phylo(tree.mamos))
max(vcv.phylo(tree.reptiles))
cocos.consensus<-read.nexus('CrocsTree_consensus.nex')
max(vcv.phylo(cocos.consensus))

plot( tree.fossils.base.calib , y.lim =  c(-5, 72.907 ))
tiplabels()
axisPhylo()
axis(3)
install.packages("remotes")
remotes::install_github("fmichonneau/phyloch")
library("phyloch")
data(strat2012)
yt <- 47.8
plot.phylo(tree.fossils.base.calib, edge.width = 2, label.offset = 0.5)
axisGeo(strat2012, tip.time = yt, unit = c("stage", "epoch"))
tree.amniota.fossils.filled <- bind.tree(tree.fossils.base.calib,tree.aves, interactive = TRUE) #con esto localizamos la rama donde tenemos que incrustar el arbol
tree.amniota.fossils.filled <- bind.tree(tree.fossils.base.calib,tree.aves, where = 2, position = 108.5918) #ahora metemos el arbol diciendo directamente la rama (2) y la edad del nodo (la hemos sacado antes con el max(vcv.phylo(tree.aves)))
tree.amniota.fossils.filled2 <- bind.tree(tree.amniota.fossils.filled,cocos.consensus, where = 26, position = 96.17918)
tree.amniota.fossils.filled3 <- bind.tree(tree.amniota.fossils.filled2,tree.reptiles, where = 20, position = 243.2984)
tree.amniota.fossils.filled4 <- bind.tree(tree.amniota.fossils.filled3,tree.mamos, where = 23, position = 188.3652)

plot(tree.amniota.fossils.filled4, show.node.label = TRUE, font = 0.1, root.edge = TRUE)

tree.amniotafossils.full <- drop.tip(tree.amniota.fossils.filled4,c("Aves", "Mammalia", "Neosuchia", "Lepidosauria"))  #ahora cortamos las ramas que sobran, las de los grandes clados que hemos utilizado para anclar las filos calibradas
tree.amniotafossils.full<-ladderize(tree.amniotafossils.full)
plot(tree.amniotafossils.full, type = "fan", show.node.label = TRUE, font = 0.1, root.edge = TRUE)
writeNexus(tree.amniotafossils.full, 'TreeAmnioteFossils.full.nex')

#ancestral reconstruction of continuous variables (angles)

angles.extantxfossils <- read.csv("angles.extantxfossils.csv",header=T,sep=";",row.names=1)
angles.extantxfossils <- as.matrix(angles.extantxfossils)
names(angles.extantxfossils) <- row.names(angles.extantxfossils)

angles.fossilsprunned <- read.csv("angles.fossils.prunned2.csv",header=T,sep=";",row.names=1)
angles.fossilsprunned <- as.matrix(angles.fossilsprunned)
names(angles.fossilsprunned) <- row.names(angles.fossilsprunned)
tree.fossils.base.calib <- read.nexus("amniota_fossils2_consensus.nex")
tree.fossils.base.calib.pru <- drop.tip(tree.fossils.base.calib, c("Thecodontosaurus","Gnathovorax","Velociraptor","Shuvuuia","Enant_Brasil","Ixalerpeton"))
plot(tree.fossils.base.calib.pru, type = "fan", show.node.label = TRUE, cex = 1, root.edge = TRUE)
plotTree.wBars(tree.fossils.base.calib.pru,angles.fossilsprunned[,10], scale = NULL, tip.labels = TRUE)

x<-angles.fossilsprunned[,1] #foramen
names(x)<-row.names(angles.fossilsprunned)
ace.foramen.withfossils <- ace(x , tree.fossils.base.calib , type = "continuous" , method = "pic" , CI = TRUE )
phenogram(tree.fossils.base.calib,x,ftype="off",spread.labels=FALSE)

x<-angles.fossilsprunned[,2] #base
ace.base.withfossils <- ace(x , tree.amniotafossils.full , type = "continuous" , method = "pic" , CI = TRUE )
phenogram(tree.amniotafossils.full,x,ftype="off",spread.labels=FALSE)
x<-angles.fossilsprunned[,3] #flexion
ace.flexion.withfossils <- ace(x , tree.amniotafossils.full , type = "continuous" , method = "pic" , CI = TRUE )
phenogram(tree.amniotafossils.full,x,ftype="off",spread.labels=FALSE)
x<-angles.fossilsprunned[,4] #clivus
ace.clivus.withfossils <- ace(x , tree.amniotafossils.full , type = "continuous" , method = "pic" , CI = TRUE )
phenogram(tree.amniotafossils.full,x,ftype="off",spread.labels=FALSE)
x<-angles.fossilsprunned[,5] #face
names(x)<-row.names(angles.extantxfossils)
x <- na.omit(x)
tree.amniotafossils.pru <- drop.tip(tree.amniotafossils.full, c("Thecodontosaurus","Gnathovorax","Velociraptor","Shuvuuia","Enant_Brasil","Ixalerpeton"))
plot(tree.amniotafossils.pru, show.node.label = TRUE, cex = 0.1, root.edge = TRUE)
nodelabels(text = tree.amniotafossils.pru$node.label,
           frame = "n", cex=0.1, col= "blue")
ace.face.withfossils <- ace(x , tree.amniotafossils.pru , type = "continuous" , method = "pic" , CI = TRUE )
phenogram(tree.amniotafossils.pru,x,ftype="off",spread.labels=FALSE)





##############################
#### ANALISIS CON FOSILES ####
##############################

LMs.fossils.new<-read.csv("RawCoords93_fossils.csv",header=T,sep=";",row.names=1)  
coords<-arrayspecs(LMs.fossils.new[,1:93],31,3)
plot3d(coords[,,1]); aspect3d("iso")

## LOOP PARA REFLEJAR LANDMARKS ##
mirrored.coords <- list() #haces una lista vacía donde vas a meter los cráneos mirroreados, puede ser una lista un array una matriz lo que sea.
for (i in 1: length(coords)) { #para probar cada línea, empezar con i = 1 e ir corriendo cada línea. Si funcionan todas, quitar la i = 1 y correr la función
  craneo.individual <- coords[,,i] #coge el elemento del número que vaya tocando en cada iteración, 1 para la primera, 2 para la segunda hasta la 250 vez.
  craneo.mirrored <- mirror2plane(craneo.individual, v1 = craneo.individual[12,], v2 = craneo.individual[15,], v3 = craneo.individual[30,])
  mirrored.coords[[i]]<-craneo.mirrored #metes el craneo mirroreado en una lista nueva
}
plot3d(mirrored.coords[[1]]); aspect3d("iso") #visualización 3D sin deformación de ejes

## CAMINO ALTERNATIVO PARA COMBINAR EL ARRAY DEL LADO LANDMARKEADO Y EL ARRAY DEL LADO REFLEJADO ##
array.mirrored.coords <- simplify2array(mirrored.coords) #transformamos la lista de cráneos reflejados creada en el loop en un 3D array
pruned.array.mirrored.coords <- array.mirrored.coords[-c(1,2,12:17,25:26,30:31),,] #borramos las filas de los landmarks repetidos (los mediales)
concatenated.arrays.with.mirrored.coords <- bindArr(coords, pruned.array.mirrored.coords, along = 1) #concatenamos los dos 3D arrays
LMs.fossils2 <- two.d.array(concatenated.arrays.with.mirrored.coords) #transformamos el 3D array combinado en una matriz 2D como las que tenemos en excel con las coordenadas de los landmarks para guardarlo así
write.csv(LMs.fossils2, 'Rawcoords150_fossils_new.csv')
plot3d(concatenated.arrays.with.mirrored.coords[,,1]); aspect3d("iso") #ploteamos para comprobar que han salido bien

## LOOP PARA CALCULAR VOLUMENES A PARTIR DE CONFIGURACIONES GEOMETRICAS ##
#BRAIN
Rawcoords_brain.paleo<-read.csv("Rawcoordsbrain_confosiles409.csv",header=T,sep=";",row.names=1)
coords.brain.paleo<-arrayspecs(Rawcoords_brain.paleo[,1:60],20,3)
gpa.brain.paleo <- gpagen(coords.brain.paleo)

brain.alphashape <- list() #haces una lista vacía donde vas a meter las alpha shape de cada cerebro, puede ser una lista un array una matriz lo que sea.
for (i in 1: length(coords.brain.paleo)) { #para probar cada línea, empezar con i = 1 e ir corriendo cada línea. Si funcionan todas, quitar la i = 1 y correr la función
  coords.individual <- coords.brain.paleo[,,i] #coge el elemento del número que vaya tocando en cada iteración, 1 para la primera, 2 para la segunda hasta la 250 vez.
  alphashape.individual <- ashape3d(coords.individual, alpha = 100, pert = FALSE, eps = 1e-09) #función simple
  vol.ashape.individual <- volume_ashape3d(alphashape.individual, byComponents = FALSE, indexAlpha = 1)
  brain.alphashape[[i]]<-vol.ashape.individual #metes la alpha-shape de cada cerebro en una lista nueva
}
plot(alphashape.individual, indexAlpha = "all")
write.csv(brain.alphashape, 'brain.convhull.paleo.csv')


## tree construction binding molecular phylogenies and fossil-based calibrated trees of extinct species  
library(paleotree)
setwd("C:/Users/Sergio Mtnez Nebreda/OneDrive - UAM/_TESIS/CAPÍTULOS/2.VENTRALIZATION")

tree.fossils.base<-read.nexus("TreeFossils_base2.nex")
plot(tree.fossils.base)
is.ultrametric(tree.fossils.base)

intervalTimes.fossils <- read.table('IntervalTimes.txt',header=TRUE,row.names=1)
taxonTimes.fossils <- read.table('TaxonInterval_TreeFossils.txt',header=TRUE,row.names=1)
timeList_fossils <- list(intervalTimes.fossils , taxonTimes.fossils)
setwd("C:/Users/Sergio Mtnez Nebreda/OneDrive - UAM/_TESIS/CAPÍTULOS/1.MACROEVOLUTION")

bintrees.fossils <- bin_timePaleoPhy(tree.fossils.base, timeList = timeList_fossils, type = "mbl", vartime = 1, ntrees = 100, nonstoch.bin = FALSE, randres = TRUE, timeres = F,sites = NULL, point.occur = FALSE, add.term = T, inc.term.adj = T, dateTreatment = "firstLast", node.mins = NULL, noisyDrop = TRUE, plot = TRUE)
writeNexus(bintrees.fossils, 'amniota_fossils_100_mlb.nex')

tree.fossils.base.calib <- read.nexus("amniota_fossils2_consensus.nex")
plot( tree.fossils.base.calib , show.tip.label = T, no.margin = T, y.lim =  c(-5, 72.907 ))

##Para combinar filogenias seleccionando la rama y la posición temporal de los nodos, creando una filo compuesta calibrada
plot( tree.fossils.base.calib , y.lim =  c(-5, 72.907 ))
tiplabels()
plot.phylo(tree.fossils.base.calib, edge.width = 2, label.offset = 0.5)
tree.amniota.fossils.filled <- bind.tree(tree.fossils.base.calib,tree.aves, where = 2, position = 108.5918) #ahora metemos el arbol diciendo directamente la rama (2) y la edad del nodo (la hemos sacado antes con el max(vcv.phylo(tree.aves)))
tree.amniota.fossils.filled2 <- bind.tree(tree.amniota.fossils.filled,cocos.consensus, where = 29, position = 96.17918)
tree.amniota.fossils.filled3 <- bind.tree(tree.amniota.fossils.filled2,tree.reptiles, where = 22, position = 243.2984)
tree.amniota.fossils.filled4 <- bind.tree(tree.amniota.fossils.filled3,tree.mamos, where = 25, position = 188.3652)
plot(tree.amniota.fossils.filled4, show.node.label = TRUE, font = 0.1, root.edge = TRUE)

tree.amniotafossils.full <- drop.tip(tree.amniota.fossils.filled4,c("Aves", "Mammalia", "Neosuchia", "Lepidosauria"))  #ahora cortamos las ramas que sobran, las de los grandes clados que hemos utilizado para anclar las filos calibradas
tree.amniotafossils.full<-ladderize(tree.amniotafossils.full)
plot(tree.amniotafossils.full, type = "fan", show.node.label = TRUE, cex = 0.1, root.edge = TRUE)
nodelabels(text = tree.amniotafossils.full$node.label,
           frame = "n", cex=0.1, col= "blue")
writeNexus(tree.amniotafossils.full, 'TreeAmnioteFossils_final.full.nex') 

tree.ages <- as.matrix(tree.age(tree.paleo.LMs))
nodes.age.paleo <- tree.ages[410:817,]
write.csv(nodes.age.paleo, "nodes.age.paleo.csv")


tree.fossils.base.calib <- read.nexus("amniota_fossils2_consensus.nex")
tree.fossils.base.calib<-ladderize(tree.fossils.base.calib)
plot( tree.fossils.base.calib , show.tip.label = T, no.margin = T, y.lim =  c(-5, 72.907 ))


##INDICES A PARTIR DE RESIDUOS FILOGENETICOS##  
variables.paleo <- read.csv("variables.paleo.csv",header=T,sep=";",row.names=1)

x <- log(variables.paleo$skull.CS)
names(x) <- row.names(variables.paleo)
y <- log(variables.paleo$brain.convhull)
names(y) <- row.names(variables.paleo)
phylresid.headbrainres <- phyl.resid(tree.paleo.LMs, x, y, method = "BM") #### ESTA PARECE LA BUENA, LA ENCEFALIZACION ES REAL

x <- log(variables.paleo$skull.CS.sinpmx)
names(x) <- row.names(variables.paleo)
y <- log(variables.paleo$brain.convhull)
names(y) <- row.names(variables.paleo)
phylresid.headbrainres.sinpmx <- phyl.resid(tree.paleo.LMs, x, y, method = "BM")

##ancestral reconstruction of continuous variables##
variables.paleo <- read.csv("variables.paleo.csv",header=T,sep=";",row.names=1)
variables.paleo <- as.matrix(variables.paleo)
names(variables.paleo) <- row.names(variables.paleo)

x<-variables.paleo[,9] #headbrain.res.conpmx
ace.paleo.headbrainres.conpmx <- ace(x , tree.paleo.LMs , type = "continuous" , method = "pic" , CI = TRUE )
phenogram(tree.paleo.LMs,x,ftype="off",spread.labels=FALSE)
x<-variables.paleo[,10] #headbrain.res.sinpmx
ace.paleo.headbrainres.sinpmx <- ace(x , tree.paleo.LMs , type = "continuous" , method = "pic" , CI = TRUE )
node.headbrain.paleo <- ace.paleo.headbrainres.sinpmx$ace
write.csv(node.headbrain.paleo, 'node.headbrain.paleo.csv')
phenogram(tree.paleo.LMs,x,ftype="off",spread.labels=FALSE)


## GPA ##
rawcoords.paleo<-read.csv("Rawcoords150_confosiles409.csv",header=T,sep=";",row.names=1)  
coords.paleo<-arrayspecs(rawcoords.paleo[,1:150],50,3)
coords.paleo49<-arrayspecs(rawcoords.paleo[,4:150],49,3)

gpa.paleo<-gpagen(coords.paleo)
gpa.paleo49<-gpagen(coords.paleo49)
plotAllSpecimens(gpa.paleo49$coords, mean = TRUE, links = NULL, label = FALSE, 
                 plot.param = list(pt.cex = 0.1, mean.cex = 5)) #para plotear todo viendo la media

## PCA ##
drop.taxa <- tree.amniotafossils.full$tip.label[ ! tree.amniotafossils.full$tip.label %in% rownames(rawcoords.paleo)]
drop.taxa # check te taxa are the expected to be drawn out
tree.paleo.LMs <- drop.tip(tree.amniotafossils.full, drop.taxa)
plot(tree.paleo.LMs, type = "fan", show.node.label = TRUE, cex = 0.1, root.edge = TRUE)
nodelabels(text = tree.paleo.LMs$node.label,
           frame = "n", cex=0.1, col= "blue")

phyloPCA.paleo<-gm.prcomp(gpa.paleo$coords, phy = tree.paleo.LMs)
summary(phyloPCA.paleo)
plot(phyloPCA.paleo, phylo = TRUE, main = "phyloPCA", phylo.par = list(node.labels = FALSE, anc.states = FALSE, node.cex = 0.2, edge.color = "grey", edge.width = 0.5, tip.txt.cex = 0.5))
phyloPCs.paleo<-phyloPCA.paleo$x


phyloPCA.paleo49<-gm.prcomp(gpa.paleo49$coords, phy = tree.paleo.LMs)
summary(phyloPCA.paleo49)
plot(phyloPCA.paleo49, phylo = TRUE, main = "phyloPCA", phylo.par = list(node.labels = TRUE, anc.states = FALSE, node.cex = 0.1, edge.color = "blue", edge.width = 0.1, tip.txt.cex = 0.1))
phyloPCs.paleo49<-phyloPCA.paleo49$x 
write.csv(phyloPCA.paleo49[["ancestors"]], 'ancestral.shapes.paleo.csv')
write.csv(phyloPCA.paleo49[["anc.x"]], 'ancestral.occupation.paleo.csv')
nodes.info.paleo <- read.csv("nodes.info.paleo.csv",header=T,sep=";",row.names=1)
Plot1 <- ggplot(nodes.info.paleo, aes(Comp1, Comp2, label = row.names(nodes.info.paleo), colour = temp.bin, size = 3)) + geom_point(alpha = 1)
Plot2 <- print(Plot1 + geom_text(color = "black", size = 2, check_overlap = T))


phyloPCA.paleobrain<-gm.prcomp(gpa.brain.paleo$coords, phy = tree.paleo.LMs)
summary(phyloPCA.paleobrain)
plot(phyloPCA.paleobrain, phylo = TRUE, main = "phyloPCA", phylo.par = list(node.labels = FALSE, anc.states = FALSE, node.cex = 0.2, edge.color = "grey", edge.width = 0.5, tip.txt.cex = 0.5))
phyloPCs.paleobrain<-phyloPCA.paleobrain$x

Rawcoords_brain3.paleo<-read.csv("Rawcoordsbrain3_confosiles409.csv",header=T,sep=";",row.names=1)
coords.brain3.paleo<-arrayspecs(Rawcoords_brain3.paleo[,1:48],16,3)
gpa.brain3.paleo <- gpagen(coords.brain3.paleo)
phyloPCA.paleobrain3<-gm.prcomp(gpa.brain3.paleo$coords, phy = tree.paleo.LMs)
summary(phyloPCA.paleobrain3)
plot(phyloPCA.paleobrain3, phylo = TRUE, main = "phyloPCA", phylo.par = list(node.labels = TRUE, anc.states = TRUE, node.cex = 0.1, edge.color = "blue", edge.width = 0.1, tip.txt.cex = 0.1))
write.csv(phyloPCA.paleobrain3[["ancestors"]], 'ancestral.shapes.brain.paleo.csv')
write.csv(phyloPCA.paleobrain3[["anc.x"]], 'ancestral.occupation.brain.paleo.csv')
nodes.info.paleo <- read.csv("nodes.info.paleo.csv",header=T,sep=";",row.names=1)
Plot1 <- ggplot(nodes.info.paleo, aes(Comp1.brain, Comp2.brain, label = row.names(nodes.info.paleo), colour = temp.bin, size = 3)) + geom_point(alpha = 1)
Plot2 <- print(Plot1 + geom_text(color = "black", size = 2, check_overlap = T))
coords.anc.brain.shapes<-arrayspecs(phyloPCA.paleobrain3[["ancestors"]][,1:48],16,3)
phyloPCs.paleobrain<-phyloPCA.paleobrain3$x
write.csv(phyloPCs.paleobrain, 'pcs.brain.paleo.csv')
PC.temp <- phyloPCA.paleobrain3$x[,2] 
preds<-shape.predictor(gpa.brain3.paleo$coords, x = as.numeric( phyloPCA.paleobrain3$x[,2] ), Intercept = F, 
                       pred1 = quantile(PC.temp, probs = 0.05), pred2 = quantile(PC.temp, probs = 0.95)) 
#modificar los valores del probs para exagerar o suavizar los cambios de forma
#calculates per landmarks variance 
library(tibble)
per_lm_variance <- function (shape.data) #para modificar los colores de los LMs, meter los codigos que gusten en cols1
{
  variances <- rowSums(apply(shape.data, c(1, 2), var))
  cols1 <- colorRampPalette(c("gray88", "yellow", "red"))
  cols <- cols1(100)
  x = (log10(variances))
  xlims <- NULL
  tol <- 1e-06
  xlims <- range(x) + c(-tol, tol)
  nbin = 100
  breaks <- 0:nbin/nbin * (xlims[2] - xlims[1]) + xlims[1]
  whichColor <- function(p, cols, breaks) {
    i <- 1
    while (p >= breaks[i] && p > breaks[i + 1]) i <- i +
        1
    cols[i]
  }
  variancecolors <- sapply(x, whichColor, cols = cols, breaks = breaks)
  variance_table <- tibble(Per_Lm_Variance = variances, Log_Variance = x,Variance_Colors = variancecolors)
  return(variance_table)
}
my.variances <- per_lm_variance(shape.data = gpa.brain3.paleo$coords)
per_lm_distance <- function (shape.data) #para modificar los colores de los LMs, meter los codigos que gusten en cols1
{
  shape.data <- preds
  df <- as.data.frame(cbind(shape.data$pred1, shape.data$pred2))
  distances <- apply(df, 1, function(x) dist(matrix(x, nrow = 2, byrow = TRUE)))
  cols1 <- colorRampPalette(c("gray88", "yellow", "red"))
  cols <- cols1(100)
  x = (log10(distances))
  xlims <- NULL
  tol <- 1e-06
  xlims <- range(x) + c(-tol, tol)
  nbin = 100
  breaks <- 0:nbin/nbin * (xlims[2] - xlims[1]) + xlims[1]
  whichColor <- function(p, cols, breaks) {
    i <- 1
    while (p >= breaks[i] && p > breaks[i + 1]) i <- i +
        1
    cols[i]
  }
  distancecolors <- sapply(x, whichColor, cols = cols, breaks = breaks)
  distance_table <- tibble(Per_Lm_Distance = distances, Log_Distance = x,Distance_Colors = distancecolors)
  return(distance_table)
}
my.distances <- per_lm_distance(shape.data = preds)
#now we need to put all shapes in the same space to plot warped mesh and hot.dots on top of it
open3d()
par3d(windowRect = c(0,0,500,500)) #para dimensiones del cuadro RGL
Sys.sleep(1)
mfrow3d(nr = 1, nc = 2, byrow = TRUE, sharedMouse = TRUE)
spheres3d(preds$pred1, col = my.distances$Distance_Colors, radius = 0.01) #radius cambia el tamaño del LM
next3d()
spheres3d(preds$pred2, col = my.distances$Distance_Colors, radius = 0.01)


angles.extantfossils5 <- read.csv("angles.extantxfossils.csv",header=T,sep=";",row.names=1)
angles.paleo <- angles.extantxfossils[,5:9]
names(angles.paleo) <- row.names(angles.extantxfossils)
angles.paleo <- angles.paleo[-c(391,395,404,406,408,410),]
names.angles.paleo <- row.names(angles.paleo)
angles.paleo <- matrix(as.numeric(angles.paleo), ncol = 5)
row.names(angles.paleo) <- names.angles.paleo

tree.amniotafossils.pru <- drop.tip(tree.amniotafossils.full, c("Thecodontosaurus","Gnathovorax","Velociraptor","Shuvuuia","Enant_Brasil","Ixalerpeton"))
plot(tree.amniotafossils.pru, show.node.label = TRUE, cex = 0.1, root.edge = TRUE)

phyloPCA.angles <- gm.prcomp(angles.paleo, phy = tree.amniotafossils.pru)   
summary(phyloPCA.angles)
plot(phyloPCA.angles, axis1 = 1, axis2 = 2)
plot(phyloPCA.angles, phylo = TRUE, main = "phyloPCA", loadings = TRUE , phylo.par = list(node.labels = FALSE, anc.states = FALSE, node.cex = 0.2, edge.color = "grey", edge.width = 0.5, tip.txt.cex = 0.5))
arrows(0, 0, phyloPCA.angles$rotation[,1]*100, phyloPCA.angles$rotation[,2]*100, length = 0.2)



#Plotting
angles.extantfossils5 <- as.matrix(angles.extantfossils5)
df<-data.frame(clades.a = angles.extantfossils5$clades.a, clades.b = angles.extantfossils5$clades.b,
               clades.c = angles.extantfossils5$clades.comp, clades.d = angles.extantfossils5$groups.temporal, 
               log.skullCS = variables.paleo$skull.CS.sinpmx,
               headbrain.res.sinpmx = variables.paleo$headbrain.res.sinpmx, headbrain.res.conpmx = variables.paleo$headbrain.res.conpmx,
               PC1skull = phyloPCs.paleo49[,1], PC2skull = phyloPCs.paleo49[,2], PC3skull = phyloPCs.paleo49[,3], 
               PC1brain = phyloPCs.paleobrain[,1], PC2brain = phyloPCs.paleobrain[,2], PC3brain = phyloPCs.paleobrain[,3], 
               PC1angles = phyloPCA.angles[,1], PC2angles = phyloPCA.angles[,2],
               row.names = row.names(variables.paleo))
df<-data.frame(clades.a = angles.extantfossils5[,1], clades.b = angles.extantfossils5[,2],
               clades.c = angles.extantfossils5[,3], clades.d = angles.extantfossils5[,4], 
               PC1angles = phyloPCA.angles$x[,1], PC2angles = phyloPCA.angles$x[,2],
               row.names = row.names(variables.paleo))

Plot1 <- ggplot(df, aes(PC1brain, PC2brain, label = row.names(df), colour = clades.d, size = log.skullCS)) + geom_point(alpha = 1)
Plot2 <- print(Plot1 + geom_text(color = "black", size = 2, check_overlap = T))


#para mapear una variable continua sobre la filogenia del morfoespacio y ver su distribución en filo y morfo
x<-variables.paleo$headbrain.res.sinpmx
y<-variables.paleo$headbrain.res.conpmx
names(y)<-row.names(variables.paleo)
tree<-tree.paleo.LMs
phylo.headbrain.paleo<-contMap(tree,y,res=100,lwd=3,legend = 40,outline=FALSE,sig=3,
                               type="phylogram",direction="rightwards",plot=FALSE)
obj<-setMap(phylo.headbrain.paleo,colors=c("#64B3E0","#64B3E0","#64B3E0","#E8E37F","#E8E37F","#D83E39","#D83E39"))
plot(obj, outline = FALSE)
phylomorph_cont_headbrain<-phylomorphospace(obj$tree, phyloPCs.paleobrain[,1:2], colors = obj$cols, lwd = 4, xlab = "PC1", ylab = "PC2")


#FILOMORFOESPACIO 3D (cronofilomorfoespacio) del PCA PLOTEANDO TIEMPO
cronophymorph.paleo <- plot(phyloPCA.paleo49, axis1 = 1, axis2 = 2, phylo = TRUE, time.plot = TRUE)
#NO FUNCIONA AUN  


##ANGLES FROM ANCESTRAL SHAPES##
anc.shapes<-read.csv("ancestral.key.shapes.paleo.csv",header=T,sep=";",row.names=1)  
coords.anc.shapes<-arrayspecs(anc.shapes[,1:147],49,3)

angulos.orientation<-list()
for (i in 1: length(coords.anc.shapes)) { #para probar cada línea, empezar con i = 1 e ir corriendo cada línea. Si funcionan todas, quitar la i = 1 y correr la función
  craneo.individual <- coords.anc.shapes[,,i] #coge el elemento del número que vaya tocando en cada iteración, 1 para la primera, 2 para la segunda hasta la 250 vez.
  vect1 <- matrix(c((coords.anc.shapes[24,,i])-(coords.anc.shapes[14,,i])), nrow = 3) #tener en cuenta que la numeración de los LMs es sin pmx
  vect2 <- matrix(c((coords.anc.shapes[15,,i])-(coords.anc.shapes[14,,i])), nrow = 3)
  angulo.individual <- acos(sum(vect1*vect2)/(sqrt(sum(vect1*vect1))*sqrt(sum(vect2*vect2))))
  angulos.orientation[[i]]<-angulo.individual #metes el craneo mirroreado en una lista nueva
}

angulos.orientation <- as.matrix(angulos.orientation)
write.csv2(angulos.orientation, 'clivus.anc.shapes.csv')

#ancestral reconstruction of continuous variables (angles)

variables.paleo <- read.csv("variables.paleo.csv",header=T,sep=";",row.names=1)
variables.paleo <- as.matrix(variables.paleo)
names(variables.paleo) <- row.names(variables.paleo)

x<-variables.paleo[,11] #foramen
ace.foramen.paleo <- ace(x , tree.paleo.LMs , type = "continuous" , method = "pic" , CI = TRUE )
x<-variables.paleo[,12] #base
ace.base.paleo <- ace(x , tree.paleo.LMs , type = "continuous" , method = "pic" , CI = TRUE )
x<-variables.paleo[,13] #flexion
ace.flexion.paleo <- ace(x , tree.paleo.LMs , type = "continuous" , method = "pic" , CI = TRUE )
x<-variables.paleo[,14] #clivus
ace.clivus.paleo <- ace(x , tree.paleo.LMs , type = "continuous" , method = "pic" , CI = TRUE )
x<-variables.paleo[,15] #face
ace.face.paleo <- ace(x , tree.paleo.LMs , type = "continuous" , method = "pic" , CI = TRUE )

x<-variables.paleo[,15] #face
ace.face.paleo <- ace(x , tree.paleo.LMs , type = "continuous" , method = "pic" , CI = TRUE )

## RECONSTRUCCION DE MORFOLOGIAS ANCESTRALES Y DIFERENCIAS ENTRE LM
#cargar meshes y configuraciones de referencia y transformar
#mesh.base.amniota <- read.ply("Pristerodon_base.ply", ShowSpecimen = F, addNormals = F )
#coords.original <- coords.paleo49[,,"Pristerodon"]
#M <- mshape(gpa.paleo49$coords)
#mesh.mshape.amniota <- tps3d(mesh.base.amniota, refmat = coords.original, tarmat = M)
anc.amniota <- coords.anc.brain.shapes[,,"410"]
#mesh.anc.amniota <- tps3d(mesh.mshape.amniota, refmat = M, tarmat = anc.amniota)
#mesh2ply(mesh.anc.amniota, filename = "mesh.anc.amniota", col = NULL, writeNormals = FALSE)

#transformaciones desde la basal
anc.mammalia <- coords.anc.brain.shapes[,,"668"]
tarmat <- anc.mammalia
#mesh.anc.aves <- tps3d(mesh.anc.theropoda, refmat = anc.theropoda, tarmat = tarmat)
#mesh2ply(mesh.anc.aves, filename = "mesh.anc.aves", col = NULL, writeNormals = FALSE)
reftar.shapes <- as.data.frame(cbind(anc.amniota, tarmat))
colnames(reftar.shapes) <- c("X","Y","Z","X","Y","Z")
library(tibble)
per_lm_variance <- function (shape.data) #para modificar los colores de los LMs, meter los codigos que gusten en cols1
{
  variances <- rowSums(apply(shape.data, c(1, 2), var))
  cols1 <- colorRampPalette(c("gray88", "yellow", "red"))
  cols <- cols1(100)
  x = (log10(variances))
  xlims <- NULL
  tol <- 1e-06
  xlims <- range(x) + c(-tol, tol)
  nbin = 100
  breaks <- 0:nbin/nbin * (xlims[2] - xlims[1]) + xlims[1]
  whichColor <- function(p, cols, breaks) {
    i <- 1
    while (p >= breaks[i] && p > breaks[i + 1]) i <- i +
        1
    cols[i]
  }
  variancecolors <- sapply(x, whichColor, cols = cols, breaks = breaks)
  variance_table <- tibble(Per_Lm_Variance = variances, Log_Variance = x,Variance_Colors = variancecolors)
  return(variance_table)
}
my.variances <- per_lm_variance(shape.data = gpa.brain3.paleo$coords)
per_lm_distance <- function (shape.data) #para modificar los colores de los LMs, meter los codigos que gusten en cols1
{
  df <- shape.data
  distances <- apply(df, 1, function(x) dist(matrix(x, nrow = 2, byrow = TRUE)))
  cols1 <- colorRampPalette(c("#1a2a6c", "#7B90E3", "#fdbb2d", "#b21f1f"))
  cols <- cols1(100)
  x = (log10(distances))
  xlims <- NULL
  tol <- 1e-06
  xlims <- range(x) + c(-tol, tol)
  nbin = 100
  breaks <- 0:nbin/nbin * (xlims[2] - xlims[1]) + xlims[1]
  whichColor <- function(p, cols, breaks) {
    i <- 1
    while (p >= breaks[i] && p > breaks[i + 1]) i <- i +
        1
    cols[i]
  }
  distancecolors <- sapply(x, whichColor, cols = cols, breaks = breaks)
  distance_table <- tibble(Per_Lm_Distance = distances, Log_Distance = x,Distance_Colors = distancecolors)
  return(distance_table)
}
my.distances <- per_lm_distance(shape.data = reftar.shapes)
#now we need to put all shapes in the same space to plot warped mesh and hot.dots on top of it
open3d()
par3d(windowRect = c(0,0,500,500)) #para dimensiones del cuadro RGL
Sys.sleep(1)
mfrow3d(nr = 1, nc = 2, byrow = TRUE, sharedMouse = TRUE)
spheres3d(anc.amniota, col = my.variances$Variance_Colors, radius = 0.01)
next3d()
spheres3d(tarmat, col = my.distances$Distance_Colors, radius = 0.01)


fin



#calculates per landmarks variance 
library(tibble)
per_lm_variance <- function (shape.data) #para modificar los colores de los LMs, meter los codigos que gusten en cols1
{
  variances <- rowSums(apply(shape.data, c(1, 2), var))
  cols1 <- colorRampPalette(c("gray88", "yellow", "red"))
  cols <- cols1(100)
  x = (log10(variances))
  xlims <- NULL
  tol <- 1e-06
  xlims <- range(x) + c(-tol, tol)
  nbin = 100
  breaks <- 0:nbin/nbin * (xlims[2] - xlims[1]) + xlims[1]
  whichColor <- function(p, cols, breaks) {
    i <- 1
    while (p >= breaks[i] && p > breaks[i + 1]) i <- i +
        1
    cols[i]
  }
  variancecolors <- sapply(x, whichColor, cols = cols, breaks = breaks)
  variance_table <- tibble(Per_Lm_Variance = variances, Log_Variance = x,Variance_Colors = variancecolors)
  return(variance_table)
}
my.variances <- per_lm_variance(shape.data = GPA.fit$coords)
per_lm_distance <- function (shape.data) #para modificar los colores de los LMs, meter los codigos que gusten en cols1
{
  df <- as.data.frame(cbind(shape.data$pred1, shape.data$pred2))
  distances <- apply(df, 1, function(x) dist(matrix(x, nrow = 2, byrow = TRUE)))
  cols1 <- colorRampPalette(c("gray88", "yellow", "red"))
  cols <- cols1(100)
  x = (log10(distances))
  xlims <- NULL
  tol <- 1e-06
  xlims <- range(x) + c(-tol, tol)
  nbin = 100
  breaks <- 0:nbin/nbin * (xlims[2] - xlims[1]) + xlims[1]
  whichColor <- function(p, cols, breaks) {
    i <- 1
    while (p >= breaks[i] && p > breaks[i + 1]) i <- i +
        1
    cols[i]
  }
  distancecolors <- sapply(x, whichColor, cols = cols, breaks = breaks)
  distance_table <- tibble(Per_Lm_Distance = distances, Log_Distance = x,Distance_Colors = distancecolors)
  return(distance_table)
}
my.distances <- per_lm_distance(shape.data = preds)

#now we need to put all shapes in the same space to plot warped mesh and hot.dots on top of it
open3d()
par3d(windowRect = c(0,0,500,500)) #para dimensiones del cuadro RGL
Sys.sleep(1)
mfrow3d(nr = 1, nc = 2, byrow = TRUE, sharedMouse = TRUE)

mesh.min <- tps3d(mesh.original, refmat = coords.original, tarmat = preds$pred1) #maps landmarks or a triangular mesh via thin plate spline based on a reference and a target configuration in 2D and 3D
shade3d(mesh.min, col= 8, alpha = 0.3) #dibuja la mesh con la forma predicha en la línea anterior #alpha cambia la transparencia de la surface
spheres3d(preds$pred1, col = my.distances$Distance_Colors, radius = 0.01) #radius cambia el tamaño del LM

next3d()
mesh.max <- tps3d(mesh.original, refmat = coords.original, tarmat = preds$pred2, lambda = 0.5)
shade3d(mesh.max, col= 8, alpha = 0.3) 
spheres3d(preds$pred2, col = my.distances$Distance_Colors, radius = 0.01)


mesh2ply(mesh.max, filename = "mesh.max", col = NULL, writeNormals = FALSE)


## RIDGES & VIOLINS FOR THE EVOLUTIONARY HISTORY OF BRAIN SIZE AND FLEXION##
# RIDGELINE PLOTS
library(ggridges)
library(ggplot2)
library(viridis)
library(hrbrthemes)

variables.paleo <- read.csv("variables.paleo.csv",header=T,sep=";",row.names=1)
df<-data.frame(clades.a = variables.paleo$clades.a, clades.b = variables.paleo$clades.b,
               clades.c = variables.paleo$clades.comp, clades.d = variables.paleo$groups.temporal, 
               headbrain.res.sinpmx = variables.paleo$headbrain.res.sinpmx, headbrain.res.conpmx = variables.paleo$headbrain.res.conpmx,
               PC1 = phyloPCs.paleobrain[,1], PC2 = phyloPCs.paleobrain[,2], PC3 = phyloPCs.paleobrain[,3], 
               row.names = row.names(variables.paleo))

plot1 <- ggplot(df, aes(x = PC2,y = clades.c, fill = stat(x))) + 
  geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.01) +
  scale_fill_gradientn(#colours = c("#1a2a6c","#1a2a6c", "#838bae", "#f5f5f5","#b21f1f","#b21f1f","#fdbb2d"),
    colours = c("#64B3E0","#64B3E0","#E8E37F","#E8E37F","#D83E39"),
    values = NULL, space = "Lab", guide = "colourbar", aesthetics = "fill") +
  labs(title = 'PC2 brain shape') +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )
plot2 <- print(plot1)



brain.fossilsprunned <- read.csv("headbrain.fossils.prunned.csv",header=T,sep=";",row.names=1)
brain.fossilsprunned <- as.matrix(brain.fossilsprunned)
names(brain.fossilsprunned) <- row.names(brain.fossilsprunned)
drop.taxa <- tree.fossils.base.calib$tip.label[ ! tree.fossils.base.calib$tip.label %in% rownames(brain.fossilsprunned)]
drop.taxa # check te taxa are the expected to be drawn out
tree.brainfossils.base.calib <- drop.tip(tree.fossils.base.calib, drop.taxa)
tree.brainfossils.base.calib<-ladderize(tree.brainfossils.base.calib)
plot(tree.brainfossils.base.calib, type = "fan", show.node.label = TRUE, cex = 1, root.edge = TRUE)
nodelabels(text = tree.brainfossils.base.calib$node.label,
           frame = "n", cex=0.1, col= "blue")
plotTree.wBars(tree.brainfossils.base.calib,brain.fossilsprunned[,3], scale = NULL, tip.labels = TRUE)



## DISPARIDAD ##
variables.paleo <- read.csv("variables.paleo.csv",header=T,sep=";",row.names=1)
gdf.paleo.disp <- geomorph.data.frame(gpa.paleo, grouping.a = variables.paleo$clades.a, grouping.b = variables.paleo$clades.b, 
                                      grouping.c = variables.paleo$clades.comp, grouping.temp = variables.paleo$groups.temporal,
                                      grouping.temp2 = variables.paleo$groups.temporal2)
gdf.paleo.disp.sinpmx <- geomorph.data.frame(gpa.paleo49, grouping.a = variables.paleo$clades.a, grouping.b = variables.paleo$clades.b, 
                                             grouping.c = variables.paleo$clades.comp, grouping.temp = variables.paleo$groups.temporal,
                                             grouping.temp2 = variables.paleo$groups.temporal2)

gdf.paleo.disp.brain <- geomorph.data.frame(gpa.brain3.paleo, grouping.a = variables.paleo$clades.a, grouping.b = variables.paleo$clades.b, 
                                            grouping.c = variables.paleo$clades.comp, grouping.temp = variables.paleo$groups.temporal,
                                            grouping.temp2 = variables.paleo$groups.temporal2)

disparity.paleo <- morphol.disparity(coords ~ grouping.temp2, groups = ~ grouping.temp2, data = gdf.paleo.disp, iter=999)
summary(disparity.paleo)

disparity.paleo.sinpmx <- morphol.disparity(coords ~ grouping.temp2, groups = ~ grouping.temp2, data = gdf.paleo.disp.sinpmx, iter=999)
summary(disparity.paleo.sinpmx)

disparity.paleo.brain <- morphol.disparity(coords ~ grouping.temp2, groups = ~ grouping.temp2, data = gdf.paleo.disp.brain, iter=999)
summary(disparity.paleo.brain)


#DISPARITY THROUGH TIME
#skull
matrix.coords <- two.d.array(gpa.paleo$coords)
DTT.paleo <- dtt(phy = tree.paleo.LMs, data = matrix.coords, index = c("avg.sq"), mdi.range = c(325,0), nsim = 1000, CI = 0.95, plot = TRUE)      
matrix.coordssin <- two.d.array(gpa.paleo49$coords)
DTT.paleo49 <- dtt(phy = tree.paleo.LMs, data = matrix.coordssin, index = c("avg.sq"), mdi.range = c(325,0), nsim = 1000, CI = 0.95, plot = TRUE)      
#brain
matrix.coordsbrain <- two.d.array(gpa.brain.paleo$coords) 
DTT.paleobrain <- dtt(phy = tree.paleo.LMs, data = matrix.coordsbrain, index = c("avg.sq"), mdi.range = c(325,0), nsim = 1000, CI = 0.95, plot = TRUE)      
#PCs skull    
PCs.con <- phyloPCs.paleo[,1:16]
DTT.paleo.PCs <- dtt(phy = tree.paleo.LMs, data = PCs.con, index = c("avg.sq"), mdi.range = c(325,0), nsim = 1000, CI = 0.95, plot = TRUE)
PCs.sin <- phyloPCs.paleo49[,1:19]
DTT.paleo49.PCs <- dtt(phy = tree.paleo.LMs, data = PCs.sin, index = c("avg.sq"), mdi.range = c(325,0), nsim = 1000, CI = 0.95, plot = TRUE)
PCs.brain <- phyloPCs.paleobrain[,1:12]
DTT.paleobrain.PCs <- dtt(phy = tree.paleo.LMs, data = PCs.brain, index = c("avg.sq"), mdi.range = c(325,0), nsim = 1000, CI = 0.95, plot = TRUE)
PCs.brain.filtrado <- phyloPCs.paleobrain[,1:3]
DTT.paleobrain.PCs.filt <- dtt(phy = tree.paleo.LMs, data = PCs.brain.filtrado, index = c("avg.sq"), mdi.range = c(325,0), nsim = 1000, CI = 0.95, plot = TRUE)
#headbrain.res
headbrain.res.sinpmx <- variables.paleo[,10]  
names(headbrain.res.sinpmx) <- row.names(variables.paleo)
DTT.paleo.headbrainsin <- dtt(phy = tree.paleo.LMs, data = headbrain.res.sinpmx, index = c("avg.sq"), mdi.range = c(325,0), nsim = 1000, CI = 0.95, plot = TRUE)
