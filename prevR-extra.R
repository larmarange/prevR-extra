require(prevR)
require(gstat)

Noptim <- function(object) {
  clusters               = slot(object,"clusters")
  clustersNumber         = nrow(clusters)
  ObservationNumber      = sum(clusters$n)
  PositiveCases          = sum(clusters$pos)
  isWeightedData         = !any(is.na(match(c("wn","wpos"),names(clusters))))
  nationalPrev           = 100*sum(clusters$pos,na.rm=T)/sum(clusters$n,na.rm=T)
  if(isWeightedData){
    nationalPrev = 100*sum(clusters$wpos,na.rm=T)/sum(clusters$wn,na.rm=T)
  }
  if (nationalPrev>50) {
    nationalPrev <- (100 - nationalPrev)
  } 
  round(14.172*ObservationNumber^0.419*nationalPrev^-0.361*clustersNumber^0.037-91.011)
}

# Source: R Graphics Cookbook pp. 317-18
theme_prevR <- function (base_size = 12) {
  require(grid)
  theme_grey(base_size) %+replace%
    theme(
      axis.title        = element_blank(),
      axis.text         = element_blank(),
      #panel.background  = element_blank(),
      #panel.grid        = element_blank(),
      axis.ticks.length = unit(0, "cm"),
      axis.ticks.margin = unit(0, "lines"),
      plot.margin       = unit(c(0, 0, 0, 0), "lines"),
      complete          = TRUE
    )
}

theme_prevR_light <- function (base_size = 12) {
  require(grid)
  theme_grey(base_size) %+replace%
    theme(
      axis.title        = element_blank(),
      axis.text         = element_blank(),
      panel.background  = element_blank(),
      panel.grid        = element_blank(),
      axis.ticks.length = unit(0, "cm"),
      axis.ticks.margin = unit(0, "lines"),
      plot.margin       = unit(c(0, 0, 0, 0), "lines"),
      complete          = TRUE
    )
}

prevR <- function(object, N=Noptim(object), nb.cells=100, cell.size=NULL, weighted=NULL, plot.results=TRUE, return.results=FALSE, legend.title="%", cex=0.7, progression=TRUE) {
  if (!require(prevR))
    stop("Package prevR is required.")
  if (!require(gstat))
    stop("Package gstat is required.")
  if (!is.prevR(object))
    stop("object must be of class prevR.")
  if (plot.results & !require(ggplot2))
    stop("Package ggplot2 is required for plotting.")
  if (plot.results & !require(directlabels))
    stop("Package directlabels is required for plotting.")
  if(is.null(weighted)) # If not precised, weighted data if available
    weighted <- "wn" %in% colnames(object@clusters)
  if (weighted & !("wn" %in% colnames(object@clusters)))
    stop("No weighted data found in object.")
  
  object <- rings(object, N=N, R=Inf, progression=progression)
  prev <- kde(object, N=N, R=Inf, weighted=weighted, nb.cells=nb.cells, cell.size=cell.size, progression=progression)
  radius <- krige(r.radius~1, object, N=N, R=Inf, nb.cells=nb.cells, cell.size=cell.size, debug.level=if (progression) -1 else 0)
  
  dimnames(prev@coords)[[2]] <- c("x","y") # Fixing name of coordinates
  dimnames(radius@coords)[[2]] <- c("x","y") # Fixing name of coordinates
  
  if(plot.results) {
    prev.df <- na.omit(as.data.frame(prev))
    radius.df <- na.omit(as.data.frame(radius))
    
    prev.df <- prev.df[order(prev.df$x,prev.df$y),]
    radius.df <- radius.df[order(radius.df$x, radius.df$y),]
    
    r <- data.frame()
    for (n in N) {
      if (weighted) {
        p.name <- paste0("k.wprev.N",n,".RInf")
      } else {
        p.name <- paste0("k.prev.N",n,".RInf")
      }
      r.name <- paste0("r.radius.N",n,".RInf")
      temp <- data.frame(N=n,x=prev.df[["x"]],y=prev.df[["y"]],prev=prev.df[[p.name]],radius=radius.df[[r.name]])
      r <- rbind(r,temp)
    }
    
    r$N <- factor(r$N, levels=N, labels=paste0("N=",N))
    
    p <- ggplot(r, aes(x=x, y=y, fill=prev, z=radius)) +
      geom_raster() +
      scale_fill_gradientn(colours=prevR.colors.red(20)) +
      stat_contour(aes(colour = ..level..)) +
      scale_colour_continuous(low = "grey50", high = "grey50") +
      facet_wrap(~ N) +
      coord_fixed() + xlab("") + ylab("") +
      labs(fill=legend.title) +
      theme_prevR()
    print(direct.label.prevR(p,list("top.pieces",cex=cex)))
  }
  
  if (return.results)
    return(list(prev=prev,radius=radius))
}

direct.label.prevR <- function (p, method = NULL, debug = FALSE) 
{
  require(ggplot2)
  getData <- function(colour.or.fill) {
    for (L in p$layers) {
      m <- p$mapping
      m[names(L$mapping)] <- L$mapping
      colvar <- m[[colour.or.fill]]
      if (!is.null(colvar)) {
        return(list(layer = L, colvar = as.character(colvar)))
      }
    }
  }
  dl.info <- getData("colour")
  if (is.null(dl.info)) {
    dl.info <- getData("fill")
  }
  if (is.null(dl.info)) {
    stop("Need colour or fill aesthetic to infer default direct labels.")
  }
  L <- dl.info$layer
  colvar <- dl.info$colvar
  geom <- L$geom$objname
  if (is.null(method)) 
    method <- default.picker("ggplot")
  data <- if ((!is.null(L$data)) && (length(L$data) > 0)) {
    L$data
  }
  else {
    NULL
  }
  a <- aes_string(label = colvar, colour = colvar)
  a2 <- structure(c(L$mapping, a), class = "uneval")
  dlgeom <- geom_dl(a2, method, stat = L$stat, debug = debug, 
                    data = data)
  dlgeom$stat_params <- L$stat_params
  # leg.info <- legends2hide(p)
  leg.info <- NULL # We want to keep the legend
  guide.args <- as.list(rep("none", length(leg.info$hide)))
  names(guide.args) <- leg.info$hide
  guide.args$colour <- "none"
  guide <- do.call(guides, guide.args)
  p + dlgeom + guide
}

prevR.comp <- function(A, B, N = Noptim(A), weighted=NULL, nb.cells=100, cell.size=NULL, plot.results=TRUE, return.results=FALSE, labA="A", labB="B", legend.title="%", cex=0.7, progression=TRUE) {
  if (!is.prevR(A))
    stop("A must be of class prevR.")
  if (!is.prevR(B))
    stop("B must be of class prevR.")
  if (length(N)>1)
    stop("N should be a single number.")
  if (plot.results & !require(ggplot2))
    stop("Package ggplot2 is required for plotting.")
  if (plot.results & !require(directlabels))
    stop("Package directlabels is required for plotting.")
  if (plot.results & !require(gridExtra))
    stop("Package gridExtra is required for plotting.")
  if (A@proj@projargs != B@proj@projargs)
    B <- changeproj(B, A@proj)
  A <- rings(A, N, progression=progression)
  
  boundary <- slot(A,"boundary")
  proj <- slot(slot(A,"proj"),"projargs")
  longlat <- F
  if(regexpr("longlat",proj) != -1 || regexpr("latlong",proj) != -1){
    longlat <- T
  }
  
  if(is.null(weighted)) # If not precised, weighted data if available
    weighted <- ("wn" %in% colnames(A@clusters) & "wn" %in% colnames(B@clusters))
  if (weighted & !("wn" %in% colnames(A@clusters) & "wn" %in% colnames(B@clusters)))
    stop("No weighted data found in A and/or B.")
  
  # Calcul de la grille
  coord = coordinates(as.SpatialGrid(A, nb.cells=nb.cells, cell.size=cell.size))
  range.x = unique(coord[,1])
  range.y = unique(coord[,2])
  
  ringName = paste0("N",N,".RInf")
  ring = A@rings[[ringName]]
  
  dataA <- merge(A@clusters, ring[["estimates"]], by="id")
  sp.dataA <- dataA
  coordinates(sp.dataA) = ~x+y
  sp.dataA@proj4string = A@proj
  
  dataB <- B@clusters
  sp.dataB <- dataB
  coordinates(sp.dataB) = ~x+y
  sp.dataB@proj4string = B@proj
  
  # Rayons pour B
  sample.vario   = variogram(r.radius~1, data=sp.dataA)
  param1Vgm      = max(sample.vario[,"gamma"])*0.66
  param2Vgm      = max(sample.vario[,"dist"])*0.5
  param          = .init.exp.model.variogram(sample.vario[,"dist"],sample.vario[,"gamma"])
  if(is.null(param)){
    warning(gettextf("problem to fit the variogram"))
  }
  one.model      = try(fit.variogram(sample.vario, model = vgm(param[1],'Exp',param[2])),silent =T ) 
  if(attr(one.model,"class")=="try-error" || attr(one.model,"singular")) one.model = vgm(param[1],'Exp',param[2])
  
  B.radius <- krige(r.radius~1, sp.dataA, sp.dataB, model = one.model, debug.level=if (progression) -1 else 0)
  dataB$r.radius <- B.radius@data[[1]]
  
  # KDE
  coord = coordinates(as.SpatialGrid(A, nb.cells=nb.cells, cell.size=cell.size))
  
  prevA <- prev.kde(dataA, coord, weighted=weighted, longlat=longlat, progression=progression)
  prevB <- prev.kde(dataB, coord, weighted=weighted, longlat=longlat, progression=progression)
  
  if(weighted)
    prev <- list(x=prevA$x, y=prevA$y, prevA=prevA$k.wprev, prevB=prevB$k.wprev, diff=prevB$k.wprev-prevA$k.wprev)
  else
    prev <- list(x=prevA$x, y=prevA$y, prevA=prevA$k.prev, prevB=prevB$k.prev, diff=prevB$k.prev-prevA$k.prev)
  
  # Passage en data.frame
  prev = xyz2dataframe(prev,'x','y',names(prev)[c(-1,-2)])
  # puis en SpatialPixelsDataFrame
  prev = SpatialPixelsDataFrame(points=prev[,c('x','y')], data=prev[,c(-1,-2)])
  prev@proj4string = A@proj
  
  # Suppression des points hors de la zone d'etude
  if (attr(boundary,"valid")) {
    prev = NA.outside.SpatialPolygons(prev, boundary)
  }
  
  # Calcul des radius
  radius <- krige(r.radius~1, A, N=N, R=Inf, nb.cells=nb.cells, cell.size=cell.size, debug.level=if (progression) -1 else 0)
  
  dimnames(prev@coords)[[2]] <- c("x","y") # Fixing name of coordinates
  dimnames(radius@coords)[[2]] <- c("x","y") # Fixing name of coordinates
  
  if(plot.results) {
    prev.df <- na.omit(as.data.frame(prev))
    radius.df <- na.omit(as.data.frame(radius))
    
    prev.df <- prev.df[order(prev.df$x,prev.df$y),]
    radius.df <- radius.df[order(radius.df$x, radius.df$y),]
    
    r <- data.frame()
    r.name <- paste0("r.radius.N",N,".RInf")
    temp <- data.frame(group="A",x=prev.df[["x"]],y=prev.df[["y"]],prev=prev.df[["prevA"]],radius=radius.df[[r.name]])
    r <- rbind(r,temp)
    temp <- data.frame(group="B",x=prev.df[["x"]],y=prev.df[["y"]],prev=prev.df[["prevB"]],radius=radius.df[[r.name]])
    r <- rbind(r,temp)
    
    diff <- data.frame(x=prev.df[["x"]], y=prev.df[["y"]], diff=prev.df[["diff"]], radius=radius.df[[r.name]])
    
    r$group <- factor(r$group, levels=c("A","B"), labels=c(labA,labB))
    
    p <- ggplot(r, aes(x=x, y=y, fill=prev, z=radius)) +
      geom_raster() +
      scale_fill_gradientn(colours=prevR.colors.red(20)) +
      stat_contour(aes(colour = ..level..)) +
      scale_colour_continuous(low = "grey50", high = "grey50") +
      facet_wrap(~ group) +
      coord_fixed() + xlab("") + ylab("") +
      labs(fill=legend.title) +
      theme_prevR()
    p.prev <- direct.label.prevR(p,list("top.pieces",cex=cex))
    
    p <- ggplot(diff, aes(x=x, y=y, fill=diff, z=radius)) +
      geom_raster() +
      scale_fill_gradient2(limits=c(-1*max(abs(diff$diff)),max(abs(diff$diff))), low="blue", high="red") +
      stat_contour(aes(colour = ..level..)) +
      scale_colour_continuous(low = "grey50", high = "grey50") +
      coord_fixed() + xlab("") + ylab("") +
      labs(fill=paste("Evolution",labA,labB)) +
      theme_prevR_light()
    p.diff <- direct.label.prevR(p,list("top.pieces",cex=cex))
    
    grid.arrange(p.prev, p.diff)
  }
  
  if (return.results)
    return(list(prev=prev,radius=radius))
}

.init.exp.model.variogram   = function(dist,gamma){
  ###############################################################################################
  # cette fonction estime par moindres carres les parametres de  lissage d'un semi variogram
  # Le modele choisi est fixe (modele Exp)
  # Ces parametres sont utilises comme valeurs initiales du programme d'ajustement fit.variogram appele  
  #       par la fonction krige (quand on est en mode auto)
  ###############################################################################################
  gammaFunc = function(h, A, a){
    A*(1 -exp(-3*h/a)) 
    A*(1 -exp(-h/a)) 
  }
  objectif = function(par,dist,gamma){
    A = par[1]
    a = par[2]
    obj = sum((gammaFunc(dist,A,a) - gamma)^2)
    obj
  }
  par = c(mean(gamma), max(dist)/2)
  opt = optim(par,objectif,dist =  dist, gamma =gamma,lower = c(-Inf,0.02),method = "L-BFGS-B")
  if(opt$convergence!=0) return(NULL)
  c(psill = opt$par[1],range = opt$par[2])
}

prev.kde <- function(data, coord, weighted, longlat, risk.ratio = FALSE, keep.details = FALSE, progression=TRUE) {
  if (progression) {
    message("Progress of calculations:",domain="R-prevR")
    barre = txtProgressBar(min=0, max=nrow(data), initial=0, style=3)
  }
  
  range.x = unique(coord[,1])
  range.y = unique(coord[,2])
  
  bw   = data[["r.radius"]]
  bwx  = bw
  bwy  = bw
  x    = data[["x"]]
  y    = data[["y"]]
  
  # Attention si x et y sont en degres il faut calculer la largeur de bande  en Km 
  # 6378.388 correspond au rayon terrestre moyen
  if(longlat){
    bwx = bw/(6378.388*cos(pi*y/180))
    bwx = bwx*180/pi
    bwy = bw/6378.388
    bwy = bwy*180/pi
  }
  k.pos = 0
  k.obs = 0
  k.wpos = 0
  k.wobs = 0
  for (i in 1:length(bw)) {
    temp   = KernSur(x=x[i],y=y[i],xbandwidth=bwx[i]/2,ybandwidth=bwy[i]/2,range.x=range.x,range.y=range.y)
    if(weighted == F || weighted == 2){
      k.pos  = k.pos + temp$zden * data[["pos"]][i]
      k.obs  = k.obs + temp$zden * data[["n"]][i]
    }
    if(weighted == T || weighted == 2){
      k.wpos = k.wpos + temp$zden * data[["wpos"]][i]
      k.wobs = k.wobs + temp$zden * data[["wn"]][i]
    }
    if (progression) {
      setTxtProgressBar(barre,value=i)
    }
  }
  if(weighted == F || weighted == 2){
    k.case     = k.pos  / sum(data[["pos"]])
    k.control  = k.obs  / sum(data[["n"]])
    k.prev     = 100*k.pos/k.obs
    k.rr       = 100*k.case/k.control
  }
  if(weighted == T || weighted == 2){
    k.wcase    = k.wpos / sum(data[["wpos"]])
    k.wcontrol = k.wobs / sum(data[["wn"]])
    k.wprev    = 100*k.wpos/k.wobs
    k.wrr      = 100*k.wcase/k.wcontrol
  }
  
  result.one = NULL
  
  if(risk.ratio == F || risk.ratio == 2){
    if(weighted == F || weighted == 2){
      result.one = c(result.one, list(k.pos=k.pos, k.obs=k.obs, k.prev=k.prev))
    }
    if(weighted == T || weighted == 2){
      result.one = c(result.one, list(k.wpos=k.wpos, k.wobs=k.wobs, k.wprev=k.wprev))
    }
  }
  if(risk.ratio == T || risk.ratio == 2){
    if(weighted == F || weighted == 2){
      result.one = c(result.one, list(k.case=k.case, k.control=k.control, k.rr=k.rr))
    }
    if(weighted == T || weighted == 2){
      result.one = c(result.one, list(k.wcase=k.wcase, k.wcontrol=k.wcontrol, k.wrr=k.wrr))
    }
  }
  
  varNames = names(result.one)
  if(!keep.details){
    ind = match(c("k.prev","k.wprev","k.rr","k.wrr"), varNames,nomatch=0)
    result.one = result.one[ind]
  }
  
  result = c(list(x=temp$xords, y=temp$yords),result.one)
  
  if (progression) { close(barre) }
  return(result)
}

loocv <- function(object, progression=TRUE) {
  clusters <- object@clusters
  result <- object@clusters
  result$k.wprev.full <- as.numeric(NA)
  result$k.wprev.loo <- as.numeric(NA)
  result$Noptim <- as.numeric(NA)
  
  if (progression) {
    message("Progress of calculations:",domain="R-prevR")
    barre = txtProgressBar(min=0, max=nrow(clusters), initial=0, style=3)
  }
  
  # Calculated estimated values (full model) for all clusters
  tmp <- rings(object,N=Noptim(object), progression=F)
  r.name <- paste0("N",Noptim(object),".RInf")
  data <- merge(clusters,tmp@rings[[r.name]]$estimates, by="id")
  
  coord <- as.matrix(clusters[,c("x","y")])
  extra.points <- matrix(c(max(data$x)+1,max(data$x)+2,max(data$y)+1,max(data$y)+2),nrow=2)
  
  col <- names(clusters)
  names(col) <- names(clusters)
  
  for (i in 1:nrow(clusters)) {
    points <- rbind(coord[i,],extra.points)
    tmp <- prev.kde(data, points, weighted=T, longlat=T, progression=F)
    result[i,"k.wprev.full"] <- tmp$k.wprev[1,1]
    
    tmp.cl <- clusters[-1*i,]
    tmp.obj <- as.prevR(tmp.cl, col)
    tmp.Nopt <- Noptim(tmp.obj)
    result[i,"Noptim"] <- tmp.Nopt
    tmp.obj <- rings(tmp.obj,N=tmp.Nopt, progression=F)
    r.name <- paste0("N",tmp.Nopt,".RInf")
    tmp.data <- merge(clusters,tmp.obj@rings[[r.name]]$estimates, by="id")
    tmp <- prev.kde(tmp.data, points, weighted=T, longlat=T, progression=F) # A grid of 100x100
    result[i,"k.wprev.loo"] <- tmp$k.wprev[1,1]
    
    if (progression) {setTxtProgressBar(barre,value=i)}
  }
  
  if (progression) {close(barre)}
  return(result)
}

rmse <- function(data, col1, col2) {
  return(sqrt(mean((data[[col1]]-data[[col2]])^2)))
}

# Partitioned Data Hold-Back
pdhb <- function(object, groups, progression=TRUE){
  names(groups) <- c("id","group")
  clusters <- object@clusters
  col <- names(clusters)
  names(col) <- names(clusters)
  clusters <- merge(clusters, groups, by="id")
  coord <- as.matrix(clusters[,c("x","y")])
  
  result <- clusters
  result$k.wprev.full <- as.numeric(NA)
  result$k.wprev.pdhb <- as.numeric(NA)
  result$k.wprev.pdhb.g1 <- as.numeric(NA)
  result$k.wprev.pdhb.g2 <- as.numeric(NA)
  result$k.wprev.pdhb.g3 <- as.numeric(NA)
  result$k.wprev.pdhb.g4 <- as.numeric(NA)
  result$k.wprev.pdhb.g5 <- as.numeric(NA)
  result$k.wprev.pdhb.g6 <- as.numeric(NA)
  result$k.wprev.pdhb.g7 <- as.numeric(NA)
  result$k.wprev.pdhb.g8 <- as.numeric(NA)
  result$k.wprev.pdhb.g9 <- as.numeric(NA)
  result$k.wprev.pdhb.g10 <- as.numeric(NA)
  
  if(progression) message("Rings - ALL")
  object <- rings(object, N=Noptim(object), progression=progression)
  r.name <- paste0("N",Noptim(object),".RInf")
  data.all <- merge(clusters,object@rings[[r.name]]$estimates, by="id")
  
  data.groups <- list()
  for(i in 1:10) {
    if(progression) message("Rings - Group ",i)
    tmp <- as.prevR(clusters[clusters$group!=i,],col,boundary=object@boundary)
    tmp <- rings(tmp, N=Noptim(tmp), progression=progression)
    r.name <- paste0("N",Noptim(tmp),".RInf")
    data.tmp <- list(merge(tmp@clusters,tmp@rings[[r.name]]$estimates, by="id"))
    names(data.tmp) <- as.character(i)
    data.groups <- c(data.groups, data.tmp)
  }
  
  if(progression) message("KDE - ALL")
  tmp <- prev.kde(data.all, coord, weighted=T, longlat=T, progression=progression)
  for (i in 1:nrow(result)) {
    x <- which(tmp$x==clusters$x[i])[1]
    y <- which(tmp$y==clusters$y[i])[1]
    result[i,"k.wprev.full"] <- tmp$k.wprev[x,y]
  }
  
  for (g in 1:10) {
    if(progression) message("KDE - Group ",g)
    tmp <- prev.kde(data.groups[[g]], coord, weighted=T, longlat=T, progression=progression)
    for (i in 1:nrow(result)) {
      x <- which(tmp$x==clusters$x[i])[1]
      y <- which(tmp$y==clusters$y[i])[1]
      result[i,paste0("k.wprev.pdhb.g",g)] <- tmp$k.wprev[x,y]
    }
    result[result$group==g,"k.wprev.pdhb"] <- result[result$group==g,paste0("k.wprev.pdhb.g",g)]
  }
  
  return(result)
}

pdhb.plot <- function(object, groups, nb.cells=100, cell.size=NULL, progression=TRUE, plot.results=TRUE, return.results=FALSE, title=NULL){
  names(groups) <- c("id","group")
  clusters <- object@clusters
  col <- names(clusters)
  names(col) <- names(clusters)
  clusters <- merge(clusters, groups, by="id")
  
  if (progression) {
    message("Progress of calculations:",domain="R-prevR")
    barre = txtProgressBar(min=0, max=20, initial=0, style=3)
  }
  
  r <- list()
  n <- c()
  for(i in 1:10) {
    tmp <- as.prevR(clusters[clusters$group!=i,],col,boundary=object@boundary)
    tmp <- rings(tmp, N=Noptim(tmp), progression=F)
    if (progression) {setTxtProgressBar(barre,value=2*i-1)}
    tmp.r <- list(group=kde(tmp,N=Noptim(tmp),weighted=T,nb.cells=nb.cells,cell.size=cell.size,progression=F))
    names(tmp.r) <- as.character(i)
    r <- c(r, tmp.r)
    n <- c(n, Noptim(tmp))
    if (progression) {setTxtProgressBar(barre,value=2*i)}
  }
  
  if (plot.results) {
    require(ggplot2)
    data <- data.frame()
    for (i in 1:10) {
      tmp <- as.data.frame(r[[i]])
      colnames(tmp) <- c("prev","wprev","x","y")
      tmp$group <- i
      data <- rbind(data,tmp)
    }
    
    data$group <- factor(data$group)
    levels(data$group) <- paste0("G ",1:10," (N=",n,")")
    
    p <- ggplot(na.omit(data), aes(x=x, y=y, fill=wprev)) +
      geom_raster() +
      scale_fill_gradientn(colours=prevR.colors.red(20)) +
      facet_wrap(~ group) +
      coord_fixed() + xlab("") + ylab("") +
      labs(fill="%") +
      ggtitle(title) +
      theme_prevR()
    print(p)
  }
  
  if (progression) {close(barre)}
  
  if (return.results)
    return(r)
}

## On calcul le MSID (cf papier CyberGeo) pour différentes valeurs de N. 
## Pour cela, on a besoin d'un fichier avec les clusters du modèle qu'on va kriger pour avoir une surface de comparaison
## Pas générique du tout. Expérimental

multi.msid <- function(object, N, modele, nb.cells=100, cell.size=NULL, progression=TRUE) {
  # On idw les clusters du modele
  # pour generer la surface 
  locations.data = as.SpatialGrid(object, nb.cells=nb.cells, cell.size=cell.size)
  sp.modele <- modele
  coordinates(sp.modele) = ~x+y
  sp.modele@proj4string = object@proj
  modele.surface <- idw(prev~1, sp.modele, locations.data, debug.level=if (progression) -1 else 0)
  modele.surface <- as(modele.surface, "SpatialPixelsDataFrame")
  boundary     = slot(object,"boundary")
  if (attr(boundary,"valid")) {
    modele.surface = NA.outside.SpatialPolygons(modele.surface, boundary)
  }
  dimnames(modele.surface@coords)[[2]] <- c("x","y") # Fixing name of coordinates
  modele.df <- as.data.frame(modele.surface)
  
  # KDE pour object
  kde.surfaces <- kde(object, N=N, R=Inf, weighted=TRUE, nb.cells=nb.cells, cell.size=cell.size, progression=progression)
  dimnames(kde.surfaces@coords)[[2]] <- c("x","y") # Fixing name of coordinates
  kde.df <- as.data.frame(kde.surfaces)
  
  # On fusionne
  modele.df <- modele.df[order(modele.df$x,modele.df$y),]
  kde.df <- kde.df[order(kde.df$x,kde.df$y),]
  kde.df$modele <- modele.df$var1.pred
  kde.df <- na.omit(kde.df)
  
  # On calcule msid
  r <- data.frame()
  for (n in N) {
    r <- rbind(r,data.frame(N=n,msid=mean((kde.df[["modele"]]-kde.df[[paste0("k.wprev.N",n,".RInf")]])^2)))
  }
  
  return(r)
}