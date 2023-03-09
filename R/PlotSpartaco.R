#' Plot SpaRTaCo
#' This function returns the ggplots of the mean and the spatial signal-to-noise ratios from a SpaRTaCo model.
#'
#' @import ggplot2
#' @export
#'
#' @param x a `spartaco` object;
#' @param type the type of plot you want to show.
#' - `1` displays the co-cluster mean levels;
#' - `2` displays the co-cluster spatial signal-to-noise ratios;
#' - `3` displays the map of the spots colored with respect to the estimated column clusters;
#' - `4` if `!is.null(k)`, it displays the average gene expression in the clusters `k` and `r`, otherwise it displays the expression of the gene given in `gene.name`.
#' @param gene.name  (used only when `type == 4`); it receives the name of the gene to display. If `!is.null(k)`, it displays the average expression of the gene clusters given by `k`.
#' @param k (used only when `type == 4`) the gene clusters to plot.
#' @param r (used only when `type == 4`) the spot clusters to plot.
#' @param manual.palette a vector of colors used when `type == 3`.
#' @param display.all.spots if `TRUE` (default) and `type == 4`, it displays the entire grid of spots.
#' @param remove.outliers (used only when `type == 4)` if TRUE remove spots with extreme gene expression value (lower-whisker or upper-whisker)
#' @param coef (used only when `type == 4` and if `remove.outliers == TRUE`) in the selection of outliers, spots with gene expression value grater that upper-whisker value times coef are removed.
#' @param range vector to set the minimum and maximum of the scale of value to display (the default is `NULL` and the minimum and maximum of the value is used).
#'
#'
#'
#' @return The requested plot is displayed. In addition, if assigned to an object, it will return the `ggplot` object.
#'
plot.spartaco <- function(x, type = 1, gene.name = readline("gene name: "), k = NULL, r = 1:ncol(x$mu), manual.palette = NULL, display.all.spots = T, remove.outliers = F, coef = 1, range = NULL, ...){
  if(class(x) != "spartaco") stop("the input file is not a spartaco object")
  if(length(type) > 1) type <- 1
  K <- nrow(x$mu)
  R <- ncol(x$mu)
  k.lab <- 1:K
  r.lab <- 1:R
  gr <- expand.grid(k.lab,r.lab)
  gr$Mu <- as.vector(x$mu)
  gr$Ratio <- as.vector(x$tau/x$xi)
  gr <- cbind(gr, expand.grid(as.vector(table(x$Cs))/nrow(x$x),as.vector(table(x$Ds))/ncol(x$x)))
  names(gr)[-c(3,4)] <- c("Y","X","height","width")

  prop.x <- as.vector(table(x$Ds))/ncol(x$x)
  prop.y <- as.vector(table(x$Cs))/nrow(x$x)
  xlim.left <- c(0, cumsum(prop.x)[-R])
  xlim.right <- cumsum(prop.x)
  ylim.left <- c(0, cumsum(prop.y)[-K])
  ylim.right <- cumsum(prop.y)

  # ---plot mu
  if(type == 1){
    if(is.null(range)){
      range[1] <- min(as.vector(x$mu))
      range[2] <- max(as.vector(x$mu))
    }
    Plots <- ggplot(gr, aes(xmin = as.vector(sapply(1:R, function(i) rep(xlim.left[i],K))),
                            xmax = as.vector(sapply(1:R, function(i) rep(xlim.right[i],K))),
                            ymin = rep(ylim.left,R),
                            ymax = rep(ylim.right,R),
                            fill = Mu)
    )+geom_rect()+theme_bw()+
      scale_fill_distiller(palette = "RdPu",
                           limits = c(range[1], range[2]),
                           breaks = round(seq(range[1], range[2], length = 3),2))+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text=element_text(size=18),
            axis.title=element_text(size=18),
            legend.text = element_text(size = 14),
            legend.title = element_text(size = 22),
            #legend.position = "bottom",
            plot.title = element_text(hjust = 0.5),
            title = element_text(size=18),
            plot.margin=grid::unit(c(3,2,3,2), "mm"))+
      labs(fill=expression(hat(mu)[kr]))+
      scale_x_continuous(breaks=(xlim.left+xlim.right)/2,
                         labels=paste("r =",1:R)
      )+
      scale_y_continuous(breaks=(ylim.left+ylim.right)/2,
                         labels=paste("k =",1:K)
      )


  }

  # ---plot tau/xi
  if(type == 2){
    if(is.null(range)){
      range[1] <- min(as.vector(x$tau/x$xi))
      range[2] <- max(as.vector(x$tau/x$xi))
    }
    Plots <- ggplot(gr, aes(xmin = as.vector(sapply(1:R, function(i) rep(xlim.left[i],K))),
                            xmax = as.vector(sapply(1:R, function(i) rep(xlim.right[i],K))),
                            ymin = rep(ylim.left,R),
                            ymax = rep(ylim.right,R),
                            fill = Ratio)
    )+geom_rect()+theme_bw()+
      viridis::scale_fill_viridis(discrete=FALSE,
                                  limits = c(range[1], range[2]),
                                  breaks = round(seq(range[1], range[2], length = 3),2))+

      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text=element_text(size=18),
            axis.title=element_text(size=18),
            legend.text = element_text(size = 18),
            legend.title = element_text(size = 22),
            #legend.position = "bottom",
            plot.title = element_text(hjust = 0.5),
            title = element_text(size=18),
            plot.margin=grid::unit(c(3,2,3,2), "mm"))+
      labs(fill=expression(hat(tau)[kr]/hat(xi)[kr]))+
      scale_x_continuous(breaks=(xlim.left+xlim.right)/2,
                         labels=paste("r =",1:R)
      )+
      scale_y_continuous(breaks=(ylim.left+ylim.right)/2,
                         labels=paste("k =",1:K)
      )


  }

  # ---plot spot clusters
  if(type == 3){
    # ---plot column clusters
        Coord <- data.frame(x = x$coordinates[,2], y = -x$coordinates[,1], z = as.factor(x$Ds))
        Coord$group <- as.factor(as.numeric(x$Ds %in% c(1,7:9)))
        Plots <- ggplot(Coord, aes(x, y, color = z))+
          geom_point(size = 3)+theme_bw()+
          labs(col = "")+
          labs(col = expression(D[r]))+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.x = element_blank(),#element_text(size=18),
                axis.text.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title=element_blank(),#element_text(size=18),
                legend.text = element_text(size = 18),
                legend.title = element_text(size = 22),
                plot.title = element_text(hjust = 0.5),
                #legend.position = "bottom",
                #legend.spacing.x = unit(0.3, 'cm'),
                title = element_text(size=18),
                plot.margin=grid::unit(c(3,2,3,2), "mm"))+
          geom_point(shape = 1,size = 3,colour = "black")
        #if(use.greys)
        #    Plots <- Plots + scale_color_grey(start = 0, end = .9) else
        Plots <- Plots + scale_fill_brewer( palette = "Set1")
  }

  # ---plot expressions
  if(type == 4){
    if(is.null(k)){
      if(gene.name != "all" & !(gene.name %in% row.names(x$x)))
        stop(paste("Gene",gene.name,"not found\n"))
      if(gene.name == "all") k <- 1:K}
    # ---plot sample means
    Coord <- data.frame(x = x$coordinates[,2], y = -x$coordinates[,1], z = as.factor(x$Ds))
    if(is.null(r)) r <- 1:R
    if(!display.all.spots){
      if(length(k) <= 1){
        if(length(k) == 0){
          x.bar <- x$x[which(row.names(x$x) == gene.name), which(x$Ds %in% r)]
          if(is.null(range)){
            min <- round(boxplot(x.bar)$stats[1,1])
            max <- round(coef*boxplot(x.bar)$stats[5,1])
          }
          else{
            min <- range[1]
            max <- range[2]
          }
          if(remove.outliers){
            x.bar[x.bar < min] <- rep(min+1e-1, sum(x.bar < min))
            x.bar[x.bar > max] <- rep(max-1e-1, sum(x.bar > max))
          }
        }
        if(length(k) == 1){
          x.bar <- colMeans(x$x[which(x$Cs == k), which(x$Ds %in% r)])
          if(is.null(range)){
            min <- round(boxplot(x.bar)$stats[1,1])
            max <- round(coef*boxplot(x.bar)$stats[5,1])
          }
          else{
            min <- range[1]
            max <- range[2]
          }
          if(remove.outliers){
            x.bar[x.bar < min] <- rep(min+1e-1, sum(x.bar < min))
            x.bar[x.bar > max] <- rep(max-1e-1, sum(x.bar > max))
          }
        }
        Plots <- ggplot(Coord[which(x$Ds %in% r),], aes(x, y, color = x.bar))+
          geom_point(size = 3)+theme_bw()+
          scale_fill_distiller(type = "seq", palette = "Spectral", direction = -1,
                               limits = c(min,max), breaks = round(seq(min,max, length = 4),2))+
          labs(col = "")+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.x = element_blank(),#element_text(size=18),
                axis.text.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title=element_blank(),#element_text(size=18),
                legend.text = element_text(size = 18),
                legend.title = element_text(size = 22),
                plot.title = element_text(hjust = 0.5),
                #legend.position = "bottom",
                #legend.spacing.x = unit(0.3, 'cm'),
                title = element_text(size=18),
                plot.margin=grid::unit(c(3,2,3,2), "mm"))+
          geom_point(shape = 1,size = 3,colour = "black")+
          ggtitle(label = ifelse(length(k) == 0, gene.name, paste("k =",k)))
      } else {
        Plots <- list()
        for(k.ind in 1:length(k)){
          Plots[[k.ind]] <- local({
            x.bar <- colMeans(x$x[x$Cs == k[k.ind], which(x$Ds %in% r)])
            if(is.null(range)){
              min <- round(boxplot(x.bar)$stats[1,1])
              max <- round(coef*boxplot(x.bar)$stats[5,1])
            }
            else{
              min <- range[1]
              max <- range[2]
            }
            if(remove.outliers){
              x.bar[x.bar < min] <- rep(min+1e-1, sum(x.bar < min))
              x.bar[x.bar > max] <- rep(max-1e-1, sum(x.bar > max))
            }
            ggplot(Coord[which(x$Ds %in% r),], aes(x, y, color = x.bar))+
              geom_point(size = 3)+theme_bw()+
              scale_color_distiller(type = "seq", palette = "Spectral", direction = -1,
                                   limits = c(min,max), breaks = round(seq(min,max, length = 4),2))+
              labs(col = "")+
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.text.x = element_blank(),#element_text(size=18),
                    axis.text.y = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.title=element_blank(),#element_text(size=18),
                    legend.text = element_text(size = 18),
                    legend.title = element_text(size = 22),
                    plot.title = element_text(hjust = 0.5),
                    title = element_text(size=18),
                    plot.margin=grid::unit(c(3,2,3,2), "mm"))+
              geom_point(shape = 1,size = 3,colour = "black")+
              ggtitle(label = paste("k =",k[k.ind]))})
        }
      }
    } else {  # if display.all.spots == T
      if(length(k) <= 1){
        if(length(k) == 0){
          x.bar <- x$x[which(row.names(x$x) == gene.name), which(x$Ds %in% r)]
          if(is.null(range)){
            min <- round(boxplot(x.bar)$stats[1,1])
            max <- round(coef*boxplot(x.bar)$stats[5,1])
          }
          else{
            min <- range[1]
            max <- range[2]
          }
          if(remove.outliers){
            x.bar[x.bar < min] <- rep(min+1e-1, sum(x.bar < min))
            x.bar[x.bar > max] <- rep(max-1e-1, sum(x.bar > max))
          }
        }
        if(length(k) == 1){
          x.bar <- colMeans(x$x[which(x$Cs == k), which(x$Ds %in% r)])
          if(is.null(range)){
            min <- round(boxplot(x.bar)$stats[1,1])
            max <- round(coef*boxplot(x.bar)$stats[5,1])
          }
          else{
            min <- range[1]
            max <- range[2]
          }
          if(remove.outliers){
            x.bar[x.bar < min] <- rep(min+1e-1, sum(x.bar < min))
            x.bar[x.bar > max] <- rep(max-1e-1, sum(x.bar > max))
          }
        }
        Plots <- ggplot(Coord[-which(x$Ds %in% r),], aes(x, y))+
          theme_bw()+
          geom_point(data = Coord[-which(x$Ds %in% r),], mapping = aes(x, y), size = 3, fill = "white", colour = "gray74", shape = 21)+
          geom_point(data = Coord[which(x$Ds %in% r),], mapping = aes(x, y, fill = x.bar), colour = "white", size = 4, shape = 21)+
          scale_fill_distiller(type = "seq", palette = "Spectral", direction = -1,
                               limits = c(min,max), breaks = round(seq(min,max, length = 4),2))+
          labs(fill = "")+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.x = element_blank(),#element_text(size=18),
                axis.text.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title=element_blank(),#element_text(size=18),
                legend.text = element_text(size = 18),
                legend.title = element_text(size = 22),
                plot.title = element_text(hjust = 0.5),
                #legend.position = "bottom",
                #legend.spacing.x = unit(0.3, 'cm'),
                title = element_text(size=18),
                plot.margin=grid::unit(c(3,2,3,2), "mm"))+
          #geom_point(shape = 1,size = 3,colour = "black")+
          ggtitle(label = ifelse(length(k) == 1, paste("k =",k), gene.name))
      } else {
        Plots <- list()
        for(k.ind in 1:length(k)){
          Plots[[k.ind]] <- local({
            x.bar <- colMeans(x$x[x$Cs == k[k.ind], which(x$Ds %in% r)])
            if(is.null(range)){
              min <- round(boxplot(x.bar)$stats[1,1])
              max <- round(coef*boxplot(x.bar)$stats[5,1])
            }
            else{
              min <- range[1]
              max <- range[2]
            }
            if(remove.outliers){
              x.bar[x.bar < min] <- rep(min+1e-1, sum(x.bar < min))
              x.bar[x.bar > max] <- rep(max-1e-1, sum(x.bar > max))
            }
            ggplot(Coord[-which(x$Ds %in% r),], aes(x, y))+
              theme_bw()+
              geom_point(data = Coord[-which(x$Ds %in% r),], mapping = aes(x, y), size = 3, fill = "white", colour = "gray74", shape = 21)+
              geom_point(data = Coord[which(x$Ds %in% r),], mapping = aes(x, y, fill = x.bar), colour = "white", size = 4, shape = 21)+
              scale_fill_distiller(type = "seq", palette = "Spectral", direction = -1,
                                   limits = c(min,max), breaks = round(seq(min,max, length = 4),2))+
              labs(col = "")+
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.text.x = element_blank(),#element_text(size=18),
                    axis.text.y = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.title=element_blank(),#element_text(size=18),
                    legend.text = element_text(size = 18),
                    legend.title = element_text(size = 22),
                    plot.title = element_text(hjust = 0.5),
                    title = element_text(size=18),
                    plot.margin=grid::unit(c(3,2,3,2), "mm"))+
              #geom_point(shape = 1,size = 3,colour = "black")+
              ggtitle(label = paste("k =",k[k.ind]))})
        }
      }
    }
  }

  Plots
}
