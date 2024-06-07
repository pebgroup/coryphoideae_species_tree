plotwe <- function (x, type = "phylogram", use.edge.length = TRUE, node.pos = NULL, 
          show.tip.label = TRUE, show.node.label = FALSE, edge.color = "black", 
          edge.width = 1, edge.lty = 1, font = 3, cex = par("cex"), 
          adj = NULL, srt = 0, no.margin = FALSE, root.edge = FALSE, 
          label.offset = 0, underscore = FALSE, x.lim = NULL, y.lim = NULL, 
          direction = "rightwards", lab4ut = NULL, tip.color = "black", 
          plot = TRUE, rotate.tree = 0, open.angle = 0, node.depth = 1, 
          align.tip.label = FALSE, link.color = "black", ...) 
{
  Ntip <- length(x$tip.label)
  if (Ntip < 2) {
    warning("found less than 2 tips in the tree")
    return(NULL)
  }
  .nodeHeight <- function(edge, Nedge, yy) .C(node_height, 
                                              as.integer(edge[, 1]), as.integer(edge[, 2]), as.integer(Nedge), 
                                              as.double(yy))[[4]]
  .nodeDepth <- function(Ntip, Nnode, edge, Nedge, node.depth) .C(node_depth, 
                                                                  as.integer(Ntip), as.integer(edge[, 1]), as.integer(edge[, 
                                                                                                                           2]), as.integer(Nedge), double(Ntip + Nnode), as.integer(node.depth))[[5]]
  .nodeDepthEdgelength <- function(Ntip, Nnode, edge, Nedge, 
                                   edge.length) .C(node_depth_edgelength, as.integer(edge[, 
                                                                                          1]), as.integer(edge[, 2]), as.integer(Nedge), as.double(edge.length), 
                                                   double(Ntip + Nnode))[[5]]
  Nedge <- dim(x$edge)[1]
  Nnode <- x$Nnode
  if (any(x$edge < 1) || any(x$edge > Ntip + Nnode)) 
    stop("tree badly conformed; cannot plot. Check the edge matrix.")
  ROOT <- Ntip + 1
  type <- match.arg(type, c("phylogram", "cladogram", "fan", 
                            "unrooted", "radial"))
  direction <- match.arg(direction, c("rightwards", "leftwards", 
                                      "upwards", "downwards"))
  if (is.null(x$edge.length)) {
    use.edge.length <- FALSE
  }
  else {
    if (use.edge.length && type != "radial") {
      tmp <- sum(is.na(x$edge.length))
      if (tmp) {
        warning(paste(tmp, "branch length(s) NA(s): branch lengths ignored in the plot"))
        use.edge.length <- FALSE
      }
    }
  }
  if (is.numeric(align.tip.label)) {
    align.tip.label.lty <- align.tip.label
    align.tip.label <- TRUE
  }
  else {
    if (align.tip.label) 
      align.tip.label.lty <- 3
  }
  if (align.tip.label) {
    if (type %in% c("unrooted", "radial") || !use.edge.length || 
        is.ultrametric(x)) 
      align.tip.label <- FALSE
  }
  if (type %in% c("unrooted", "radial") || !use.edge.length || 
      is.null(x$root.edge) || !x$root.edge) 
    root.edge <- FALSE
  phyloORclado <- type %in% c("phylogram", "cladogram")
  horizontal <- direction %in% c("rightwards", "leftwards")
  xe <- x$edge
  if (phyloORclado) {
    phyOrder <- attr(x, "order")
    if (is.null(phyOrder) || phyOrder != "cladewise") {
      x <- reorder(x)
      if (!identical(x$edge, xe)) {
        ereorder <- match(x$edge[, 2], xe[, 2])
        if (length(edge.color) > 1) {
          edge.color <- rep(edge.color, length.out = Nedge)
          edge.color <- edge.color[ereorder]
        }
        if (length(edge.width) > 1) {
          edge.width <- rep(edge.width, length.out = Nedge)
          edge.width <- edge.width[ereorder]
        }
        if (length(edge.lty) > 1) {
          edge.lty <- rep(edge.lty, length.out = Nedge)
          edge.lty <- edge.lty[ereorder]
        }
      }
    }
    yy <- numeric(Ntip + Nnode)
    TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
    yy[TIPS] <- 1:Ntip
  }
  z <- reorder(x, order = "postorder")
  if (phyloORclado) {
    if (is.null(node.pos)) 
      node.pos <- if (type == "cladogram" && !use.edge.length) 
        2
    else 1
    if (node.pos == 1) 
      yy <- .nodeHeight(z$edge, Nedge, yy)
    else {
      ans <- .C(node_height_clado, as.integer(Ntip), as.integer(z$edge[, 
                                                                       1]), as.integer(z$edge[, 2]), as.integer(Nedge), 
                double(Ntip + Nnode), as.double(yy))
      xx <- ans[[5]] - 1
      yy <- ans[[6]]
    }
    if (!use.edge.length) {
      if (node.pos != 2) 
        xx <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, 
                         node.depth) - 1
      xx <- max(xx) - xx
    }
    else {
      xx <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge, 
                                 z$edge.length)
    }
  }
  else {
    twopi <- 2 * pi
    rotate.tree <- twopi * rotate.tree/360
    if (type != "unrooted") {
      TIPS <- x$edge[which(x$edge[, 2] <= Ntip), 2]
      xx <- seq(0, twopi * (1 - 1/Ntip) - twopi * open.angle/360, 
                length.out = Ntip)
      theta <- double(Ntip)
      theta[TIPS] <- xx
      theta <- c(theta, numeric(Nnode))
    }
    switch(type, fan = {
      theta <- .nodeHeight(z$edge, Nedge, theta)
      if (use.edge.length) {
        r <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, 
                                  Nedge, z$edge.length)
      } else {
        r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
        r <- 1/r
      }
      theta <- theta + rotate.tree
      if (root.edge) r <- r + x$root.edge
      xx <- r * cos(theta)
      yy <- r * sin(theta)
    }, unrooted = {
      nb.sp <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
      XY <- if (use.edge.length) unrooted.xy(Ntip, Nnode, 
                                             z$edge, z$edge.length, nb.sp, rotate.tree) else unrooted.xy(Ntip, 
                                                                                                         Nnode, z$edge, rep(1, Nedge), nb.sp, rotate.tree)
      xx <- XY$M[, 1] - min(XY$M[, 1])
      yy <- XY$M[, 2] - min(XY$M[, 2])
    }, radial = {
      r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
      r[r == 1] <- 0
      r <- 1 - r/Ntip
      theta <- .nodeHeight(z$edge, Nedge, theta) + rotate.tree
      xx <- r * cos(theta)
      yy <- r * sin(theta)
    })
  }
  if (phyloORclado) {
    if (!horizontal) {
      tmp <- yy
      yy <- xx
      xx <- tmp - min(tmp) + 1
    }
    if (root.edge) {
      if (direction == "rightwards") 
        xx <- xx + x$root.edge
      if (direction == "upwards") 
        yy <- yy + x$root.edge
    }
  }
  if (no.margin) 
    par(mai = rep(0, 4))
  if (show.tip.label) 
    nchar.tip.label <- nchar(x$tip.label)
  max.yy <- max(yy)
  getLimit <- function(x, lab, sin, cex) {
    s <- strwidth(lab, "inches", cex = cex)
    if (any(s > sin)) 
      return(1.5 * max(x))
    Limit <- 0
    while (any(x > Limit)) {
      i <- which.max(x)
      alp <- x[i]/(sin - s[i])
      Limit <- x[i] + alp * s[i]
      x <- x + alp * s
    }
    Limit
  }
  if (is.null(x.lim)) {
    if (phyloORclado) {
      if (horizontal) {
        xx.tips <- xx[1:Ntip]
        if (show.tip.label) {
          pin1 <- par("pin")[1]
          tmp <- getLimit(xx.tips, x$tip.label, pin1, 
                          cex)
          tmp <- tmp + label.offset
        }
        else tmp <- max(xx.tips)
        x.lim <- c(0, tmp)
      }
      else x.lim <- c(1, Ntip)
    }
    else switch(type, fan = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
                        cex)
        x.lim <- range(xx) + c(-offset, offset)
      } else x.lim <- range(xx)
    }, unrooted = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
                        cex)
        x.lim <- c(0 - offset, max(xx) + offset)
      } else x.lim <- c(0, max(xx))
    }, radial = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.03 * cex)
        x.lim <- c(-1 - offset, 1 + offset)
      } else x.lim <- c(-1, 1)
    })
  }
  else if (length(x.lim) == 1) {
    x.lim <- c(0, x.lim)
    if (phyloORclado && !horizontal) 
      x.lim[1] <- 1
    if (type %in% c("fan", "unrooted") && show.tip.label) 
      x.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy * 
                         cex)
    if (type == "radial") 
      x.lim[1] <- if (show.tip.label) 
        -1 - max(nchar.tip.label * 0.03 * cex)
    else -1
  }
  if (phyloORclado && direction == "leftwards") 
    xx <- x.lim[2] - xx
  if (is.null(y.lim)) {
    if (phyloORclado) {
      if (horizontal) 
        y.lim <- c(1, Ntip)
      else {
        pin2 <- par("pin")[2]
        yy.tips <- yy[1:Ntip]
        if (show.tip.label) {
          tmp <- getLimit(yy.tips, x$tip.label, pin2, 
                          cex)
          tmp <- tmp + label.offset
        }
        else tmp <- max(yy.tips)
        y.lim <- c(0, tmp)
      }
    }
    else switch(type, fan = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
                        cex)
        y.lim <- c(min(yy) - offset, max.yy + offset)
      } else y.lim <- c(min(yy), max.yy)
    }, unrooted = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
                        cex)
        y.lim <- c(0 - offset, max.yy + offset)
      } else y.lim <- c(0, max.yy)
    }, radial = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.03 * cex)
        y.lim <- c(-1 - offset, 1 + offset)
      } else y.lim <- c(-1, 1)
    })
  }
  else if (length(y.lim) == 1) {
    y.lim <- c(0, y.lim)
    if (phyloORclado && horizontal) 
      y.lim[1] <- 1
    if (type %in% c("fan", "unrooted") && show.tip.label) 
      y.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy * 
                         cex)
    if (type == "radial") 
      y.lim[1] <- if (show.tip.label) 
        -1 - max(nchar.tip.label * 0.018 * max.yy * cex)
    else -1
  }
  if (phyloORclado && direction == "downwards") 
    yy <- y.lim[2] - yy
  if (phyloORclado && root.edge) {
    if (direction == "leftwards") 
      x.lim[2] <- x.lim[2] + x$root.edge
    if (direction == "downwards") 
      y.lim[2] <- y.lim[2] + x$root.edge
  }
  asp <- if (type %in% c("fan", "radial", "unrooted")) 
    1
  else NA
  plot.default(0, type = "n", xlim = x.lim, ylim = y.lim, xlab = "", 
               ylab = "", axes = FALSE, asp = asp, ...)
  if (plot) {
    if (is.null(adj)) 
      adj <- if (phyloORclado && direction == "leftwards") 
        1
    else 0
    if (phyloORclado && show.tip.label) {
      MAXSTRING <- max(strwidth(x$tip.label, cex = cex))
      loy <- 0
      if (direction == "rightwards") {
        lox <- label.offset + MAXSTRING * 1.05 * adj
      }
      if (direction == "leftwards") {
        lox <- -label.offset - MAXSTRING * 1.05 * (1 - 
                                                     adj)
      }
      if (!horizontal) {
        psr <- par("usr")
        MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3])/(psr[2] - 
                                                             psr[1])
        loy <- label.offset + MAXSTRING * 1.05 * adj
        lox <- 0
        srt <- 90 + srt
        if (direction == "downwards") {
          loy <- -loy
          srt <- 180 + srt
        }
      }
    }
    if (type == "phylogram") {
      phylogram.plot(x$edge, Ntip, Nnode, xx, yy, horizontal, 
                     edge.color, edge.width, edge.lty)
    }
    else {
      if (type == "fan") {
        ereorder <- match(z$edge[, 2], x$edge[, 2])
        if (length(edge.color) > 1) {
          edge.color <- rep(edge.color, length.out = Nedge)
          edge.color <- edge.color[ereorder]
        }
        if (length(edge.width) > 1) {
          edge.width <- rep(edge.width, length.out = Nedge)
          edge.width <- edge.width[ereorder]
        }
        if (length(edge.lty) > 1) {
          edge.lty <- rep(edge.lty, length.out = Nedge)
          edge.lty <- edge.lty[ereorder]
        }
        circular.plot(z$edge, Ntip, Nnode, xx, yy, theta, 
                      r, edge.color, edge.width, edge.lty)
      }
      else cladogram.plot(x$edge, xx, yy, edge.color, edge.width, 
                          edge.lty)
    }
    if (root.edge) {
      rootcol <- if (length(edge.color) == 1) 
        edge.color
      else "black"
      rootw <- if (length(edge.width) == 1) 
        edge.width
      else 1
      rootlty <- if (length(edge.lty) == 1) 
        edge.lty
      else 1
      if (type == "fan") {
        tmp <- polar2rect(x$root.edge, theta[ROOT])
        segments(0, 0, tmp$x, tmp$y, col = rootcol, lwd = rootw, 
                 lty = rootlty)
      }
      else {
        switch(direction, rightwards = segments(0, yy[ROOT], 
                                                x$root.edge, yy[ROOT], col = rootcol, lwd = rootw, 
                                                lty = rootlty), leftwards = segments(xx[ROOT], 
                                                                                     yy[ROOT], xx[ROOT] + x$root.edge, yy[ROOT], 
                                                                                     col = rootcol, lwd = rootw, lty = rootlty), 
               upwards = segments(xx[ROOT], 0, xx[ROOT], x$root.edge, 
                                  col = rootcol, lwd = rootw, lty = rootlty), 
               downwards = segments(xx[ROOT], yy[ROOT], xx[ROOT], 
                                    yy[ROOT] + x$root.edge, col = rootcol, lwd = rootw, 
                                    lty = rootlty))
      }
    }
    if (show.tip.label) {
      if (is.expression(x$tip.label)) 
        underscore <- TRUE
      if (!underscore) 
        x$tip.label <- gsub("_", " ", x$tip.label)
      if (phyloORclado) {
        if (align.tip.label) {
          xx.tmp <- switch(direction, rightwards = max(xx[1:Ntip]), 
                           leftwards = min(xx[1:Ntip]), upwards = xx[1:Ntip], 
                           downwards = xx[1:Ntip])
          yy.tmp <- switch(direction, rightwards = yy[1:Ntip], 
                           leftwards = yy[1:Ntip], upwards = max(yy[1:Ntip]), 
                           downwards = min(yy[1:Ntip]))
          segments(xx[1:Ntip], yy[1:Ntip], xx.tmp, yy.tmp, 
                   lty = align.tip.label.lty, col=link.color)
        }
        else {
          xx.tmp <- xx[1:Ntip]
          yy.tmp <- yy[1:Ntip]
        }
        text(xx.tmp + lox, yy.tmp + loy, x$tip.label, 
             adj = adj, font = font, srt = srt, cex = cex, 
             col = tip.color)
      }
      else {
        angle <- if (type == "unrooted") 
          XY$axe
        else atan2(yy[1:Ntip], xx[1:Ntip])
        lab4ut <- if (is.null(lab4ut)) {
          if (type == "unrooted") 
            "horizontal"
          else "axial"
        }
        else match.arg(lab4ut, c("horizontal", "axial"))
        xx.tips <- xx[1:Ntip]
        yy.tips <- yy[1:Ntip]
        if (label.offset) {
          xx.tips <- xx.tips + label.offset * cos(angle)
          yy.tips <- yy.tips + label.offset * sin(angle)
        }
        if (lab4ut == "horizontal") {
          y.adj <- x.adj <- numeric(Ntip)
          sel <- abs(angle) > 0.75 * pi
          x.adj[sel] <- -strwidth(x$tip.label)[sel] * 
            1.05
          sel <- abs(angle) > pi/4 & abs(angle) < 0.75 * 
            pi
          x.adj[sel] <- -strwidth(x$tip.label)[sel] * 
            (2 * abs(angle)[sel]/pi - 0.5)
          sel <- angle > pi/4 & angle < 0.75 * pi
          y.adj[sel] <- strheight(x$tip.label)[sel]/2
          sel <- angle < -pi/4 & angle > -0.75 * pi
          y.adj[sel] <- -strheight(x$tip.label)[sel] * 
            0.75
          text(xx.tips + x.adj * cex, yy.tips + y.adj * 
                 cex, x$tip.label, adj = c(adj, 0), font = font, 
               srt = srt, cex = cex, col = tip.color)
        }
        else {
          if (align.tip.label) {
            POL <- rect2polar(xx.tips, yy.tips)
            POL$r[] <- max(POL$r)
            REC <- polar2rect(POL$r, POL$angle)
            xx.tips <- REC$x
            yy.tips <- REC$y
            segments(xx[1:Ntip], yy[1:Ntip], xx.tips, 
                     yy.tips, lty = align.tip.label.lty)
          }
          if (type == "unrooted") {
            adj <- abs(angle) > pi/2
            angle <- angle * 180/pi
            angle[adj] <- angle[adj] - 180
            adj <- as.numeric(adj)
          }
          else {
            s <- xx.tips < 0
            angle <- angle * 180/pi
            angle[s] <- angle[s] + 180
            adj <- as.numeric(s)
          }
          font <- rep(font, length.out = Ntip)
          tip.color <- rep(tip.color, length.out = Ntip)
          cex <- rep(cex, length.out = Ntip)
          for (i in 1:Ntip) text(xx.tips[i], yy.tips[i], 
                                 x$tip.label[i], font = font[i], cex = cex[i], 
                                 srt = angle[i], adj = adj[i], col = tip.color[i])
        }
      }
    }
    if (show.node.label) 
      text(xx[ROOT:length(xx)] + label.offset, yy[ROOT:length(yy)], 
           x$node.label, adj = adj, font = font, srt = srt, 
           cex = cex)
  }
  L <- list(type = type, use.edge.length = use.edge.length, 
            node.pos = node.pos, node.depth = node.depth, show.tip.label = show.tip.label, 
            show.node.label = show.node.label, font = font, cex = cex, 
            adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset, 
            x.lim = x.lim, y.lim = y.lim, direction = direction, 
            tip.color = tip.color, Ntip = Ntip, Nnode = Nnode, root.time = x$root.time, 
            align.tip.label = align.tip.label)
  assign("last_plot.phylo", c(L, list(edge = xe, xx = xx, yy = yy)), 
         envir = .PlotPhyloEnv)
  invisible(L)
}

#########################################################################################################
#########################################################################################################
#########################################################################################################
process_phylogeny <- function(MRCAnode, phylogeny, phylogeny_en, tipcols_orthologs, gene_counts, colors_scheme, jenks_value) {
  
  # Extract tree using MRCA node
  astral_tree <- extract.clade(phylogeny, MRCAnode)
  astral_tree_EN <- extract.clade(phylogeny_en, MRCAnode)
  
  # Get support values as edge labels
  astral_tree$node.label -> supportvals
  names(supportvals) <- 1:length(astral_tree$node.label) + Ntip(astral_tree)
  
  # Translate node labels into edge labels
  edgelabs <- c()
  for(i in 1:nrow(astral_tree$edge)){ # loops through all edges
    dscndnt <- as.character(astral_tree$edge[i,2]) # selects the crown node of each edge
    if(dscndnt %in% names(supportvals)){ # if this is an internal branch.... 
      edgelabs <- c(edgelabs, supportvals[dscndnt]) # assign the support value of the crown node
    } else {
      edgelabs <- c(edgelabs, "") # if not, assign nothing. 
    }
  }
  
  edgelabs[edgelabs=="1"] <- "" # suppress support values == 1
  
  # Get edge colours
  EN <- as.numeric(astral_tree$node.label)
  names(EN) <- 1:length(astral_tree$node.label) + Ntip(astral_tree)

  # Create intervals for the number of genes
  nodeclass <- classIntervals(EN[!is.na(EN)], jenks_value, style = "jenks") 
  pal <- c(met.brewer(colors_scheme, jenks_value, type = "discrete"))
  as.vector(findColours(nodeclass, pal)) -> ENcolor
  
  names(ENcolor) <- names(EN[!is.na(EN)]) 
  edgecols <- rep("black", nrow(astral_tree$edge))
  
  # Loop through internal nodes and colour branches
  for(nl in (1:length(astral_tree$node.label) + Ntip(astral_tree))){
    if(nl %in% names(ENcolor)){
      edgecols[astral_tree$edge[,2] == nl] <- ENcolor[as.character(nl)]
    }
  }
  
  # Defining tip colours
  tipres <- table(gene_counts$tipname) 
  tipclass <- nodeclass
  tipclass$var <- tipres
  tipclass$brks[jenks_value + 1] <- max(tipres)
  as.vector(findColours(tipclass, pal)) -> tipcols
  names(tipcols) <- names(tipres)
  print(names(tipcols))
  
  
  # Find tip colors by comparing tip labels
  tip_labels <- astral_tree$tip.label
  tipcols <- rep("black", length(tip_labels))
  for (i in 1:length(tip_labels)) {
    tip_label <- tip_labels[i]
    if (tip_label %in% names(tipcols_orthologs)) {
      tipcols[i] <- tipcols_orthologs[[tip_label]]
      names(tipcols)[i] <- tip_label
    }
  }

  # Loop through terminal nodes and colour branches
  print(names(tipcols))
  i = 1
  for (j in astral_tree_EN$tip.label) {
  cat("Is ", j, " in names(tipcols): " , j %in% names(tipcols), "\n")
    if (j %in% names(tipcols)) {
      edgecols[astral_tree_EN$edge[,2] == i] <- tipcols[j]
      i = i+1
    }
  }

  # Return the modified tree, edgecols, edgelabs, and tipcols
  return(list(astral_tree, edgecols, edgelabs, tipcols))
}

#######################################################################################################
#######################################################################################################
#######################################################################################################

split.plotTree <- function(tree, splits = NULL, file = NULL, fn = NULL, edge.color = "black", tip.color = "black", ...) {
  ef <- 0.037037037037
  if (!is.null(file)) pdf(file, width = 8.5, height = 11)
  if (is.null(splits)) splits <- (floor(0.5 * Ntip(tree)) + 0.5) / Ntip(tree)
  S <- matrix(c(0, splits, splits, 1 + 1/Ntip(tree)), length(splits) + 1, 2)
  S <- cbind(S[, 1] + ef * (S[, 2] - S[, 1]), S[, 2] - ef * (S[, 2] - S[, 1]))
  for (i in nrow(S):1) {
    if (is.null(file) && i < nrow(S)) par(ask = TRUE)
    plotwe(tree, plot = FALSE, edge.color = edge.color, tip.color = tip.color, ...)
    plotwe(tree, plot = TRUE, start = (S[i, 1] * Ntip(tree) + 1), end = (S[i, 2] * Ntip(tree)), edge.color = edge.color, tip.color = tip.color, ...)
    if (!is.null(fn)) fn()  # Apply custom function for additional annotations or adjustments
  }
  if (!is.null(file)) dev.off()
}

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

# split.plotTree <- function(tree, splits = NULL, file = NULL, cex = 1.0, x_lim = NULL, tipcols = NULL,clades = NULL){
# 	ef <- 0.037037037037 # percentage to remove plottingarea which causes problems when printing
# 	if (!is.null(file)) pdf(file, width = 8.5, height = 11) # Checking if a file is specified
# 	if (is.null(splits)) splits <- (floor(0.5 * Ntip(tree)) + 0.5) / Ntip(tree) # Checking if the number of splits is specified
# 	S <- matrix(c(0, splits, splits, 1 + 1/Ntip(tree)), length(splits) + 1, 2) # Creating the splits
# 	S <- cbind(S[, 1] + ef * (S[, 2] - S[, 1]), S[, 2] - ef * (S[, 2] - S[, 1])) # Adjusting the splits based on paper edges
# 	for (i in nrow(S):1) { # Looping through the splits
# 		if (is.null(file) && i < nrow(S)) par(ask = TRUE) # Checking if a file is specified
# 		par(fg = "transparent") # Making the tips transparent
# 		plotSimmap(tree, ylim = Ntip(tree) * S[i, ], xlim = x_lim, split.vertical = TRUE, ...) # Plotting the tree
# 		obj <- get("last_plot.phylo", envir = .PlotPhyloEnv) # Get the phylogeny object
# 		par(fg = "black") # Making the Tips Black
# 		text(rep(max(obj$xx[1:Ntip(tree)]), Ntip(tree)), obj$yy[1:Ntip(tree)],labels = tree$tip.label, font = 3, pos = 4, cex = cex) # Add the tip labels
# 		for (i in 1:Ntip(tree)) { # looping through the tips in the tree to add the dots
# 			colour_name <- names(test_tree$tip.label[i]) # Getting the tip name to get the colour
# 			colour_tip <- tipcols[[colour_name]] # Getting the colour of the dots
# 			#cat("Colour: ", colour_tip, " Tip: ", tree$tip.label[i] ,"\n") # writing the tip and the colour
# 			lines(c(obj$xx[i], max(obj$xx[1:Ntip(tree)])), rep(obj$yy[i], 2), lty = "dashed", col = colour_tip) # Add the tip lines
# 		}
# 		edgelabels(text = edgelabs_orthologs, frame = "none", cex = 0.5, adj = c(-0.5, -0.5)) # adding the edge labels
# 	}
# 	if (!is.null(file)) oo <- dev.off()
# }




# # for(i in 1:Ntip(test_tree)){
# # 	colour_name <- names(test_tree$tip.label[i])
# # 	colour_tip <- tipcols_orthologs[[colour_name]]
# # 	lines(c(obj$xx[i],max(obj$xx[1:Ntip(tree)])),rep(obj$yy[i],2),lty="dotted",colour_tip )
# # }

# # tipcols_orthologs

# # names <- test_tree$tip.label
# # names
# # write.csv2(names, "names_in_tree_05_09_2024.csv")



# ##############################################################################################################################

# # Abandoned attempt on adding clade labels using R

# #Sabaleae
#   #Sabal 
#   # Sabal Yapa Sabal Etonio
# MRCA_sabal <- MRCA(test_tree, c("Sabal yapa", "Sabal etonia")) #502

# # Cryosophileae
#   # Schippia, Itaya, Cryosophila, Trithrinax, Chelyocarpus, Zombia, Coccothrinax, Thrinax, Hemithrinax, Leucothrinax
#   # Cryosophila nana, Coccothrinax macroglossa
# MRCA_cryosophileae <- MRCA(test_tree, c("Sabinaria magnifica", "Coccothrinax macroglossa")) #495

# #Phoeniceae
#   # Phoenix
#   # Phoenix rupicola, Phoenix dactylifera
# MRCA_phoeniceae <- MRCA(test_tree, c("Phoenix rupicola", "Phoenix dactylifera")) #447

# # Brahea
#   # Brahea
#   # Brahea aculeata , Brahea calcarea
# MRCA_brahea <- MRCA(test_tree,c("Brahea aculeata","Brahea calcarea"))

# #Trachycarpeae
# # Subtribe Rhapidinae
#   # Chamaerops, Guihaia, Trachycarpus, Rhapidophyllum, Maxburretia, Rhapis, (Brahea, Colpothrinax)
#   # Colpothrinax aphanopetala, Rhapis puhuongensis
# MRCA_rhapidinae <- MRCA(test_tree, c("Rhapidophyllum hystrix", "Rhapis humilis")) #636

# #Subtribe livistoninae 
#   #Livistona, Licuala, johannesteijsmannia, Pholidocarpus, Saribus, Acoelorrhaphe, Serenoa
#   # Serenoa repens, Licuala simplex
# MRCA_livistoninae <-MRCA(test_tree, c("Livistona carinensis", "Licuala simplex")) #636

# # Chuniophoeniceae
#   # Chuniophoenix, Kerriodoxa, Nannorrhops, Tahina
#   # Nannorrhops ritchieana, chuniophoenix hainanensis
# MRCA_chuniophoeniceae <- MRCA(test_tree, c("Nannorrhops ritchieana", "Chuniophoenix hainanensis")) #409

# # Caryoteae
#   # Caryota, Wallichia, Arenga
#   # Caryota obtusa, Wallichia gracilis
# MRCA_caryoteae <- MRCA(test_tree, c("Caryota obtusa", "Wallichia gracilis")) #436

# #Corypheae
#   # Corypha
#   # Corypha lecomtei, Corypha taliera
# MRCA_corypheae <- MRCA(test_tree, c("Corypha lecomtei 1", "Corypha taliera")) #415

# # Borasseae
#   # Subtribe Hyphaeninae
#   # Bismarckia, Satranala, Hyphaene, Medemia
#   # Bismarkia nobilis, Hyphaene thebaica
# MRCA_borasseae <- MRCA(test_tree, c("Bismarckia nobilis", "Hyphaene thebaica")) #419

#   # Subtribe Lataniinae
#   # Borassodendron, Latania, Borassus, Lodoicea
#   # Borassus aethiopum, Latania lontaroides
# MRCA_lataniinae <- MRCA(test_tree, c("Borassus aethiopum", "Latania lontaroides")) #427



# dat_subfam <- data.frame(
#            node = c(MRCA_sabal,MRCA_cryosophileae,MRCA_phoeniceae,MRCA_rhapidinae,MRCA_livistoninae,MRCA_chuniophoeniceae,MRCA_caryoteae,MRCA_corypheae,MRCA_borasseae,MRCA_lataniinae),
#            name = c("Sabaleae","Cryosophileae","Phoeniceae","Rhapidineae","Livistoninae","Chuniophoeniceae","Caryoteae","Corypheae", "Hyphaeninae","Lataniinae")
#        )

split.plotTree_clades <- function(tree, splits = NULL, file = NULL, cex = 1.0, x_lim = NULL, tipcols = NULL, clades = NULL, ...) {
    ef <- 0.037037037037 # percentage to remove plotting area which causes problems when printing
    if (!is.null(file)) pdf(file, width = 8.5, height = 11) # Checking if a file is specified
    if (is.null(splits)) splits <- (floor(0.5 * Ntip(tree)) + 0.5) / Ntip(tree) # Checking if the number of splits is specified
    S <- matrix(c(0, splits, splits, 1 + 1/Ntip(tree)), length(splits) + 1, 2) # Creating the splits
    S <- cbind(S[, 1] + ef * (S[, 2] - S[, 1]), S[, 2] - ef * (S[, 2] - S[, 1])) # Adjusting the splits based on paper edges
    
    for (i in nrow(S):1) { # Looping through the splits
        if (is.null(file) && i < nrow(S)) par(ask = TRUE) # Checking if a file is specified
        par(fg = "transparent") # Making the tips transparent
        plotSimmap(tree, ylim = Ntip(tree) * S[i, ], xlim = x_lim, split.vertical = TRUE, ...) # Plotting the tree
        obj <- get("last_plot.phylo", envir = .PlotPhyloEnv) # Get the phylogeny object
        par(fg = "black") # Making the Tips Black
        text(rep(max(obj$xx[1:Ntip(tree)]), Ntip(tree)), obj$yy[1:Ntip(tree)], labels = tree$tip.label, font = 3, pos = 4, cex = cex) # Add the tip labels

        for (j in 1:Ntip(tree)) { # Looping through the tips in the tree to add the dots
            colour_name <- names(tree$tip.label[j]) # Getting the tip name to get the colour
            colour_tip <- tipcols[[colour_name]] # Getting the colour of the dots
            lines(c(obj$xx[j], max(obj$xx[1:Ntip(tree)])), rep(obj$yy[j], 2), lty = "dashed", col = colour_tip) # Add the tip lines
        }
        
        edgelabels(text = edgelabs_orthologs, frame = "none", cex = 0.5, adj = c(-0.5, -0.5)) # Adding the edge labels
        
        # Adding clade labels
        if (!is.null(clades)) {
            for (k in 1:nrow(clades)) {
                clade_name <- clades[k, 2]
                node <- as.numeric(clades[k, 1])
                if (node %in% obj$edge[, 2]) { # Check if the node is in the current plot
                  cladelabels( text = clade_name, node = node ,orientation = "vertical", cex = 0.5, offset = 10)
                }
            }
        }
    }
    if (!is.null(file)) dev.off()
}



