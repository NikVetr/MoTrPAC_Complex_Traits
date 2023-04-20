
#### define functions ####
polar2cart <- function(t, r){
  if(length(t) == 1 & length(r) == 1){
    return(c(r*cos(t), r * sin(t)))  
  } else if(length(t) == 1){
    return(do.call(rbind, lapply(r, function(ri) polar2cart(t,ri))))
  } else if(length(r) == 1){
    return(do.call(rbind, lapply(t, function(ti) polar2cart(ti,r))))
  } else {
    if(length(t) > length(r)){
      rep(r, times = ceiling(length(t) / length(r)))[1:length(t)]
    } else {
      rep(t, times = ceiling(length(r) / length(t)))[1:length(r)]
    }
    return(do.call(rbind, lapply(1:length(t), function(i) polar2cart(t[i],r[i]))))
  }
  
}

polarp <- function (t, r, col = 1, pch = 1, cex = 1, center = c(0,0)) {
  n <- length(t)
  z <- cbind(t, r)
  xy <- pol2cart(z)
  if (n == 1) 
    dim(xy) <- c(1, 2)
  hy <- hypot(xy[, 1], xy[, 2])
  points(xy[, 1] + center[1], xy[, 2] + center[2], cex = cex, col = col, pch = pch)
}

polarl <- function (t, r, col = 1, lwd = 1, center = c(0,0), adjx = 1) {
  n <- length(t)
  z <- cbind(t, r)
  xy <- pol2cart(z)
  if (n == 1) 
    dim(xy) <- c(1, 2)
  hy <- hypot(xy[, 1], xy[, 2])
  lines(xy[, 1] * adjx + center[1], xy[, 2] + center[2], lwd = lwd, col = col)
}

arc <- function(t1,t2,r1,r2,res=50,lwd = 1,col=1, mindist = T, self_adjust = 2 * pi / 45, random_selfing = T, clockwise_selfing = T, pointy_selfing = F,
                center = c(0,0), adjx = 1){
  if(mindist){
    if(abs(t1-t2) > pi){
      if(t1 > t2){
        t1 <- t1-2*pi
      } else {
        t2 <- t2-2*pi
      }
    }
  }
  if(abs(t1-t2) < 1E-6){
    
    if(random_selfing){
      lor <- sample(c(-1,1), 1)
    } else {
      lor <- ifelse(clockwise_selfing, -1, 1)
    }
    
    if(pointy_selfing){
      ts <- c(seq(t1, t1 + lor * self_adjust, length.out = res/2), seq(t1 + lor * self_adjust, t1, length.out = res/2))
      rs <- seq(r1, r2, length.out = res)
    } else {
      center <- center + pol2cart(cbind(t2, mean(c(r1, r2))))
      if(random_selfing){
        ts <- seq(t1, t1 + sample(c(-pi, pi), 1), length.out = res)
      } else {
        ts <- seq(t1, t1 + ifelse(clockwise_selfing, -pi, pi), length.out = res)
      }
      rs <- seq(max(c(r1,r2)) - mean(c(r1,r2)), mean(c(r1,r2)) - min(c(r1,r2)), length.out = res)
    }
    
  } else {
    ts <- seq(t1, t2, length.out = res)
    rs <- seq(r1, r2, length.out = res)
  }
  polarl(ts, rs, lwd = lwd,col=col, center = center, adjx = adjx)
}

plotMatrix <- function(mobject, size, location, lwd = 2, lwd_inner = 1.5, grid = T, font = 1, cex = 1, rownames = T, colnames = T, title = T, title.label = "Matrix Object"){
  lines(rbind(location, location + c(0,size[2])), lwd = lwd)
  lines(rbind(location, location + c(size[1]/8,0)), lwd = lwd)
  lines(rbind(location + c(0, size[2]), location + c(size[1]/8,size[2])), lwd = lwd)
  lines(rbind(location + c(size[1],0), location + size), lwd = lwd)
  lines(rbind(location + size, location + size - c(size[1]/8,0)), lwd = lwd)
  lines(rbind(location + c(size[1],0), location + c(size[1],0) - c(size[1]/8,0)), lwd = lwd)
  if(grid == T){
    for(i in 1:(dim(mobject)[1]-1)){
      lines(rbind(location + c(0,i*size[2]/dim(mobject)[1]), location + c(size[1], i*size[2]/dim(mobject)[1])), lwd = lwd_inner)
    }
    for(j in 1:(dim(mobject)[2]-1)){
      lines(rbind(location + c(j*size[1]/dim(mobject)[2],0), location + c(j*size[1]/dim(mobject)[2], size[2])), lwd = lwd_inner)
    }
  }
  if(class(mobject[1,1]) != "expression" & class(mobject[1,1]) != "character"){mobject <- matrix(as.character(mobject), nrow = dim(mobject)[1], ncol = dim(mobject)[2])}
  for(i in 1:(dim(mobject)[1])){
    for(j in 1:dim(mobject)[2]){
      text(labels = mobject[i,j], x = location[1] + (j-1/2)*size[1]/dim(mobject)[2], y = location[2] + size[2] - (i-1/2)*size[2]/dim(mobject)[1], font = font, cex = cex)
    }
  }
  if(title){
    text(title.label, x = location[1] + size[1]/2, y = location[2] + size[2] + strheight(title.label, font = 2, cex = 1.5)/1.5, cex = 1.5, font = 2)
  }
  if(rownames){
    for(i in 1:dim(mobject)[1]){
      text(rownames(mobject)[i], x = location[1] - strwidth(rownames(mobject)[i])/2 - size[1]/(ncol(mobject)*6), y = location[2] + size[2] - (i-1/2)*size[2]/dim(mobject)[2])
    }
  }
  if(colnames){
    for(i in 1:dim(mobject)[1]){
      text(colnames(mobject)[i], x = location[1] + (i-1/2)*size[1]/dim(mobject)[1], y = location[2] - strheight(colnames(mobject)[i])/2- size[2]/(nrow(mobject)*6))
    }
  }
}

draw.contour<-function(a,alpha=0.95,plot.dens=FALSE, line.width=2, line.type=1, limits=NULL, density.res=300,spline.smooth=-1,...){
  ##a is a list or matrix of x and y coordinates (e.g., a=list("x"=rnorm(100),"y"=rnorm(100)))
  ## if a is a list or dataframe, the components must be labeled "x" and "y"
  ## if a is a matrix, the first column is assumed to be x, the second y
  ##alpha is the contour level desired
  ##if plot.dens==TRUE, then the joint density of x and y are plotted,
  ##   otherwise the contour is added to the current plot.
  ##density.res controls the resolution of the density plot
  
  ##A key assumption of this function is that very little probability mass lies outside the limits of
  ## the x and y values in "a". This is likely reasonable if the number of observations in a is large.
  
  require(MASS)
  require(ks)
  if(length(line.width)!=length(alpha)){
    line.width <- rep(line.width[1],length(alpha))
  }
  
  if(length(line.type)!=length(alpha)){
    line.type <- rep(line.type[1],length(alpha))
  }
  
  if(is.matrix(a)){
    a=list("x"=a[,1],"y"=a[,2])
  }
  ##generate approximate density values
  if(is.null(limits)){
    limits=c(range(a$x),range(a$y))
  }
  f1<-kde2d(a$x,a$y,n=density.res,lims=limits)
  
  ##plot empirical density
  if(plot.dens) image(f1,...)
  
  if(is.null(dev.list())){
    ##ensure that there is a window in which to draw the contour
    plot(a,type="n",xlim=limits[1:2],ylim=limits[3:4],...)
  }
  
  ##estimate critical contour value
  ## assume that density outside of plot is very small
  
  zdens <- rev(sort(f1$z))
  Czdens <- cumsum(zdens)
  Czdens <- (Czdens/Czdens[length(zdens)])
  for(cont.level in 1:length(alpha)){
    ##This loop allows for multiple contour levels
    crit.val <- zdens[max(which(Czdens<=alpha[cont.level]))]
    
    ##determine coordinates of critical contour
    b.full=contourLines(f1,levels=crit.val)
    for(c in 1:length(b.full)){
      ##This loop is used in case the density is multimodal or if the desired contour
      ##  extends outside the plotting region
      b=list("x"=as.vector(unlist(b.full[[c]][2])),"y"=as.vector(unlist(b.full[[c]][3])))
      
      ##plot desired contour
      line.dat<-xspline(b,shape=spline.smooth,open=TRUE,draw=FALSE)
      lines(line.dat,lty=line.type[cont.level],lwd=line.width[cont.level])
    }
  }
}

minwhich <- function(x){min(which(x))}

maxwhich <- function(x){max(which(x))}

shadowtext <- function(x, y=NULL, labels, col='white', bg='black',
                       theta= seq(pi/4, 2*pi, length.out=8), r=0.1, ... ) {
  
  
  
  xy <- xy.coords(x,y)
  
  xo <- r*strwidth('A')
  
  yo <- r*strheight('A')
  
  for (i in theta) {
    
    text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, 
          
          labels, col=bg, ... )
    
  }
  
  text(xy$x, xy$y, labels, col=col, ... )
  
}

overlaps_with_zero <- function(qi){qi[1,] < 0 & qi[2,] > 0}

addImg <- function(
  obj, # an image file imported as an array (e.g. png::readPNG, jpeg::readJPEG)
  x = NULL, # mid x coordinate for image
  y = NULL, # mid y coordinate for image
  width = NULL, # width of image (in x coordinate units)
  interpolate = TRUE # (passed to graphics::rasterImage) A logical vector (or scalar) indicating whether to apply linear interpolation to the image when drawing. 
){
  if(is.null(x) | is.null(y) | is.null(width)){stop("Must provide args 'x', 'y', and 'width'")}
  USR <- par()$usr # A vector of the form c(x1, x2, y1, y2) giving the extremes of the user coordinates of the plotting region
  PIN <- par()$pin # The current plot dimensions, (width, height), in inches
  DIM <- dim(obj) # number of x-y pixels for the image
  ARp <- DIM[1]/DIM[2] # pixel aspect ratio (y/x)
  WIDi <- width/(USR[2]-USR[1])*PIN[1] # convert width units to inches
  HEIi <- WIDi * ARp # height in inches
  HEIu <- HEIi/PIN[2]*(USR[4]-USR[3]) # height in units
  rasterImage(image = obj, 
              xleft = x-(width/2), xright = x+(width/2),
              ybottom = y-(HEIu/2), ytop = y+(HEIu/2), 
              interpolate = interpolate)
}

subset_samps <- function(include = "", exclude = "", samps){
  incl_inds <- unique(unlist(lapply(include, function(i) grep(i, colnames(samps)))))
  excl_inds <- unique(unlist(lapply(exclude, function(i) grep(i, colnames(samps)))))
  return_inds <- setdiff(incl_inds, excl_inds)
  return(samps[,return_inds])
}


prop_greater_than_0 <- function(x) mean(x>0)

logit <- function(p) log(p/(1-p))
logit_prime <- function(p) (1/p + 1 / (1-p))
invlogit <- function(x) {exp(x) / (1+exp(x))}

intersect_rec <- function(x) {
  if(length(x) <= 1){
    return(x)
  }else if(length(x) > 2){
    x[[1]] <- intersect(x[[1]], x[[2]])
    intersect_rec(x[-2])
  } else {
    intersect(x[[1]], x[[2]])
  }
}

logit <- function(p) log(p / (1-p))
invlogit <- function(x) exp(x)/(1+exp(x))
squish_middle_p <- function(p,f) invlogit(logit(p)*f)
unsquish_middle_p <- function(p,f) invlogit(logit(p)/f)
# squish_middle_x <- function(x,f) log(abs(x)+1)/log(f)*sign(x)
# unsquish_middle_x <- function(x,f) (f^(abs(x))-1)*sign(x)
squish_middle_x <- function(x,f) asinh(x*f)
unsquish_middle_x <- function(x,f) sinh(x)/f
redistribute <- function(x, incr){
  new_x <- seq(from = min(x), length.out = length(x), by = incr)[rank(x)]
  new_x - max(new_x) + max(x)
}


fig_label <- function(text, region="figure", pos="topleft", cex=NULL, shrinkX = 0.95, shrinkY = 0.95, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1*shrinkX, y1*shrinkY, text, cex=cex, ...)
  return(invisible(c(x,y)))
}

line <- function(t,r1,r2,lwd = 1,col=1, center = c(0,0)){
  polarl(rep(t, 2), c(r1,r2), lwd = lwd,col=col, center = center)
}
logit <- function(p) log(p / (1-p))
invlogit <- function(x) exp(x)/(1+exp(x))
squish_middle_p <- function(p,f) invlogit(logit(p)*f)
unsquish_middle_p <- function(p,f) invlogit(logit(p)/f)
# squish_middle_x <- function(x,f) log(abs(x)+1)/log(f)*sign(x)
# unsquish_middle_x <- function(x,f) (f^(abs(x))-1)*sign(x)
squish_middle_x <- function(x,f) asinh(x*f)
unsquish_middle_x <- function(x,f) sinh(x)/f
redistribute <- function(x, incr){
  new_x <- seq(from = min(x), length.out = length(x), by = incr)[rank(x)]
  new_x - max(new_x) + max(x)
}
ifelse2 <- function(bool, opt1, opt2){if(bool){return(opt1)}else{return(opt2)}}

text_indivcolor <- function(labels_mat, colors_mat, xloc, cex, ...){
  text(labels = labels_mat[,1], col = colors_mat[,1], x = xloc, cex = cex, ...)
  for(i in 2:ncol(labels_mat)){
    prior_string_widths <- apply(do.call(cbind, lapply(1:(i-1), function(li) strwidth(labels_mat[,li], units = "user"))), 1, sum) * cex
    text(labels = labels_mat[,i], col = colors_mat[,i], x = xloc + prior_string_widths, cex = cex, ...)
  }
}

ifelse3 <- function(bool, opt1, opt2){sapply(bool, function(bool_i) ifelse(bool_i, opt1, opt2))}
min2 <- function(vec, opt2){sapply(vec, function(veci) min(veci, opt2))}

#functions from transcr vs. prot stuff
text_cols <- function(string, cols, x, y, cex = 1, ...){
  for(char_i in 1:nchar(string)){
    txt_exploded <- c(substr(string, 1, char_i-1), substr(string, char_i, char_i), substr(string, char_i+1, nchar(string)))
    text(x = x, y = y, labels = bquote(phantom(.(txt_exploded[1])) * .(txt_exploded[2]) * phantom(.(txt_exploded[3]))), col = cols[char_i], cex = cex, ...)
  }
}

find_optimal_cex_and_lines <- function(txt, rect_coords, rect_rescaling_ratio = 0.95, srt_height_rescaling_ratio = 1.4){
  
  strwidths <- strwidth(txt)
  strheight <- strheight(txt[1]) * srt_height_rescaling_ratio 
  space_width_min <- strwidth(" ")
  rectwidth <- abs(rect_coords$x1 - rect_coords$x0) * rect_rescaling_ratio
  rectheight <- abs(rect_coords$y1 - rect_coords$y0) * rect_rescaling_ratio
  
  # ceiling(cumsum(strwidths) / rectwidth)
  data <- list(strwidths = strwidths, strheight = strheight, space_width_min = space_width_min, rectwidth = rectwidth, rectheight = rectheight)
  par <- log((rectwidth * rectheight) / ((sum(strwidths) + space_width_min * length(strwidths)) * strheight) * 0.5) #initialize cex
  while(compute_space_remaining(data, par) == Inf){
    par <- par + log(0.5)
  }
  # plot(1:120/100, sapply(log(1:120/100), function(cex) compute_space_remaining(data, cex)), type = "l")
  opt_cex <- suppressWarnings(optimx::optimx(par = par, fn = compute_space_remaining, data = data, hess = NULL, 
                                             method = c('Nelder-Mead'), hessian=FALSE, #can't compute hessian bc of sharp jumps when new line is formed? or maybe not?
                                             control = list(maxit = 1E4, trace = 0, kkt=FALSE)))
  # compute_space_remaining(data = data, par = opt_cex$p1)
  return(list(cex = exp(opt_cex$p1), 
              words_on_lines = put_words_on_lines(data = data, par = exp(opt_cex$p1)),
              space_width_min = space_width_min * exp(opt_cex$p1),
              vertical_space = strheight * exp(opt_cex$p1))
  )
  
}

compute_space_remaining <- function(data, par){
  
  #clarify par-cex relationship
  cex <- exp(par)
  
  #unwrap data
  strwidths_mod <- data$strwidths * cex
  strheight_mod <- data$strheight * cex
  space_width_min_mod <- data$space_width_min * cex
  rectwidth_mod <- data$rectwidth
  rectheight_mod <- data$rectheight
  
  #check that no words are wider than a line
  if(any(strwidths_mod > rectwidth_mod)){
    return(Inf)
  }
  
  txi <- 1
  linei <- 1
  current_width <- strwidths_mod[txi]
  txt_lines <- list()
  txt_lines[[linei]] <- txi
  
  while(txi < length(txt)){
    txi <- txi + 1
    txt_lines[[linei]] <- c(txt_lines[[linei]], txi)
    current_width <- current_width + strwidths_mod[txi] + space_width_min_mod
    if(current_width > rectwidth_mod){
      txt_lines[[linei]] <- txt_lines[[linei]][-length(txt_lines[[linei]])]
      linei <- linei + 1
      txt_lines[[linei]] <- txi
      current_width <- strwidths_mod[txi]
    }
  }
  
  last_line_width_remaining <- rectwidth_mod - sum(strwidths_mod[txt_lines[[linei]]])
  current_height <- linei * strheight_mod
  space_remaining <- (rectheight_mod - current_height) * rectwidth_mod + (last_line_width_remaining * strheight_mod)
  
  if(space_remaining < 0){return(Inf)} else {return(space_remaining)}
  
}

put_words_on_lines <- function(data, par){
  
  #clarify par-cex relationship
  cex <- par
  
  #unwrap data
  strwidths_mod <- data$strwidths * cex
  strheight_mod <- data$strheight * cex
  space_width_min_mod <- data$space_width_min * cex
  rectwidth_mod <- data$rectwidth
  rectheight_mod <- data$rectheight
  
  txi <- 1
  linei <- 1
  current_width <- strwidths_mod[txi]
  txt_lines <- list()
  txt_lines[[linei]] <- txi
  
  while(txi < length(txt)){
    txi <- txi + 1
    txt_lines[[linei]] <- c(txt_lines[[linei]], txi)
    current_width <- current_width + strwidths_mod[txi] + space_width_min_mod
    if(current_width > rectwidth_mod){
      txt_lines[[linei]] <- txt_lines[[linei]][-length(txt_lines[[linei]])]
      linei <- linei + 1
      txt_lines[[linei]] <- txi
      current_width <- strwidths_mod[txi]
    }
  }
  
  return(txt_lines)
  
}

text_wrapped_words <- function(txt, rect_coords, optimal_word_placement_inf, justified = F, str_height_lower_start_ratio = 0.75, 
                               str_width_lefter_start_ratio = 0.01, rect_rescaling_ratio = 0.95, col = "black", multicolor_words = F, cols_list, ...){
  
  curr_x <- rect_coords$x0 - abs(rect_coords$x0 - rect_coords$x1) * str_width_lefter_start_ratio
  curr_y <- rect_coords$y0 - optimal_word_placement_inf$vertical_space * str_height_lower_start_ratio
  nlines <- length(optimal_word_placement_inf$words_on_lines)
  strwidths_plotting <- strwidth(txt) * optimal_word_placement_inf$cex
  space_left_on_lines <- sapply(1:length(optimal_word_placement_inf$words_on_lines), function(linei)
    abs(rect_coords$x0 - rect_coords$x1) * rect_rescaling_ratio - sum(strwidths_plotting[optimal_word_placement_inf$words_on_lines[[linei]]]))
  justified_space_between_words <- sapply(1:length(optimal_word_placement_inf$words_on_lines), function(linei)
    space_left_on_lines[linei] / (length(optimal_word_placement_inf$words_on_lines[[linei]]) - 1))
  
  words_written <- 0
  for(linei in 1:nlines){
    for(wordi in 1:length(optimal_word_placement_inf$words_on_lines[[linei]])){
      words_written <- words_written + 1
      word_to_write <- txt[optimal_word_placement_inf$words_on_lines[[linei]][wordi]]
      if(multicolor_words){
        text_cols(x = curr_x, y = curr_y, cex = optimal_word_placement_inf$cex,
                  string = word_to_write, pos = 4, cols = cols_list[[words_written]])
      } else {
        text(x = curr_x, y = curr_y, cex = optimal_word_placement_inf$cex,
             labels = word_to_write, pos = 4, col = col)
      }
      if(justified){
        curr_x <- curr_x + strwidth(word_to_write) * optimal_word_placement_inf$cex + justified_space_between_words[linei]
      } else {
        curr_x <- curr_x + strwidth(word_to_write) * optimal_word_placement_inf$cex + optimal_word_placement_inf$space_width_min  
      }
      
    }
    curr_x <- rect_coords$x0 - abs(rect_coords$x0 - rect_coords$x1) * str_width_lefter_start_ratio
    curr_y <- rect_coords$y0 - optimal_word_placement_inf$vertical_space * str_height_lower_start_ratio - optimal_word_placement_inf$vertical_space * linei
  }
  
}


seq_incl_0 <- function(from, to, length.out){
  incr <- (to - from) / length.out
  testseq <- c(rev(seq(from = 0, to = from - incr, by = -incr)), seq(from = 0, to = to + incr, by = incr)[-1])
  foo <- 1
  while(length(testseq) > length.out){
    incr <- (to - from) / (length.out-foo)
    testseq <- c(rev(seq(from = 0, to = from - incr, by = -incr)), seq(from = 0, to = to + incr, by = incr)[-1])
    foo <- foo + 1
  }
  foo <- -1
  while(length(testseq) < length.out){
    incr <- (to - from) / (length.out-foo)
    testseq <- c(rev(seq(from = 0, to = from - incr, by = -incr)), seq(from = 0, to = to + incr, by = incr)[-1])
    foo <- foo - 1
  }
  return(testseq)
}

dexp_smooth <- function(x, y, r, reweight_trunc_tail = F, reweight_n_pts = F, fix_endpoints = F, interpolate_at = NA){
  
  nx <- length(x)
  if(is.na(interpolate_at[1])){
    n <- nx
    w <- t(sapply(1:(n-1), function(i) c(dexp((x[i] - x[0:(i-1)]), rate = r), dexp(0, rate = r), dexp(-(x[i] - x[(i+1):n]), rate = r))))
    w <- rbind(w, dexp(x[n] - x, rate = r))
  } else {
    n <- length(interpolate_at)
    w <- t(sapply(1:(n-1), function(i) c(dexp((interpolate_at[i] - x[x<=interpolate_at[i]]), rate = r), 
                                         dexp((x[x>interpolate_at[i]] - interpolate_at[i]), rate = r))))
    w <- rbind(w, dexp(x[nx] - x, rate = r))
  }
  # dim(w)
  # plot(w[2300,])
  # plot(w[2400,])
  # plot(w[2500,])
  # plot(w[2600,])
  
  if(reweight_trunc_tail){
    if(is.na(interpolate_at[1])){
      tw <- t(sapply(1:(n-1), function(i) c(pexp((x[i] - x[1]), rate = r), pexp(-(x[i] - x[n]), rate = r))))
    } else {
      tw <- t(sapply(1:(n-1), function(i) c(pexp((interpolate_at[i] - x[1]), rate = r), pexp(-(interpolate_at[i] - x[nx]), rate = r))))
    }
    tw[1,] <- tw[2,]
    tw <- rbind(tw, tw[n-1,])
    tw <- 1/tw
    if(is.na(interpolate_at[1])){
      tw <- t(sapply(1:n, function(ri) c(rep(tw[ri,1], ri), rep(1, n-ri)) * c(rep(1, ri), rep(tw[ri,2], n-ri)))) #eh let's give it to the nth pt too
    } else {
      tw <- t(sapply(1:n, function(ri) c(rep(tw[ri,1], sum(interpolate_at[ri] >= x)), rep(1, nx-sum(interpolate_at[ri] >= x))) * 
                       c(rep(1, sum(interpolate_at[ri] >= x)), rep(tw[ri,2], nx-sum(interpolate_at[ri] >= x)))))
    }
    w <- t(sapply(1:n, function(ri) w[ri,] * tw[ri,]))
    
  }
  # dim(w)
  # plot(w[2300,])
  # plot(w[2400,])
  # plot(w[2500,])
  # plot(w[2600,])
  
  if(reweight_n_pts){
    if(is.na(interpolate_at[1])){
      tw <- cbind(0:(n-1), (n-1):0)
      tw <- 1/tw
    } else {
      tw <- sapply(1:n, function(i) sum(interpolate_at[i] >= x))
      tw <- cbind(tw, nx - tw)
    }
    mintw <- apply(tw, 1, min)
    tw[,1] <- tw[,1] / mintw 
    tw[,2] <- tw[,2] / mintw 
    
    tw[1,] <- tw[2,]; tw[n,] <- tw[n-1,]
    
    if(is.na(interpolate_at[1])){
      tw <- t(sapply(1:n, function(ri) c(rep(tw[ri,1], ri), rep(1, n-ri)) * c(rep(1, ri), rep(tw[ri,2], n-ri)))) #eh let's give it to the nth pt too
    } else {
      tw <- t(sapply(1:n, function(ri) c(rep(tw[ri,1], sum(interpolate_at[ri] >= x)), rep(1, nx-sum(interpolate_at[ri] >= x))) * 
                       c(rep(1, sum(interpolate_at[ri] >= x)), rep(tw[ri,2], nx-sum(interpolate_at[ri] >= x)))))
    }
    
    w <- t(sapply(1:n, function(ri) w[ri,] * tw[ri,]))
  }
  # dim(w)
  # plot(w[2300,])
  # plot(w[2400,])
  # plot(w[2500,])
  # plot(w[2600,])
  
  if(fix_endpoints){
    w[1,] <- c(1, rep(0,nx-1))
    w[n,] <- c(rep(0,nx-1), 1)
  }
  
  wy <- c(w %*% t(t(y)))
  wy <- wy / apply(w,1,sum)
  return(wy)
}

trim_n <- function(x, n){
  ox <- order(x)
  xi_to_trim <- c(head(ox, n = n / 2), tail(ox, n = n / 2))
  return(x[-xi_to_trim])
}

grad_arrow <- function(arrow_locs, prop_head_length = 0.2, prop_shaft_width = 0.5,
                       cols = c("white", "black"), nslices = 100){
  colsgrad <- colorRampPalette(colors = cols)(nslices)
  # split the base of the arrow in nslices and fill each progressively
  ys <- seq(arrow_locs[3], 
            arrow_locs[4] - diff(arrow_locs[3:4]) * prop_head_length, 
            len = nslices + 1)
  for (i in 1:nslices) {
    polygon(c(arrow_locs[1] + diff(arrow_locs[1:2]) * prop_shaft_width / 2, 
              arrow_locs[1] + diff(arrow_locs[1:2]) * prop_shaft_width / 2, 
              arrow_locs[2] - diff(arrow_locs[1:2]) * prop_shaft_width / 2, 
              arrow_locs[2] - diff(arrow_locs[1:2]) * prop_shaft_width / 2), 
            c(ys[i], ys[i+1], ys[i+1], ys[i]), col = colsgrad[i], border = NA)
  }
  # add a filled arrowhead
  polygon(c(arrow_locs[1], mean(arrow_locs[1:2]), arrow_locs[2], arrow_locs[1]), 
          c(arrow_locs[4] - diff(arrow_locs[3:4]) * prop_head_length, 
            arrow_locs[4], 
            arrow_locs[4] - diff(arrow_locs[3:4]) * prop_head_length, 
            arrow_locs[4] - diff(arrow_locs[3:4]) * prop_head_length), 
          col = cols[2], border = cols[2])
  
}

mcprint <- function(...){
  system(sprintf('echo "%s"', paste0(..., collapse="")))
}

jaccard <- function(x1, x2) length(intersect(x1, x2)) / length(union(x1, x2))

#functions for estimation of probit correlations

probs_quadrants <- function(mus, rij){
  probs <- matrix(0, 2, 2)
  probs[2,1] <- mvtnorm::pmvnorm(lower = c(-Inf, -Inf), upper = c(0, 0), corr = matrix(c(1,rij,rij,1), 2, 2), mean = mus)[1]
  probs[1,1] <- mvtnorm::pmvnorm(lower = c(0, -Inf), upper = c(Inf, 0), corr = matrix(c(1,rij,rij,1), 2, 2), mean = mus)[1]
  probs[1,2] <- mvtnorm::pmvnorm(lower = c(0, 0), upper = c(Inf, Inf), corr = matrix(c(1,rij,rij,1), 2, 2), mean = mus)[1]
  probs[2,2] <- 1 - sum(probs)
  return(probs)
}

ll_biv_probit <- function(data, par) {
  
  #undo params
  rij <- par[1]
  
  #undo data
  ct2x2 <- data$ct2x2
  mus <- data$mus
  
  #compute log likelihood
  probs <- probs_quadrants(mus, rij)
  log_probs <- log(probs)
  ll <- sum(log_probs * ct2x2)
  
  #regularize?
  ll <- ll + dbeta((par + 1) / 2, shape1 = 5, shape2 = 5, log = T)
  
  #return log likelihood
  return(-ll)
  
}

ll_biv_probit_mus <- function(data, par) {
  
  #undo params
  rij <- par[1]
  mus <- par[2:3]
  
  #undo data
  ct2x2 <- data$ct2x2
  
  #compute log likelihood
  probs <- probs_quadrants(mus, rij)
  log_probs <- log(probs)
  ll <- sum(log_probs * ct2x2)
  
  #regularize?
  ll <- ll + dbeta((rij + 1) / 2, shape1 = 5, shape2 = 5, log = T)
  
  #return log likelihood
  return(-ll)
  
}

estimate_correlations <- function(y, t, optimize_jointly = F, ncores = 12, nearPDadj = T, print_progress = F) {
  
  #retrieve metadata
  ncols <- length(y)
  
  #reorder 't's to match 'y's
  if(!all(names(y) %in% names(t))){return(NA)}
  t <- t[names(y)]
  
  #iterate through all pairs of dimensions, optimizing bivariate probit probability
  estimated_corrs <- mclapply(1:(ncols-1), function(ci1) sapply((ci1+1):ncols, function(ci2){
    
    #find compatible intersect and subset
    all_poss <- intersect(t[[ci1]], t[[ci2]])
    n <- length(all_poss)
    y1 <- y[[ci1]][y[[ci1]] %in% all_poss]
    y2 <- y[[ci2]][y[[ci2]] %in% all_poss]
    
    #return 0 if there are no obs in y1 or y2
    if(length(y1) == 0 | length(y2) == 0){
      return(0)
    }
    
    #snag 2x2 contingency table
    ct2x2 <- matrix(c(length(setdiff(y1,y2)), 
                      n - length(union(y1,y2)), 
                      length(intersect(y1,y2)), 
                      length(setdiff(y2,y1))), 
                    2, 2)
    
    #can optimize jointly or marginally, knowing the location of the mus
    if(optimize_jointly){
      optim_out <- optimx::optimx(par = rep(0,3), 
                                  fn = ll_biv_probit_mus, 
                                  data = list(ct2x2 = ct2x2), method = "nlm", lower = c(-1,-Inf,-Inf), upper = c(1,Inf,Inf),
                                  control = list(maxit = 1E3, trace = 0, dowarn = F))
    } else {
      mus <- qnorm(c(sum(ct2x2[1,]), sum(ct2x2[,2])) / sum(ct2x2))
      optim_out <- optimx::optimx(par = 0, 
                                  fn = ll_biv_probit, 
                                  data = list(ct2x2 = ct2x2, mus = mus), method = "nlm", lower = -1, upper = 1,
                                  control = list(maxit = 1E3, trace = 0, dowarn = F))
    }
    
    #print current location
    if(print_progress){
      mcprint(paste0(ci1, ", ", ci2, ", ", optim_out$p1))  
    }
    
    #return estimate
    return(optim_out$p1)
    
  }), mc.cores = ncores)
  
  #did any traits fail?
  failed_traits <- which(sapply(estimated_corrs, class) == "try-error")
  
  #reconstruct and (optionally) adjust estimated matrix
  estimated_corr_mat <- matrix(0, nrow = ncols, ncol = ncols)
  for(ri in 1:(ncols-1)){
    corrs_to_insert <- estimated_corrs[[ri]]
    estimated_corr_mat[ri, (ri+1):ncols] <- corrs_to_insert
  }
  estimated_corr_mat <- estimated_corr_mat + t(estimated_corr_mat) + diag(ncols)
  if(nearPDadj){
    estimated_corr_mat <- as.matrix(Matrix::nearPD(estimated_corr_mat, corr = T)$mat  )
  }
  
  rownames(estimated_corr_mat) <- colnames(estimated_corr_mat) <- names(y)
  
  return(estimated_corr_mat)
  
}

ll_biv_probit_multi <- function(data, par) {
  
  #undo params
  rij <- par[1]
  
  #undo data
  cts <- data$cts
  mus <- data$mus
  rown <- setNames(names(mus), names(mus))
  
  #compute log likelihood
  ll <- sum(unlist(lapply(rown, function(rown_i){
    probs <- probs_quadrants(mus[[rown_i]], rij)
    log_probs <- log(probs)
    ll <- sum(log_probs * cts[[rown_i]])
    ll
  })))
  
  
  #regularize?
  ll <- ll + dbeta((par + 1) / 2, shape1 = 10, shape2 = 10, log = T)
  
  #return log likelihood
  return(-ll)
  
}

estimate_correlations_multi <- function(y, t, ncores = 12, nearPDadj = T, print_progress = F) {
  
  #retrieve metadata
  ncols <- length(y)
  coln <- setNames(names(y), names(y))
  nrows <- length(y[[1]])
  rown <- setNames(names(y[[1]]), names(y[[1]]))
  
  #reorder 't's to match 'y's
  if(!all(names(y) %in% names(t))){return(NA)}
  t <- t[names(y)]
  
  #iterate through all pairs of dimensions, optimizing bivariate probit probability
  estimated_corrs <- mclapply(1:(ncols-1), function(ci1) sapply((ci1+1):ncols, function(ci2){
    
    #get CTs per tissue
    cts <- lapply(rown, function(ri){
      #find compatible intersect and subset
      all_poss <- intersect(t[[ci1]][[ri]], t[[ci2]][[ri]])
      n <- length(all_poss)
      y1 <- y[[ci1]][[ri]][ y[[ci1]][[ri]] %in% all_poss ]
      y2 <- y[[ci2]][[ri]][ y[[ci2]][[ri]] %in% all_poss ]
      
      #return 0 if there are no obs in y1 or y2
      if(length(y1) == 0 | length(y2) == 0){
        return(NULL)
      }
      
      #snag 2x2 contingency table
      ct2x2 <- matrix(c(length(setdiff(y1,y2)), 
                        n - length(union(y1,y2)), 
                        length(intersect(y1,y2)), 
                        length(setdiff(y2,y1))), 
                      2, 2)
    })
    
    #get rid of null values
    cts <- cts[!sapply(cts, is.null)]
    if(length(cts) == 0){return(0)}
    
    if(print_progress){
      mcprint(paste0("starting (", ci1, ", ", ci2, "): ", length(cts)))  
    }
    
    #can optimize jointly or marginally, knowing the location of the mus
    mus <- lapply(cts, function(ct2x2) qnorm(c(sum(ct2x2[1,]), sum(ct2x2[,2])) / sum(ct2x2)))
    optim_out <- optimx::optimx(par = 0, 
                                fn = ll_biv_probit_multi, 
                                data = list(cts = cts, mus = mus), method = "nlm", lower = -1, upper = 1,
                                control = list(maxit = 1E3, trace = 0, dowarn = F))
  
    
    #print current location
    if(print_progress){
      mcprint(paste0("(", ci1, ", ", ci2, "): ", optim_out$p1))  
    }
    
    #return estimate
    return(optim_out$p1)
    
  }), mc.cores = ncores)
  
  #did any traits fail?
  failed_traits <- which(sapply(estimated_corrs, class) == "try-error")
  
  #reconstruct and (optionally) adjust estimated matrix
  estimated_corr_mat <- matrix(0, nrow = ncols, ncol = ncols)
  for(ri in 1:(ncols-1)){
    corrs_to_insert <- estimated_corrs[[ri]]
    estimated_corr_mat[ri, (ri+1):ncols] <- corrs_to_insert
  }
  estimated_corr_mat <- estimated_corr_mat + t(estimated_corr_mat) + diag(ncols)
  if(nearPDadj){
    estimated_corr_mat <- as.matrix(Matrix::nearPD(estimated_corr_mat, corr = T)$mat  )
  }
  
  rownames(estimated_corr_mat) <- colnames(estimated_corr_mat) <- names(y)
  
  return(estimated_corr_mat)
  
}

seq3 <- function(xlims, lout, contains = NA, err_up = F){
  if(!is.na(contains)){
    nlims <- c(min(contains, xlims), max(contains, xlims))
  } else {
    nlims <- xlims
  }
  incr <- abs(diff(nlims)) / lout
  mag <- 10^floor(log10(incr))
  incr <- ifelse(err_up, floor(incr / mag) * mag, ceiling(incr / mag) * mag)
  if(!is.na(contains)){
    down <- seq(contains, nlims[1], by = -incr)
    up <- seq(contains, nlims[2], by = incr)
    return(sort(unique(c(down, up))))
  } else {
    return(seq(nlims[1] + abs(nlims[1]) %% incr, nlims[2], by = incr))
  }
}

xyrat <- function(){
  prop <- c(diff(par("usr")[1:2]), diff(par("usr")[3:4])) / par("pin")
  prop[1] / prop[2]
}

rrect <- function(loc, w, h, pe = 0.25, npts = 50, rot = 0, hat_prop = 0.15, bold_border = 0,
                  col = 0, lwd = 1, border = 1, background_col = 0, ...){
  
  #get arc where yval has height 1, xval has height xyrat
  xyr <- xyrat()
  tr <- polar2cart(t = seq(0, pi/2, length.out = round(npts/4)), r = 1) %*% diag(c(xyr, 1))
  
  #scale to smaller edge
  tallboi <- w < (h * xyr)
  tr <- tr * ifelse(tallboi, pe * w / 2 / xyr, pe * h / 2)
  
  #shift to circumscribed centered rectangle quadrant
  tr[,1] <- (tr[,1] + w / 2 - tr[1,1])
  tr[,2] <- (tr[,2] + h / 2 - tr[nrow(tr),2])
  # rect(xleft = loc[1] - w/2, ybottom = loc[2] - h/2, xright = loc[1] + w/2, ytop = loc[2] + h/2)
  
  #find polygon points and recenter
  outer_coords_nc <- rbind(tr, cbind(-rev(tr[,1]), rev(tr[,2])), cbind(-tr[,1], -tr[,2]), cbind(rev(tr[,1]), -rev(tr[,2])))
  outer_coords <- outer_coords_nc + rep(loc, each = nrow(outer_coords_nc))
  
  #first draw outer border if desired & inner polygon
  if(bold_border > 1E-6){
    
    #draw solid outer polygon
    polygon(outer_coords, col = border, border = border, lwd = 1E-1)
    
    #adjust border to bounds
    bold_border <- max(c(min(c(bold_border, 1)), 0))
    
    if(tallboi){
      nw <- (1-bold_border) * w
      nh <- h - bold_border * w / xyr
    } else { #shortboi
      nw <- w - bold_border * h * xyr
      nh <- (1-bold_border) * h
    }
    
    tr <- polar2cart(t = seq(0, pi/2, length.out = round(npts/4)), r = 1) %*% diag(c(xyr, 1))
    tr <- tr * ifelse(tallboi, pe * nw / 2 / xyr, pe * nh / 2)
    tr[,1] <- (tr[,1] + nw / 2 - tr[1,1])
    tr[,2] <- (tr[,2] + nh / 2 - tr[nrow(tr),2])
    
    inner_coords_nc <- rbind(tr, cbind(-rev(tr[,1]), rev(tr[,2])), cbind(-tr[,1], -tr[,2]), cbind(rev(tr[,1]), -rev(tr[,2])))
    inner_coords <- inner_coords_nc + rep(loc, each = nrow(inner_coords_nc))
    # rect(xleft = loc[1] - nw/2, ybottom = loc[2] - nh/2, xright = loc[1] + nw/2, ytop = loc[2] + nh/2)
    
    # inner_coords <- outer_coords_nc %*% diag(ifelse2(tallboi,
    #                                                  c(1-bold_border, 1 - bold_border * w * xyr / h),
    #                                                  c(1 - bold_border * h / xyr / w, 1-bold_border))
    #                                          ) +
    #   rep(loc, each = nrow(inner_coords))
    
    polygon(inner_coords, col = background_col, lwd = 1E-1, border = background_col)
    polygon(inner_coords, col = col, lwd = 1E-1, border = border)
  } else {
    polygon(outer_coords, col = col, lwd = lwd, border = border)
  }
  
  #add a hat if desired
  if(hat_prop > 1E-6){
    sub_coords <- outer_coords[outer_coords[,2] > loc[2] - h/2 + h*(1-hat_prop),]
    if(bold_border > 1E-6){
      polygon(sub_coords, col = background_col, lwd = lwd, border = background_col)
      polygon(sub_coords, col = col, lwd = 1E-6, border = border)
    } else {
      polygon(sub_coords, col = col, lwd = lwd, border = border)  
    }
  }
  
}

grad_arrow_curve <- function(arrow_locs, prop_head_width = 2, prop_shaft_length = 0.1,
                             cols = c("white", "black"), nslices = 500, direc = "h", w = 0.25,
                             raster = T, xpd = NA, raster_res = 51, interp_raster = T,
                             col_alpha = 1, col_pow = 0.5, taper_ratio = 0.45, taper_pow = 1.5,
                             outline = T, outline_lwd = 2.5, outline_col = "black"){
  
  if(direc == "h"){
    dx <- diff(arrow_locs[1:2]) * (1-prop_shaft_length)
    dy <- diff(arrow_locs[3:4])
  } else if(direc == "v"){
    dx <- diff(arrow_locs[1:2])
    dy <- diff(arrow_locs[3:4]) * (1-prop_shaft_length)
  }
  
  colsgrad <- colorRampPalette(colors = cols)(nslices * 10)
  colsgrad <- colsgrad[round((seq(1, (nslices*10)^col_pow, length.out=nslices))^(1/col_pow))]
  colsgrad <- adjustcolor(colsgrad, col_alpha)
  
  # split the base of the arrow in nslices and fill each progressively
  piseq <- seq(-pi/2,pi/2,length.out=nslices+1)
  
  if(direc == "h"){
    xs <- seq(arrow_locs[1],arrow_locs[1]+dx, length.out=nslices+1)
    ys <- dy*(sin(piseq)+1)/2 + arrow_locs[3]
    m <- cos(piseq)/2 * sign(dy) * sign(dx)
    t <- atan(m)
    dispx <- w * sin(t) / 2
    dispy <- w * cos(t) / 2
  } else if (direc == "v"){
    ys <- seq(arrow_locs[3],arrow_locs[3]+dy, length.out=nslices+1)
    xs <- dx*(sin(piseq)+1)/2 + arrow_locs[1]
    m <- cos(piseq)/2 * sign(dy) * sign(dx)
    t <- atan(m)
    dispx <- w * cos(t) / 2
    dispy <- w * sin(t) / 2
  }
  
  #taper the arrow if desired
  taper <- seq(1, taper_ratio^taper_pow, length.out=nslices+1)^(1/taper_pow)
  dispx <- dispx * taper
  dispy <- dispy * taper
  
  #final coords
  coords <- data.frame(x1 = xs - dispx, y1 = ys + dispy,
                       x2 = xs + dispx, y2 = ys - dispy
  )
  if(direc == "h"){
    head_coords <- data.frame(x = rev(c(xs[nslices], arrow_locs[2], xs[nslices])), 
                              y = rev(c(ys[nslices] - prop_head_width * w / 2 * taper_ratio, ys[nslices], ys[nslices] + prop_head_width * w / 2 * taper_ratio)))
  } else if(direc == "v"){
    head_coords <- data.frame(x = c(xs[nslices] - prop_head_width * w / 2 * taper_ratio, xs[nslices], xs[nslices] + prop_head_width * w / 2 * taper_ratio), 
                              y = c(ys[nslices], arrow_locs[4], ys[nslices]))
  }
  
  
  if(raster){
    
    #get plotting params
    usr <- par("usr")
    upct <- par("plt")
    gr_usr <- usr + c(diff(usr[1:2]) / diff(upct[1:2]) * (c(0,1) - upct[1:2]),
                      diff(usr[3:4]) / diff(upct[3:4]) * (c(0,1) - upct[3:4]))
    
    #write to temporary png
    tmp <- tempfile()
    
    png(tmp, width = par("din")[1], height = par("din")[2], units = "in", res = raster_res, bg = "transparent", type="cairo")
    par(mar = c(0,0,0,0), xpd = NA)
    plot.new(); plot.window(xlim=gr_usr[1:2], ylim=gr_usr[3:4], xaxs = "i", yaxs = "i")
    
    #plot the arrows
    if(length(cols) == 1){
      polygon(x = c(coords$x1, rev(coords$x2)), y = c(coords$y1, rev(coords$y2)), cols = cols, border = NA)  
    } else {
      for(i in 1:nslices){
        polygon(x = c(coords$x1[i], coords$x1[i+1], coords$x2[i+1], coords$x2[i]), 
                y = c(coords$y1[i], coords$y1[i+1], coords$y2[i+1], coords$y2[i]), 
                col = colsgrad[i], border = NA) 
      }
    }
    if(direc == "h"){
      polygon(x = head_coords$x, y = head_coords$y, col = colsgrad[nslices], border = NA)
    } else if(direc == "v"){
      polygon(x = head_coords$x, y = head_coords$y, col = colsgrad[nslices], border = NA)
    }
    dev.off()
    
    #draw to file
    rasterImage(png::readPNG(tmp), gr_usr[1], gr_usr[3], gr_usr[2], gr_usr[4], interpolate = interp_raster, xpd = xpd)
    
    #delete temp file
    rm(tmp)
    
  } else {
    
    #or alternatively draw it in full vector graphics
    if(length(cols) == 1){
      polygon(x = c(coords$x1, rev(coords$x2)), y = c(coords$y1, rev(coords$y2)), cols = cols, border = NA, xpd = xpd)  
    } else {
      for(i in 1:nslices){
        polygon(x = c(coords$x1[i], coords$x1[i+1], coords$x2[i+1], coords$x2[i]), 
                y = c(coords$y1[i], coords$y1[i+1], coords$y2[i+1], coords$y2[i]), 
                col = colsgrad[i], border = NA, xpd = xpd) 
      }
    }
    if(direc == "h"){
      polygon(x = head_coords$x, y = head_coords$y, col = colsgrad[nslices], border = NA)
    } else if(direc == "v"){
      polygon(x = head_coords$x, y = head_coords$y, col = colsgrad[nslices], border = NA)
    }
  }
  
  if(outline){
    polygon(x = c(coords$x1, head_coords$x, rev(coords$x2)), 
            y = c(coords$y1, head_coords$y, rev(coords$y2)), 
            border = outline_col, xpd = xpd, lwd = outline_lwd)
  }
  
}

colblend <- function(c1, c2, brighten = 0){
  nc <- (col2rgb(c1) + col2rgb(c2)) / 2
  nc <- nc + (255-max(nc))*brighten
  rgb(nc[1], nc[2], nc[3],m=255)
}

#### redux of text-wrapping functions ####

latext <- function(labs, x, y, cex = 1, boxh = NA, first_line_col = 1, first_line_hadj = NA, col = 1, pos = NULL, first_line_center_loc = NA, ...){
  new_labs <- strsplit(labs, split = "\n")[[1]]
  new_labs_no_uflow <- lapply(new_labs, function(l) latex2exp::TeX(gsub("g|j|p|q|y|,|_|\\(|\\)|Q|u", "a", l)))
  new_labs_no_oflow <- lapply(new_labs, function(l) latex2exp::TeX(gsub("\\^", "a", l)))
  new_labs <- lapply(new_labs, function(l) latex2exp::TeX(l))
  
  wsh <- strheight("M\nM", cex = cex) - strheight("M", cex = cex) * 2
  lineh <- sapply(new_labs, strheight, cex = cex)  
  lineh_no_uflow <- sapply(new_labs_no_uflow, strheight, cex = cex)
  lineh_no_oflow <- sapply(new_labs_no_oflow, strheight, cex = cex)
  ebot <- lineh_no_uflow - lineh
  etop <- lineh_no_oflow - lineh
  uflow_adj <- (lineh_no_uflow - lineh) / 2 - (lineh - lineh_no_oflow) / 2
  uflow_adj <- (lineh_no_uflow - lineh)
  flow_adj <- ebot
  flow_adj[etop < -1E-6] <- flow_adj[etop < -1E-6] / 2 + etop[etop < -1E-6] / 2
  
  charh <- rep(strheight("A", cex = cex), length(new_labs)-1)
  yadj <- -cumsum(c(0, charh + wsh)) + lineh/2
  
  if(!is.na(boxh)){
    yadj[-1] <- yadj[-1] - (boxh + tail(yadj,1) + yadj[2]) / 5
  }
  
  for(i in 1:length(new_labs)){
    # print(ifelse(i==1 & !is.na(first_line_hadj), paste0(new_labs[[i]], ": ", (first_line_hadj - strwidth(new_labs[[i]], cex = cex)) / 2), 0))
    if(!is.na(first_line_center_loc) & i == 1){
      text(x = first_line_center_loc, 
           y = y + yadj[i] + flow_adj[i], labels = new_labs[[i]], cex = cex, col = first_line_col,
           pos = NULL, ...)
    } else {
      text(x = x + ifelse(i==1 & !is.na(first_line_hadj), (first_line_hadj - strwidth(new_labs[[i]], cex = cex)) / 2, 0), 
           y = y + yadj[i] + flow_adj[i], labels = new_labs[[i]], cex = cex, col = ifelse(i==1, first_line_col, col),
           pos = pos, ...)  
    }
    
    # points(x = x, y = y+yadj[i]-lineh[i]/2)
    # segments(x0 = x, x1 = x + strwidth(new_labs[[i]], cex = cex) * 1.5, y0 = y+yadj[i]-lineh[i]/2+0.0, y1 = y+yadj[i]-lineh[i]/2+0.0)
    # rect(xleft = x, xright = x + strwidth(new_labs[[i]], cex = cex), ybottom = y+yadj[i]-lineh[i]/2, ytop = y+yadj[i]+lineh[i]/2)
  }
}

text2 <- function(x, y, pos = NULL, cex = 1, labels = NULL, drect = F, ...){
  adj_x <- x + ifelse(any(pos %in% c(2,4)), ifelse(any(pos == 2), -1, 1) * strwidth(labels, cex = cex) / 2, 0)
  adj_y <- y + ifelse(any(pos %in% c(1,3)), ifelse(any(pos == 1), -1, 1) * strheight(labels, cex = cex) / 2, 0)
  text(x = adj_x, y = adj_y, labels = labels, pos = NULL, cex = cex, ...)
  if(drect){
    rect(xleft = adj_x - strwidth(labels, cex = cex) / 2, 
         xright = adj_x + strwidth(labels, cex = cex) / 2, 
         ybottom = adj_y - strheight(labels, cex = cex) / 2, 
         ytop = adj_y + strheight(labels, cex = cex) / 2)
    # abline(h = adj_y - strheight(labels, cex = cex) / 2, lwd = 0.5)
  }
}


remove_bottom <- function(x, replacement){
  nobot <- gsub("g|j|p|q|y|,|\\(|\\)|Q", replacement, x)
  nobot <- gsub("\\_s*\\{[^\\)]+\\}", replacement, nobot) #underscore in brackets
  nobot <- gsub("_[a-z|0-9|A-Z]{1}", replacement, nobot) #underscore w/ one letter following
  nobot
}

remove_top <- function(x, replacement){
  notop <- gsub("\\^s*\\{[^\\)]+\\}", replacement, x)
  notop <- gsub("\\^[a-z|0-9|A-Z]{1}", replacement, notop)
  notop
}

remove_tb <- function(x, replacement){
  remove_top(remove_bottom(x, replacement), replacement)
}

text3 <- function(x, y, pos = NULL, cex = 1, labels = NULL, drect = F, col = 1, replacement = "a", ...){
  
  #convert text label to expression
  if(all(class(labels) == "character")){
    word_expression <- correct_l2xp_vec(labels)  
  } else {
    word_expression <- labels
  }
  
  #find general params
  strw <- strwidth(word_expression, cex = cex)
  strh <- strheight(word_expression, cex = cex)
  base_strh <- strheight(correct_l2xp_vec("GIs"), cex = cex)
  
  #adjust base location
  adj_x <- x + ifelse(any(pos %in% c(2,4)), ifelse(any(pos == 2), -1, 1) * strw / 2, 0)
  adj_y <- y + ifelse(any(pos %in% c(1,3)), ifelse(any(pos == 1), -1, 1) * strh / 2, 0)
  
  #adjust in case of ledding
  nobot <- remove_bottom(labels, replacement)
  ebot <- strheight(latex2exp::TeX(nobot), cex = cex) - strh
  
  notop <- remove_top(labels, replacement)
  etop <- strheight(latex2exp::TeX(notop), cex = cex) - strh
  
  nobottop <- remove_tb(labels, replacement)
  ebottop <- strheight(latex2exp::TeX(nobottop), cex = cex) - strh
  
  #ugh this was obnoxious to figure out
  ebt_delta <- ebottop - (ebot + etop)
  adj_ledding <- ifelse(abs(ebt_delta) > 1E-6, 
                        ebot / 2 - (ebottop - ebot) / 2, 
                        ebot / 2 - etop / 2)
  adj_ledding <- adj_ledding - ifelse(base_strh > strh, (base_strh - strh) / 2, 0)
  
  #print the text itself
  text(x = adj_x, y = adj_y + adj_ledding, labels = word_expression, pos = NULL, cex = cex, col = col, ...)
  
  #draw a box around it if desired
  if(drect){
    rect(xleft = adj_x - strw / 2, 
         xright = adj_x + strw / 2, 
         ybottom = adj_y - strh / 2 + adj_ledding, 
         ytop = adj_y + strh / 2 + adj_ledding, border = col)
    abline(h=y - strheight(latex2exp::TeX("GIs"), cex = cex) / 2, lwd = 0.5)
  }
}


text_cols <- function(string, cols, x, y, cex = 1, ...){
  for(char_i in 1:nchar(string)){
    txt_exploded <- c(substr(string, 1, char_i-1), substr(string, char_i, char_i), substr(string, char_i+1, nchar(string)))
    text2(x = x, y = y, labels = bquote(phantom(.(txt_exploded[1])) * .(txt_exploded[2]) * phantom(.(txt_exploded[3]))), col = cols[char_i], cex = cex, ...)
  }
}

find_optimal_cex_and_lines <- function(txt, rect_coords, rect_rescaling_ratio = 1, str_height_rescaling_ratio = 1, fixed_cex = NA,
                                       newlines = NA){
  
  strwidths <- strwidth(correct_l2xp_vec(txt))
  strheight <- (strheight("M\nM") - strheight("M")) * str_height_rescaling_ratio
  space_width_min <- strwidth(" ")
  rectwidth <- abs(rect_coords$x1 - rect_coords$x0) * rect_rescaling_ratio
  rectwidth <- abs(rect_coords$x1 - rect_coords$x0) * rect_rescaling_ratio
  rectheight <- abs(rect_coords$y1 - rect_coords$y0) * rect_rescaling_ratio
  
  if(all(is.na(newlines))){
    newlines <- rep(F, length(txt))
  }
  
  # ceiling(cumsum(strwidths) / rectwidth)
  data <- list(strwidths = strwidths, strheight = strheight, space_width_min = space_width_min, rectwidth = rectwidth, rectheight = rectheight,
               newlines = newlines)
  par <- log((rectwidth * rectheight) / ((sum(strwidths) + space_width_min * length(strwidths)) * strheight) * 0.5) #initialize cex
  while(compute_space_remaining(data, par) == Inf){
    par <- par + log(0.95)
  }
  # plot(1:120/100, sapply(log(1:120/100), function(cex) compute_space_remaining(data, cex)), type = "l")
  if(is.na(fixed_cex)){
    opt_cex <- suppressWarnings(optimx::optimx(par = par, fn = compute_space_remaining, data = data, hess = NULL, 
                                               method = c('nlminb', 'nlm', 'BFGS', 'Nelder-Mead')[4], hessian=FALSE, #can't compute hessian bc of sharp jumps when new line is formed? or maybe not?
                                               control = list(maxit = 1E4, trace = 0, kkt=FALSE)))
    
    #check if it worked
    compute_space_remaining(data = data, par = opt_cex$p1)
    compute_space_remaining(data = data, par = opt_cex$p1 + 0.0001)
    
    return(list(cex = exp(opt_cex$p1), 
                words_on_lines = put_words_on_lines(data = data, par = opt_cex$p1),
                space_width_min = space_width_min * exp(opt_cex$p1),
                vertical_space = strheight * exp(opt_cex$p1),
                vertical_space_noWS = strheight("M") * exp(opt_cex$p1))
    )
  } else {
    return(list(cex = fixed_cex, 
                words_on_lines = put_words_on_lines(data = data, par = fixed_cex),
                space_width_min = space_width_min * fixed_cex,
                vertical_space = strheight * fixed_cex,
                vertical_space_noWS = strheight("M") * fixed_cex)
    )
  }
  
  
}

compute_space_remaining <- function(data, par){
  put_words_on_lines(data, par, return_space_remaining = T)
}

put_words_on_lines <- function(data, par, return_space_remaining = F){
  
  #clarify par-cex relationship
  cex <- exp(par)
  
  #unwrap data
  strwidths_mod <- data$strwidths * cex
  strheight_mod <- data$strheight * cex
  space_width_min_mod <- data$space_width_min * cex
  rectwidth_mod <- data$rectwidth
  rectheight_mod <- data$rectheight
  newlines <- data$newlines
  
  #check that no words are wider than a line
  if(any(strwidths_mod > rectwidth_mod) & return_space_remaining){
    return(Inf)
  }
  
  txi <- 1
  linei <- 1
  current_width <- strwidths_mod[txi] + space_width_min_mod #if first line has a line break, need to clear current_width
  txt_lines <- list()
  txt_lines[[1]] <- txi
  txt_lines[[2]] <- integer(0)
  
  #check if first word has a line break
  linei <- linei + newlines[txi]
  current_width <- ifelse(newlines[txi], 0, current_width)
  
  while(txi < length(txt)){
    
    txi <- txi + 1
    
    prop_current_width <- current_width + strwidths_mod[txi]
    add_to_line <- (prop_current_width < rectwidth_mod) | (length(txt_lines[[linei]]) == 0)
    
    if(add_to_line){
      txt_lines[[linei]] <- c(txt_lines[[linei]], txi)
      current_width <- prop_current_width + space_width_min_mod
    } else {
      current_width <- strwidths_mod[txi] + space_width_min_mod
      linei <- linei + 1
      txt_lines[[linei]] <- txi
    }
    
    if(newlines[txi]){
      linei <- linei + 1
      current_width <- 0
      txt_lines[[linei]] <- integer(0)
    }
    
    #print for debugging
    # print(paste0("line: ", linei, ", txt_i: ", txi, ", txt: ", txt[txi], ", nl: ", newlines[txi], ", curr_w: ", round(current_width, 4)))
    
  }
  
  if(return_space_remaining){
    last_line_width_remaining <- rectwidth_mod - sum(strwidths_mod[txt_lines[[linei]]])
    current_height <- linei * strheight_mod
    vspace_remaining <- (rectheight_mod - current_height) * rectwidth_mod
    hspace_remaining <- (last_line_width_remaining * strheight_mod)
    
    if(vspace_remaining < 0 | hspace_remaining < 0){
      return(Inf)
    } else {
      total_space_remaining <- vspace_remaining + hspace_remaining
      return(total_space_remaining)
    }
  } else {
    return(txt_lines)
  }
  
}

text_wrapped_words <- function(txt, rect_coords, optimal_word_placement_inf, justified = F, str_height_lower_start_ratio = 0, 
                               str_width_lefter_start_ratio = 0, rect_rescaling_ratio = 1, col = "black", multicolor_words = F, cols_list,
                               vertically_justified = F,...){
  
  ws_height <- optimal_word_placement_inf$vertical_space - optimal_word_placement_inf$vertical_space_noWS
  cex <- optimal_word_placement_inf$cex
  curr_x <- rect_coords$x0 - abs(rect_coords$x0 - rect_coords$x1) * str_width_lefter_start_ratio
  curr_y <- rect_coords$y0 - optimal_word_placement_inf$vertical_space_noWS * (1 + str_height_lower_start_ratio) / 2
  nlines <- length(optimal_word_placement_inf$words_on_lines)
  strwidths_plotting <- strwidth(correct_l2xp_vec(txt), cex = cex)
  
  space_left_on_lines <- sapply(1:length(optimal_word_placement_inf$words_on_lines), function(linei)
    abs(rect_coords$x0 - rect_coords$x1) * rect_rescaling_ratio - sum(strwidths_plotting[optimal_word_placement_inf$words_on_lines[[linei]]]))
  justified_space_between_words <- sapply(1:length(optimal_word_placement_inf$words_on_lines), function(linei)
    space_left_on_lines[linei] / (length(optimal_word_placement_inf$words_on_lines[[linei]]) - 1))
  
  if(vertically_justified){
    space_taken_vertically <- optimal_word_placement_inf$vertical_space_noWS * nlines + 
      (optimal_word_placement_inf$vertical_space - optimal_word_placement_inf$vertical_space_noWS) * (nlines - 1)
    space_left_vertically <- abs(rect_coords$y0 - rect_coords$y1) - space_taken_vertically
    extra_leading <- space_left_vertically / (nlines - 1)
  } else {
    extra_leading <- 0
  }
  
  words_written <- 0
  for(linei in 1:nlines){
    for(wordi in 1:length(optimal_word_placement_inf$words_on_lines[[linei]])){
      
      words_written <- words_written + 1
      word_to_write <- txt[optimal_word_placement_inf$words_on_lines[[linei]][wordi]]
      
      #adjust for sticky-outy bits
      word_expression <- correct_l2xp_vec(word_to_write)
      
      
      if(multicolor_words){
        text_cols(x = curr_x, y = curr_y, cex = cex,
                  string = word_expression, pos = 4, cols = cols_list[[words_written]])
      } else {
        # text2(x = curr_x, y = curr_y, cex = cex,
        #       labels = word_expression, pos = 4, col = col, drect = T)
        text3(x = curr_x, y = curr_y, cex = cex,
              labels = word_to_write, pos = 4, col = col, drect = T)
        abline(h=curr_y - strheight(correct_l2xp_vec("GIs"), cex = cex) / 2, lwd = 0.5)
      }
      if(justified & (linei != nlines)){
        curr_x <- curr_x + strwidth(word_expression) * cex + justified_space_between_words[linei]
      } else {
        curr_x <- curr_x + strwidth(word_expression) * cex + optimal_word_placement_inf$space_width_min  
      }
      
    }
    curr_x <- rect_coords$x0 - abs(rect_coords$x0 - rect_coords$x1) * str_width_lefter_start_ratio
    curr_y <- curr_y - optimal_word_placement_inf$vertical_space * (1 + str_height_lower_start_ratio) - extra_leading
  }
  
}

swap_one <- function(x, y_i){
  list(append(x[[1]], values = y_i[[1]], after = y_i[[2]] + x[[2]])[-(y_i[[2]] + x[[2]])], x[[2]] + length(y_i[[1]]) - 1)
}

swap <- function(x = vector(), y = list(), inds = vector()){
  #swaps y items into vector x for elements of x at locations given by inds
  x <- list(x, 0)
  y <- y[order(inds)]
  inds <- sort(inds)
  out <- Reduce(f = swap_one, x = c(list(x), lapply(seq_along(y), function(i) list(y[[i]], inds[i]))))
  out[[1]]
}


string_to_tokens <- function(txt_string){
  
  #get basic string
  txt <- strsplit(txt_string, split = " ")[[1]]
  
  #first adjust for new characters
  has_newlines <- grep(pattern = "\n", x = txt)
  newline_swap <- lapply(txt[has_newlines], function(x){
    y <- strsplit(x, "\n")[[1]]
    num_nl <- lengths(regmatches(x, gregexpr("\n", x)))
    y[1:num_nl] <- paste0(y[1:num_nl], "\n")
    y
  })
  txt <- swap(txt, newline_swap, has_newlines)
  
  #find relevant LaTeX notation & modify tokens accordingly
  math_pairs <- do.call(rbind, lapply(seq_along(txt), function(i){
    mathi = which(strsplit(txt[i], "")[[1]] == "$")
    if(length(mathi) != 0){
      return(data.frame(token = i, math_i = mathi))
    } else {
      return(integer(0))
    }
  }))
  if(length(math_pairs) != 0){
    math_pairs <- data.frame(cbind(matrix(math_pairs$token, ncol = 2, byrow = T), 
                                   matrix(math_pairs$math_i, ncol = 2, byrow = T)))
    colnames(math_pairs) <- c("open_i", "close_i", "char_i_open", "char_i_close")
    for(i in 1:nrow(math_pairs)){
      if(math_pairs$open_i[i] == math_pairs$close_i[i]){
        next()
      } else {
        txt[math_pairs$open_i[i]] <- paste0(txt[math_pairs$open_i[i]:math_pairs$close_i[i]], collapse = " ")
        txt <- txt[-((math_pairs$open_i[i]+1):math_pairs$close_i[i])]
      }
    }
  }
  
  
  bracket_pairs <- list(
    do.call(rbind, lapply(seq_along(txt), function(i){
      openi <- which(strsplit(txt[i], "")[[1]] == "{")
      if(length(openi) != 0){
        slashi <- which(sapply(strsplit(txt[i], "")[[1]], grepl, pattern = '\t'))
        if(length(slashi) == 0){
          return(data.frame(token_open = i, 
                            open_i = openi,
                            slash_i = "",
                            style = ""))
        } else {
          styles <- sapply(seq_along(slashi), function(ci) substr(txt[i], slashi[ci], openi[ci]))
          return(data.frame(token_open = i, 
                            open_i = openi,
                            slash_i = slashi,
                            style = styles))  
        }
      } else{
        return(integer(0))
      }
    })
    ),
    
    do.call(rbind, lapply(seq_along(txt), function(i){
      closei <- which(strsplit(txt[i], "")[[1]] == "}")
      if(length(closei) != 0){
        return(data.frame(token_close = i, 
                          close_i = closei))
      } else{
        return(integer(0))
      }
    })
    )
  )
  
  bracket_pos <- rbind(data.frame(type = "o", 
                                  pos = bracket_pairs[[1]]$token_open + 
                                    bracket_pairs[[1]]$open_i / 
                                    nchar(txt[bracket_pairs[[1]]$token_open]),
                                  ind = 1:nrow(bracket_pairs[[1]])),
                       data.frame(type = "c", 
                                  pos = bracket_pairs[[2]]$token_close + 
                                    bracket_pairs[[2]]$close_i / 
                                    nchar(txt[bracket_pairs[[2]]$token_close]),
                                  ind = 1:nrow(bracket_pairs[[2]]))
  )
  bracket_pos <- bracket_pos[order(bracket_pos$pos), c("type", "ind")]
  bracket_pos$score <- as.numeric(bracket_pos$type == "o") - as.numeric(bracket_pos$type == "c")
  bracket_pos$cumscore <- cumsum(bracket_pos$score)
  
  bracket_match <- data.frame(do.call(rbind, lapply(1:max(bracket_pos$cumscore), function(mcs){
    cbind(o = bracket_pos$ind[which(bracket_pos$type == "o" & bracket_pos$cumscore == mcs)], 
          c = bracket_pos$ind[which(bracket_pos$type == "c" & bracket_pos$cumscore == (mcs - 1))]) 
  })))
  bracket_match <- bracket_match[order(bracket_match$o),]
  
  bracket_pairs <- do.call(rbind, apply(bracket_match, 1, function(i) cbind(bracket_pairs[[1]][i[1],], bracket_pairs[[2]][i[2],])))
  
  #adjust tokens to reflect text modifications
  for(i in 1:nrow(bracket_pairs)){
    if(bracket_pairs$token_open[i] == bracket_pairs$token_close[i]){
      next()
    } else {
      inds <- bracket_pairs$token_open[i]:bracket_pairs$token_close[i]
      txt[inds] <- sapply(inds, function(j){
        new_txt <- txt[j]
        if(j != bracket_pairs$token_open[i]){
          new_txt <- paste0(bracket_pairs$style[i], new_txt)
        }
        if(j != bracket_pairs$token_close[i]){
          new_txt <- paste0(new_txt, "}")
        }
        new_txt
      })
    }
  }
  
  txt <- gsub(pattern = "\t", replacement = "\\t", x = txt, fixed = T)
  
  list(tokens = gsub(x = txt, pattern = "\n", replacement = ""),
       newlines = grepl(txt, pattern = "\n"))
  
}

correct_l2xp_vec <- function(x){
  # latex2exp::TeX(x) <- this doesn't handle bolditalic correctly
  empties <- which(x == "")
  x[empties] <- "placeholder"
  out <- lapply(x, latex2exp::TeX)
  out_str <- sapply(seq_along(out), function(i) as.character(out[[i]]))
  bis <- intersect(grep(pattern = "bold", out_str), grep(pattern = "italic", out_str))
  new_out_str <- sapply(out_str[bis], function(i) paste0("bolditalic(\"", strsplit(i, "\"")[[1]][2], "\")"))
  new_out_str <- swap(out_str, new_out_str, bis)
  new_out <- sapply(new_out_str, function(i) parse(text = i))
  new_out[empties] <- ""
  new_out
}

NULL


