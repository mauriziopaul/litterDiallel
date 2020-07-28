#' litterDiallel: A package for reproducing the analysis used in the litter diallel manuscript.
#'
#' A package for reproducing the analysis used in the litter diallel manuscript.
#' 
#' @docType package

#' @name litterDiallel
NULL

#' @section require namespaces:
requireNamespace("MCMCglmm", quietly=TRUE)

#' @title incidence.matrix: Make an incidence matrix
#' @description Convert a factor into an incidence matrix. (From the WVmisc package, Will Valdar).
#' @param fact factor
#' @param ... additional arguments
#' @return generates an incidence matrix
#' @examples
#' ## not run
#' @export
incidence.matrix <- function(fact, ...){
    m = diag(nlevels(fact))[fact, ]
    colnames(m) = levels(fact)
    m
}

#' @title plot_hpd: Plot highest posterior density intervals
#' @description Plot HPD intervals. (From the BayesDiallel package, Will Valdar and Alan Lenarcic).
#' @param coda.object coda object
#' @param wanted variable names for coda object
#' @param prob.wide outer width of posterior probability
#' @param prob.narrow inner width of posterior probability
#' @param xlab x-axis label
#' @param names names
#' @param type type of plot
#' @param name.margin margin for name
#' @param plt.left left of plot margin
#' @param plt.right right of plot margin
#' @param plt.bottom bottom of plot margin
#' @param plt.title title of plot margin
#' @param ylab y-axis label
#' @param name.line line where names should be printed
#' @param main text of main title
#' @param main.line line where main title should be printed
#' @param ... additional arguments
#' @return returns HPD plot
#' @examples
#' ## not run
#' @export plot_hpd
plot_hpd <- function(coda.object,
    wanted=coda::varnames(coda.object),
    prob.wide=0.95,
    prob.narrow=0.50,
    xlab="HPD interval",
    names=NULL,
    type="p",
    name.margin = 6.1,
    plt.left=NULL, plt.right=NULL, plt.bottom=NULL, plt.title=NULL,
    ylab="",  name.line = 3.9, main="", main.line=2,
    ...)
{
  which.wanted=ifow(is.integer(wanted), wanted, match(wanted, coda::varnames(coda.object)))
  num.wanted=length(which.wanted)
  if(!exists("name.margin") || is.null(name.margin)) { name.margin = 6.1; }
  chain <- mcmc.stack(coda.object)
  mu    <- colMeans(chain[,which.wanted])
  med   <- apply(coda::HPDinterval(chain, prob=0.01)[which.wanted,],
      1, mean)
  hpd.wide    <- coda::HPDinterval(chain, prob=prob.wide)[which.wanted,]
  hpd.narrow  <- coda::HPDinterval(chain, prob=prob.narrow)[which.wanted,]

    mid.vals <- med;
  if (is.null(names)) names <- coda::varnames(chain)[which.wanted]
  else names <- rep(names, length.out=length(wanted))
  ypos <- plot_ci(med, hpd.narrow, hpd.wide, names=names, xlab=xlab, col.midvals="white", pch.midvals="|", type=type,
    name.margin=name.margin,plt.left=plt.left, plt.right=plt.right,
    plt.bottom=plt.bottom, plt.title=plt.title, ylab=ylab,
    name.line = name.line, main=main, main.line=main.line, ...)
  if ("p"==type)
  {
    graphics::points(mu, ypos, pch="|")
  }
  invisible(ypos)
}

# Source: WVmisc package
#' @title ifow
#' @description if, otherwise
#' @param test test
#' @param yes specify yes
#' @param no specify no
#' @return returns yes if test is true, no if false
#' @examples
#' ## not run
#' @export ifow
ifow <- function(test, yes, no)
{
  if (test)
  {
    return (yes)
  }
  no
}

# Source: WVgraphics package
#' @title sides
#' @description define sides of a plot
#' @param default specify overall default
#' @param bottom specify bottom
#' @param left specify left
#' @param top specify top
#' @param right specify default
#' @return returns named side arguments
#' @examples
#' ## not run
#' @export sides
sides <- function(default=NA, bottom=default, left=default, top=default, right=default)
{
  x=c(bottom, left, top, right)
  names(x)=c("bottom", "left", "top", "right")
  x
}

# Source: WVgraphics package
#' @title update_sides
#' @description update sides of a plot
#' @param old defaults to par("mar")
#' @param new defaults to 4 NA's
#' @return returns updated sides
#' @examples
#' ## not run
#' @export update_sides
update_sides=function(old=graphics::par("mar"), new=rep(NA, 4))
{
  old[!is.na(new)]=new[!is.na(new)]
  old
}

# Source: BayesDiallel package
#' @title mcmc.stack
#' @description stack multiple mcmc chains into one mcmc object
#' @param coda.object coda object with multiple chains
#' @param ... additional arguments
#' @return returns stacked mcmc
#' @examples
#' ## not run
#' @export mcmc.stack
mcmc.stack <- function (coda.object, ...){
  ## This function is from Will; also part of BayesDiallel
    if (inherits(coda.object, "mcmc")) {
        return(coda.object)
    }
    if (!inherits(coda.object, "mcmc.list")) {
        stop("Non-mcmc object passed to function\n")
    }
    chain <- coda.object[[1]]
    for (i in 2:coda::nchain(coda.object)) {
        chain <- rbind(chain, coda.object[[i]])
    }
    coda::as.mcmc(chain)
}

# Source: BayesDiallel package
#' @title plot_ci
#' @description stack multiple mcmc chains into one mcmc object
#' @param midvals midvals
#' @param narrow.intervals narrow.intervals
#' @param wide.intervals wide.intervals
#' @param names names
#' @param add add
#' @param main main
#' @param main.line main.line
#' @param xlab xlab
#' @param xlab.line xlab.line
#' @param xlim xlim
#' @param ylab ylab
#' @param yaxis yaxis
#' @param ylim ylim
#' @param name.line name.line
#' @param pch.midvals pch.midvals
#' @param col col
#' @param col.midvals col.midvals
#' @param cex.labels cex.labels
#' @param type type
#' @param name.margin name.margin
#' @param title.margin title.margin
#' @param title.line title.line
#' @param bottom.margin bottom.margin
#' @param bottom.line bottom.line
#' @param right.margin right.margin
#' @param right.line right.line
#' @param mar mar
#' @param mar.update mar.update
#' @param before.data before.data
#' @param plt.left plt.left
#' @param plt.right plt.right
#' @param plt.bottom plt.bottom
#' @param plt.title plt.title
#' @param ... additional arguments
#' @return returns plotted caterpillar plot with confidence intervals
#' @examples
#' ## not run
#' @export plot_ci
plot_ci <- function(midvals, narrow.intervals, wide.intervals,
    names=1:length(midvals),
    add=FALSE,
    main="", main.line=2,
    xlab="Estimate",
    xlab.line=2.5,
    xlim=NULL,
    ylab="",
    yaxis=TRUE,
    ylim=c(0, length(midvals)),
    name.line=4,
    pch.midvals=19,
    col="black",
    col.midvals=col,
    cex.labels=1,
    type="p",
    name.margin=6.1,
    title.margin=4.1, title.line = 3.5,
    bottom.margin=5.1, bottom.line=4.5,
    right.margin=2.1, right.line=1.5,
    mar=sides(left=name.margin, bottom=bottom.margin, top=title.margin, right=right.margin),
    mar.update=sides(),
    before.data=function(){},
    plt.left=NULL, plt.right=NULL, plt.bottom=NULL, plt.title=NULL,
    ...)
# Example: plot_ci( c(0,10), narrow.intervals=rbind(c(-1,1), c(8,12)), wide.intervals=rbind(c(-3,4), c(5,15)), names=c("Fred", "Barney"))
{
  nvals <- length(midvals)
  col.midvals <- rep(col.midvals, length.out=nvals)
  y.pos <- (1:nvals)-0.5
  if (!add)
  {
    if (is.null(xlim))
    {
      xlim <- range(c(wide.intervals,narrow.intervals,midvals), na.rm=TRUE)
      xlim <- c(-1,1) * diff(xlim)*0.1 + xlim
    }

    if (name.margin == 6.1 && !is.null(plt.left)  && is.numeric(plt.left) &&
      plt.left >= 0 && plt.left <= 1.0 ) {
      name.margin=6.1 * plt.left / .2;
    }
    if (right.margin == 2.1 && !is.null(plt.right)  && is.numeric(plt.right) &&
      plt.right >= 0 && plt.right <= 1.0 ) {
      right.margin=2.1 * (1.0-plt.right) / .05;
    }
    if (title.margin == 4.1 && !is.null(plt.title)  && is.numeric(plt.title) &&
      plt.title >= 0 && plt.title <= 1.0 ) {
      title.margin=4.1* (1.0-plt.title) / .12;
    }
    if (bottom.margin == 5.1 && !is.null(plt.bottom)  && is.numeric(plt.bottom) &&
      plt.bottom >= 0 && plt.bottom <= 1.0 ) {
      bottom.margin=5.1* plt.bottom / .16;
    }
    mar <- c(bottom.margin, name.margin, title.margin, right.margin)+0.1
    mar=update_sides(mar, mar.update)
    oldmar <- graphics::par(mar=mar); on.exit(graphics::par(mar=oldmar))

    MyD = FALSE;
    AT = "plot(x = xlim, y=ylim, type=\"n\", axes=FALSE, ylab=ylab, xlim=xlim, ylim=ylim, xlab=\"\", main=\"\",
       ...); MyD = TRUE";
    try(eval(parse(text=AT)), silent=TRUE);
    if (MyD == FALSE) {
      try(plot(x=xlim, y=ylim, type="n", axes=FALSE,ylab=ylab, ylim=ylim, xlim=xlim, xlab="", main="", ));
    }

    if (!is.null(main) && is.character(main) && main[1] != "") {
      try(graphics::title(main=main, line=main.line, cex.main = 1.5));
    }
    graphics::title(xlab=xlab, line=xlab.line)
    graphics::axis(1)
    graphics::axis(3, line=-.8)
    if (yaxis)
    {
      graphics::axis(2, at=y.pos, labels=rev(names), las=1, lty=0, hadj=0, line=name.line, cex.axis=cex.labels)
    }
  }
  before.data()
  if ("p"==type)
  {
    for (i in 1:nvals)
    {
      pos <- nvals-i + 0.5
      graphics::lines(wide.intervals[i,], rep(pos,2))
      graphics::lines(narrow.intervals[i,], rep(pos,2), lwd=3)
      graphics::points(midvals[i], pos, pch=pch.midvals, col=col.midvals[i])
    }
  }
  invisible(rev(y.pos))
}

#' @title makeRotationMatrix: Make a rotation matrix
#' @description Turn an n-column design matrix into an n-1, sum to 0, design matrix
#'              while maintaining independence.
#'              See Appendix (Lenarcic, 2012, Genetics; Crowley, 2014, Genetics)
#' @param X original design matrix
#' @param n number of original columns of design matrix
#' @param ... additional arguments
#' @return generates rotation matrix to reduce dimensions for better sampling
#' @examples
#' ## not run
#' @export
makeRotationMatrix <- function(X, n, ...)
{
 j <- ncol(X)
 k <- (-1 + sqrt(j))*(j - 1)^(-3/2)
 m <- 1/sqrt(j - 1)
 c <- (n - 2)*k + m
 M <- diag(j - 1)
 M[0 == M] <- -k
 M <- rbind(M, rep(-m, j - 1))
 return(M)
}

#' @title diallelMatrixMaker
#' @description Make design matrices for diallel
#' @param data data frame
#' @param dam.col.name dam column name
#' @param sire.col.name sire column name
#' @param batch.col.name name of batch/random effect column
#' @param batch.1.col.name name of additional batch/random effect column
#' @param strains character vector with names of inbred parental strains
#' @param ... additional arguments
#' @return returns diallel incidence matrices
#' @examples
#' ## not run
#' @export
diallelMatrixMaker <- function(data, dam.col.name, sire.col.name, batch.col.name = NULL, batch.1.col.name = NULL,
    strains=c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB"), ...){
  dam.mat <- incidence.matrix(as.factor(data[, as.character(dam.col.name)]))[,strains]
  sire.mat <- incidence.matrix(as.factor(data[, as.character(sire.col.name)]))[,strains]
  add.mat <- dam.mat + sire.mat
  mat.mat <- dam.mat - sire.mat
  colnames(dam.mat) <- paste0("dam:", strains)
  colnames(sire.mat) <- paste0("sire:", strains)
  colnames(add.mat) <- paste0("additive:", strains)
  colnames(mat.mat) <- paste0("maternal:", strains)
  # names according to lower triangular matrix, down first row, down second, etc.
  cross.n <- length(strains)
  jk.names <- NULL

  for(k in 1:cross.n){
    for(j in 1:cross.n){
      if(j > k){
        jk.names <- c(jk.names, paste(strains[j], strains[k], sep=";"))
      }
    }
  }

  jk.bind <- cbind.data.frame(data[, as.character(dam.col.name)],
                               data[, as.character(sire.col.name)])
  jk <- apply(X=jk.bind, MARGIN=1, FUN=function(x){ifelse(paste(x, collapse=";") %in% jk.names, paste(x, collapse=";"), paste(rev(x), collapse=";"))})
  jk.asymm <- apply(X=jk.bind, MARGIN=1, FUN=function(x){paste(rev(x), collapse=";")})
  asymm <- apply(cbind.data.frame(jk, jk.asymm), 1, function(x) {
    ifelse(x[1] == x[2], -1, 1)})
  jk.mat <- incidence.matrix(as.factor(jk))
  drops <- paste(strains, strains, sep=";")
  data$inbred <- ifelse(data[, as.character(dam.col.name)] ==
                          data[, as.character(sire.col.name)], 1, 0)
  inbred.mat <- diag(data$inbred) %*% jk.mat
  jk.mat <- jk.mat[, !(colnames(jk.mat) %in% drops)]
  inbred.mat <- inbred.mat[, colnames(inbred.mat) %in% drops]
  asymm.mat <- diag(asymm) %*% jk.mat
  colnames(jk.mat) <- paste0("v:", colnames(jk.mat))
  colnames(asymm.mat) <- paste0("w:", colnames(asymm.mat))
  colnames(inbred.mat) <- paste0("inbred:", colnames(inbred.mat))
  if (is.null(batch.col.name) & is.null(batch.1.col.name)) {
    return(list(dam.mat = dam.mat, sire.mat = sire.mat, add.mat = add.mat,
                mat.mat = mat.mat, inbred.mat = inbred.mat, jk.mat = jk.mat,
                asymm.mat = asymm.mat))
  }
  else {
    if(1==length(batch.col.name)){
      batch.mat <- incidence.matrix(as.factor(data[, as.character(batch.col.name)]))
      batch.1.mat <- NULL
      try(batch.1.mat <- incidence.matrix(as.factor(data[,
                                                         as.character(batch.1.col.name)])))
      return(list(dam.mat = dam.mat, sire.mat = sire.mat,
                  add.mat = add.mat, mat.mat = mat.mat, inbred.mat = inbred.mat,
                  jk.mat = jk.mat, asymm.mat = asymm.mat, batch.mat = batch.mat,
                  batch.1.mat = batch.1.mat))
    }
    else{
      stop("Not implemented for length(batch.col.name)>1")
    }
  }
}

#' @title diallelMatrixMaker_av
#' @description Make design matrices for diallel, model (a, v)
#' @param data data frame
#' @param dam.col.name dam column name
#' @param sire.col.name sire column name
#' @param batch.col.name name of batch/random effect column
#' @param batch.1.col.name name of additional batch/random effect column
#' @param strains character vector with names of inbred parental strains
#' @param ... additional arguments
#' @return returns diallel incidence matrices
#' @examples
#' ## not run
#' @export
diallelMatrixMaker_av <- function(data, dam.col.name, sire.col.name, batch.col.name = NULL, batch.1.col.name = NULL,
                               strains=c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB"), ...){
  dam.mat <- incidence.matrix(as.factor(data[, as.character(dam.col.name)]))[,strains]
  sire.mat <- incidence.matrix(as.factor(data[, as.character(sire.col.name)]))[,strains]
  add.mat <- dam.mat + sire.mat
  colnames(add.mat) <- paste0("additive:", strains)
  # names according to lower triangular matrix, down first row, down second, etc.
  cross.n <- length(strains)
  jk.names <- NULL

  for(k in 1:cross.n){
    for(j in 1:cross.n){
      if(j > k){
        jk.names <- c(jk.names, paste(strains[j], strains[k], sep=";"))
      }
    }
  }

  jk.bind <- cbind.data.frame(data[, as.character(dam.col.name)],
                              data[, as.character(sire.col.name)])
  jk <- apply(X=jk.bind, MARGIN=1, FUN=function(x){ifelse(paste(x, collapse=";") %in% jk.names, paste(x, collapse=";"), paste(rev(x), collapse=";"))})
  jk.asymm <- apply(X=jk.bind, MARGIN=1, FUN=function(x){paste(rev(x), collapse=";")})
  jk.mat <- incidence.matrix(as.factor(jk))
  drops <- paste(strains, strains, sep=";")
  jk.mat <- jk.mat[, !(colnames(jk.mat) %in% drops)]
  colnames(jk.mat) <- paste0("v:", colnames(jk.mat))
  if (is.null(batch.col.name) & is.null(batch.1.col.name)) {
    return(list(add.mat = add.mat, jk.mat = jk.mat))
  }
  else {
    if(1==length(batch.col.name)){
      batch.mat <- incidence.matrix(as.factor(data[, as.character(batch.col.name)]))
      batch.1.mat <- NULL
      try(batch.1.mat <- incidence.matrix(as.factor(data[,
                                                         as.character(batch.1.col.name)])))
      return(list(add.mat = add.mat, jk.mat = jk.mat, batch.mat = batch.mat,
                  batch.1.mat = batch.1.mat))
    }
    else{
      stop("Not implemented for length(batch.col.name)>1")
    }
  }
}

#' @title diallelMatrixMakeAndRotate
#' @description Make design matrices for diallel, rotate to n-1 space.
#' @param data data frame
#' @param dam.col.name dam column name
#' @param sire.col.name sire column name
#' @param batch.col.name name of batch/random effect column
#' @param batch.1.col.name name of additional batch/random effect column
#' @param strains character vector with names of inbred parental strains
#' @param ... additional arguments
#' @return returns diallel incidence matrices, rotated
#' @examples
#' ## not run
#' @export
diallelMatrixMakeAndRotate <- function(data, dam.col.name, sire.col.name,
                                       batch.col.name=NULL, batch.1.col.name=NULL,
                                       strains=c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB"), ...){
  matrices <- diallelMatrixMaker(data=data, dam.col.name=dam.col.name, sire.col.name=sire.col.name, batch.col.name=batch.col.name,
                                 batch.1.col.name=batch.1.col.name, strains=strains)
  n.strains <- length(strains)
  dam.mat <- matrices$dam.mat
  sire.mat <- matrices$sire.mat
  add.mat <- matrices$add.mat
  mat.mat <- matrices$mat.mat
  inbred.mat <- matrices$inbred.mat
  jk.mat <- matrices$jk.mat
  asymm.mat <- matrices$asymm.mat
  n.jk <- n.strains*(n.strains-1)/2
  M.dam <- makeRotationMatrix(X = dam.mat, n = n.strains)
  M.sire <- makeRotationMatrix(X = sire.mat, n = n.strains)
  M.add <- makeRotationMatrix(X = add.mat, n = n.strains)
  M.mat <- makeRotationMatrix(X = mat.mat, n = n.strains)
  M.inbred <- makeRotationMatrix(X = inbred.mat, n = n.strains)
  M.jk <- makeRotationMatrix(X = jk.mat, n = n.jk)
  M.asymm <- makeRotationMatrix(X = asymm.mat, n = n.jk)
  t.dam.mat <- dam.mat %*% M.dam
  t.sire.mat <- sire.mat %*% M.sire
  t.add.mat <- add.mat %*% M.add
  t.mat.mat <- mat.mat %*% M.mat
  t.inbred.mat <- inbred.mat %*% M.inbred
  t.jk.mat <- jk.mat %*% M.jk
  t.asymm.mat <- asymm.mat %*% M.asymm
  if(!is.null(batch.col.name)){
    batch.mat <- matrices$batch.mat
    M.batch <- makeRotationMatrix(X = batch.mat, n = ncol(batch.mat))
    t.batch.mat <- batch.mat %*% M.batch
    if(!is.null(batch.1.col.name)){
      batch.1.mat <- matrices$batch.1.mat
      M.batch.1 <- makeRotationMatrix(X = batch.1.mat, n = ncol(batch.1.mat))
      t.batch.1.mat <- batch.1.mat %*% M.batch.1
      return(list(RotMat=list(	M.dam=M.dam, M.sire=M.sire,
                               M.add=M.add, M.mat=M.mat, M.inbred=M.inbred,
                               M.jk=M.jk, M.asymm=M.asymm,
                               M.batch=M.batch, M.batch.1=M.batch.1),
                  DesignMat=list(t.dam.mat=t.dam.mat, t.sire.mat=t.sire.mat,
                                 t.add.mat = t.add.mat, t.mat.mat = t.mat.mat, t.inbred.mat = t.inbred.mat,
                                 t.jk.mat = t.jk.mat, t.asymm.mat = t.asymm.mat,
                                 t.batch.mat = t.batch.mat, t.batch.1.mat = t.batch.1.mat)))
    }else{
      return(list(RotMat=list(	M.dam=M.dam, M.sire=M.sire,
                               M.add=M.add, M.mat=M.mat, M.inbred=M.inbred,
                               M.jk=M.jk, M.asymm=M.asymm,
                               M.batch=M.batch),
                  DesignMat=list(t.dam.mat=t.dam.mat, t.sire.mat=t.sire.mat,
                                 t.add.mat = t.add.mat, t.mat.mat = t.mat.mat, t.inbred.mat = t.inbred.mat,
                                 t.jk.mat = t.jk.mat, t.asymm.mat = t.asymm.mat,
                                 t.batch.mat = t.batch.mat)))
    }
  }else{
    return(list(RotMat=list(	M.dam=M.dam, M.sire=M.sire,
                             M.add=M.add, M.mat=M.mat, M.inbred=M.inbred,
                             M.jk=M.jk, M.asymm=M.asymm),
                DesignMat=list(t.dam.mat=t.dam.mat, t.sire.mat=t.sire.mat,
                               t.add.mat = t.add.mat, t.mat.mat = t.mat.mat, t.inbred.mat = t.inbred.mat,
                               t.jk.mat = t.jk.mat, t.asymm.mat = t.asymm.mat)))
  }
}

#' @title diallelMatrixMakeAndRotate_av
#' @description Make design matrices for diallel, rotate to n-1 space. (a, v) model
#' @param data data frame
#' @param dam.col.name dam column name
#' @param sire.col.name sire column name
#' @param batch.col.name name of batch/random effect column
#' @param batch.1.col.name name of additional batch/random effect column
#' @param strains character vector with names of inbred parental strains
#' @param ... additional arguments
#' @return returns diallel incidence matrices, rotated
#' @examples
#' ## not run
#' @export
diallelMatrixMakeAndRotate_av <- function(data, dam.col.name, sire.col.name,
                                       batch.col.name=NULL, batch.1.col.name=NULL,
                                       strains=c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB"), ...){
  matrices <- diallelMatrixMaker_av(data=data, dam.col.name=dam.col.name, sire.col.name=sire.col.name, batch.col.name=batch.col.name,
                                 batch.1.col.name=batch.1.col.name, strains=strains)
  n.strains <- length(strains)
  add.mat <- matrices$add.mat
  jk.mat <- matrices$jk.mat
  n.jk <- n.strains*(n.strains-1)/2
  M.add <- makeRotationMatrix(X = add.mat, n = n.strains)
  M.jk <- makeRotationMatrix(X = jk.mat, n = n.jk)
  t.add.mat <- add.mat %*% M.add
  t.jk.mat <- jk.mat %*% M.jk
  if(!is.null(batch.col.name)){
    batch.mat <- matrices$batch.mat
    M.batch <- makeRotationMatrix(X = batch.mat, n = ncol(batch.mat))
    t.batch.mat <- batch.mat %*% M.batch
    if(!is.null(batch.1.col.name)){
      batch.1.mat <- matrices$batch.1.mat
      M.batch.1 <- makeRotationMatrix(X = batch.1.mat, n = ncol(batch.1.mat))
      t.batch.1.mat <- batch.1.mat %*% M.batch.1
      return(list(RotMat=list(M.add=M.add,
                               M.jk=M.jk,
                               M.batch=M.batch, M.batch.1=M.batch.1),
                  DesignMat=list(t.add.mat = t.add.mat,
                                 t.jk.mat = t.jk.mat,
                                 t.batch.mat = t.batch.mat, t.batch.1.mat = t.batch.1.mat)))
    }else{
      return(list(RotMat=list(M.add=M.add,
                               M.jk=M.jk,
                               M.batch=M.batch),
                  DesignMat=list(t.add.mat = t.add.mat,
                                 t.jk.mat = t.jk.mat,
                                 t.batch.mat = t.batch.mat)))
    }
  }else{
    return(list(RotMat=list(M.add=M.add, M.jk=M.jk),
                DesignMat=list(t.add.mat = t.add.mat, t.jk.mat = t.jk.mat)))
  }
}

#' @title diallelMatrixMakerSexed
#' @description Make design matrices for diallel, fulls model
#' @param data data frame
#' @param dam.col.name dam column name
#' @param sire.col.name sire column name
#' @param batch.col.name name of batch/random effect column
#' @param batch.1.col.name name of additional batch/random effect column
#' @param strains character vector indicating expected strain names
#' @param sex.col.name the column name for sex
#' @param female.sex.label the character string indicating female sex
#' @param ... additional arguments
#' @return returns diallel incidence matrices
#' @examples
#' ## not run
#' @export
diallelMatrixMakerSexed <- function(data, dam.col.name, sire.col.name, batch.col.name = NULL, batch.1.col.name = NULL,
                               strains=c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB"), sex.col.name, female.sex.label,
                               ...){
  dam.mat <- incidence.matrix(data[, as.character(dam.col.name)])[,strains]
  sire.mat <- incidence.matrix(data[, as.character(sire.col.name)])[,strains]
  add.mat <- dam.mat + sire.mat
  mat.mat <- dam.mat - sire.mat
  colnames(dam.mat) <- paste0("dam:", strains)
  colnames(sire.mat) <- paste0("sire:", strains)
  colnames(add.mat) <- paste0("additive:", strains)
  colnames(mat.mat) <- paste0("maternal:", strains)

  # names according to lower triangular matrix, down first row, down second, etc.
  cross.n <- length(strains)
  jk.names <- f.jk.names <- NULL

  for(k in 1:cross.n){
    for(j in 1:cross.n){
      if(j > k){
        jk.names <- c(jk.names, paste(strains[j], strains[k], sep=";"))
        f.jk.names <- c(f.jk.names, paste("f:", strains[j], strains[k], sep=";"))
      }
    }
  }

  jk.bind <- cbind.data.frame(data[, as.character(dam.col.name)],
                              data[, as.character(sire.col.name)])
  jk <- apply(X=jk.bind, MARGIN=1, FUN=function(x){ifelse(paste(x, collapse=";") %in% jk.names, paste(x, collapse=";"), paste(rev(x), collapse=";"))})
  jk.asymm <- apply(X=jk.bind, MARGIN=1, FUN=function(x){paste(rev(x), collapse=";")})
  asymm <- apply(cbind.data.frame(jk, jk.asymm), 1, function(x) {ifelse(x[1] == x[2], -1, 1)})
  jk.mat <- incidence.matrix(as.factor(jk))
  drops <- paste(strains, strains, sep=";")
  data$inbred <- ifelse(data[, as.character(dam.col.name)]==data[, as.character(sire.col.name)], 1, 0)
  inbred.mat <- diag(data$inbred) %*% jk.mat
  jk.mat <- jk.mat[, !(colnames(jk.mat) %in% drops)]
  inbred.mat <- inbred.mat[, colnames(inbred.mat) %in% drops]
  asymm.mat <- diag(asymm) %*% jk.mat
  colnames(jk.mat) <- paste0("v:", colnames(jk.mat))
  colnames(asymm.mat) <- paste0("w:", colnames(asymm.mat))
  colnames(inbred.mat) <- paste0("inbred:", colnames(inbred.mat))

  ## assumes indicator-variable for sex
  is.female <- ifelse(data[,sex.col.name]==female.sex.label,1,0)
  f.add.mat <- is.female * add.mat
  f.mat.mat <- is.female * mat.mat
  f.dam.mat <- is.female * dam.mat
  f.sire.mat <- is.female * sire.mat
  f.inbred.mat <- is.female * inbred.mat
  f.jk.mat <- is.female * jk.mat
  f.asymm.mat <- is.female * asymm.mat

  colnames(f.add.mat) <- paste0("f:", colnames(f.add.mat))
  colnames(f.mat.mat) <- paste0("f:", colnames(f.mat.mat))
  colnames(f.dam.mat) <- paste0("f:", colnames(f.dam.mat))
  colnames(f.sire.mat) <- paste0("f:", colnames(f.sire.mat))
  colnames(f.inbred.mat) <- paste0("f:", colnames(f.inbred.mat))
  colnames(f.jk.mat) <- paste0("f:", colnames(f.jk.mat))
  colnames(f.asymm.mat) <- paste0("f:", colnames(f.asymm.mat))

  if (is.null(batch.col.name) & is.null(batch.1.col.name)) {
    return(list(dam.mat = dam.mat, sire.mat = sire.mat, add.mat = add.mat,
                mat.mat = mat.mat, inbred.mat = inbred.mat, jk.mat = jk.mat,
                asymm.mat = asymm.mat,
                f.dam.mat = f.dam.mat, f.sire.mat = f.sire.mat, f.add.mat = f.add.mat,
                f.mat.mat = f.mat.mat, f.inbred.mat = f.inbred.mat, f.jk.mat = f.jk.mat,
                f.asymm.mat = f.asymm.mat))
  }
  else {
    if(1==length(batch.col.name)){
      batch.mat <- incidence.matrix(as.factor(data[, as.character(batch.col.name)]))
      batch.1.mat <- NULL
      try(batch.1.mat <- incidence.matrix(as.factor(data[,
                                                         as.character(batch.1.col.name)])))
      return(list(dam.mat = dam.mat, sire.mat = sire.mat, add.mat = add.mat,
                  mat.mat = mat.mat, inbred.mat = inbred.mat, jk.mat = jk.mat,
                  asymm.mat = asymm.mat,
                  f.dam.mat = f.dam.mat, f.sire.mat = f.sire.mat, f.add.mat = f.add.mat,
                  f.mat.mat = f.mat.mat, f.inbred.mat = f.inbred.mat, f.jk.mat = f.jk.mat,
                  f.asymm.mat = f.asymm.mat, batch.mat = batch.mat,
                  batch.1.mat = batch.1.mat))
    }
    else{
      stop("Not implemented for length(batch.col.name)>1")
    }
  }
}


#' @title diallelMatrixMakeAndRotateSexed
#' @description Make design matrices for diallel, rotate to n-1 space.
#' @param data data frame
#' @param dam.col.name dam column name
#' @param sire.col.name sire column name
#' @param batch.col.name name of batch/random effect column
#' @param batch.1.col.name name of additional batch/random effect column
#' @param strains character vector indicating expected strain names
#' @param sex.col.name the column name for sex
#' @param female.sex.label the character string indicating female sex
#' @param ... additional arguments
#' @return returns diallel incidence matrices, rotated
#' @examples
#' ## not run
#' @export
diallelMatrixMakeAndRotateSexed <- function(data, dam.col.name, sire.col.name,
                                       batch.col.name=NULL, batch.1.col.name=NULL,
                                       strains=c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB"),
                                       sex.col.name, female.sex.label, ...){
  n.strains <- length(strains)
  matrices <- diallelMatrixMakerSexed(data=data, dam.col.name=dam.col.name, sire.col.name=sire.col.name, batch.col.name=batch.col.name,
                                 batch.1.col.name=batch.1.col.name, strains=strains, sex.col.name=sex.col.name,
                                 female.sex.label=female.sex.label)
  dam.mat <- matrices$dam.mat
  sire.mat <- matrices$sire.mat
  add.mat <- matrices$add.mat
  mat.mat <- matrices$mat.mat
  inbred.mat <- matrices$inbred.mat
  jk.mat <- matrices$jk.mat
  asymm.mat <- matrices$asymm.mat
  f.dam.mat <- matrices$f.dam.mat
  f.sire.mat <- matrices$f.sire.mat
  f.add.mat <- matrices$f.add.mat
  f.mat.mat <- matrices$f.mat.mat
  f.inbred.mat <- matrices$f.inbred.mat
  f.jk.mat <- matrices$f.jk.mat
  f.asymm.mat <- matrices$f.asymm.mat

  n.jk <- n.strains*(n.strains-1)/2
  M.dam <- makeRotationMatrix(X = dam.mat, n = n.strains)
  M.sire <- makeRotationMatrix(X = sire.mat, n = n.strains)
  M.add <- makeRotationMatrix(X = add.mat, n = n.strains)
  M.mat <- makeRotationMatrix(X = mat.mat, n = n.strains)
  M.inbred <- makeRotationMatrix(X = inbred.mat, n = n.strains)
  M.jk <- makeRotationMatrix(X = jk.mat, n = n.jk)
  M.asymm <- makeRotationMatrix(X = asymm.mat, n = n.jk)

  M.f.dam <- makeRotationMatrix(X = f.dam.mat, n = n.strains)
  M.f.sire <- makeRotationMatrix(X = f.sire.mat, n = n.strains)
  M.f.add <- makeRotationMatrix(X = f.add.mat, n = n.strains)
  M.f.mat <- makeRotationMatrix(X = f.mat.mat, n = n.strains)
  M.f.inbred <- makeRotationMatrix(X = f.inbred.mat, n = n.strains)
  M.f.jk <- makeRotationMatrix(X = f.jk.mat, n = n.jk)
  M.f.asymm <- makeRotationMatrix(X = f.asymm.mat, n = n.jk)

  t.dam.mat <- dam.mat %*% M.dam
  t.sire.mat <- sire.mat %*% M.sire
  t.add.mat <- add.mat %*% M.add
  t.mat.mat <- mat.mat %*% M.mat
  t.inbred.mat <- inbred.mat %*% M.inbred
  t.jk.mat <- jk.mat %*% M.jk
  t.asymm.mat <- asymm.mat %*% M.asymm
  t.f.dam.mat <- f.dam.mat %*% M.f.dam
  t.f.sire.mat <- f.sire.mat %*% M.f.sire
  t.f.add.mat <- f.add.mat %*% M.f.add
  t.f.mat.mat <- f.mat.mat %*% M.f.mat
  t.f.inbred.mat <- f.inbred.mat %*% M.f.inbred
  t.f.jk.mat <- f.jk.mat %*% M.f.jk
  t.f.asymm.mat <- f.asymm.mat %*% M.f.asymm
  if(!is.null(batch.col.name)){
    batch.mat <- matrices$batch.mat
    M.batch <- makeRotationMatrix(X = batch.mat, n = ncol(batch.mat))
    t.batch.mat <- batch.mat %*% M.batch
    if(!is.null(batch.1.col.name)){
      batch.1.mat <- matrices$batch.1.mat
      M.batch.1 <- makeRotationMatrix(X = batch.1.mat, n = ncol(batch.1.mat))
      t.batch.1.mat <- batch.1.mat %*% M.batch.1
      return(list(RotMat=list( M.dam=M.dam, M.sire=M.sire,
                               M.add=M.add, M.mat=M.mat, M.inbred=M.inbred,
                               M.jk=M.jk, M.asymm=M.asymm,
                               M.batch=M.batch, M.batch.1=M.batch.1,
                               M.f.dam=M.f.dam, M.f.sire=M.f.sire,
                               M.f.add=M.f.add, M.f.mat=M.f.mat, M.f.inbred=M.f.inbred,
                               M.f.jk=M.f.jk, M.f.asymm=M.f.asymm),
                  DesignMat=list(t.dam.mat=t.dam.mat, t.sire.mat=t.sire.mat,
                                 t.add.mat = t.add.mat, t.mat.mat = t.mat.mat, t.inbred.mat = t.inbred.mat,
                                 t.jk.mat = t.jk.mat, t.asymm.mat = t.asymm.mat,
                                 t.f.dam.mat=t.f.dam.mat, t.f.sire.mat=t.f.sire.mat,
                                 t.f.add.mat = t.f.add.mat, t.f.mat.mat = t.f.mat.mat, t.f.inbred.mat = t.f.inbred.mat,
                                 t.f.jk.mat = t.f.jk.mat, t.f.asymm.mat = t.f.asymm.mat,
                                 t.batch.mat = t.batch.mat, t.batch.1.mat = t.batch.1.mat)))
    }else{
      return(list(RotMat=list( M.dam=M.dam, M.sire=M.sire,
                               M.add=M.add, M.mat=M.mat, M.inbred=M.inbred,
                               M.jk=M.jk, M.asymm=M.asymm,
                               M.batch=M.batch,
                               M.f.dam=M.f.dam, M.f.sire=M.f.sire,
                               M.f.add=M.f.add, M.f.mat=M.f.mat, M.f.inbred=M.f.inbred,
                               M.f.jk=M.f.jk, M.f.asymm=M.f.asymm),
                  DesignMat=list(t.dam.mat=t.dam.mat, t.sire.mat=t.sire.mat,
                                 t.add.mat = t.add.mat, t.mat.mat = t.mat.mat, t.inbred.mat = t.inbred.mat,
                                 t.jk.mat = t.jk.mat, t.asymm.mat = t.asymm.mat,
                                 t.f.dam.mat=t.f.dam.mat, t.f.sire.mat=t.f.sire.mat,
                                 t.f.add.mat = t.f.add.mat, t.f.mat.mat = t.f.mat.mat, t.f.inbred.mat = t.f.inbred.mat,
                                 t.f.jk.mat = t.f.jk.mat, t.f.asymm.mat = t.f.asymm.mat,
                                 t.batch.mat = t.batch.mat)))
    }
  }else{
    return(list(RotMat=list( M.dam=M.dam, M.sire=M.sire,
                             M.add=M.add, M.mat=M.mat, M.inbred=M.inbred,
                             M.jk=M.jk, M.asymm=M.asymm,
                             M.f.dam=M.f.dam, M.f.sire=M.f.sire,
                             M.f.add=M.f.add, M.f.mat=M.f.mat, M.f.inbred=M.f.inbred,
                             M.f.jk=M.f.jk, M.f.asymm=M.f.asymm),
                DesignMat=list(t.dam.mat=t.dam.mat, t.sire.mat=t.sire.mat,
                               t.add.mat = t.add.mat, t.mat.mat = t.mat.mat, t.inbred.mat = t.inbred.mat,
                               t.jk.mat = t.jk.mat, t.asymm.mat = t.asymm.mat,
                               t.f.dam.mat=t.f.dam.mat, t.f.sire.mat=t.f.sire.mat,
                               t.f.add.mat = t.f.add.mat, t.f.mat.mat = t.f.mat.mat, t.f.inbred.mat = t.f.inbred.mat,
                               t.f.jk.mat = t.f.jk.mat, t.f.asymm.mat = t.f.asymm.mat)))
  }
}
