
################################################################
# Basic funtions
################################################################

##' @param df Input the data frame to be normalized.
##'
##' @export

quantile_normalization <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)

  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }

  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

#' Title
#'
#' @param m the matrix need to zscore transform
#' @param min the minimum value for the zscore
#' @param max the maximum value for the zscore
#' @param limit Whether to limit the maximum and minimum values of zscore
#'
#' @return matrix after the zscore transform
#' @export
#'
#' @examples

rowZscores <- function(m = NULL, min = -2, max = 2, limit = FALSE){
  z <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m),`/`)
  if(limit){
    z[z > max] <- max
    z[z < min] <- min
  }
  return(z)
}


#' Calculate the overlap interval
#'
#' @param inter_result Result tibble from valr::bed_intersect, contains (start.x, end.x) and (start.y, end.y)
#'
#' @return A tibble contains the interval(start, end) of overlap region
#' @export
#' @importFrom dplyr mutate
#'
#' @examples
#' get_overlap(valr::bed_intersect(a, b))
get_overlap <- function(inter_result) {
  result <- dplyr::mutate(inter_result,
                          start = pmax(start.x, start.y),
                          end = pmin(end.x, end.y),
  )[c('chrom','start','end')]

  result <- unique(result)

  return(result)
}


#' Get distinctive records of A
#'
#' @param inter_result Result tibble from valr::bed_intersect
#'
#' @return A tibble contains distinctive intervals of x
#' @export
#'
#' @examples
#' get_unique(valr::bed_intersect(a, b))
get_unique <- function(inter_result) {
  return(unique(inter_result[c('chrom','start.x','end.x')]))
}


#' Change identity of the beds$rec based on intersection result
#'
#' @param bed_list A list contains records of all bed files(Read by function get_beds())
#' @param id The ID of bed in bed_list
#' @param rec Result tibble from valr::intersection, contains (chrom, start, end, ident)
#'
#' @return A list which has been updated bed$res$ident based on intersection result
#' @importFrom dplyr left_join
change_res_ident <- function(bed_list, id, rec){
  source <- bed_list[[id]]$res
  newdata <- dplyr::left_join(source, rec, by=c('chrom','start', 'end'))
  newdata$ident.y[is.na(newdata$ident.y)] <- newdata$ident.x[is.na(newdata$ident.y)]
  newdata$ident.x <- NULL
  colnames(newdata)[5] <- 'ident'
  bed_list[[id]]$res <- newdata
  return(bed_list)
}


#' Update ident of bed_list based on intersection result
#'
#' @param combine Combination of beds ident
#' @param bed_list A list contains records of all bed files(Read by function get_beds())
#' @param entrance Mark the entrance of recursion
#' @param fparam Minimum overlap required as a fraction of A. Default is 1E-9 (i.e. 1bp)
#' @param Fparam Minimum overlap required as a fraction of B. Default is 1E-9 (i.e. 1bp)
#' @param rparam equire that the fraction of overlap be reciprocal for A and B
#'
#' @return A list which has been updated ident
#' @importFrom valr bed_intersect
get_intersect <- function(combine, bed_list, entrance=T, fparam=1e-9, Fparam=1e-9, rparam=F){
  if(length(combine)==2){
    a <- bed_list[[combine[1]]]$interval
    b <- bed_list[[combine[2]]]$interval

    inter_result <- valr::bed_intersect(a, b)
    inter_result <- dplyr::filter(inter_result, .overlap/len.x >= fparam)
    if(rparam){
      inter_result <- dplyr::filter(inter_result, .overlap/len.y >= fparam)
    }else{
      inter_result <- dplyr::filter(inter_result, .overlap/len.y >= Fparam)
    }


    if(entrance){
      res_a  <- unique(
        dplyr::mutate(inter_result,
                      ident.x=paste(ident.x, ident.y, sep="")
        )[c('chrom','start.x', 'end.x','ident.x')])
      colnames(res_a) <- c('chrom','start', 'end','ident')
      bed_list <- change_res_ident(bed_list, combine[1], res_a)
      res_b  <- unique(
        dplyr::mutate(inter_result,
                      ident.y=paste(ident.x, ident.y, sep="")
        )[c('chrom','start.y', 'end.y','ident.y')])
      colnames(res_b) <- c('chrom','start', 'end','ident')
      bed_list <- change_res_ident(bed_list, combine[2], res_b)

      return(bed_list)
    }else{
      return(inter_result)
    }
  }else{
    f_inter_result <- get_overlap(get_intersect(combine[-1], bed_list, entrance=F))

    a <- bed_list[[combine[1]]]$interval
    b <- f_inter_result

    inter_result <- valr::bed_intersect(a, b)
    if(entrance){
      inter_result <- get_overlap(inter_result)
      for(ident in combine){
        a <- bed_list[[ident]]$interval
        b <- inter_result

        tmp_inter <- valr::bed_intersect(a,b)
        tmp_inter <- dplyr::filter(tmp_inter, .overlap/len.x >= fparam)

        tmp_res <- unique(
          dplyr::mutate(tmp_inter,
                        ident.x=paste(combine, collapse="")
          )[c('chrom','start.x', 'end.x','ident.x')])
        colnames(tmp_res) <- c('chrom','start', 'end','ident')
        bed_list <- change_res_ident(bed_list, ident, tmp_res)
      }

      return(bed_list)
    }else{
      return(inter_result)
    }
  }
}


#' Init dataframe which would be used in intersection results based on bed ident
#'
#' @param n_combine The number of samples
#' @param num_i Combination number
#'
#' @return Initiated dataframe with columns key and value
#'
#' @examples
#' init_result_df(4, 2)
init_result_df <- function(n_combine, num_i){
  ident_df <- data.frame(t(combn(c(1:n_combine), num_i)))
  key <- apply(ident_df,1,paste, collapse="")
  result_df <- data.frame(key=key, value=rep(0, length(key)))
  return(result_df)
}


#' Get intersection result which can be used as input of venn/upset plot
#'
#' @param bed_list A list contains records of all bed files(Read by function get_beds())
#'
#' @return A dataframe with columns key and value, where key represents combination of each bed and value represents the number of intervals
#' @export
#' @importFrom dplyr left_join
#'
#' @examples
#' # Get beds using get_beds
#' beds <- get_beds(c('path1', 'path2', ... ))
#'
#' # Use get_combine_result() to get intersection result
#' # The processing time may vary depending on the number of sets involved
#' inputs <- get_combine_result(beds)
get_combine_result <- function(bed_list){
  n_combine <- length(bed_list)

  if(n_combine>=10){
    message("ERROR: Number of datasets can not be greater than 9")
    return()
  }
  if(n_combine<2){
    message("ERROR: At least 2 datasets are needed")
    return()
  }

  result_df <- data.frame(key=c(1:n_combine), value=rep(0, length(n_combine)))

  for(i in c(2:n_combine)){
    # Get result structure
    result_df <- rbind(result_df, init_result_df(n_combine, i))
    # Get Intersection
    comb_df <- data.frame(combn(c(1:n_combine), i))
    for(row in comb_df){
      bed_list <- get_intersect(row, bed_list)
    }
  }

  for(beds_rec in bed_list){
    res <- data.frame(table(beds_rec$res$ident))
    colnames(res) <- c('key', 'value')

    tmp_result <- dplyr::left_join(result_df, res, by='key')
    tmp_result$value.x <- ifelse(
      tmp_result$value.x == 0 & !is.na(tmp_result$value.y),
      tmp_result$value.y,
      tmp_result$value.x
    )
    tmp_result$value.y <- NULL
    colnames(tmp_result)[2] <- 'value'
    result_df <- tmp_result
  }

  return(result_df)

}


#' Get keys which can be used in venn plot
#'
#' @param source_key The combination of number: 1,2,3,12,13,23,123
#' @param key_len The num of different numbers
#'
#' @return Binary keys such as 001,010,100...
#'
#' @examples
#' get_venn_key(c(1,2,3,12,13,23,123), 3)
get_venn_key <- function(source_key, key_len){
  tmp_key <- rep('0', key_len)
  for(pos in strsplit(source_key, "")[[1]]){
    tmp_key[as.integer(pos)] <- 1
  }
  return(paste0(tmp_key, collapse = ""))
}


#' Get expression input which can be used in upset plot
#'
#' @param source_df A dataframe contains key and value(From get_combine_result())
#' @param names A vector includes the names of each bed
#'
#' @return A vector can be used in UpSetR fromExpression
#'
#' @examples
#' # Get beds using get_beds
#' beds <- get_beds(c('path1', 'path2', ... ))
#'
#' # Use get_combine_result() to get intersection result
#' # The processing time may vary depending on the number of sets involved
#' inputs <- get_combine_result(beds)
#'
#' get_upset_input(inputs, names)
get_upset_input <- function(source_df, names) {
  string <- apply(source_df,1, function(rec, names){
    key <- rec['key']
    value <- rec['value']
    result <- names[sapply(strsplit(key,"")[[1]], as.integer)]
    result <- paste(result, collapse = "&")
    result <- paste0("`", result, "`=", value)
  }, names=names)

  result <- eval(parse(text=paste0('c(',paste(string, collapse = ","),')')))
  return(result)
}


#' Draw flower plot
#'
#' @param sample A vector contains the name of each petal
#' @param value A vector contains the value of each petal
#' @param labels The label of the core
#' @param start Start deg, default is 90
#' @param a Width of the petal, default is 0.5
#' @param b Length of the petal, default is 2
#' @param r Radius of the core, default is 1.5
#' @param width Width of the plot, default is 10
#' @param height Heigth of the plot, default is 10
#' @param ellipse_col Color of the petal
#' @param circle_col Color of the core
#' @param circle_text_cex Text size of label around the petal, default is 1
#'
#' @return Flower plot
#' @export
#'
#' @examples
#' flower_plot(c("WSM419", "A321", "M1", "M2", "M22", "M58", "M102", "M161",
#' "KH36b", "KH36c", "KH36d", "KH53a", "KH53b"),c(519, 556, 83, 62, 415, 425,
#' 357, 441, 22, 41, 33, 44, 43), 90, 0.5, 2, labels="core")
flower_plot <- function(sample, value, labels, start=90, a=0.5, b=2, r=1.5, width=10, height=10,
                        ellipse_col = rgb(135, 206, 235, 150, max = 255),
                        circle_col = rgb(0, 162, 214, max = 255),
                        circle_text_cex = 1){
  par( bty = "n", ann = F, xaxt = "n", yaxt = "n", mar = c(1,1,1,1))
  plot(c(0,width),c(0,height),type="n")
  n <- length(sample)
  deg <- 360 / n
  res <- lapply(1:n, function(t){
    plotrix::draw.ellipse(x = width/2 + cos((start + deg * (t - 1)) * pi / 180),
                          y = height/2 + sin((start + deg * (t - 1)) * pi / 180),
                          col = ellipse_col,
                          border = ellipse_col,
                          a = a, b = b, angle = deg * (t - 1))
    text(x = width/2 + width/4 * cos((start + deg * (t - 1)) * pi / 180),
         y = height/2 + height/4 * sin((start + deg * (t - 1)) * pi / 180),
         value[t])

    if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
      text(x = width/2 + width/3.03 * cos((start + deg * (t - 1)) * pi / 180),
           y = height/2 + height/3.03 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) - start,
           adj = 1,
           cex = circle_text_cex
      )

    } else {
      text(x = width/2 + width/3.03 * cos((start + deg * (t - 1)) * pi / 180),
           y = height/2 + height/3.03 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) + start,
           adj = 0,
           cex = circle_text_cex
      )
    }
  })
  plotrix::draw.circle(x = width/2, y = height/2, r = r, col = circle_col, border = circle_col)
  text(x = width/2, y = height/2, labels=labels)
}


## https://www.r-bloggers.com/bigcor-large-correlation-matrices-in-r/
#' Title
#'
#' @param x
#' @param nblocks
#' @param verbose
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
bigcor <- function(x, nblocks = 10, verbose = TRUE, ...){
  library(ff, quietly = TRUE)
  NCOL <- ncol(x)
  ## test if ncol(x) %% nblocks gives remainder 0
  if (NCOL %% nblocks != 0) stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")
  ## preallocate square matrix of dimension
  ## ncol(x) in 'ff' single format
  corMAT <- ff(vmode = "single", dim = c(NCOL, NCOL))
  ## split column numbers into 'nblocks' groups
  SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))
  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)
  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  for (i in 1:nrow(COMBS)) {
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
    flush.console()
    COR <- cor(x[, G1], x[, G2], ...)
    corMAT[G1, G2] <- COR
    corMAT[G2, G1] <- t(COR)
    COR <- NULL
  }
  gc()
  return(corMAT)
}

## R code - plot LD decay and heatmap
## Simple R code to plot LD decay and heatmap
################################################# https://fabiomarroni.wordpress.com/2011/08/09/estimate-decay-of-linkage-disequilibrium-with-distance/
Er <- function(C_, d, n){
  length(d)
  res <- ((10+C_*d)/((2+C_*d)*(11+C_*d)))*(1+((3+C_*d)*(12+12*C_*d+(C_*d)^2)/(n*(2+C_*d)*(11+C_*d))))
  return(res)
}

LDit <- function(x, n, ...){
  d <- x$Dist;
  y <- x[order(d), ]
  ld <- y$rsq

  ## quartz()
  # plot(d[order(d)], ld, cex=.5, pch=19, col="grey", xlab="Distance", ylab=expression(r^2), ...)
  plot(0, type='n', xlim=c(0, max(d)), ylim=c(0, 1), xlab="Distance (Kb)", ylab=expression(r^2), ...)
  nlm <- nls(ld~Er(C_, d[order(d)], n), start=list(C_=0.01))
  C_ <- summary(nlm)$coefficients[1]
  lines(d[order(d)], Er(C_, d[order(d)], n), col="black", lwd=2)
  abline(h=0.1, col='red', lty=5)
  ## quartz()
  require(lattice)
  row_bp <- unique(x$Si)
  col_bp <- unique(x$Sj)
  segsites <- length(row_bp)+1
  vlen <- length(x$rsq)
  #cm <- matrix(nrow=segsites, ncol=segsites, dimnames=list(c(row_bp, col_bp[segsites-1]), col_bp))#
  cm <- matrix(nrow=segsites, ncol=segsites, dimnames=list(c(row_bp, col_bp[segsites-1]), c(row_bp[1], col_bp)))
  count <- 0;
  nrows <- segsites-1;
  for(r in 1:nrows){
    start <- r+1
    for(c in start:segsites){
      count <- count+1
      cm[r, c] <- x$rsq[count]
    }
  }
  levelplot(cm, at=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), col.regions=rev(grey.colors(10)),
            scales=list(cex=.4, x=list(cex=.4), tck=c(1, 0), alternating=c(1)), ylab="", xlab="")
}

LDplot <- function(x, n, ...){
  d <- x$Dist;
  od <- order(d)
  y <- x[od, ]
  ld <- y$rsq

  plot(0, type='n', xlim=c(0, max(d)), ylim=c(0, 1), ...)
  nlm <- nls(ld~Er(C_, d[od], n), start=list(C_=0.01))
  C_ <- summary(nlm)$coefficients[1]
  fit <- Er(C_, d[od], n)
  lines(d[od], fit, col="black", lwd=2)
  abline(h=0.2, col='red', lty=5)
  r0.1 <- d[od][order(abs(fit-0.1))[1]]
  r0.2 <- d[od][order(abs(fit-0.2))[1]]
  r0.3 <- d[od][order(abs(fit-0.3))[1]]
  r0.4 <- d[od][order(abs(fit-0.4))[1]]
  r0.5 <- d[od][order(abs(fit-0.5))[1]]
  res <- c(r0.1, r0.2, r0.3, r0.4, r0.5)
  text(res, seq(0.1, 0.5, by=0.1), res, col='red', cex=0.8)
  segments(res, rep(0, 5), res, seq(0.1, 0.5, by=0.1), col='red', lty=5)
  return(res)
}


#' Title  Message a timestamp.
#'
#' @param msg The message.
#'
#' @return
#' @export
#'
#' @examples
logfile <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(paste(timestamp, msg))
}

checkFTAnno <- function(){
  if (!exists("FTAnno")) {
    stop("Footprint annotation not exist. Please **addFootprint** first.")
  }
}

checkGeAnno  <- function(){
  if (!exists("CATAnno")) {
    stop("Necessary Annotation not exist. Please **addAnnotation** first.")
  }
}

