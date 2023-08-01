#' Title
#'
#' @param quant_matrix   The normalized matrix obtained by **quantification**.
#' @param max_dis   The max distance of two peaks to calculate the correlations.
#' @param window   The window size to shift.
#' @param peak_anno The merged peaks annotation.
#'
#' @return
#' @export
#'
#' @examples   getCoaccessible(peak_anno="peak_anno.txt", quant_matrix="quant_mat.tsv", max_dis=20000, window=20000)
getCoaccessible <- function(peak_anno, quant_matrix, max_dis, window){
  checkGeAnno()
  anno <- read.table(peak_anno, head=T, sep="\t", stringsAsFactors=FALSE)
  summit <- anno[,-c(4)]
  colnames(summit) <- c('chr', 'bp1', 'bp2', 'summit')
  #rownames(summit) <- sprintf("%s_%s_%s", summit$chr, summit$bp1, summit$bp2)
  #category <- read.table(peak_type, head=FALSE, sep="\t", stringsAsFactors=FALSE)
  category <- anno[,-c(5)]
  colnames(category) <- c('chr', 'bp1', 'bp2', 'category')
  #rownames(category) <- sprintf("%s_%s_%s", category$chr, category$bp1, category$bp2)

  # processing the matrix
  mat1 <- mat1 <- read.table(quant_matrix,row.names=1,head=T)
  peak <- data.frame(peak = rownames(mat1))
  peak$chr <- sapply(strsplit(peak$peak,":"), `[`, 1 )
  peak$bed <- sapply(strsplit(peak$peak,":"), `[`, 2 )
  mat2 <- peak %>% tidyr::separate(bed, c("bp1", "bp2"), "-")
  mat2 <- mat2[,c(2:4)]
  mat<- cbind(mat2, mat1)
  rownames(mat) <- NULL
  rownames(mat) <- peaks <- sprintf("%s_%s_%s", mat$chr, mat$bp1, mat$bp2)

  ## raw correlation
  # maxn <- floor(nrow(mat)/100) * 10
  # cormat <- bigcor(t(head(mat[,-c(1:3)], n=maxn)), nblocks=nblocks)
  # rawcor <- matrix(0, nrow=nrow(mat), ncol=nrow(mat))
  # rownames(rawcor) <- colnames(rawcor) <- rownames(mat)
  # rawcor[1:maxn, 1:maxn] <- cormat[1:maxn, 1:maxn]
  # if(nrow(mat) > maxn){
  #   xcor <- cor(t(mat[(maxn+1):nrow(mat),-c(1:3)]))
  #   rawcor[rownames(xcor), colnames(xcor)] <- xcor
  #   xcor <- cor(t(mat[(maxn+1):nrow(mat),-c(1:3)]), t(head(mat[,-c(1:3)], n=maxn)))
  #   rawcor[rownames(xcor), colnames(xcor)] <- xcor
  #   rawcor[colnames(xcor), rownames(xcor)] <- t(xcor)
  # }

  ## regularized correlation
  window <- window
  genomic_coords <- data.frame(CATAnno$genome)
  regcon <- run_cicero(mat, window=window, genomic_coords)
  regcor <- matrix(NA, ncol=length(peaks), nrow=length(peaks))
  rownames(regcor) <- colnames(regcor) <- peaks
  regcor[cbind(as.character(regcon$Peak1),as.character(regcon$Peak2))] <- regcon$coaccess
  ## save.image('connections.RData')
  #dim(regcor)
  #dim(rawcor)
  rawcor[is.na(regcor)] <- NA

  rawcon <- cbind(regcon[,c('Peak1','Peak2')],coaccess=rawcor[cbind(as.character(regcon$Peak1),as.character(regcon$Peak2))])

  colist <- na.omit(regcon)
  colist <- colist[summit[as.character(colist$Peak2),'summit'] > summit[as.character(colist$Peak1),'summit'], ]
  colist$dist <- summit[as.character(colist$Peak2),'summit'] - summit[as.character(colist$Peak1),'summit']
  colist$Center1 <- summit[as.character(colist$Peak1),'summit']
  colist$Center2 <- summit[as.character(colist$Peak2),'summit']
  colist$dist10k <- round(colist$dist / 10000)
  colist$dist2k <- round(colist$dist / 2000)*2000
  qlts <- c(-1, quantile(colist$coaccess, probs=1:2/3), 1)
  colist$colevel <- factor(cut(colist$coaccess, qlts), labels=c('Low','Medium','High'))
  colist$group <- ifelse(category[as.character(colist$Peak1),'category']=="Distal" & category[as.character(colist$Peak2),'category']=="Distal", "EE", ifelse(category[as.character(colist$Peak1),'category']=="Proximal" & category[as.character(colist$Peak2),'category']=="Proximal", "PP", "PE"))
  return(colist)
}
