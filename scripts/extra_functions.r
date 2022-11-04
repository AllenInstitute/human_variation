
findFromGroups <- function (datExpr, groupVector, fn = "mean"){
    groups = names(table(groupVector))
    fn = match.fun(fn)
    datMeans = matrix(0, nrow = dim(datExpr)[1], ncol = length(groups))
    for (i in 1:length(groups)) {
        datIn = datExpr[, groupVector == groups[i] ]
        if (is.null(dim(datIn)[1])) {
            datMeans[, i] = as.numeric(datIn)
        }
        else {
            datMeans[, i] = as.numeric(apply(datIn, 1, fn))
        }
    }
    colnames(datMeans) = groups
    rownames(datMeans) = rownames(datExpr)
    return(datMeans)
}

getBetaScore <- function (propExpr, returnScore = TRUE, spec.exp = 2) 
{
    calc_beta <- function(y, spec.exp = 2) {
        d1 <- as.matrix(dist(y))
        eps1 <- 1e-10
        score1 <- sum(d1^spec.exp)/(sum(d1) + eps1)
        score1
    }
    betaScore <- apply(propExpr, 1, calc_beta)
    betaScore[is.na(betaScore)] <- 0
    if (returnScore) 
        return(betaScore)
    scoreRank <- rank(-betaScore)
    scoreRank
}

subsampleCells <- function (clusters, subSamp = 25, seed = 5) 
{
    if (length(subSamp) == 1) 
        subSamp = rep(subSamp, length(unique(as.character(clusters))))
    if (is.null(names(subSamp))) 
        names(subSamp) <- unique(as.character(clusters))
    kpSamp <- rep(FALSE, length(clusters))
    for (cli in unique(as.character(clusters))) {
        val = subSamp[cli]
        if (!is.na(val)[1]) {
            set.seed(seed)
            seed <- seed + 1
            kp <- which(clusters == cli)
            kpSamp[kp[sample(1:length(kp), min(length(kp), val))]] <- TRUE
        }
    }
    kpSamp
}

logCPM <- function (counts) 
{
    norm.dat <- cpm(counts)
    if (is.matrix(norm.dat)) {
        norm.dat <- log2(norm.dat + 1)
    }
    else {
        norm.dat@x <- log2(norm.dat@x + 1)
    }
    norm.dat
}

cpm <- function (counts) 
{
    sf <- Matrix::colSums(counts)/1e+06
    if (is.matrix(counts)) {
        return(t(t(counts)/sf))
    }
    else if (class(counts) == "dgCMatrix") {
        sep <- counts@p
        sep <- sep[-1] - sep[-length(sep)]
        j <- S4Vectors::Rle(1:length(sep), sep)
        counts@x <- counts@x/sf[as.integer(j)]
    }
    else if (class(counts) == "dgTMatrix") {
        j = counts@j
        counts@x = counts@x/sf[j + 1]
    }
    else {
        stop(paste("cpm function for", class(counts)[1], "not supported"))
    }
    return(counts)
}

tmean      <- function(x) mean(x,trim=0.25)

groupSamples <- function (groups, nGroups = 4, seed = 5){
  kpSamp <- rep(0, length(groups))
  for (cli in unique(as.character(groups))) {
    set.seed(seed)
    seed <- seed + 1
    n <- sum(groups==cli)
    val = (sample(1:n,n)+seed)%%nGroups + 1
    kpSamp[groups==cli] = val
  }
  kpSamp
}

rf_prediction_from_seurat <- function(INPUT, pcs=30, seed=1, scvi_normalized=FALSE){
  ## Calculate PCs
  if(scvi_normalized){
    INPUT@assays$RNA@data <- logCPM(INPUT@assays$RNA@scale.data)
  }
  INPUT  <- FindVariableFeatures(INPUT, selection.method = "vst",  nfeatures = 2000, verbose = FALSE)
  INPUT  <- ScaleData(INPUT, verbose = FALSE)
  INPUT  <- RunPCA(INPUT, features = VariableFeatures(object = INPUT), verbose = FALSE, npcs=pcs)
  datPCA <- FetchData(INPUT,vars = c("donor_name",paste0("PC_",1:pcs)))
  datPCA$donor_name <- droplevels(as.factor(datPCA$donor_name))

  ## Divide data into four groups and run random forest prediction 
  sample <- groupSamples(datPCA$donor_name, nGroups=4, seed=seed)
  predictions <- NULL
  for (j in 1:4){
    train  <- subset(datPCA, sample != j)
    test   <- subset(datPCA, sample == j)
    rf     <- randomForest(donor_name ~ ., data=train)
    pred   <- predict(rf, newdata=test[-1])  # the first column is donor name
    predictions <- rbind(predictions, cbind(as.character(test[,1]),as.character(pred)))
  }
  colnames(predictions) <- c("donor","predicted_donor")
  predictions 
} 

rf_prediction_from_logCPM <- function(INPUT, genes=NULL, nfeatures = 2000, seed=5){
  if(is.null(genes))
    genes <- VariableFeatures(FindVariableFeatures(INPUT, selection.method = "vst",  nfeatures = nfeatures, verbose = FALSE))
  datExpr <- data.frame(donor_name=FetchData(INPUT,vars = "donor_name"))
  datExpr <- cbind(datExpr,t(INPUT@assays$RNA@data[genes,colnames(INPUT)]))
  datExpr$donor_name <- droplevels(as.factor(datExpr$donor_name))
  genes <- setNames(colnames(datExpr),make.names(colnames(datExpr)))
  colnames(datExpr) <- names(genes)
  
  ## Divide data into four groups and run random forest prediction 
  sample <- groupSamples(datExpr$donor_name, nGroups=4, seed=seed)
  predictions <- NULL
  importance  <- list()
  for (j in 1:4){
    train  <- subset(datExpr, sample != j)
    test   <- subset(datExpr, sample == j)
    rf     <- randomForest(donor_name ~ ., data=train, importance=TRUE)
    pred   <- predict(rf, newdata=test[-1])  # the first column is donor name
    predictions <- rbind(predictions, cbind(as.character(test[,1]),as.character(pred)))
    importance[[j]] <- rf$importance
    rownames(importance[[j]]) <- as.character(genes[rownames(importance[[j]])])
  }
  colnames(predictions) <- c("donor","predicted_donor")
  list(predictions=predictions,importance=importance)
} 


run_dmvhyper <- function(or1,or2,n){
  or2 <- factor(or2,levels = names(table(or1)))
  or1 <- setNames(as.numeric(table(or1)),names(table(or1)))
  or2 <- setNames(as.numeric(table(or2)),names(table(or2)))
  signif(dmvhyper(x = or2, n = or1, k = sum(or2)),3)
}


getStats <- function(x,l=round(length(x)/2)){
  perm <- x[1:l]
  real <- x[(l+1):length(x)]
  tt   <- t.test(perm,real,paired=FALSE,na.action=na.omit)
	c(mean(perm,na.rm=TRUE),mean(real,na.rm=TRUE),
	  sd(perm,na.rm=TRUE),sd(real,na.rm=TRUE),tt$statistic,tt$p.val)
}

getStats_initial <- function(x){
  l    <- round(length(x)/2)
  perm <- x[1:l]
  real <- x[(l+1):length(x)]
  tt   <- t.test(perm,real,paired=TRUE,na.action=na.omit)
  c(mean(perm,na.rm=TRUE),mean(real,na.rm=TRUE),
    sd(perm,na.rm=TRUE),sd(real,na.rm=TRUE),tt$statistic,tt$p.val)
}

