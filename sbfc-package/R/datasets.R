#' @docType data
#' @keywords datasets
#' @name corral_augmented
#' @usage data(corral_augmented)
#' @title Augmented corral data set: synthetic data with correlated attributes augmented with noise features
#' @description This is an artificial domain where the target concept is (X1^X2) V (X3^X4). \cr
#' Data set from John et al (1994). Training and test splits from SGI. \cr
#' The first 6 features are the real features from the original corral data set.
#' The rest are noise features added by V. Krakovna by shuffling copies of real features.\cr
#' The SBFC paper uses subsets of this data set with the first 100 and 1000 features. 
#' @format \describe{
#' \item{\code{TrainX}}{A matrix with 128 rows and 10000 columns.}
#' \item{\code{TrainY}}{A vector with 128 rows.}
#' }
#' @examples corral_result = sbfc(data=list(TrainX=corral_augmented$TrainX[,1:6],
#'                                          TrainY = corral_augmented$TrainY))
#' corral100_result = sbfc(data=list(TrainX=corral_augmented$TrainX[,1:100], 
#'                                   TrainY = corral_augmented$TrainY))
#' @references \href{https://ai.stanford.edu/~ronnyk/ml94.pdf}{John et al (1994) paper introducing the corral data set}
#' @references \href{https://arxiv.org/abs/1506.02371}{SBFC paper describing augmentation of corral data set}
NULL

#' @docType data
#' @keywords datasets
#' @name heart
#' @usage data(heart)
#' @title Heart disease data set: disease outcomes given health attributes
#' @description Data set from UCI repository, discretized using the \code{mdlp} package.
#' @format \describe{
#' \item{\code{TrainX}}{A matrix with 270 rows and 13 columns.}
#' \item{\code{TrainY}}{A vector with 270 rows.}
#' }
#' @references \href{https://archive.ics.uci.edu/ml/datasets/Statlog+(Heart)}{UCI heart data set}
NULL

#' @docType data
#' @keywords datasets
#' @name madelon
#' @usage data(madelon)
#' @title Madelon data set: synthetic data from NIPS 2003 feature selection challenge
#' @description This is a two-class classification problem. 
#' The difficulty is that the problem is multivariate and highly non-linear. 
#' Of the 500 features, 20 are real features, 480 are noise features. \cr
#' Data set from UCI repository, discretized using median cutoffs.
#' @format \describe{
#' \item{\code{TrainX}}{A matrix with 2000 rows and 500 columns.}
#' \item{\code{TrainY}}{A vector with 2000 rows.}
#' \item{\code{TestX}}{A matrix with 600 rows and 500 columns.}
#' \item{\code{TestY}}{A vector with 600 rows.}
#' }
#' @references \href{https://archive.ics.uci.edu/ml/datasets/Madelon}{UCI madelon data set}
NULL
