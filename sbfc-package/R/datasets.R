#' @docType data
#' @keywords datasets
#' @name chess
#' @usage data(chess)
#' @title Chess End-Game - King+Rook versus King+Pawn on a7
#' @description Outcomes of chess games given the board descriptions. \cr
#' Data set from UCI repository. Training and test splits from SGI. 
#' @format \describe{
#' \item{\code{TrainX}}{A matrix with 2130 rows and 36 columns.}
#' \item{\code{TrainY}}{A vector with 2130 rows.}
#' \item{\code{TestX}}{A matrix with 1066 rows and 36 columns.}
#' \item{\code{TestY}}{A vector with 1066 rows.}
#' }
#' @references \href{https://archive.ics.uci.edu/ml/datasets/Chess+(King-Rook+vs.+King-Pawn)}{UCI chess data set}
#' @references \href{http://www.sgi.com/tech/mlc/db/chess.names}{SGI listing for chess data set}
NULL

#' @docType data
#' @keywords datasets
#' @name corral
#' @usage data(corral)
#' @title Corral: synthetic data with correlated attributes
#' @description This is an artificial domain where the target concept is (X1^X2) V (X3^X4). \cr
#' Data set by R. Kohavi. Training and test splits from SGI. 
#' @format \describe{
#' \item{\code{TrainX}}{A matrix with 128 rows and 6 columns.}
#' \item{\code{TrainY}}{A vector with 128 rows.}
#' }
#' @references \href{http://www.sgi.com/tech/mlc/db/corral.names}{SGI listing for corral data set}
NULL

#' @docType data
#' @keywords datasets
#' @name corral_augmented
#' @usage data(corral_augmented)
#' @title Augmented corral: synthetic data with correlated attributes augmented with noise features
#' @description This is an artificial domain where the target concept is (X1^X2) V (X3^X4). \cr
#' Data set by R. Kohavi. Training and test splits from SGI. \cr
#' Noise features added by V. Krakovna by shuffling copies of real features.\cr
#' The SBFC paper uses subsets of this data set with the first 100 and 1000 features. 
#' @format \describe{
#' \item{\code{TrainX}}{A matrix with 128 rows and 10000 columns.}
#' \item{\code{TrainY}}{A vector with 128 rows.}
#' }
#' @references \href{http://www.sgi.com/tech/mlc/db/corral.names}{SGI listing for corral data set}
#' @references \href{arxiv.org/abs/1506.02371}{SBFC paper describing augmentation of corral data set}
NULL

#' @docType data
#' @keywords datasets
#' @name heart
#' @usage data(heart)
#' @title Heart disease outcomes given health attributes
#' @description Data set from UCI repository, discretized using the \code{mdlp} package.
#' @format \describe{
#' \item{\code{TrainX}}{A matrix with 270 rows and 13 columns.}
#' \item{\code{TrainY}}{A vector with 270 rows.}
#' }
#' @references \href{https://archive.ics.uci.edu/ml/datasets/Statlog+(Heart)}{UCI heart data set}
#' @references \href{http://www.sgi.com/tech/mlc/db/heart.names}{SGI listing for heart data set}
NULL

#' @docType data
#' @keywords datasets
#' @name madelon
#' @usage data(madelon)
#' @title Madelon: synthetic data set from NIPS 2003 feature selection challenge
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