#' Area Under the Curve of the Reciever Operating Curve
#' 
#' \code{auc} estimates the AUC of the ROC using a Mann-Whitney U statistic.
#' \cr \cr \bold{Note:} this method will exclude any missing data.
#' 
#' 
#' @param obs a vector of observed values which must be 0 for absences and 1
#' for occurrences
#' @param pred a vector of the same length as \code{obs} representing the
#' predicted values. Values must be between 0 & 1 representing a likelihood.
#' @return Returns a single value represting the AUC value.
#' @author Jeremy VanDerWal (extracted originally from package SDMTools) \email{jjvanderwal@@gmail.com}
#' 
auc <- function (obs, pred) 
{
    if (length(obs) != length(pred)) 
        stop("this requires the same number of observed & predicted values")
    if (length(which(is.na(c(obs, pred)))) > 0) {
        na = union(which(is.na(obs)), which(is.na(pred)))
        warning(length(na), " data points removed due to missing data")
        obs = obs[-na]
        pred = pred[-na]
    }
    n = length(obs)
    if (length(which(obs %in% c(0, 1))) != n) 
        stop("observed values must be 0 or 1")
    n1 = as.double(length(which(obs == 1)))
    n0 = as.double(length(which(obs == 0)))
    if (n1 == 0 || n1 == n) 
        return(NaN)
    pred0 = pred[which(obs == 0)]
    pred1 = pred[which(obs == 1)]
    ranks = rank(pred, ties.method = "average")
    ranks0 = ranks[which(obs == 0)]
    ranks1 = ranks[which(obs == 1)]
    U = n0 * n1 + (n0 * (n0 + 1))/2 - sum(ranks0)
    AUC = U/(n0 * n1)
    if (AUC < 0.5) 
        AUC = 1 - AUC
    return(AUC)
}
