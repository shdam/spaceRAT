#' Helper function to compute ranks
#' @noRd
ranking <- function(mat, rank_scale = FALSE) {
    ranked_mat <- apply(mat, 2, function(x) {
        r <- rank(x)
        if (rank_scale) r <- r * (1 + (min(r) / length(r))) #- min(r) #(1 + (sum(r == min(r)) / length(r)))
        return(r)
        }
    )
    return(ranked_mat)
}
