#Back and Forward Filling without Selecting
fill_na <- function(vector, n){
  if (n == 0) {
    vector
  } else {
    fill_na(
      vector = dplyr::coalesce(vector, dplyr::lag(vector)),
      n = n - 1
    )
  }
}
annotatedchasedf$ChaseID <- fill_na(annotatedchasedf$ChaseID, n=24)
annotatedchasedf_tmp <- annotatedchasedf
annotatedchasedf_tmp <- annotatedchasedf_tmp %>%
  group_by(track) %>%
  lead(annotatedchasedf_tmp$ChaseID, n=25)