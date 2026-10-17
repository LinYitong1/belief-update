block_mv_slope <- function(dt, item_cols, response_col = "response",
                           response_scale = 1, expected_repeats = 2L) {
  x <- copy(dt)[, response_mv := get(response_col) / response_scale]
  stopifnot(all(x$response_mv >= 0 & x$response_mv <= 1))

  item_mv <- x[, .(n_repeat = .N,
                   v = var(response_mv),
                   mv_term = mean(response_mv) * (1 - mean(response_mv))),
               by = item_cols]
  stopifnot(all(item_mv$n_repeat == expected_repeats))

  coef(lm(v ~ mv_term, data = item_mv))[["mv_term"]]
}
