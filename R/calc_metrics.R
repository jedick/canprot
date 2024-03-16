# canprot/calc_metrics.R
# Calculate one or more chemical metrics for proteins
# 20191027 initial version as canprot/metrics.R
# 20230704 adapted for chem16S/calc_metrics.R
# 20240302 moved to canprot
calc_metrics <- function(AAcomp, metrics = c("Zc", "nO2", "nH2O"), ...) {

  # Replace shortcuts with function names for metrics 20240305
  metrics.orig <- metrics
  metrics <- tolower(metrics)
  metrics[metrics == "length"] <- "plength"
  metrics[metrics %in% c("h/c", "h_c")] <- "hc"
  metrics[metrics %in% c("n/c", "n_c")] <- "nc"
  metrics[metrics %in% c("o/c", "o_c")] <- "oc"
  metrics[metrics %in% c("s/c", "s_c")] <- "sc"
  imetric <- match(metrics, tolower(names(cplab)))
  ina <- is.na(imetric)
  if(any(ina)) stop(paste("metric(s) not available:", paste(metrics.orig[ina], collapse = ", ")))
  metrics <- names(cplab)[imetric]

  values <- lapply(metrics, function(metric) {
    getFromNamespace(metric, "canprot")(AAcomp, ...)
  })

  values <- do.call(cbind, values)
  colnames(values) <- metrics
  as.data.frame(values)

}
