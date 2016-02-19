blm.likelihood <- function(type = c("norm", "t"), df = 3)
{
  type <- match.arg(type)
  list(type = type, errorCov = NULL, df = df)
}


