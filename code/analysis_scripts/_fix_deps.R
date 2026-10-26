.libPaths(c('E:/workspace/rlibs', .libPaths()))
options(repos=c(CRAN='https://cloud.r-project.org'))
for (i in 1:20) {
  res <- tryCatch({ suppressMessages(suppressWarnings(library(SoupX))); "OK" },
                  error=function(e) conditionMessage(e))
  if (identical(res, "OK")) { cat("SoupX 加载成功\n"); break }
  m <- regmatches(res, regexpr("there is no package called '[^']+'", res))
  if (length(m) == 0) { cat("其他错误:", res, "\n"); break }
  p <- sub(".*'([^']+)'.*", "\\1", m)
  cat(sprintf("缺 %s -> 安装\n", p))
  install.packages(p, lib='E:/workspace/rlibs', quiet=TRUE)
}
cat("最终:", tryCatch({suppressMessages(library(SoupX)); paste("OK", as.character(packageVersion('SoupX')))}, error=function(e) conditionMessage(e)), "\n")
