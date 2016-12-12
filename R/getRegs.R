getTransRegistries = function(chrs=c(6, 17)) {
 tags = lapply(chrs, function(x) gsub("CC", x, "registries/chrCC-export"))
 suppressMessages({
  suppressWarnings({
    ans = lapply(tags, function(x) {
     tmp = loadRegistry(system.file(x, package="gQTLstats"))
     availJobs = as.numeric(dir(file.path(
           system.file(x, package="gQTLstats", "jobs"))))
     tmp$files.dir = system.file(x, package="gQTLstats")
     tmp$work.dir = system.file(x, package="gQTLstats")
     tmp$availJobs = availJobs
     tmp
     })
    })
   })
  ans
}
