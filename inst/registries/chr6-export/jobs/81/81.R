Sys.sleep(0.000000)
options(BatchJobs.on.slave = TRUE, BatchJobs.resources.path = '/udd/stvjc/VM/TRANS_2016/PROCODE/PERM_3/procoMult2016_chr6-files/resources/resources_1481032968.RData')
library(checkmate)
library(BatchJobs)
res = BatchJobs:::doJob(
	reg = loadRegistry('/udd/stvjc/VM/TRANS_2016/PROCODE/PERM_3/procoMult2016_chr6-files'),
	ids = c(81L),
	multiple.result.files = FALSE,
	disable.mail = FALSE,
	first = 4L,
	last = 100L,
	array.id = NA)
BatchJobs:::setOnSlave(FALSE)