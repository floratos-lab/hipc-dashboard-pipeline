library(HGNChelper)

new.hgnc.table <- HGNChelper::getCurrentHumanMap()
save(new.hgnc.table, file = "../data/reference_files/hgnc_complete_set.RData", compress="bzip2")
