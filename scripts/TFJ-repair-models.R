# This script renames the "Single talker" condition to "Single-talker" in all previously run models.
# It doesn't change any of the results. It only relables the condition to facilitate consistent plotting
# and tablurizing.
file_list <- list.files(path="models/", pattern = "^exp1.*.rds", full.names = T)

for (i in file_list) {
  print(i)
  m = readRDS(i) 
  
  if (is.brmsfit(m)) {
    m$data %<>%
      mutate_all(.funs = function(x) { 
        if (length(which(x == "Single talker")) == 0) return(x) 
        dimnames(contrasts(x))[[1]] = rev(levels.exposure)
        c = contrasts(x)
        x = ifelse(x == "Single talker", "Single-talker", if(is.factor(x)) as.character(x) else x)
        x = factor(x, levels = rev(levels.exposure))
        contrasts(x) = c
        return(x)
      })
  }
  
  saveRDS(m, i)
}
