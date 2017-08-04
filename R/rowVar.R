rowVar <-
function(x){rowSums((x-rowMeans(x))^2)/(dim(x)[2]-1)}
