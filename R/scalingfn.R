scalingfn <-
function(Y,type){
  tY=t(Y)
  if(type=="pareto"){
    sdY=rowVar(tY)^(.25)
    tY=tY/sdY}
  if(type=="unit"){
    sdY=rowVar(tY)^(.5)
    tY=tY/sdY}
  if(type=="none"){tY=tY}
  return(Y=t(tY))}
