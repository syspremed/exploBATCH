iniA <-
function(i,u,Covars,L){
  if(L<2){dat=data.frame(cbind(u[i,],as.matrix(Covars[2:(L+1),])))}
  else{dat=data.frame(cbind(u[i,],t(Covars[2:(L+1),])))}
glm(dat,family=gaussian)$coef}
