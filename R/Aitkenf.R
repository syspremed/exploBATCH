Aitkenf <-
function(ll,lla,v,eps){la=0
tol=eps+1
if(v>2){av=(ll[v]-ll[v-1])/(ll[v-1]-ll[v-2])
lla[v]=ll[v-1]+(ll[v]-ll[v-1])/(1-av)
la=lla[v]
tol=abs(lla[v]-lla[v-1])}
list(tol=tol,la=la)}
