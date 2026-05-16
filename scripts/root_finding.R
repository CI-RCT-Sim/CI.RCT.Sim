#Use while waiting for the simulations to finish.
x<-round(runif(1)*10)+1
stop<-FALSE
while(!stop) {
  zz<-readline("An welche Zahl denke ich? ")
  z<-as.numeric(zz)
  if(z<x) print("Zu klein!")
  if(z>x) print("Zu groß!")
  if(x==z) {
    print("Super, du hast die Zahl erraten!")
    stop<-TRUE
  }
  if(zz=="stop") stop<-TRUE
}

