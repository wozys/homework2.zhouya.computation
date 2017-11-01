library(jpeg)
img = readJPEG("./灰度图.jpg")
img = img[,,1]
par(mar = c(0,0,0,0))
plot(0,0,type = "n")
rasterImage(img, -1, -1, 1, 1)

set.seed(100)
a<-matrix(,642,500)
for(i in 1:642){
  a[i,]=rep(i,500)
}
b<-matrix(,642,500)
for(i in 1:500){
  b[,i]=rep(i,642)
}
a1<-sample(a,642*500*0.6)
b1<-sample(b,642*500*0.6)
c1<-matrix(c(a1,b1),642*500*0.6,2)
Xomega<-matrix(0,642,500)

X=img
for(i in 1:192600){
  Xomega[c1[i,1],c1[i,2]]=X[c1[i,1],c1[i,2]]
}

rasterImage(Xomega, -1, -1, 1, 1)


img2<-homework2(Xomega,c1,0.1)

rasterImage(img2, -1, -1, 1, 1)
