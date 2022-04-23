# load('n_lam_out.RData')
# load('Nt_out_c.RData')
# 
# eq.d=N.t-N




required.packages=c("coda",
                    "fBasics",
                    "fields",
                    "ggmap",
                    "ggplot2",
                    "gridExtra",
                    "gstat",
                    "inline",
                    "maptools",
                    "msm",
                    "raster",
                    "rasterVis",
                    "RColorBrewer",
                    "RcppArmadillo",
                    "rgdal",
                    "rgeos")
lapply(required.packages,library,character.only=TRUE)

# library(foreach)
# library(doParallel)
# registerDoParallel(cores=20)

##
## C++ function for propagation
##

## modified with so H gets converted to sparse matrix for multiplication -- way fast!
# ###
# library(devtools)  # not sure this is needed
# Sys.setenv(PATH = paste("C:/Users/jeisaguirre/Documents/rtools40/usr/bin", Sys.getenv("PATH"), sep=";"))
# Sys.setenv(BINPREF = "C:/Users/jeisaguirre/Documents/rtools40/mingw$(WIN)/bin/")
# ###
code <- '
arma::mat Hmat = Rcpp::as<arma::mat>(H);
arma::vec c0vec = Rcpp::as<arma::vec>(c0);
arma::vec gvec = Rcpp::as<arma::vec>(gamma);
arma::vec kvec = Rcpp::as<arma::vec>(ktilde);
double dt = Rcpp::as<double>(deltat);
int n = Rcpp::as<int>(timesteps);
int k = Hmat.n_rows;
arma::sp_mat temp(Hmat);
arma::mat call(k, n);
call.col(0) = temp*c0vec;
for (int i = 1; i < n; ++i) {
    call.col(i) = temp*call.col(i - 1) + dt*gvec%call.col(i - 1) -
                  dt*gvec%(pow(call.col(i - 1), 2)/kvec);
}
return Rcpp::wrap(call);
'
calcc=cxxfunction(signature(H="numeric",
                            c0="numeric",
                            gamma="numeric",
                            ktilde="numeric",
                            timesteps="numeric",
                            deltat="numeric"),
                  body=code,
                  plugin="RcppArmadillo")


## Neighborhood matrix calculator
neighborhood=function(raster,boundary){
  nn=matrix(,length(raster[]),4)
  for(i in 1:dim(nn)[1]){
    if(raster::extract(boundary,i)==1){
      loc=adjacent(raster,i)[,2]
      if(length(loc[which((loc+1)==i)])>0){
        if(raster::extract(boundary,loc[which((loc+1)==i)])==1){
          ln=loc[which((loc+1)==i)]
        }else{ln=0}
      }else{ln=NA}
      
      
      if(length(loc[which((loc-1)==i)])>0){
        if(raster::extract(boundary,loc[which((loc-1)==i)])==1){
          rn=loc[which((loc-1)==i)]
        }else{rn=0}
      }else{rn=NA}
      
      if(length(loc[which((loc-dim(raster)[2])==i)])>0){
        if(raster::extract(
          boundary,
          loc[which((loc-dim(raster)[2])==i)])==1){
          bn=loc[which((loc-dim(raster)[2])==i)]
        }else{bn=0}
      }else{bn=NA}
      
      if(length(loc[which((loc+dim(raster)[2])==i)])>0){
        if(raster::extract(
          boundary,
          loc[which((loc+dim(raster)[2])==i)])==1){
          tn=loc[which((loc+dim(raster)[2])==i)]
        }else{tn=0}
      }else{tn=NA}
      nn[i,1]=ln
      nn[i,2]=rn
      nn[i,3]=bn
      nn[i,4]=tn
    }else{nn[i,]=0}
  }
  nn
}

##
## Propagator matrix
##

propagator.matrix=function(NN,delta,gamma,dx,dy,dt){
  H=matrix(0,dim(NN)[1],dim(NN)[1])
  for(i in 1:dim(H)[1]){
    if(length(which(NN[i,]>0))>0){
      ind.tmp=ifelse(NN[i,]>0,1,0)
      H[i,i]=sum(c(1-2*delta[i]*(dt/dx^2 + dt/dy^2), #+
                   #dt*gamma[i],
                   dt/dx^2*delta[i]*(1-ind.tmp[1]),
                   dt/dx^2*delta[i]*(1-ind.tmp[2]),
                   dt/dy^2*delta[i]*(1-ind.tmp[3]),
                   dt/dy^2*delta[i]*(1-ind.tmp[4])
      ),na.rm=TRUE)
      H[i,NN[i,1]]=dt/dx^2*delta[i]*ind.tmp[1]
      H[i,NN[i,2]]=dt/dx^2*delta[i]*ind.tmp[2]
      H[i,NN[i,3]]=dt/dy^2*delta[i]*ind.tmp[3]
      H[i,NN[i,4]]=dt/dy^2*delta[i]*ind.tmp[4]
    }else{H[i,i]=0}
  }
  H
}



# load(paste0(getwd(),'/seak_out/','beta_out.RData'))
# load(paste0(getwd(),'/seak_out/','alpha_out.RData'))
# load(paste0(getwd(),'/seak_out/','tau_out.RData'))
# load(paste0(getwd(),'/seak_out/','p_out.RData'))
# load(paste0(getwd(),'/seak_out/','kappa_out.RData'))
# load(paste0(getwd(),'/seak_out/','theta_out.RData'))
# load(paste0(getwd(),'/seak_out/','gamma_out.RData'))

load(paste0(getwd(),'/beta_out1.RData'))
load(paste0(getwd(),'/alpha_out1.RData'))
load(paste0(getwd(),'/tau_out1.RData'))
load(paste0(getwd(),'/p_out1.RData'))
load(paste0(getwd(),'/kappa_out1.RData'))
load(paste0(getwd(),'/theta_out1.RData'))
load(paste0(getwd(),'/gamma_out1.RData'))

beta.out=beta
alpha.out=alpha
tau.out=tau
p.out=p
kappa.out=kappa
theta.out=theta
gamma.out=gamma

#load('SEAK_boundary.RData')
load('pde_input_ws.RData')

load('cuc_ind_which.RData')


##
## seed
##
#seed=data$seed
set.seed(23)

##
## Priors
##

gamma.mean=priors$gamma.prior[1]
gamma.var=priors$gamma.prior[2]
beta.mean=priors$beta.prior[1]
beta.var=priors$beta.prior[2]
theta.mean=priors$theta.prior[1:7]
theta.var=priors$theta.prior[8:14]
kappa.mean=priors$kappa.prior[1:7]
kappa.var=priors$kappa.prior[8:14]
q.p=priors$p.prior[1]
r.p=priors$p.prior[2]
p.alpha=priors$p.2017[1]
p.beta=priors$p.2017[2]
q.tau=priors$tau.prior[1]
r.tau=priors$tau.prior[2]
alpha.mean=priors$alpha.prior[1]
alpha.var=priors$alpha.prior[2]


## Boundary layers
##

Boundary=Background$Boundary
BoundaryNA=Background$BoundaryNA
BoundaryInf=Background$BoundaryInf
#BoundaryInf.v=Background$BoundaryInf.v
Boundary.us=Background$Boundary.us

##
## Dimensions
##

dt=st.info$dt
time.frame=st.info$time.frame
us.fact=st.info$us.fact
smooth.fact=st.info$smooth.fact
res=st.info$res
d=st.info$d
time.steps=1/dt*length(time.frame)
keep=seq(1,time.steps,time.steps/length(time.frame))
dx=res*us.fact
dy=res*us.fact
xmin=st.info$extent[1]
xmax=st.info$extent[2]
ymin=st.info$extent[3]
ymax=st.info$extent[4]
y = dim(Boundary)[1]  ## changed to get these straight from raster
x = dim(Boundary)[2]
q=x*y

##
## Years to store lambda
##
years.keep=st.info$years.keep


##
## Data
##

#Y=data$Y
X=data$X
CISU=data$CISU
ISU.i=data$ISU.i
N.vec=data$N.vec
NISU=data$NISU

p.l=list()
for(i in 1:length(attributes(data$ISU)$names)){
  p.l[[i]]=data$ISU[[i]]
}


# Y.2017=data$Y.aerial$Y.2017
# Y.2018=data$Y.aerial$Y.2018
# Y.2019=data$Y.aerial$Y.2019
# 
# if(length(Y.2017)>0){
#   max.Y.2017=suppressWarnings(ifelse(apply(Y.2017,1,max,
#                                            na.rm=TRUE)>0,
#                                      apply(Y.2017,1,max,
#                                            na.rm=TRUE),
#                                      NA))
# }
# if(length(Y.2018)>0){
#   max.Y.2018=suppressWarnings(ifelse(apply(Y.2018,1,max,
#                                            na.rm=TRUE)>0,
#                                      apply(Y.2018,1,max,
#                                      NA))
# }
# if(length(Y.2019)>0){
#   max.Y.2019=suppressWarnings(ifelse(apply(Y.2019,1,max,
#                                            na.rm=TRUE)>0,
#                                      apply(Y.2019,1,max,
#                                            na.rm=TRUE),
#                                      NA))
# }




# polygon that encompasses GLBA
g.p = Polygon(matrix(c(-150000,240000,-40000,270000,
                       -31000,323000,-150000,323000),
                     ncol=2,byrow=T))
glba.ind = cellFromPolygon(Background$bath,
                           SpatialPolygons(list(Polygons(list(g.p),ID='a')),
                                           proj4string = crs(Background$bath)))
glba.ind = unlist(glba.ind)
glba.ind = glba.ind[!is.na(Background$bath[][glba.ind])]
# gdb with fish closures (sea cucs)
fgdb <- "SEAK_SeaCucumber.gdb"
fc <- readOGR(dsn=fgdb)
cuc.clos=fc[fc$Status_Rea=='Closed - officially by regulation (Alaska Board of Fish)' | fc$Status_Rea=='Closed - officially by regulation (Federal government)' & fc$OBJECTID==114
            |fc$OBJECTID==115 |fc$OBJECTID==116 |fc$OBJECTID==117,] # numbered ones are native land S of Ketchikan
cuc.clos.t=spTransform(cuc.clos,crs(raster('SEAK_bath_rot_50.tif')))
fgdb <- "SEAK_Geoduck_RedUrchin.gdb"
fc <- readOGR(dsn=fgdb,layer="SEAK_RedSeaUrchin")
urch.clos=fc[fc$CurrentS_1=='Closed - officially by regulation (Alaska Board of Fish)' | fc$CurrentS_1=='Closed - officially by regulation (Federal Government)',]
urch.clos.t=spTransform(urch.clos,crs(raster('SEAK_bath_rot_50.tif')))
fc <- readOGR(dsn=fgdb,layer="SEAK_Geoduck")
geo.clos=fc[fc$Rotation=='Control Site' |fc$Rotation2=='Control Site' ,]
geo.clos.t=spTransform(geo.clos,crs(raster('SEAK_bath_rot_50.tif')))
fish.ind1=cellFromPolygon(Background$bath,cuc.clos.t)
fish.ind2=cellFromPolygon(Background$bath,urch.clos.t)
fish.ind3=cellFromPolygon(Background$bath,geo.clos.t)
fish.ind1=unlist(fish.ind1)
fish.ind2=unlist(fish.ind2)
fish.ind3=unlist(fish.ind3)
fish.ind=c(fish.ind1,fish.ind2,fish.ind3)
fish.ind=unique(fish.ind)
fish.ind = fish.ind[!is.na(Background$bath[][fish.ind])]

cell=raster(,nrows=y,ncols=x,xmn=xmin,xmx=xmax,
            ymn=ymin,ymx=ymax,crs=NA)
cell[]=1:q


### containers
N.cuc=array(,dim=c(length(time.frame),length(cuc.which),dim(beta.out)[1]))
Nt.cuc=matrix(,length(cuc.which),dim(beta.out)[1])
eqd.cuc=array(,dim=c(length(time.frame),length(cuc.which),dim(beta.out)[1]))

# load('eqd_cuc.RData')
# load('N_cuc.RData')
# load('Nt_cuc.RData')

n.iter=dim(beta.out)[1]

### computre pde and derived quantities
tock=Sys.time()
#for(k in 1:n.iter){
for(k in 1:3300){

#for(k in 1:10){
  
  gamma=gamma.out[k,]
  beta=beta.out[k,]
  theta=theta.out[k,]
  kappa=kappa.out[k,]
  p.ind=p.out[k,]
  p=rep(p.ind[1],q)
  for(i in 2:length(p.ind)){
    p.tmp=rep(p.ind[i],q)
    p=c(p,p.tmp)
  }
  tau=tau.out[k,]
  alpha=alpha.out[k,]
  delta=raster(,nrows=y,ncols=x,xmn=xmin,xmx=xmax,
               ymn=ymin,ymx=ymax,crs=NA)
  delta[]=exp(X%*%beta)
  delta.inv=1/delta
  delta.inv=focal(delta.inv, 
                  w=matrix(rep(1/smooth.fact^2, smooth.fact^2), 
                           nrow=smooth.fact), 
                  na.rm=T, pad=T, padValue=0)
  delta.inv2=1/delta^2
  delta.inv2=focal(delta.inv2, 
                   w=matrix(rep(1/smooth.fact^2, smooth.fact^2), 
                            nrow=smooth.fact, ncol=smooth.fact), 
                   na.rm=T, pad=T, padValue=0)
  delta.star=delta
  delta.bar=aggregate(1/delta.inv, fact=us.fact, 
                      na.rm=TRUE, fun=function(x, na.rm) x[length(x)/2]) 
  
  
  
  
  
  gamma.r=delta
  gamma.r[]=gamma
  gamma.r.star=gamma.r
  gamma.bar=aggregate(gamma.r,
                      fact=us.fact,
                      na.rm=TRUE,
                      fun=mean)
  K.bar.tmp=exp(alpha[1])*delta.inv/delta.inv2
  K.bar.tmp[glba.ind]=exp(alpha[1]+alpha[2])*delta.inv[glba.ind]/delta.inv2[glba.ind]
  K.bar.tmp[fish.ind]=exp(alpha[1]+alpha[3])*delta.inv[fish.ind]/delta.inv2[fish.ind]
  K.bar <- aggregate(K.bar.tmp,
                     fact=us.fact, na.rm=TRUE,
                     fun=function(x, na.rm) x[length(x)/2])
  
  if(!exists('nt.sum')){NN=neighborhood(delta.bar,Boundary.us)}
  H=propagator.matrix(NN=NN,
                      delta=delta.bar[],
                      gamma=gamma.bar[],
                      dx=dx,dy=dy,dt=dt)
  
  lambda0=raster(,nrows=y,ncols=x,xmn=xmin,xmx=xmax,
                 ymn=ymin,ymx=ymax,crs=NA)
  c0=raster(,nrows=y/us.fact,ncols=x/us.fact,
            xmn=xmin,xmx=xmax,ymn=ymin,ymx=ymax,crs=NA)
  D = list()
  lam0 = list()
  D[[1]] = rdist(data.frame(SpatialPoints(lambda0)),matrix(d[[1]],1,2,byrow=T))/1000
  lam0[[1]] = (exp(-D[[1]]^2/kappa[1]^2)/sum(exp(-D[[1]]^2/kappa[1]^2))*theta[1])
  for(i in 2:length(d)){
    D[[i]] = rdist(data.frame(SpatialPoints(lambda0)),matrix(d[[i]],1,2,byrow=T))/1000
    lam0[[i]] = (exp(-D[[i]]^2/kappa[i]^2)/sum(exp(-D[[i]]^2/kappa[i]^2))*theta[i])
  }
  us.cells=aggregate(cell,fact=us.fact,fun=mean)
  us.cells[]=1:length(us.cells[])
  lambda0[]=Reduce('+',lam0)
  
  #lambda0.star=lambda0
  c0[]=raster::extract(delta*lambda0,SpatialPoints(us.cells))
  c0.star=c0
  
  ## c.all and lambda.all and N.all
  c.all=brick(nrows=dim(c0)[1],ncols=dim(c0)[2],xmn=xmin,xmx=xmax,
              ymn=ymin,ymx=ymax)
  c.all=setValues(c.all,
                  calcc(H,
                        vec(c0[]),
                        vec(gamma.bar[]),
                        vec(K.bar[]),
                        time.steps,
                        dt)[,keep])
  lambda=vec(disaggregate(c.all,us.fact)/delta)
  N=rnbinom(length(lambda),size=tau,,mu=lambda)
  N.t=suppressWarnings(rnbinom(n=q,  # equilibrium abundance
                               size=tau,
                               ,mu=vec(disaggregate(K.bar,us.fact)*
                                         delta.inv*Background$BoundaryInf)))
  
  ## current abundance
  N.21=N[(length(lambda)-q+1):length(lambda)]
  if(k==1){N.21.m=N.21/n.iter}  ## running mean
  if(k>1){N.21.m=N.21.m+N.21/n.iter}
  
  ## equilibrium abundance
  if(k==1){N.t.m=N.t/n.iter}  ## running mean
  if(k>1){N.t.m=N.t.m+N.t/n.iter}
  
  ## equilibrium differential across study area
  eqd=N.t-N.21
  if(k==1){eqd.m=eqd/n.iter}  ## running mean
  if(k>1){eqd.m=eqd.m+eqd/n.iter}
  
  ## fisheries units derived quantities
  for(kk in 1:length(cuc.which)){
    Nt.cuc[kk,k]=sum(N.t[cuc.which[[kk]]],na.rm=TRUE)
    for(tt in 1:length(time.frame)){
      N.cuc[tt,kk,k]=sum(N[cuc.which[[kk]]+(tt-1)*q],na.rm=TRUE)
      eqd.cuc[tt,kk,k]=sum((Nt.cuc[kk,k]-N.cuc[tt,kk,k]),na.rm=TRUE)
    }
    
  }
  if(k%%100==0){print(k)}
  if(k%%100==0){
    save(Nt.cuc,file='Nt_cuc1.RData')
    save(N.cuc,file='N_cuc1.RData')
    save(eqd.cuc,file='eqd_cuc1.RData')
    save(N.t.m,file='Ntm_seak1.RData')
    save(eqd.m,file='eqdm_seak1.RData')
    save(N.21.m,file='N21m_seak1.RData')
    }
}
tick=Sys.time()
tick-tock
