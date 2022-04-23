#########################################################################
#########################################################################
###
### Function to fit the Glacier Bay model to data using MCMC
###
#########################################################################
#########################################################################

MCMC=function(data,
              priors,
              inits,
              parameters,
              st.info,
              Background,
              n.iter,
              checkpoint,
              savepoint,
              output.location){
    
    ##
    ##  Subroutines and Packages
    ##
    
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
    
    ##
    ## C++ function for propagation
    ##
    
    ## modified with so H gets converted to sparse matrix for multiplication -- way fast!
    ## updated with logistic version (Perry gave me exp version)
    
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

##
## seed
##
seed=data$seed
set.seed(seed)

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
# K.alpha=priors$K.prior[1]
# K.beta=priors$K.prior[2]
# K.g.alpha=priors$K.g.prior[1]
# K.g.beta=priors$K.g.prior[2]
# K.f.alpha=priors$K.f.prior[1]
# K.f.beta=priors$K.f.prior[2]


##

## Boundary layers
##

Boundary=Background$Boundary
BoundaryNA=Background$BoundaryNA
BoundaryInf=Background$BoundaryInf
BoundaryInf.v=Background$BoundaryInf.v
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
# y=(ymax-ymin)/res
# x=(xmax-xmin)/res
# q=x*y
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

Y=data$Y
X=data$X
CISU=data$CISU
ISU.i=data$ISU.i
N.vec=data$N.vec
NISU=data$NISU

p.l=list()
for(i in 1:length(attributes(data$ISU)$names)){
    p.l[[i]]=data$ISU[[i]]
}


Y.2017=data$Y.aerial$Y.2017
Y.2018=data$Y.aerial$Y.2018
Y.2019=data$Y.aerial$Y.2019

if(length(Y.2017)>0){
    max.Y.2017=suppressWarnings(ifelse(apply(Y.2017,1,max,
                                             na.rm=TRUE)>0,
                                       apply(Y.2017,1,max,
                                             na.rm=TRUE),
                                       NA))
}
if(length(Y.2018)>0){
    max.Y.2018=suppressWarnings(ifelse(apply(Y.2018,1,max,
                                             na.rm=TRUE)>0,
                                       apply(Y.2018,1,max,
                                             na.rm=TRUE),
                                       NA))
}
if(length(Y.2019)>0){
    max.Y.2019=suppressWarnings(ifelse(apply(Y.2019,1,max,
                                             na.rm=TRUE)>0,
                                       apply(Y.2019,1,max,
                                             na.rm=TRUE),
                                       NA))
}

##
## Starting values
##

gamma=inits$gamma
beta=inits$beta
theta=inits$theta
kappa=inits$kappa
p.ind=inits$p.ind
p=rep(p.ind[1],q)
for(i in 2:length(p.ind)){
    p.tmp=rep(p.ind[i],q)
    p=c(p,p.tmp)
}
tau=inits$tau
# K=inits$K
# K.g=inits$K.g
# K.f=inits$K.f
alpha=inits$alpha
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
# delta.bar=aggregate(delta,fact=us.fact,
#                     fun=function(x,na.rm){
#                         (1/mean(1/x,na.rm=TRUE))
#                     })
delta.bar=aggregate(1/delta.inv, fact=us.fact, 
                    na.rm=TRUE, fun=function(x, na.rm) x[length(x)/2]) 
## make this different for glba and elsewhere
# gamma.bar=delta.bar*aggregate(gamma/delta,
#                               fact=us.fact,
#                               fun=mean,
#                               na.rm=TRUE
# )
cell=raster(,nrows=y,ncols=x,xmn=xmin,xmx=xmax,
            ymn=ymin,ymx=ymax,crs=NA)
cell[]=1:q

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

NN=neighborhood(delta.bar,Boundary.us)
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
# for(i in 1:length(d)){
#     lam0[[i]][X[,3]>-0.3400708 & X[,3]!=0] = 0  # restricts intensity to within 5km of shore
# }
#D=rdist(data.frame(SpatialPoints(lambda0)),matrix(d,1,2))/1000
us.cells=aggregate(cell,fact=us.fact,fun=mean)
us.cells[]=1:length(us.cells[])
lambda0[]=Reduce('+',lam0)
#lambda0[]=lambda0[]*(sum(theta)/sum(lambda0[]))  # this scales intensity surface back up so it integrates to ~413 after multiplying by Boundary (inherent when delta*lambda0)  
#(quick and dirty truncation; make better later); probably doesn't make much of a difference

#lambda0[]=(exp(-D^2/kappa^2)/sum(exp(-D^2/kappa^2))*theta)
lambda0.star=lambda0
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
Y=Y*rep(vec(Background$BoundaryNA),20)  # remove otters on land (because bathymetry layer sucks)

N[(23*q+1):(23*q+q*20)]=ifelse(N[(23*q+1):(23*q+q*20)]<Y,
                               Y,
                               N[(23*q+1):(23*q+q*20)])
#N[(23*q+1):(23*q+q*20)][Y==0]=0
if(length(Y.2017)>0){
    N[(47*q+1):(48*q)][!is.na(max.Y.2017)]=
        ifelse(N[(47*q+1):(48*q)][!is.na(max.Y.2017)]<
                   max.Y.2017[!is.na(max.Y.2017)],
               max.Y.2017[!is.na(max.Y.2017)]+1,
               N[(47*q+1):(48*q)][!is.na(max.Y.2017)]
        )
}
if(length(Y.2018)>0){
    N[(48*q+1):(49*q)][!is.na(max.Y.2018)]=
        ifelse(N[(48*q+1):(49*q)][!is.na(max.Y.2018)]<
                   max.Y.2018[!is.na(max.Y.2018)],
               max.Y.2018[!is.na(max.Y.2018)]+1,
               N[(48*q+1):(49*q)][!is.na(max.Y.2018)]
        )
}
if(length(Y.2019)>0){
    N[(49*q+1):(50*q)][!is.na(max.Y.2019)]=
        ifelse(N[(49*q+1):(50*q)][!is.na(max.Y.2019)]<
                   max.Y.2019[!is.na(max.Y.2019)],
               max.Y.2019[!is.na(max.Y.2019)]+1,
               N[(49*q+1):(50*q)][!is.na(max.Y.2019)]
        )
}

# ##
# ## index for where Y observed
# ##
# 
# ind.y=rep(FALSE,q*length(keep))
# ind.y[(23*q+1):(23*q+q*20)]=!is.na(Y)
# ind.y[(47*q+1):(48*q)]=apply(Y.2017,1,function(x,xx=dim(Y.2017)[2]) sum(is.na((x)))<xx)
# ind.y[(48*q+1):(49*q)]=apply(Y.2018,1,function(x,xx=dim(Y.2018)[2]) sum(is.na((x)))<xx)
# ind.y[(49*q+1):(50*q)]=apply(Y.2019,1,function(x,xx=dim(Y.2019)[2]) sum(is.na((x)))<xx)
# 
# n.tot.v=suppressWarnings(rnbinom(n=length(lambda),
#                                  size=tau,
#                                  ,mu=lambda*BoundaryInf.v
# ))
# N[!ind.y]=n.tot.v[!ind.y]*BoundaryInf.v[!ind.y]  # update N where no Y with pp


## Create indicator vector for ISU cells
ind=1:length(N)
# ISU.ind=matrix(0,nr=length(N),nc=1)
# ISU.ind[ISU]=1

## Identify which cells were: 1) observed to be occupied
##                            2) not in ISU sample
# N.occ.tmp=ind[Y==1&ISU.ind==0]
# N.occ=N.occ.tmp[!is.na(N.occ.tmp)]

## Identify which cells were: 1) observed to be unoccupied
##                            2) not in ISU sample
# N.abb.tmp = ind[Y==0&ISU.ind==0]
# N.abb=N.abb.tmp[!is.na(N.abb.tmp)]

## Identify which sites were surveyed 
SurvSites=ind[(23*q+1):(23*q+q*20)][!is.na(Y)]
SurvSites=c(SurvSites,ind[(47*q+1):(48*q)][!(apply(Y.2017,1,function(x)all(is.na(x))))])
SurvSites=c(SurvSites,ind[(48*q+1):(49*q)][!(apply(Y.2018,1,function(x)all(is.na(x))))])
SurvSites=c(SurvSites,ind[(49*q+1):(50*q)][!(apply(Y.2019,1,function(x)all(is.na(x))))])

## ISU site index
ISU.i=ISU.i+23*q

## fix N at ISU sites
N[ISU.i]=NISU

## index for sites to update with pp
pp.i=ind[-SurvSites]
pp.i=pp.i[!(pp.i%in%ISU.i)]


##
## Tuning parameters
##

gamma.tune=0.0005097501
gamma.tune=.001
#gamma.tune=0.01
# beta.tune=c(0.02125764,
#             0.04374,
#             0.045,
#             0.045,
#             0.01945226,
#             0.002658881,
#             0.05,
#             0.05,
#             0.05)
# beta.tune=c(.001,
#             .001,
#             .001,
#             .001,
#             .001,
#             .001,
#             .001,
#             .001)
# alpha.tune=c(.01,
#              .1,
#              .1)
# theta.tune=c(15,
#              1,
#              1,
#              15,
#              1,
#              15,
#              15)
# kappa.tune=c(1.5,
#              1.5,
#              1.5,
#              .15,
#              1.5,
#              1.5,
#              1.5)
# tau.tune=0.0012
beta.tune=c(.005,
            .005,
            .005,
            .005,
            .005,
            .005,
            .005,
            .005)
alpha.tune=c(.04,
             .3,
             .4)
theta.tune=c(50,
             3,
             3,
             65,
             3,
             65,
             35)
kappa.tune=c(4,
             1.5,
             6,
             .25,
             1.5,
             5,
             2)
tau.tune=0.0013
# K.tune=0.005
# K.g.tune=.85
# K.f.tune=.5
accept.gamma=0
accept.beta=rep(0,length(beta))
accept.alpha=rep(0,length(alpha))
accept.theta=rep(0,length(theta))
accept.kappa=rep(0,length(kappa))
accept.N=rep(0,length(Y))
accept.tau=0
accept.K=0
accept.K.g=0
accept.K.f=0
H.check=0
lam=list()
N.t=0

##
## Containers
##

MCMC.Chains=vector('list',length(parameters))
names(MCMC.Chains)=parameters
#ind.est <- seq(floor(n.iter/2), n.iter, by=10)
ind.est <- seq(100, n.iter, by=100)
if ('K'%in%parameters){
    MCMC.Chains$K=matrix(,n.iter,1)
}
if ('K.g'%in%parameters){
    MCMC.Chains$K.g=matrix(,n.iter,1)
}
if ('K.f'%in%parameters){
    MCMC.Chains$K.f=matrix(,n.iter,1)
}
if('gamma'%in%parameters){
    MCMC.Chains$gamma=matrix(,n.iter,1)
}
if('beta'%in%parameters){
    MCMC.Chains$beta=matrix(,n.iter,length(beta))
}
if('alpha'%in%parameters){
    MCMC.Chains$alpha=matrix(,n.iter,length(alpha))
}
if('theta'%in%parameters){
    MCMC.Chains$theta=matrix(,n.iter,length(theta))
}
if('kappa'%in%parameters){
    MCMC.Chains$kappa=matrix(,n.iter,length(kappa))
}
if('p'%in%parameters){
    MCMC.Chains$p=matrix(,n.iter,sum(!is.na(p.ind)))
}
if('tau'%in%parameters){
    MCMC.Chains$tau=matrix(,n.iter,1)
}
if('n.tot'%in%parameters){
    MCMC.Chains$n.tot=matrix(,n.iter,length(keep))
}
if('n.tot.v'%in%parameters){
    MCMC.Chains$n.tot.v=matrix(,q,n.iter)
}
if('lambda'%in%parameters){
    MCMC.Chains$lambda=matrix(,q,n.iter)
}
if('lambda11'%in%parameters){
    MCMC.Chains$lambda11=matrix(,q,n.iter)
}
if('lambda.all'%in%parameters){
    MCMC.Chains$lambda.all=matrix(,q*length(years.keep),n.iter)
}
if('lambda.tot'%in%parameters){
    MCMC.Chains$lambda.tot=matrix(NA,nrow=n.iter, ncol=length(keep))
}
if('score'%in%parameters){
    MCMC.Chains$score=matrix(,n.iter,1)
}
if('tuners'%in%parameters){
    MCMC.Chains$tuners=matrix(NA,nrow=n.iter,ncol=length(c(beta.tune,
                                                           gamma.tune,
                                                           alpha.tune,
                                                           # K.tune,
                                                           # K.g.tune,
                                                           # K.f.tune,
                                                           theta.tune,
                                                           kappa.tune,
                                                           tau.tune)))
}
if('N'%in%parameters){
    MCMC.Chains$N=matrix(NA,nrow=q,ncol=length(time.frame))
    MCMC.Chains$N=lambda/length(ind.est)
}
if('pp'%in%parameters){
    MCMC.Chains$pp=matrix(,length(SurvSites),n.iter)
}
if('pp'%in%parameters){
    MCMC.Chains$E.y=matrix(,length(SurvSites),n.iter)
}
if('N.t'%in%parameters){
    MCMC.Chains$N.t=matrix(,q,n.iter)
}
if('H.check'%in%parameters){
    MCMC.Chains$H.check=matrix(,n.iter,1)
}
if('K.bar'%in%parameters){
    MCMC.Chains$K.bar=matrix(,ncell(K.bar),n.iter)
}
if('delta.bar'%in%parameters){
    MCMC.Chains$delta.bar=matrix(,ncell(delta.bar),n.iter)
}

##
## Gibbs loop
##

# seed=data$seed
# set.seed(seed)
#sys_t = Sys.time()
for(k in 1:n.iter){
    
    ##
    ## Sample gamma
    ##
    
    gamma.star=rnorm(1,gamma,gamma.tune)
    #if(gamma.star> q.gamma & gamma.star< r.gamma){
    # delta.gamma = gamma.star/delta
    # gamma.bar.star=delta.bar*aggregate(delta.gamma,
    #                                    fact=us.fact,
    #                                    fun=mean,
    #                                    na.rm=TRUE)
    gamma.r.star[]=gamma.star
    gamma.bar.star=aggregate(gamma.r.star,
                             fact=us.fact,
                             na.rm=TRUE,
                             fun=mean)
    H.star=propagator.matrix(NN,delta.bar[],gamma.bar.star[],
                             dx=dx,dy=dy,dt=dt)
    #if(min(range(H.star,na.rm=TRUE))>=0){
    c.all.star=setValues(c.all,
                         calcc(H.star,vec(c0[]),
                               vec(gamma.bar.star[]),
                               vec(K.bar[]),time.steps,dt)[,keep]
    )
    if(min(c.all.star[[50]][],na.rm=TRUE)>=0){
        lambda.star=vec(disaggregate(c.all.star,us.fact)/
                            delta)
        mh1=sum(foo=dnbinom(x=N[SurvSites],size=tau,,mu=lambda.star[SurvSites],log=TRUE),  # only evaluate at points where there's observations
                na.rm=TRUE)+
            dnorm(gamma.star,gamma.mean,gamma.var^.5,log=TRUE)
        mh2=sum(dnbinom(x=N[SurvSites],size=tau,,mu=lambda[SurvSites],log=TRUE),
                na.rm=TRUE)+
            dnorm(gamma,gamma.mean,gamma.var^.5,log=TRUE)
        mh=exp(mh1-mh2)
        if(mh>runif(1)){
            gamma=gamma.star
            gamma.bar=gamma.bar.star
            H=H.star
            c.all=c.all.star
            lambda=lambda.star
            accept.gamma=accept.gamma+1
        }
    }
    #}
    
    
    
    
    ##
    ## Sample beta
    ##
    
    for(i in 1:length(beta)){
        #for(i in 1:1){
        beta.tmp=rnorm(1,mean=beta[i],sd=beta.tune[i])
        beta.star=beta
        beta.star[i]=beta.tmp
        delta.star[]=exp(X%*%beta.star)
        delta.inv.star=1/delta.star
        delta.inv.star=focal(delta.inv.star,
                             w=matrix(rep(1/smooth.fact^2, smooth.fact^2),
                                      nrow=smooth.fact),
                             na.rm=T, pad=T, padValue=0)
        delta.inv2.star=1/delta.star^2
        delta.inv2.star=focal(delta.inv2.star,
                              w=matrix(rep(1/smooth.fact^2, smooth.fact^2),
                                       nrow=smooth.fact, ncol=smooth.fact),
                              na.rm=T, pad=T, padValue=0)
        # K.bar.star.tmp=K*delta.inv.star/delta.inv2.star
        # K.bar.star.tmp[glba.ind]=K.g*delta.inv.star[glba.ind]/delta.inv2.star[glba.ind]
        # K.bar.star.tmp[fish.ind]=K.f*delta.inv.star[fish.ind]/delta.inv2.star[fish.ind]
        K.bar.star.tmp=exp(alpha[1])*delta.inv/delta.inv2
        K.bar.star.tmp[glba.ind]=exp(alpha[1]+alpha[2])*delta.inv[glba.ind]/delta.inv2[glba.ind]
        K.bar.star.tmp[fish.ind]=exp(alpha[1]+alpha[3])*delta.inv[fish.ind]/delta.inv2[fish.ind]
        K.bar.star <- aggregate(K.bar.star.tmp,
                                fact=us.fact, na.rm=TRUE,
                                fun=function(x, na.rm) x[length(x)/2])
        delta.bar.star=aggregate(1/delta.inv.star,
                                 fact=us.fact,
                                 fun=function(x,na.rm){
                                     x[length(x)/2]
                                 }
        )
        gamma.r[]=gamma
        gamma.bar=aggregate(gamma.r,
                            fact=us.fact,
                            na.rm=TRUE,
                            fun=mean)
        H.star=propagator.matrix(NN,delta.bar.star[],
                                 gamma.bar[],
                                 dx=dx,dy=dy,dt=dt)
        #if(min(H.star,na.rm=TRUE)<0){H.check=H.check+1}
        if(min(H.star,na.rm=TRUE)>=0){
            c0.star[]=raster::extract(delta.star*lambda0,
                                      SpatialPoints(us.cells))
            c.all.star=setValues(c.all,
                                 calcc(H.star,vec(c0.star[]),
                                       vec(gamma.bar[]),
                                       vec(K.bar.star[]),
                                       time.steps,dt)[,keep]
            )
            if(min(c.all.star[[50]][],na.rm=TRUE)>=0){
                lambda.star=vec(disaggregate(c.all.star,us.fact)/
                                    delta.star)
                mh1=sum(dnbinom(x=N[SurvSites],size=tau,,mu=lambda.star[SurvSites],log=TRUE),
                        na.rm=TRUE)+
                    sum(dnorm(beta.star,beta.mean,beta.var^0.5,log=TRUE))
                mh2=sum(dnbinom(x=N[SurvSites],size=tau,,mu=lambda[SurvSites],log=TRUE),
                        na.rm=TRUE)+
                    sum(dnorm(beta,beta.mean,beta.var^0.5,log=TRUE))
                mh=exp(mh1-mh2)
                if(mh>runif(1)){
                    beta=beta.star
                    delta=delta.star
                    delta.bar=delta.bar.star
                    gamma.bar=gamma.bar.star
                    K.bar=K.bar.star
                    delta.inv=delta.inv.star
                    delta.inv2=delta.inv2.star
                    H=H.star
                    c0=c0.star
                    c.all=c.all.star
                    lambda=lambda.star
                    accept.beta[i]=accept.beta[i]+1
                }
            }
        }
    }
    
    ##
    ## Sample theta
    ##
    
    for(i in 1:length(theta)){
        theta.star.tmp=rnorm(1,theta[i],theta.tune[i])
        theta.star=theta
        theta.star[i]=theta.star.tmp
        if(theta.star[i]>1){  # for computation
            for(ii in 1:length(d)){
                lam0[[ii]] = (exp(-D[[ii]]^2/kappa[ii]^2)/sum(exp(-D[[ii]]^2/kappa[ii]^2))*theta.star[ii])
            }
            lambda0.star[]=Reduce('+',lam0)
            c0.star=raster::extract(delta*lambda0.star,
                                    SpatialPoints(us.cells))
            c.all.star=setValues(c.all,
                                 calcc(H,vec(c0.star[]),
                                       vec(gamma.bar[]),vec(K.bar[]),
                                       time.steps,dt)[,keep]
            )
            if(min(c.all.star[[50]][],na.rm=TRUE)>=0){
                lambda.star=vec(disaggregate(c.all.star,us.fact)/
                                    delta)
                mh1=sum(dnbinom(x=N[SurvSites],size=tau,,mu=lambda.star[SurvSites],log=TRUE),
                        na.rm=TRUE)+
                    sum(dnorm(theta.star,theta.mean,theta.var^0.5,log=TRUE))
                mh2=sum(dnbinom(x=N[SurvSites],size=tau,,mu=lambda[SurvSites],log=TRUE),
                        na.rm=TRUE)+
                    sum(dnorm(theta,theta.mean,theta.var^0.5,log=TRUE))
                mh=exp(mh1-mh2)
                if(mh>runif(1)){
                    theta=theta.star
                    lambda0=lambda0.star
                    c0=c0.star
                    c.all=c.all.star
                    lambda=lambda.star
                    accept.theta[i]=accept.theta[i]+1
                }
            }
        }
    }
    
    ##
    ## Sample kappa
    ##
    
    for(i in 1:length(kappa)){
        kappa.star.tmp=rnorm(1,kappa[i],kappa.tune[i])
        kappa.star=kappa
        kappa.star[i]=kappa.star.tmp
        if(kappa.star[i]>.5){  # .5 instaed of 0 for computational stability
            for(ii in 1:length(d)){
                lam0[[ii]] = (exp(-D[[ii]]^2/kappa.star[ii]^2)/sum(exp(-D[[ii]]^2/kappa.star[ii]^2))*theta[ii])
            }
            lambda0.star[]=Reduce('+',lam0)
            c0.star=raster::extract(delta*lambda0.star,
                                    SpatialPoints(us.cells))
            c.all.star=setValues(c.all,
                                 calcc(H,vec(c0.star[]),
                                       vec(gamma.bar[]),vec(K.bar[]),
                                       time.steps,dt)[,keep]
            )
            if(min(c.all.star[[50]][],na.rm=TRUE)>=0){
                lambda.star=vec(disaggregate(c.all.star,us.fact)/
                                    delta)
                mh1=sum(dnbinom(x=N[SurvSites],size=tau,,mu=lambda.star[SurvSites],log=TRUE),
                        na.rm=TRUE)+
                    sum(dnorm(kappa.star,kappa.mean,kappa.var^0.5,log=TRUE))
                mh2=sum(dnbinom(x=N[SurvSites],size=tau,,mu=lambda[SurvSites],log=TRUE),
                        na.rm=TRUE)+
                    sum(dnorm(kappa,kappa.mean,kappa.var^0.5,log=TRUE))
                mh=exp(mh1-mh2)
                if(mh>runif(1)){
                    kappa=kappa.star
                    lambda0=lambda0.star
                    c0=c0.star
                    c.all=c.all.star
                    lambda=lambda.star
                    accept.kappa[i]=accept.kappa[i]+1
                }
            }
        }
    }
    
    ##
    ## Sample N (1993-2012)
    ## 24 years where we don't have data (note q*20 for 20 years of data)
    
    N.star=sample(-1:1,length(Y),replace=TRUE)+N[(23*q+1):(23*q+q*20)]
    N.star=ifelse(N.star<0,0,N.star)
    mh1=dbinom(Y,N.star,p[(23*q+1):(23*q+q*20)],log=TRUE)+   # NA's are in Y where's no data so only evaluated where there's data
        dnbinom(N.star,tau,,lambda[(23*q+1):(23*q+q*20)],log=TRUE)  # added to vector with NA's where no Y, so all good
    mh2=dbinom(Y,N[(23*q+1):(23*q+q*20)],prob=p[(23*q+1):(23*q+q*20)],log=TRUE)+
        dnbinom(N[(23*q+1):(23*q+q*20)],tau,,lambda[(23*q+1):(23*q+q*20)],log=TRUE)
    mh=exp(mh1-mh2)
    ru.tmp=runif(length(mh))
    N[(23*q+1):(23*q+q*20)]=ifelse(mh>ru.tmp,
                                   N.star,
                                   N[(23*q+1):(23*q+q*20)])
    # accept.N[(23*q+1):(23*q+q*20)]=ifelse(mh>ru.tmp,
    #                                       accept.N[(23*q+1):(23*q+q*20)]+1,
    #                                       accept.N[(23*q+1):(23*q+q*20)])
    
    ##
    ## Sample N (2017)
    ##
    
    if(length(Y.2017)>0){
        N.star=sample(-1:1,dim(Y.2017)[1],replace=TRUE)+
            N[(47*q+1):(48*q)]
        N.star=ifelse(N.star<0,0,N.star)
        mh1=apply(dbinom(Y.2017,
                         N.star,
                         p[(47*q+1):(48*q)],
                         log=TRUE),
                  1,sum,na.rm=TRUE)+
            dnbinom(N.star,size=tau,
                    ,mu=(60*90)/(400*400)*lambda[(47*q+1):(48*q)],
                    log=TRUE)
        mh2=apply(dbinom(Y.2017,
                         N[(47*q+1):(48*q)],
                         p[(47*q+1):(48*q)],
                         log=TRUE),
                  1,sum,na.rm=TRUE)+
            dnbinom(N[(47*q+1):(48*q)],size=tau,
                    ,mu=(60*90)/(400*400)*lambda[(47*q+1):(48*q)],
                    log=TRUE)
        mh=exp(mh1-mh2)
        mh[is.na(mh)]=0
        rcut=runif(length(mh))
        N[(47*q+1):(48*q)]=ifelse(mh>rcut,N.star,N[(47*q+1):(48*q)])
        # accept.N[(47*q+1):(48*q)]=ifelse(mh>rcut,
        #                                  accept.N[(47*q+1):(48*q)]+1,
        #                                  accept.N[(47*q+1):(48*q)])
    }
    
    ##
    ## Sample N (2018)
    ##
    
    if(length(Y.2018)>0){
        N.star=sample(-1:1,dim(Y.2018)[1],replace=TRUE)+
            N[(48*q+1):(49*q)]
        N.star=ifelse(N.star<0,0,N.star)
        mh1=apply(dbinom(Y.2018,
                         N.star,
                         p[(48*q+1):(49*q)],
                         log=TRUE),
                  1,sum,na.rm=TRUE)+
            dnbinom(N.star,size=tau,
                    ,mu=(60*90)/(400*400)*lambda[(48*q+1):(49*q)],
                    log=TRUE)
        mh2=apply(dbinom(Y.2018,
                         N[(48*q+1):(49*q)],
                         p[(48*q+1):(49*q)],
                         log=TRUE),
                  1,sum,na.rm=TRUE)+
            dnbinom(N[(48*q+1):(49*q)],size=tau,
                    ,mu=(60*90)/(400*400)*lambda[(48*q+1):(49*q)],
                    log=TRUE)
        mh=exp(mh1-mh2)
        mh[is.na(mh)]=0
        rcut=runif(length(mh))
        N[(48*q+1):(49*q)]=ifelse(mh>rcut,N.star,N[(48*q+1):(49*q)])
        # accept.N[(48*q+1):(49*q)]=ifelse(mh>rcut,
        #                                  accept.N[(48*q+1):(49*q)]+1,
        #                                  accept.N[(48*q+1):(49*q)])
    }
    
    ##
    ## Sample N (2019)
    ##
    
    if(length(Y.2019)>0){
        N.star=sample(-1:1,dim(Y.2019)[1],replace=TRUE)+
            N[(49*q+1):(50*q)]
        N.star=ifelse(N.star<0,0,N.star)
        mh1=apply(dbinom(Y.2019,
                         N.star,
                         p[(49*q+1):(50*q)],
                         log=TRUE),
                  1,sum,na.rm=TRUE)+
            dnbinom(N.star,size=tau,
                    ,mu=(60*90)/(400*400)*lambda[(49*q+1):(50*q)],
                    log=TRUE)
        mh2=apply(dbinom(Y.2019,
                         N[(49*q+1):(50*q)],
                         p[(49*q+1):(50*q)],
                         log=TRUE),
                  1,sum,na.rm=TRUE)+
            dnbinom(N[(49*q+1):(50*q)],size=tau,
                    ,mu=(60*90)/(400*400)*lambda[(49*q+1):(50*q)],
                    log=TRUE)
        mh=exp(mh1-mh2)
        mh[is.na(mh)]=0
        rcut=runif(length(mh))
        N[(49*q+1):(50*q)]=ifelse(mh>rcut,N.star,N[(49*q+1):(50*q)])
        # accept.N[(48*q+1):(49*q)]=ifelse(mh>rcut,
        #                                  accept.N[(48*q+1):(49*q)]+1,
        #                                  accept.N[(48*q+1):(49*q)])
    }
    
    ## keep these fixed
    N[ISU.i]=NISU
    
    ##
    ## Sample tau
    ##
    
    tau.star=rnorm(1,tau,tau.tune)
    if(tau.star>q.tau&tau.star<r.tau){
        mh1=sum(dnbinom(x=N[SurvSites],size=tau.star,,mu=lambda[SurvSites],log=TRUE),
                na.rm=TRUE)
        mh2=sum(dnbinom(x=N[SurvSites],size=tau,,mu=lambda[SurvSites],log=TRUE),
                na.rm=TRUE)
        mh=exp(mh1-mh2)
        if(mh>runif(1)){
            tau=tau.star
            accept.tau=accept.tau+1
        }
    }
    
    ##
    ## Sample p 
    ##
    # ISU data
    for(i in 1:(sum(!is.na(p.ind))-3)){
        p.ind[!is.na(p.ind)][i]=rbeta(1,
                                      sum(p.l[[i*2-1]])+q.p,
                                      sum(p.l[[i*2]]-p.l[[i*2-1]])+
                                          r.p
        )
    }
    # camera data, just update with moment-matched prior
    p.ind[!is.na(p.ind)][12]=rbeta(1, p.alpha, p.beta)
    p.ind[!is.na(p.ind)][13]=rbeta(1, p.alpha, p.beta)
    p.ind[!is.na(p.ind)][14]=rbeta(1, p.alpha, p.beta)
    
    p=rep(p.ind[1],q)
    for(i in 2:length(p.ind)){
        p.tmp=rep(p.ind[i],q)
        p=c(p,p.tmp)
    }
    
    ##
    ## Sample alpha_0
    ##
    
    alpha0.star=rnorm(1,mean=alpha[1],sd=alpha.tune[1])
    #if(K.star>=K.alpha&K.star<=K.beta) {
    #if(K.star>=0){
    # K.bar.star.tmp=K.star*delta.inv/delta.inv2
    # K.bar.star.tmp[glba.ind]=K.g*delta.inv[glba.ind]/delta.inv2[glba.ind]
    # K.bar.star.tmp[fish.ind]=K.f*delta.inv[fish.ind]/delta.inv2[fish.ind]
    K.bar.star.tmp=exp(alpha0.star)*delta.inv/delta.inv2
    K.bar.star.tmp[glba.ind]=exp(alpha0.star+alpha[2])*delta.inv[glba.ind]/delta.inv2[glba.ind]
    K.bar.star.tmp[fish.ind]=exp(alpha0.star+alpha[3])*delta.inv[fish.ind]/delta.inv2[fish.ind]
    K.bar.star <- aggregate(K.bar.star.tmp,
                            fact=us.fact, na.rm=TRUE,
                            fun=function(x, na.rm) x[length(x)/2])
    c.all.star=setValues(c.all,
                         calcc(H,vec(c0[]),vec(gamma.bar[]),
                               vec(K.bar.star[]),
                               time.steps,dt)[,keep]
    )
    if(min(c.all.star[[50]][],na.rm=TRUE)>=0){
        lambda.star=vec(disaggregate(c.all.star,us.fact)/
                            delta)
        mh1=sum(dnbinom(x=N[SurvSites],size=tau,,mu=lambda.star[SurvSites],log=TRUE),
                na.rm=TRUE)+
            dnorm(alpha0.star,alpha.mean,alpha.var^.5,log=TRUE)
        mh2=sum(dnbinom(x=N[SurvSites],size=tau,,mu=lambda[SurvSites],log=TRUE),
                na.rm=TRUE)+
            dnorm(alpha[1],alpha.mean,alpha.var^.5,log=TRUE)
        mh=exp(mh1-mh2)
        if(mh>runif(1)){
            alpha[1]=alpha0.star
            K.bar=K.bar.star
            c.all=c.all.star
            lambda=lambda.star
            accept.alpha[1]=accept.alpha[1]+1
        }
    }
    #}
    
    ##
    ## Sample alpha_1
    ##
    
    alpha1.star=rnorm(1,mean=alpha[2],sd=alpha.tune[2])
    #if(K.star>=K.alpha&K.star<=K.beta) {
    #if(K.star>=0){
    # K.bar.star.tmp=K.star*delta.inv/delta.inv2
    # K.bar.star.tmp[glba.ind]=K.g*delta.inv[glba.ind]/delta.inv2[glba.ind]
    # K.bar.star.tmp[fish.ind]=K.f*delta.inv[fish.ind]/delta.inv2[fish.ind]
    K.bar.star.tmp=exp(alpha[1])*delta.inv/delta.inv2
    K.bar.star.tmp[glba.ind]=exp(alpha[1]+alpha1.star)*delta.inv[glba.ind]/delta.inv2[glba.ind]
    K.bar.star.tmp[fish.ind]=exp(alpha[1]+alpha[3])*delta.inv[fish.ind]/delta.inv2[fish.ind]
    K.bar.star <- aggregate(K.bar.star.tmp,
                            fact=us.fact, na.rm=TRUE,
                            fun=function(x, na.rm) x[length(x)/2])
    c.all.star=setValues(c.all,
                         calcc(H,vec(c0[]),vec(gamma.bar[]),
                               vec(K.bar.star[]),
                               time.steps,dt)[,keep]
    )
    if(min(c.all.star[[50]][],na.rm=TRUE)>=0){
        lambda.star=vec(disaggregate(c.all.star,us.fact)/
                            delta)
        mh1=sum(dnbinom(x=N[SurvSites],size=tau,,mu=lambda.star[SurvSites],log=TRUE),
                na.rm=TRUE)+
            dnorm(alpha1.star,alpha.mean,alpha.var^.5,log=TRUE)
        mh2=sum(dnbinom(x=N[SurvSites],size=tau,,mu=lambda[SurvSites],log=TRUE),
                na.rm=TRUE)+
            dnorm(alpha[2],alpha.mean,alpha.var^.5,log=TRUE)
        mh=exp(mh1-mh2)
        if(mh>runif(1)){
            alpha[2]=alpha1.star
            K.bar=K.bar.star
            c.all=c.all.star
            lambda=lambda.star
            accept.alpha[2]=accept.alpha[2]+1
        }
    }
    #}
    
    ##
    ## Sample alpha_2
    ##
    
    alpha2.star=rnorm(1,mean=alpha[3],sd=alpha.tune[3])
    #if(K.star>=K.alpha&K.star<=K.beta) {
    #if(K.star>=0){
    # K.bar.star.tmp=K.star*delta.inv/delta.inv2
    # K.bar.star.tmp[glba.ind]=K.g*delta.inv[glba.ind]/delta.inv2[glba.ind]
    # K.bar.star.tmp[fish.ind]=K.f*delta.inv[fish.ind]/delta.inv2[fish.ind]
    K.bar.star.tmp=exp(alpha[1])*delta.inv/delta.inv2
    K.bar.star.tmp[glba.ind]=exp(alpha[1]+alpha[2])*delta.inv[glba.ind]/delta.inv2[glba.ind]
    K.bar.star.tmp[fish.ind]=exp(alpha[1]+alpha2.star)*delta.inv[fish.ind]/delta.inv2[fish.ind]
    K.bar.star <- aggregate(K.bar.star.tmp,
                            fact=us.fact, na.rm=TRUE,
                            fun=function(x, na.rm) x[length(x)/2])
    c.all.star=setValues(c.all,
                         calcc(H,vec(c0[]),vec(gamma.bar[]),
                               vec(K.bar.star[]),
                               time.steps,dt)[,keep]
    )
    if(min(c.all.star[[50]][],na.rm=TRUE)>=0){
        lambda.star=vec(disaggregate(c.all.star,us.fact)/
                            delta)
        mh1=sum(dnbinom(x=N[SurvSites],size=tau,,mu=lambda.star[SurvSites],log=TRUE),
                na.rm=TRUE)+
            dnorm(alpha2.star,alpha.mean,alpha.var^.5,log=TRUE)
        mh2=sum(dnbinom(x=N[SurvSites],size=tau,,mu=lambda[SurvSites],log=TRUE),
                na.rm=TRUE)+
            dnorm(alpha[3],alpha.mean,alpha.var^.5,log=TRUE)
        mh=exp(mh1-mh2)
        if(mh>runif(1)){
            alpha[3]=alpha2.star
            K.bar=K.bar.star
            c.all=c.all.star
            lambda=lambda.star
            accept.alpha[3]=accept.alpha[3]+1
        }
    }
    #}
    
    # ##
    # ## Sample K
    # ##
    # 
    # K.star=rnorm(1,mean=K,sd=K.tune)
    # #if(K.star>=K.alpha&K.star<=K.beta) {
    # if(K.star>=0){
    #     K.bar.star.tmp=K.star*delta.inv/delta.inv2
    #     K.bar.star.tmp[glba.ind]=K.g*delta.inv[glba.ind]/delta.inv2[glba.ind]
    #     K.bar.star.tmp[fish.ind]=K.f*delta.inv[fish.ind]/delta.inv2[fish.ind]
    #     K.bar.star <- aggregate(K.bar.star.tmp,
    #                             fact=us.fact, na.rm=TRUE,
    #                             fun=function(x, na.rm) x[length(x)/2])
    #     c.all.star=setValues(c.all,
    #                          calcc(H,vec(c0[]),vec(gamma.bar[]),
    #                                vec(K.bar.star[]),
    #                                time.steps,dt)[,keep]
    #     )
    #     if(min(c.all.star[[50]][],na.rm=TRUE)>=0){
    #         lambda.star=vec(disaggregate(c.all.star,us.fact)/
    #                             delta)
    #         mh1=sum(dnbinom(x=N[SurvSites],size=tau,,mu=lambda.star[SurvSites],log=TRUE),
    #                 na.rm=TRUE)+
    #             dnorm(K.star,K.alpha,K.beta^.5,log=TRUE)
    #         mh2=sum(dnbinom(x=N[SurvSites],size=tau,,mu=lambda[SurvSites],log=TRUE),
    #                 na.rm=TRUE)+
    #             dnorm(K,K.alpha,K.beta^.5,log=TRUE)
    #         mh=exp(mh1-mh2)
    #         if(mh>runif(1)){
    #             K=K.star
    #             K.bar=K.bar.star
    #             c.all=c.all.star
    #             lambda=lambda.star
    #             accept.K=accept.K+1
    #         }
    #     }
    # }
    # 
    # ##
    # ## Sample K.g
    # ##
    # 
    # K.g.star=rnorm(1,mean=K.g,sd=K.g.tune)
    # #if(K.g.star>=K.g.alpha&K.g.star<=K.g.beta) {
    # if(K.g.star>=0){
    #     K.bar.star.tmp=K*delta.inv/delta.inv2
    #     K.bar.star.tmp[glba.ind]=K.g.star*delta.inv[glba.ind]/delta.inv2[glba.ind]
    #     K.bar.star.tmp[fish.ind]=K.f*delta.inv[fish.ind]/delta.inv2[fish.ind]
    #     K.bar.star <- aggregate(K.bar.star.tmp,
    #                             fact=us.fact, na.rm=TRUE,
    #                             fun=function(x, na.rm) x[length(x)/2])
    #     c.all.star=setValues(c.all,
    #                          calcc(H,vec(c0[]),vec(gamma.bar[]),
    #                                vec(K.bar.star[]),
    #                                time.steps,dt)[,keep]
    #     )
    #     if(min(c.all.star[[50]][],na.rm=TRUE)>=0){
    #         lambda.star=vec(disaggregate(c.all.star,us.fact)/
    #                             delta)
    #         mh1=sum(dnbinom(x=N[SurvSites],size=tau,,mu=lambda.star[SurvSites],log=TRUE),
    #                 na.rm=TRUE)+
    #             dnorm(K.g.star,K.g.alpha,K.g.beta^.5,log=TRUE)
    #         mh2=sum(dnbinom(x=N[SurvSites],size=tau,,mu=lambda[SurvSites],log=TRUE),
    #                 na.rm=TRUE)+
    #             dnorm(K.g,K.g.alpha,K.g.beta^.5,log=TRUE)
    #         mh=exp(mh1-mh2)
    #         if(mh>runif(1)){
    #             K.g=K.g.star
    #             K.bar=K.bar.star
    #             c.all=c.all.star
    #             lambda=lambda.star
    #             accept.K.g=accept.K.g+1
    #         }
    #     }
    # }
    # 
    # ##
    # ## Sample K.f
    # ##
    # 
    # K.f.star=rnorm(1,mean=K.f,sd=K.f.tune)
    # #if(K.f.star>=K.f.alpha&K.f.star<=K.f.beta) {
    # if(K.f.star>=0){
    #     K.bar.star.tmp=K*delta.inv/delta.inv2
    #     K.bar.star.tmp[glba.ind]=K.g*delta.inv[glba.ind]/delta.inv2[glba.ind]
    #     K.bar.star.tmp[fish.ind]=K.f.star*delta.inv[fish.ind]/delta.inv2[fish.ind]
    #     K.bar.star <- aggregate(K.bar.star.tmp,
    #                             fact=us.fact, na.rm=TRUE,
    #                             fun=function(x, na.rm) x[length(x)/2])
    #     c.all.star=setValues(c.all,
    #                          calcc(H,vec(c0[]),vec(gamma.bar[]),
    #                                vec(K.bar.star[]),
    #                                time.steps,dt)[,keep]
    #     )
    #     if(min(c.all.star[[50]][],na.rm=TRUE)>=0){
    #         lambda.star=vec(disaggregate(c.all.star,us.fact)/
    #                             delta)
    #         mh1=sum(dnbinom(x=N,size=tau,,mu=lambda.star,log=TRUE),
    #                 na.rm=TRUE)+
    #             dnorm(K.f.star,K.f.alpha,K.f.beta^.5,log=TRUE)
    #         mh2=sum(dnbinom(x=N,size=tau,,mu=lambda,log=TRUE),
    #                 na.rm=TRUE)+
    #             dnorm(K.f,K.f.alpha,K.f.beta^.5,log=TRUE)
    #         mh=exp(mh1-mh2)
    #         if(mh>runif(1)){
    #             K.f=K.f.star
    #             K.bar=K.bar.star
    #             c.all=c.all.star
    #             lambda=lambda.star
    #             accept.K.f=accept.K.f+1
    #         }
    #     }
    # }
    # 
    ##
    ## Derived parameters
    ##
    
    n.tot.v=suppressWarnings(rnbinom(n=length(lambda),
                                     size=tau,
                                     ,mu=lambda*BoundaryInf.v
    ))
    N[pp.i]=n.tot.v[pp.i]*BoundaryInf.v[pp.i]  # update N where no Y with pp
    # n.tot=unname(
    #     tapply(n.tot.v,(seq_along(n.tot.v)-1)%/%q,sum,na.rm=TRUE)
    # )
    n.tot=unname(
        tapply(N,(seq_along(N)-1)%/%q,sum,na.rm=TRUE)
    )
    if('N.t'%in%parameters){
        N.t=suppressWarnings(rnbinom(n=q,  # equilibrium abundance
                                     size=tau,
                                     ,mu=vec(disaggregate(K.bar,us.fact)*
                                                 delta.inv*Background$BoundaryInf)))
    }
    
    ##
    ## Posterior predictive of Y where Y observed
    ##
    
    pp=suppressWarnings(rbinom(length(N[SurvSites]),N[SurvSites],p[SurvSites]))
    E.y=N[SurvSites]*p[SurvSites]  # expected N
    
    ##
    ## Store results
    ##
    
    if('K'%in%parameters){
        MCMC.Chains$K[k,]=K
    }
    if('K.g'%in%parameters){
        MCMC.Chains$K.g[k,]=K.g
    }
    if('K.f'%in%parameters){
        MCMC.Chains$K.f[k,]=K.f
    }
    if('gamma'%in%parameters){
        MCMC.Chains$gamma[k,]=gamma
    }
    if('beta'%in%parameters){
        MCMC.Chains$beta[k,]=beta
    }
    if('alpha'%in%parameters){
        MCMC.Chains$alpha[k,]=alpha
    }
    if('theta'%in%parameters){
        MCMC.Chains$theta[k,]=theta
    }
    if('kappa'%in%parameters){
        MCMC.Chains$kappa[k,]=kappa
    }
    if('p'%in%parameters){
        MCMC.Chains$p[k,]=p.ind[!is.na(p.ind)]
    }
    if('tau'%in%parameters){
        MCMC.Chains$tau[k,]=tau
    }
    if('n.tot'%in%parameters){
        MCMC.Chains$n.tot[k,]=n.tot
    }
    if('n.tot.v'%in%parameters){
        MCMC.Chains$n.tot.v[,k]=n.tot.v[((length(keep)-1)*q+1):
                                            (length(keep)*q)]
    }
    if('lambda'%in%parameters){
        MCMC.Chains$lambda[,k]=lambda[((length(keep)-1)*q+1):
                                          (length(keep)*q)]
    }
    if('lambda.all'%in%parameters){
        MCMC.Chains$lambda.all[,k]=c(lambda[1:(years.keep[1]*q)],
                                     lambda[((years.keep[2]-1)*q+1):(years.keep[2]*q)],
                                     lambda[((years.keep[3]-1)*q+1):(years.keep[3]*q)],
                                     lambda[((years.keep[4]-1)*q+1):(years.keep[4]*q)],
                                     lambda[((years.keep[5]-1)*q+1):(years.keep[5]*q)],
                                     lambda[((years.keep[6]-1)*q+1):(years.keep[6]*q)],
                                     lambda[((years.keep[7]-1)*q+1):(years.keep[7]*q)],
                                     lambda[((years.keep[8]-1)*q+1):(years.keep[8]*q)],
                                     lambda[((years.keep[9]-1)*q+1):(years.keep[9]*q)],
                                     lambda[((years.keep[10]-1)*q+1):(years.keep[10]*q)],
                                     lambda[((years.keep[11]-1)*q+1):(years.keep[11]*q)]
        )
    }
    if('lambda11'%in%parameters){
        MCMC.Chains$lambda11[,k]=lambda[((length(keep)-10)*q+1):
                                            ((length(keep)-9)*q)]
    }
    if('tuners'%in%parameters){
        MCMC.Chains$tuners[k,]=c(gamma.tune,
                                 beta.tune,
                                 alpha.tune,
                                 theta.tune,
                                 kappa.tune,
                                 tau.tune)
    }
    if('H.check'%in%parameters){
        MCMC.Chains$H.check[k,]=H.check
    }
    if('pp'%in%parameters){
        MCMC.Chains$pp[,k]=pp
    }
    if('pp'%in%parameters){
        MCMC.Chains$E.y[,k]=E.y
    }
    if('N.t'%in%parameters){
        MCMC.Chains$N.t[,k]=N.t
    }
    if (('N'%in%parameters) & (k%in%ind.est)) {
        MCMC.Chains$N <- MCMC.Chains$N + lambda/length(ind.est)
    }
    if('K.bar'%in%parameters){
        MCMC.Chains$K.bar[,k]=K.bar[]
    }
    if('delta.bar'%in%parameters){
        MCMC.Chains$delta.bar[,k]=delta.bar[]
    }
    
    ##
    ## Checkpoint
    ##
    
    # if(k%%checkpoint==0){
    #     
    #     ##
    #     ## Update tuning parameters
    #     ##
    #     
    #     if(accept.gamma/k<0.3){
    #         gamma.tune=gamma.tune*0.9
    #     }
    #     if(accept.gamma/k>0.5){
    #         gamma.tune=gamma.tune*1.1
    #     }
    #     beta.tune=ifelse(accept.beta/k<0.3,
    #                      beta.tune*0.9,
    #                      ifelse(accept.beta/k>0.5,
    #                             beta.tune*1.1,
    #                             beta.tune))
    #     kappa.tune=ifelse(accept.kappa/k<0.3,
    #                       kappa.tune*0.9,
    #                       ifelse(accept.kappa/k>0.5,
    #                              kappa.tune*1.1,
    #                              kappa.tune))
    #     theta.tune=ifelse(accept.theta/k<0.3,
    #                       theta.tune*0.9,
    #                       ifelse(accept.theta/k>0.5,
    #                              theta.tune*1.1,
    #                              theta.tune))
    #     if(accept.tau/k<0.3){
    #         tau.tune=tau.tune*0.9
    #     }
    #     if(accept.tau/k>0.5){
    #         tau.tune=tau.tune*1.1
    #     }
    #     if(accept.alpha[1]/k<0.3){
    #         alpha.tune[1]=alpha.tune[1]*0.9
    #     }
    #     if(accept.alpha[1]/k>0.5){
    #         alpha.tune[1]=alpha.tune[1]*1.1
    #     }
    #     if(accept.alpha[2]/k<0.3){
    #         alpha.tune[2]=alpha.tune[2]*0.9
    #     }
    #     if(accept.alpha[2]/k>0.5){
    #         alpha.tune[2]=alpha.tune[2]*1.1
    #     }
    #     if(accept.alpha[3]/k<0.3){
    #         alpha.tune[3]=alpha.tune[3]*0.9
    #     }
    #     if(accept.alpha[3]/k>0.5){
    #         alpha.tune[3]=alpha.tune[3]*1.1
    #     }
    #     
    # }
    
    ##
    ## Save to list
    ##
    
    if(k%%savepoint==0){
        save(MCMC.Chains,file=output.location)
        write.csv('',paste0(output.location,k,'.csv'))
    }
    
}

}

