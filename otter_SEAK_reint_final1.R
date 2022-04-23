rm(list=ls())

########################################################################
########################################################################
########################################################################
###
### Real Data Analysis
###
########################################################################
########################################################################
########################################################################

###
### Packages and dependencies
###

required.packages=c("coda",
                    "fBasics",
                    "fields",
                    "ggmap",
                    "ggplot2",
                    "gridExtra",
                    "mapview",  ## added
                    "gstat",
                    "inline",
                    "maptools",
                    "parallel",
                    "raster",
                    "rasterVis",
                    "RCurl",
                    "RColorBrewer",
                    "RcppArmadillo",
                    "rgdal",
                    "rgeos",
                    'sf',
                    'gdistance')
lapply(required.packages,library,character.only=TRUE)


###
### Function  to rotate a raster
###

# rotateProj=function(spobj,angle){
#   ## get bounding box as spatialpoints object
#   boxpts=SpatialPoints(t(bbox(spobj)),
#                        proj4string=CRS(proj4string(spobj)))
#   ## convert to lat-long
#   boxLL=bbox(spTransform(boxpts,CRS("+init=epsg:4326")))
#   ## find the center
#   llc=apply(boxLL,1,mean)
#   ## construct the proj4string
#   prj=paste0("+proj=omerc +lat_0=",llc[2]," +lonc=",llc[1],
#              " +alpha=",
#              angle,
#              paste0(" +gamma=0.0 +k=1.000000 +x_0=0.000 ",
#                     "+y_0=0.000 +ellps=WGS84 +units=m "))
#   ## return as a CRS:
#   CRS(prj)
# }


###
### Wrapper for parallel processing
###

run.chain.2pl.list=function(data=data,
                            priors=priors,
                            inits=inits,
                            parameters=parameters,
                            st.info=st.info,
                            Background=Background,
                            n.iter=n.iter,
                            checkpoint=checkpoint,
                            savepoint=savepoint,
                            out.loc=out.loc.l
){
  chain.list=mclapply(1:length(out.loc),
                      function(j){
                        this.data=data
                        this.out.loc=out.loc[[j]]
                        this.data$seed=data$seed[j]
                        ## Run the chain on this core
                        this.chain=MCMC(data=this.data,
                                        priors=priors,
                                        inits=inits,
                                        parameters=parameters,
                                        st.info=st.info,
                                        Background=Background,
                                        n.iter=n.iter,
                                        checkpoint=checkpoint,
                                        savepoint=savepoint,
                                        output.location=this.out.loc
                        )
                        return(this.chain)
                      },
                      mc.cores=min(length(data$seed),
                                   detectCores())
  )
  
  ## Save the initial random seed as the name of the chain
  names(chain.list)=data$seed
  return(chain.list)
}


########################################################################
### Load Bathymetry Data
########################################################################

#bath.raster=raster(paste0("~/Dropbox/Projects/SeaOtters/",
#                         "OriginalDataFiles/glb_bath/dblbnd.adf"))

#bath.raster = raster('Raster_400m_Bathy_10_11_18_Clip_UTM8_NAD27.tif')
#proj = crs(raster(paste0("~/Google Drive/otters/glb_bath/dblbnd.adf")))
#bath.raster = projectRaster(bath.raster, crs=proj)

bath.raster = raster('SEAK_bath_rot_50.tif')  # this raster was created from UTM projection (from Tom Diltz) and projected
# and roated 50 degrees so SEAK is horizontal


###
### Study area extent
###

#xmin=-421000-800-800-1600
xmin=-421000+800
xmax=480000
ymin=114800
ymax=390000+800
#ymin=114800+4000

###
### Initial population intensity in 1969
###

d=list()
d[[1]]=c(290000,215000)  # Maurelle Islands
d[[2]]=c(414000,224000)  # Barrier Islands
d[[3]]=c(120000,200000)  # Necker Islands
d[[4]]=c(27000,219000)   # Khaz Bay
d[[5]]=c(-287000,208100) # Yakutat Bay
d[[6]]=c(-12000,219000)  # Yakobi Island (south of Surge Bay)
d[[7]]=c(-45000,237000)  # Cape Spencer
#d[[8]]=c(-40000,280000)  # Glacier bay
#d=c(20000,200000)

#theta = c(50,85,20,110,20,110,20) # 154 released in d[1-3], 223 released d[4] and d[6]
theta = c(51,55,48,194,10,30,25)*1
#theta=c(166.503,195.8774,155.2982,605.1446,34.18157,102.7206,89.20142)

# kappa = c(20,20,20,20,20,20,20)
# kappa=c(20.10774,19.3128,20.85343,20.1137,19.78165,19.45277,19.8688)

theta=c(60,5,24,50,5,30,30)
kappa=c(15,3,6,6,1.25,3.75,3.125)
theta=c(30,2.5,12,25,2.5,15,15)
kappa=theta/4/2
kappa[4]=12
kappa=c(15,3,6,6,1.25,3.75,3.125)
kappa=rep(3,7)*2
theta=round(theta/2)
# ### changed to 1975 from tinker et al 2019.
# 
# d=list()
# d[[1]]=c(290000,215000)  # Maurelle Islands
# d[[2]]=c(414000,224000)  # Barrier Islands
# d[[3]]=c(120000,200000)  # Necker Islands
# d[[4]]=c(27000,219000)   # Khaz Bay
# d[[5]]=c(-287000,208100) # Yakutat Bay
# d[[6]]=c(-12000,219000)  # Yakobi Island (south of Surge Bay)
# d[[7]]=c(-55000,232000)  # Cape Spencer
# d[[8]]=c(245000,200000)  # coronation island
# d[[9]]=c(75000,201500)   # N04 from tinker
# d[[10]]=c(-40000,270000) # Icy straight, near mouth of GLBA
# #d[[8]]=c(-40000,280000)  # Glacier bay
# #d=c(20000,200000)
# 
# #theta = c(50,85,20,110,20,110,20) # 154 released in d[1-3], 223 released d[4] and d[6]
# theta = c(46,19,40,150,10,30,25,74,10,80)
# #kappa = c(100,100,100,100,100,100,100)
# kappa = c(10,10,10,10,10,10,10,10,10,10)*5
# 


###
### Discretization values
###

dt=1/365
#time.frame=1993:2019
time.frame=1970:2021
time.steps=1/dt*length(time.frame)
keep=seq(1,time.steps,time.steps/length(time.frame))
res=400

###
### Homogenization value
###

us.fact=10
smooth.fact=15

###
### Bathymetry data
###

bath=crop(bath.raster,extent(xmin,xmax,ymin,ymax))

###
### Add Boulder Island
###

# just use this to find boulder island
Bath.g = '+proj=utm +zone=8 +datum=NAD27 +units=m +no_defs +ellps=clrk66
+nadgrids=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat'


xbi=440626.7
ybi=6491505.7
bi = SpatialPoints(matrix(c(xbi,ybi), ncol=2), CRS(Bath.g))
bi_r = spTransform(bi, crs(bath))
bit=cellFromXY(bath,cbind(bi_r@coords[,1],bi_r@coords[,2]))
bath[bit]=NA


x=dim(bath)[2]  ## moved these
y=dim(bath)[1]
q=x*y
proj = crs(bath)

########################################################################
### Load SE AK Sea Otter Data
########################################################################

data=read.csv(paste0("SEAK_otter.csv"))
dataDistr=read.csv(paste0("DistribSurveysC.csv"))
ind.pred.tmp=1:dim(dataDistr)[1]
ind.pred=ind.pred.tmp[dataDistr$year==2004|dataDistr$year==2006]  ### cuts out years when both distr/transects flown
dataDistr=dataDistr[-ind.pred,]
ind=1993:2012

### added by JME ###
### project distr data
dist_pts = spTransform(SpatialPoints(matrix(c(dataDistr$POINT_X, dataDistr$POINT_Y),ncol=2), CRS(Bath.g)), crs(bath))
dataDistr$POINT_X = dist_pts@coords[,1]
dataDistr$POINT_Y = dist_pts@coords[,2]

#### added by JME
## remove distribution data from years with transects for SEAK and boat data, and 
# removed distr data from transect years in other script
# remove SEAK data that overlaps with GLBA
gl_poly_02 = with(data.frame(X=data$GROUP_X[data$year==2002 & data$region=='GLBA' & !is.na(data$GROUP_Y)],
                             Y=data$GROUP_Y[data$year==2002 & data$region=='GLBA' & !is.na(data$GROUP_Y)]), chull(X,Y))
gl_coords_02 = data.frame(X=data$GROUP_X[data$year==2002 & data$region=='GLBA' & !is.na(data$GROUP_Y)],
                          Y=data$GROUP_Y[data$year==2002 & data$region=='GLBA' & !is.na(data$GROUP_Y)])[gl_poly_02,]
gl_poly_03 = with(data.frame(X=data$GROUP_X[data$year==2003 & data$region=='GLBA' & !is.na(data$GROUP_Y)],
                             Y=data$GROUP_Y[data$year==2003 & data$region=='GLBA' & !is.na(data$GROUP_Y)]), chull(X,Y))
gl_coords_03 = data.frame(X=data$GROUP_X[data$year==2003 & data$region=='GLBA' & !is.na(data$GROUP_Y)],
                          Y=data$GROUP_Y[data$year==2003 & data$region=='GLBA' & !is.na(data$GROUP_Y)])[gl_poly_03,]

seak_pts_02_e = matrix(c(data$east_long[ data$year==2002 & data$region=='SEAK'], 
                         data$east_lat[ data$year==2002 & data$region=='SEAK']),ncol=2)
seak_pts_02_w = matrix(c(data$west_long[ data$year==2002 & data$region=='SEAK'], 
                         data$west_lat[ data$year==2002 & data$region=='SEAK']),ncol=2)
seak_pts_03_e = matrix(c(data$east_long[ data$year==2003 & data$region=='SEAK'], 
                         data$east_lat[ data$year==2003 & data$region=='SEAK']),ncol=2)
seak_pts_03_w = matrix(c(data$west_long[ data$year==2003 & data$region=='SEAK'], 
                         data$west_lat[ data$year==2003 & data$region=='SEAK']),ncol=2)
pts_keep_02_e = seak_pts_02_e[point.in.polygon(seak_pts_02_e[,1], seak_pts_02_e[,2], gl_coords_02$X, gl_coords_02$Y)==1,]
pts_keep_02_w = seak_pts_02_w[point.in.polygon(seak_pts_02_w[,1], seak_pts_02_w[,2], gl_coords_02$X, gl_coords_02$Y)!=1,]
pts_keep_03_e = seak_pts_03_e[point.in.polygon(seak_pts_03_e[,1], seak_pts_03_e[,2], gl_coords_03$X, gl_coords_03$Y)==1,]
pts_keep_03_w = seak_pts_03_w[point.in.polygon(seak_pts_03_w[,1], seak_pts_03_w[,2], gl_coords_03$X, gl_coords_03$Y)==1,]

data = data[!(data$east_lat%in%data$east_lat[data$year==2002 & data$region=='SEAK'][point.in.polygon(seak_pts_02_e[,1], seak_pts_02_e[,2], gl_coords_02$X, gl_coords_02$Y)==1]),]
data = data[!(data$east_lat%in%data$east_lat[data$year==2003 & data$region=='SEAK'][point.in.polygon(seak_pts_02_e[,1], seak_pts_02_e[,2], gl_coords_02$X, gl_coords_02$Y)==1]),]

#######


###
### Convert Data to Raster
###

cell=raster(,extent(bath), 
            crs=proj, 
            ncols=x, nrows=y, 
            xmx=xmax(bath), xmn=xmin(bath), ymx=ymax(bath), ymn=ymin(bath))  ## modified to match bath
CISU=N.r=
  ISUind.r=
  Y.r=Counts=
  Boundary=
  BoundaryDist=
  DistCov=
  DepthCov=
  gamma=
  delta=lambda0=SurvR=SimulatedDataR=cell
c0=cell  ## changed these
data.r=cell
data1.sub=subset(data,year=="1993")
data1.sub=data1.sub[!is.na(data1.sub$GROUP_X),]
data2.sub=subset(dataDistr,year=="1993")
x_loc=c(data1.sub$Group_X,data2.sub$POINT_X)
y_loc=c(data1.sub$Group_Y,data2.sub$POINT_Y)
Counts.tmp=c(data1.sub$animals,data2.sub$animals)
Counts.sub=rasterize(x=cbind(x_loc,y_loc),y=data.r,
                     field=Counts.tmp,fun="max",na.rm=TRUE,background=0)
Counts[]=Counts.sub[]
Counts=stack(mget(rep("Counts",25)))

years=1994:2012
ind=2
for(t in years){
  data1.sub=subset(data,year==t)
  data1.sub=data1.sub[!is.na(data1.sub$GROUP_X),]
  data2.sub=subset(dataDistr,year==t)
  
  x_loc=c(data1.sub$GROUP_X,data2.sub$POINT_X)
  y_loc=c(data1.sub$GROUP_Y,data2.sub$POINT_Y)
  Counts.tmp=c(data1.sub$animals,data2.sub$animals)
  if(length(x_loc)>0){
    Counts.sub=rasterize(x=cbind(x_loc,y_loc),
                         y=data.r,field=Counts.tmp,
                         fun="max",na.rm=TRUE,
                         background=0)
    Counts[[ind]][]=Counts.sub[]
  }else{
    Counts[[ind]][]=rep(NA,q)
  }
  ind=ind+1
}

########################################################################
### Add 2017 aerial photograph data
########################################################################

rand.2017=read.csv(paste0("SO_D_20170719_R.csv"))
rand.2017=rand.2017[!is.na(rand.2017$LATITUDE_WGS84),]
opt.2017=read.csv(paste0("SO_D_20170721_O.csv"))
opt.2017=opt.2017[!is.na(opt.2017$LATITUDE_WGS84),]
abund.2017=read.csv(paste0("SO_D_20170728_A.csv"))
abund.2017=abund.2017[!is.na(abund.2017$LATITUDE_WGS84),]
all.2017.tmp=rbind(rand.2017,opt.2017,abund.2017)
s.ind=seq(1,nrow(all.2017.tmp),2)
all.2017=all.2017.tmp[s.ind,]
xyc=cbind(all.2017$LONGITUDE_WGS84,
          all.2017$LATITUDE_WGS84)

###
### Assign projection
###

DDcoor=SpatialPoints(xyc,CRS(
  "+proj=longlat +datum=WGS84")
)

###
### Change projection
###

utmcoor=spTransform(DDcoor,
                    proj      
)

###
### Counts of sea otters
###

count.2017=all.2017$COUNT_ADULT+all.2017$COUNT_PUP

###
### Identify cell id
###

cell[]=1:q
cell.ID=raster::extract(cell,utmcoor@coords)

###
### build matrix with coords, cell id, and count values
###

mat=cbind(utmcoor@coords,cell.ID,count.2017)
mat=mat[order(mat[,3]),]

###
### Create Y matrix
###

uniq.ID=unique(mat[,3])
max.photos=rep(NA,length(uniq.ID))
for(i in uniq.ID){
  data.tmp=subset(mat,cell.ID==i)
  max.photos[i]=dim(data.tmp)[1]
}
max.dim.Y2=max(max.photos,na.rm=TRUE)

Y.tmp=matrix(,length(uniq.ID),max.dim.Y2)
for(i in 1:length(uniq.ID)){
  data.tmp=subset(mat,mat[,3]==uniq.ID[i])
  Y.tmp[i,]=c(data.tmp[,4],rep(NA,dim(Y.tmp)[2]-length(data.tmp[,4])))
}
Y=cbind(uniq.ID,Y.tmp)
Y.2017=matrix(NA,q,max.dim.Y2)
Y.2017[uniq.ID,]=Y.tmp


########################################################################
### Add 2018 aerial photograph data
########################################################################

rand.2018=read.csv(paste0("SO_D_20180720_R.csv"))
rand.2018=rand.2018[!is.na(rand.2018$LATITUDE_WGS84),]
opt.2018=read.csv(paste0("SO_D_20180719_O.csv"))
opt.2018=opt.2018[!is.na(opt.2018$LATITUDE_WGS84),]
abund.2018=read.csv(paste0("SO_D_20180718_A.csv"))
abund.2018=abund.2018[!is.na(abund.2018$LATITUDE_WGS84),]
all.2018.tmp=rbind(rand.2018,opt.2018,abund.2018)
s.ind=seq(1,nrow(all.2018.tmp),2)
all.2018=all.2018.tmp[s.ind,]


xyc=cbind(all.2018$LONGITUDE_WGS84,
          all.2018$LATITUDE_WGS84)

###
### Assign projection
###

DDcoor=SpatialPoints(xyc,CRS(
  "+proj=longlat +datum=WGS84")
)

###
### Change projection
###

utmcoor=spTransform(DDcoor,
                    proj      ### jme changed
)

###
### Counts of sea otters
###

count.2018.tmp=all.2018$COUNT_ADULT+all.2018$COUNT_PUP

###
### Identify cell id
###

cell[]=1:q
cell.ID.tmp=raster::extract(cell,utmcoor@coords)

###
### Remove samples outside study area
###

ind=sort(which(!is.na(cell.ID.tmp)))

###
### build matrix with coords, cell id, and count values
###

cell.ID=cell.ID.tmp[ind]
count.2018=count.2018.tmp[ind]
mat=cbind(utmcoor@coords[ind,],cell.ID,count.2018)
mat=mat[order(mat[,3]),]

###
### Create Y.2018 matrix
###

uniq.ID=unique(mat[,3])
max.photos=rep(NA,length(uniq.ID))
for(i in uniq.ID){
  data.tmp=subset(mat,cell.ID==i)
  max.photos[i]=dim(data.tmp)[1]
}
max.dim.Y2=max(max.photos,na.rm=TRUE)

Y.tmp=matrix(,length(uniq.ID),max.dim.Y2)
for(i in 1:length(uniq.ID)){
  data.tmp=subset(mat,mat[,3]==uniq.ID[i])
  Y.tmp[i,]=c(data.tmp[,4],rep(NA,dim(Y.tmp)[2]-length(data.tmp[,4])))
}
Y=cbind(uniq.ID,Y.tmp)
Y.2018=matrix(NA,q,max.dim.Y2)
Y.2018[uniq.ID,]=Y.tmp

########################################################################
### Add 2019 aerial photograph data
########################################################################

rand.2019=read.csv(paste0("SO_D_20190710_R.csv"))
rand.2019=rand.2019[!is.na(rand.2019$LATITUDE_WGS84),]
opt.2019=read.csv(paste0("SO_D_20190801_O.csv"))
opt.2019=opt.2019[!is.na(opt.2019$LATITUDE_WGS84),]
abund.2019=read.csv(paste0("SO_D_20190712_A.csv"))
abund.2019=abund.2019[!is.na(abund.2019$LATITUDE_WGS84),]
all.2019.tmp=rbind(rand.2019,opt.2019,abund.2019)
s.ind=seq(1,nrow(all.2019.tmp),2)
all.2019=all.2019.tmp[s.ind,]


xyc=cbind(all.2019$LONGITUDE_WGS84,
          all.2019$LATITUDE_WGS84)

###
### Assign projection
###

DDcoor=SpatialPoints(xyc,CRS(
  "+proj=longlat +datum=WGS84")
)

###
### Change projection
###

utmcoor=spTransform(DDcoor,
                    proj      ### jme changed
)

###
### Counts of sea otters
###

count.2019.tmp=all.2019$COUNT_ADULT+all.2019$COUNT_PUP

###
### Identify cell id
###

cell[]=1:q
cell.ID.tmp=raster::extract(cell,utmcoor@coords)

###
### Remove samples outside study area
###

ind=sort(which(!is.na(cell.ID.tmp)))

###
### build matrix with coords, cell id, and count values
###

cell.ID=cell.ID.tmp[ind]
count.2019=count.2019.tmp[ind]
mat=cbind(utmcoor@coords[ind,],cell.ID,count.2019)
mat=mat[order(mat[,3]),]

###
### Create Y.2019 matrix
###

uniq.ID=unique(mat[,3])
max.photos=rep(NA,length(uniq.ID))
for(i in uniq.ID){
  data.tmp=subset(mat,cell.ID==i)
  max.photos[i]=dim(data.tmp)[1]
}
max.dim.Y2=max(max.photos,na.rm=TRUE)

Y.tmp=matrix(,length(uniq.ID),max.dim.Y2)
for(i in 1:length(uniq.ID)){
  data.tmp=subset(mat,mat[,3]==uniq.ID[i])
  Y.tmp[i,]=c(data.tmp[,4],rep(NA,dim(Y.tmp)[2]-length(data.tmp[,4])))
}
Y=cbind(uniq.ID,Y.tmp)
Y.2019=matrix(NA,q,max.dim.Y2)
Y.2019[uniq.ID,]=Y.tmp


########################################################################
### Transects
########################################################################
OrigData=read.csv(paste0("SEAK_otter.csv"),
                  header=TRUE)

transect.sub=subset(OrigData,year=="1999")
#l=vector("list",nrow(transect.sub))
l.l=list()
begin.coord.tmp=data.frame(lon=transect.sub$west_long,
                           lat=transect.sub$west_lat)
end.coord.tmp=data.frame(lon=transect.sub$east_long,
                         lat=transect.sub$east_lat)
begin.coord=begin.coord.tmp[!is.na(begin.coord.tmp$lon) | !is.na(begin.coord.tmp$lon) ,]
end.coord=end.coord.tmp[!is.na(begin.coord.tmp$lon) | !is.na(begin.coord.tmp$lon) ,]
begin.coord=begin.coord.tmp[!is.na(end.coord.tmp$lon) | !is.na(end.coord.tmp$lon) ,]
end.coord=end.coord.tmp[!is.na(end.coord.tmp$lon) | !is.na(end.coord.tmp$lon) ,]
l=vector("list",nrow(end.coord))
for(i in seq_along(l)){
  l[[i]]=Lines(list(Line(
    rbind(begin.coord[i,],
          end.coord[i,]))),as.character(i))
}
l.l[[1]]=l
transect.r=raster(,nrows=y,ncols=x,xmn=xmin,
                  xmx=xmax,ymn=ymin,ymx=ymax,crs=NA)
transect.sub=subset(data,year=="1999")
transect.sub.r=rasterize(x=SpatialLines(l),y=transect.r,
                         field=1,fun="last",background=NA)
transectAll=transect.sub.r
transectAll[]=transect.sub.r[]
#transectAll[Counts[[1]][]>0]=1
transectAll=stack(mget(rep("transectAll",20)))
transectAll[[1]][]=rep(NA,q)
ind=2
for(t in years){
  transect.sub=subset(OrigData,year==t)
  #if(dim(transect.sub)[1]>0){
  if(!is.na(transect.sub$transect[1])){
    #l=vector("list",nrow(transect.sub))
    begin.coord.tmp=data.frame(lon=transect.sub$west_long,
                               lat=transect.sub$west_lat)
    end.coord.tmp=data.frame(lon=transect.sub$east_long,
                             lat=transect.sub$east_lat)
    begin.coord=begin.coord.tmp[!is.na(begin.coord.tmp$lon) | !is.na(begin.coord.tmp$lon) ,]
    end.coord=end.coord.tmp[!is.na(begin.coord.tmp$lon) | !is.na(begin.coord.tmp$lon) ,]
    begin.coord=begin.coord.tmp[!is.na(end.coord.tmp$lon) | !is.na(end.coord.tmp$lon) ,]
    end.coord=end.coord.tmp[!is.na(end.coord.tmp$lon) | !is.na(end.coord.tmp$lon) ,]
    l=vector("list",nrow(end.coord))
    for(i in seq_along(l)){
      l[[i]]=Lines(list(Line(rbind(
        begin.coord[i,],
        end.coord[i,]))),
        as.character(i))
    }
    l.l[[ind]]=l
    transect.sub.r=rasterize(x=SpatialLines(l),
                             y=transect.r,field=1,fun="last",
                             background=NA)
    transectAll[[ind]][]=transect.sub.r[]
    #transectAll[[ind]][Counts[[ind]][]>0]=1
  }else{
    transectAll[[ind]][]=rep(NA,q)
  }
  ind=ind+1
}

########################################################################
### Combine transect data and count data
########################################################################
#### this creates absence data 

Y.r[]=Counts[[7]][]*transectAll[[7]][]
Y.r=stack(mget(rep("Y.r",20)))
for (i in 1:20) {
  Y.r[[i]][]=ifelse(Counts[[i]][]>0,Counts[[i]][],
                    transectAll[[i]][]-1)
}
## plot(Y.r[[20]])

######################################################
### Load ISU data
######################################################
# 
all.ISU.data=read.csv(paste("SEAK_otter_ISU.csv"))
ISU.N=all.ISU.data$circle.adt+
  all.ISU.data$circle.pup
ISU.Y=all.ISU.data$strip.adt+
  all.ISU.data$strip.pup
ISU.N[ISU.Y>ISU.N]=ISU.Y[ISU.Y>ISU.N]
sum(ISU.N<ISU.Y)
ISU.year=all.ISU.data$year
ISU.data=data.frame(ISU.N,ISU.Y,ISU.year)
Y.1995=ISU.data$ISU.Y[ISU.data$ISU.year==1995]
N.1995=ISU.data$ISU.N[ISU.data$ISU.year==1995]
Y.1999=ISU.data$ISU.Y[ISU.data$ISU.year==1999]
N.1999=ISU.data$ISU.N[ISU.data$ISU.year==1999]
Y.2000=ISU.data$ISU.Y[ISU.data$ISU.year==2000]
N.2000=ISU.data$ISU.N[ISU.data$ISU.year==2000]
Y.2001=ISU.data$ISU.Y[ISU.data$ISU.year==2001]
N.2001=ISU.data$ISU.N[ISU.data$ISU.year==2001]
Y.2002=ISU.data$ISU.Y[ISU.data$ISU.year==2002]
N.2002=ISU.data$ISU.N[ISU.data$ISU.year==2002]
Y.2003=ISU.data$ISU.Y[ISU.data$ISU.year==2003]
N.2003=ISU.data$ISU.N[ISU.data$ISU.year==2003]
Y.2004=ISU.data$ISU.Y[ISU.data$ISU.year==2004]
N.2004=ISU.data$ISU.N[ISU.data$ISU.year==2004]
Y.2005=ISU.data$ISU.Y[ISU.data$ISU.year==2005]
N.2005=ISU.data$ISU.N[ISU.data$ISU.year==2005]
Y.2006=ISU.data$ISU.Y[ISU.data$ISU.year==2006]
N.2006=ISU.data$ISU.N[ISU.data$ISU.year==2006]
Y.2010=ISU.data$ISU.Y[ISU.data$ISU.year==2010]
N.2010=ISU.data$ISU.N[ISU.data$ISU.year==2010]
Y.2012=ISU.data$ISU.Y[ISU.data$ISU.year==2012]
N.2012=ISU.data$ISU.N[ISU.data$ISU.year==2012]
ISU.n=list(Y.1995=Y.1995,
         N.1995=N.1995,
         Y.1999=Y.1999,
         N.1999=N.1999,
         Y.2000=Y.2000,
         N.2000=N.2000,
         Y.2001=Y.2001,
         N.2001=N.2001,
         Y.2002=Y.2002,
         N.2002=N.2002,
         Y.2003=Y.2003,
         N.2003=N.2003,
         Y.2004=Y.2004,
         N.2004=N.2004,
         Y.2005=Y.2005,
         N.2005=N.2005,
         Y.2006=Y.2006,
         N.2006=N.2006,
         Y.2010=Y.2010,
         N.2010=N.2010,
         Y.2012=Y.2012,
         N.2012=N.2012)
#all.ISU.data=read.csv("All ISUs 1999_2012_02032016.csv")
## get locs of isu's
isux=0
isuy=0
#all.ISU.data=all.ISU.data[-209,]
for(i in 1:nrow(all.ISU.data)){
  tmpx=data$GROUP_X[data$ISU==all.ISU.data$isunum[i] & data$year==all.ISU.data$year[i]]
  tmpy=data$GROUP_Y[data$ISU==all.ISU.data$isunum[i] & data$year==all.ISU.data$year[i]]
  if(length(tmpx[!is.na(tmpx)])>0){
    isux[i]=tmpx[!is.na(tmpx)]
    isuy[i]=tmpy[!is.na(tmpy)]
  }else{
    isux[i]=NA
    isuy[i]=NA
  }
}
all.ISU.data$x=isux
all.ISU.data$y=isuy

ISUind.r=raster(,nrows=y,ncols=x,xmn=xmin,
                xmx=xmax,ymn=ymin,ymx=ymax,crs=NA)

ISU.df=all.ISU.data[!is.na(all.ISU.data$x),]
ISU.df$N=ISU.df$circle.adt+ISU.df$circle.pup
ISU.df$C=ISU.df$strip.adt+ISU.df$strip.pup

## remove data where C larger than N
ISU.df=ISU.df[-which(ISU.df$N<ISU.df$C),]

ISU.sub=subset(ISU.df,year==1999)
ISU.sub.r=rasterize(cbind(ISU.sub$x,ISU.sub$y),
                    ISUind.r,field=1,fun="last",background=0)
isuInd=ISU.sub.r
isuInd[]=ISU.sub.r[]
isuInd=stack(mget(rep("isuInd",20)))
ind=1
years=1993:2012
for(t in years){
  ISU.sub=subset(ISU.df,year==t)
  if(dim(ISU.sub)[1]>0){
    ISU.sub.r=rasterize(
      x=cbind(
        ISU.sub$x,
        ISU.sub$y),
      y=ISUind.r,
      field=1,
      fun="last",
      background=0)
    isuInd[[ind]][]=ISU.sub.r[]
  }else{
    isuInd[[ind]][]=rep(0,q)
  }
  ind=ind+1
}

ISU.ind=vec(isuInd[])
indic=1:length(ISU.ind)
ISU=sort(indic[ISU.ind==1])

## N data (known Ns from ISU)
N.r=raster(,nrows=y,ncols=x,xmn=xmin,
           xmx=xmax,ymn=ymin,ymx=ymax,crs=NA)



## Remove aggregation affect of Cs from ISU units
N.sub=subset(ISU.df,year==1999)
N.sub.r=rasterize(cbind(N.sub$x,N.sub$y),N.r,field=N.sub$N,
                  fun="sum",background=0)
CISU.sub.r=rasterize(cbind(N.sub$x,N.sub$y),N.r,field=N.sub$C,
                     fun="sum",background=0)

##
N=N.sub.r
CISU=CISU.sub.r
N=stack(mget(rep("N",20)))
CISU=stack(mget(rep("CISU",20)))
ind=1
for(t in years){
  N.sub=subset(ISU.df,year==t)
  if(dim(N.sub)[1]>0){
    N.sub.r=rasterize(x=cbind(N.sub$x,N.sub$y),
                      y=N.r,field=N.sub$N,fun="sum",
                      background=0)
    CISU.sub.r=rasterize(cbind(N.sub$x,N.sub$y),N.r,field=N.sub$C,
                         fun="sum",background=0)
    N[[ind]][]=N.sub.r[]
    CISU[[ind]][]=CISU.sub.r[]
  }else{
    N[[ind]][]=rep(0,q)
    CISU[[ind]][]=rep(0,q)
  }
  ind=ind+1
}


###
### Data for analysis
###

Y=c(Y.r[[1]][],
    Y.r[[2]][],
    Y.r[[3]][],
    Y.r[[4]][],
    Y.r[[5]][],
    Y.r[[6]][],
    Y.r[[7]][],
    Y.r[[8]][],
    Y.r[[9]][],
    Y.r[[10]][],
    Y.r[[11]][],
    Y.r[[12]][],
    Y.r[[13]][],
    Y.r[[14]][],
    Y.r[[15]][],
    Y.r[[16]][],
    Y.r[[17]][],
    Y.r[[18]][],
    Y.r[[19]][],
    Y.r[[20]][]
)


## Include distribution surveys in N
TransectVec.tmp=vec(transectAll[])
TransectVec=ifelse(Y>=0,1,NA)
N.vec=vec(N[])*TransectVec
CISU=vec(CISU[])
CISU=CISU[ISU]
Y[ISU]=CISU
NISU=N.vec[ISU]

###
### Create `Cell' raster
###

cell[]=1:q
## plot(cell)

###
### Create Bounday raster
###

Boundary[is.na(bath)]=0
Boundary[bath<0]=1
Boundary[bath>=0]=0
# fill in missing chunk of ocean with 1's
fill_poly = coords2Polygons(matrix(c(360000,100000,xmax(Boundary)+1000,250000,xmax(Boundary)+1000,100000),
                                   byrow=2,ncol=2),ID='p')
Boundary[fill_poly] <- 1

###
### Exposed shore
###

ExpCov=Boundary
isl.poly=coords2Polygons(matrix(c(360000,180000,390000,180000,390000,140000,360000,140000),ncol=2,byrow=T),
                         ID='p')
ExpCov[isl.poly]=1  # don't let the little island block the swell
exp.mat=as.matrix(ExpCov)
exp.mat1=matrix(0L,nrow = nrow(exp.mat),ncol=ncol(exp.mat))
for(i in 1:dim(exp.mat)[2]){
  tmp=max(which(exp.mat[,i]==0))
  exp.mat1[tmp:dim(exp.mat)[1],i]=1
}
values(ExpCov)=0
ExpCov=raster(exp.mat1,template=ExpCov)
ExpCov=ExpCov*Boundary

###
### Distance to shore covariate
###

DistCov.tmp1=Boundary
isl.poly=coords2Polygons(matrix(c(360000,180000,390000,180000,390000,140000,360000,140000),ncol=2,byrow=T),
                         ID='p')
DistCov.tmp1[isl.poly]=1  # don't let the little island block the swell
DistCov.tmp1[Boundary==1]=NA
DistCov.tmp1[DistCov.tmp1==1]=NA
DistCov.tmp2=distance(DistCov.tmp1)
DistCov.tmp2[is.na(DistCov.tmp2)]=0

###
### make boundary at 5 km offshore or up to 20 km if depth<100
###
d.r=bath*ExpCov
d.r[ExpCov==0]=0
d.r[bath>-100&bath<0]=1
d.r[d.r!=1]=0
d.r=d.r*ExpCov
b.tmp=Boundary
b.tmp[DistCov.tmp2>=20000 & ExpCov==1]=0
Boundary=Boundary*b.tmp

Boundary[d.r==0 & DistCov.tmp2>=5000 & ExpCov==1]=0
Boundary[Boundary==1 & b.tmp==0 & ExpCov==1]=0
#Boundary[ExpCov==0 & !is.na(bath) & bath<0]=1

ExpCov=ExpCov*Boundary
DistCov=scale(DistCov.tmp2*Boundary,center=TRUE,scale=TRUE)
DistCov=DistCov*Boundary
## plot(Boundary)

BoundaryNA=Boundary
BoundaryNA[Boundary==0]=NA
## plot(BoundaryNA)


###
### Larger Domain
###

BoundaryInf=BoundaryNA
#BoundaryInf[145:180,0:80]=NA
#BoundaryInf[160:180,120:150]=NA
BoundaryInf.v=rep(BoundaryInf[],length(keep))

###
### Up-scaled boundary layer surrounded by zeros
###

Boundary.us=aggregate(Boundary,fact=us.fact,na.rm=TRUE,fun=max)
Boundary.us[Boundary.us[]>1]=1
#Boundary.us[1,]=0  ## no dispersal out of top
#Boundary.us[2,5:6]=0 ## no dispersal between two arms
#Boundary.us[2:5,1]=0 ## no dispersal out of west arm
#Boundary.us[18,1:2]=1 ## yes dispersal out of bottom left
## plot(Boundary.us)

###
### Depth covariate
###

depth=40
DepthCov=bath
DepthCov[bath[]< -depth]=0
DepthCov[bath[]>= -depth]=1
DepthCov[is.na(bath[])]=0
isl.poly=coords2Polygons(matrix(c(360000,180000,390000,180000,390000,140000,360000,140000),ncol=2,byrow=T),
                         ID='p')
DepthCov[isl.poly]=0
DepthCov=DepthCov*Boundary



###
### Bottom slope
###
SlopeCov.tmp=terrain(bath,opt='slope',unit='degrees',neighbors=8)  
SlopeCov.tmp[is.na(SlopeCov.tmp)]=0
SlopeCov=scale(SlopeCov.tmp*Boundary, center=TRUE, scale=TRUE)
SlopeCov=SlopeCov*Boundary*DepthCov

###
### Shoreline complexity
###

ShoreCov.tmp=boundaries(DistCov.tmp1,type='inner')
ShoreCov.tmp[is.na(ShoreCov.tmp)]=0
ShoreCov=focal(ShoreCov.tmp,w=matrix(1,nr=11,nc=11))*Boundary
ShoreCov[is.na(ShoreCov)]=0
ShoreCov=scale(ShoreCov*Boundary,center=TRUE, scale=TRUE)
ShoreCov=ShoreCov*Boundary

###
### Dist to town
###

towns=st_read('mv_town_pt.shp')
towns=as_Spatial(towns)
towns=as(towns,'SpatialPoints')
towns=raster::intersect(towns,projectRaster(Boundary,crs=crs(towns)))
towns=spTransform(towns,crs(Boundary))
Boundary.f=Boundary
Boundary.f[]=as.factor(Boundary[])
trans=transition(Boundary.f,mean,4,symm=F)
r.coords=xyFromCell(Boundary,1:q)
distTown=list() # note towns 5 and 31 outside of study area
for(i in 1:length(towns)){
  if(i!=5&i!=31) distTown[[i]]=costDistance(trans,towns[i],r.coords)
}
distTown=distTown[-5]
distTown=distTown[-30]
TownCov=Boundary
TownCov[]=Reduce('+',distTown)
TownCov[is.na(TownCov)]=0
TownCov=scale(TownCov*Boundary,center=TRUE, scale = TRUE)
TownCov=TownCov*Boundary

###
### GLBA
###
# polygon that encompasses GLBA
g.p = Polygon(matrix(c(-150000,240000,-40000,270000,-31000,323000,-150000,323000),ncol=2,byrow=T))
glba.ind = cellFromPolygon(bath,SpatialPolygons(list(Polygons(list(g.p),ID='a')),
                                                           proj4string = crs(bath)))
glba.ind = unlist(glba.ind)
glba.ind = glba.ind[!is.na(bath[][glba.ind])]
glbaCov=Boundary
glbaCov[]=0
glbaCov[glba.ind]=1
glbaCov=glbaCov*Boundary

###
### Fisheries closures
###
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
fish.ind1=cellFromPolygon(bath,cuc.clos.t)
fish.ind2=cellFromPolygon(bath,urch.clos.t)
fish.ind3=cellFromPolygon(bath,geo.clos.t)
fish.ind1=unlist(fish.ind1)
fish.ind2=unlist(fish.ind2)
fish.ind3=unlist(fish.ind3)
fish.ind=c(fish.ind1,fish.ind2,fish.ind3)
fish.ind=unique(fish.ind)
fish.ind = fish.ind[!is.na(bath[][fish.ind])]
fishCov=Boundary
fishCov[]=0
fishCov[fish.ind]=1
fishCov=fishCov*Boundary

###
### Create co-variate matrices
###

X=cbind(1,DepthCov[],DistCov[],SlopeCov[]*DepthCov[],ShoreCov[],TownCov[],glbaCov[],fishCov[])
W=matrix(1,nr=q,nc=1)

# 
# ###
# ### Data for analysis
# ###
# 
# Y=c(Y.r[[1]][],
#     Y.r[[2]][],
#     Y.r[[3]][],
#     Y.r[[4]][],
#     Y.r[[5]][],
#     Y.r[[6]][],
#     Y.r[[7]][],
#     Y.r[[8]][],
#     Y.r[[9]][],
#     Y.r[[10]][],
#     Y.r[[11]][],
#     Y.r[[12]][],
#     Y.r[[13]][],
#     Y.r[[14]][],
#     Y.r[[15]][],
#     Y.r[[16]][],
#     Y.r[[17]][],
#     Y.r[[18]][],
#     Y.r[[19]][],
#     Y.r[[20]][]
# )

########################################################################
########################################################################
###
### Fit the model to the real data using MCMC
###
########################################################################
########################################################################

No.chains=1
seed=c(13,23,33,43,53,63,73,83,93,103,113,123)#[1:No.chains]
i=1  # use this seed

###
### Bundle data
###

data=list(Y=Y,
          N.vec=N.vec,
          CISU=CISU,
          Y.aerial=list(Y.2017=Y.2017,Y.2018=Y.2018,Y.2019=Y.2019),
          ISU=ISU.n,
          ISU.i=ISU,
          NISU=NISU,
          X=X,
          seed=seed[i])

###
### Spatio-temporal settings
###

extent=c(xmin,xmax,ymin,ymax)
st.info=list(dt=dt,
             time.frame=time.frame,
             us.fact=us.fact,
             smooth.fact=smooth.fact,
             res=res,
             d=d,
             extent=extent
)

###
### MCMC Settings
###

n.iter=1000
checkpoint=100
savepoint=100

###
### Priors
###

## beta Normal(mu.beta,var.beta)
mu.beta=0
var.beta=10^2

## gamma (normal)
q.gamma=.25
r.gamma=0.01^2



## theta  (changed by jme)
mu.theta=c(100,10,10,100,10,100,100)
var.theta=c(20^2,1^2,1^2,20^2,1^2,20^2,20^2)

# kappa
mu.kappa=c(10,2,10,10,2,10,10)
var.kappa=c(3^2,1^2,3^2,3^2,1^2,3^2,3^2)

## p (beta(1,1))
q.p=1
r.p=1

## Aerial photograph priors (moment matching from MEE paper -
##  use Appendix S2 from ms)
p.alpha=44.04937
p.beta=13.40566

## tau (unif(0,1))
q.tau=0
r.tau=1

## alpha (normal)
mu.alpha=0
var.alpha=10^2
# 
# ## K
# K.alpha=1
# K.beta=.25^2
# 
# ## K.g
# K.g.alpha=6
# K.g.beta=1^2
# 
# ## K.f
# K.f.alpha=1
# K.f.beta=.25^2

priors=list(
  beta.prior=c(mu.beta,var.beta),
  gamma.prior=c(q.gamma,r.gamma),
  theta.prior=c(mu.theta,var.theta),
  kappa.prior=c(mu.kappa,var.kappa),
  p.prior=c(q.p,r.p),
  p.2017=c(p.alpha,p.beta),
  tau.prior=c(q.tau,r.tau),
  alpha.prior=c(mu.alpha,var.alpha)
  # K.prior=c(K.alpha,K.beta),
  # K.g.prior=c(K.g.alpha,K.g.beta),
  # K.f.prior=c(K.f.alpha,K.f.beta)
)

###
### Starting values
###

gamma=0.295

beta=c(16.35, 
       -1.7,
       .3,
       .23,
       .18,
       .4,
       -.23,
       -1.37)

kappa=c(28,
        2.5,
        10,
        .65,
        2,
        8,
        10)

theta=c(140,
        10,
        10,
        100,
        10,
        100,
        110)
#theta=theta*10

alpha=c(-1.7,
        3.2,
        1.3)


#theta=7800
#kappa=10000
#theta = c()
p.1995=0.5067195
p.1999=0.7972976
p.2000=0.7614920
p.2001=0.8527629
p.2002=0.8995615
p.2003=0.8035992
p.2004=0.7751074
p.2005=0.5427057
p.2006=0.7375375
p.2010=0.9192737
p.2012=0.5799473
p.2017=0.7933035
p.2018=0.8121646
p.2019=0.8121646
p.ind=rep(NA,length(time.frame))
p.ind[26]=p.1995
p.ind[30]=p.1999
p.ind[31]=p.2000
p.ind[32]=p.2001
p.ind[33]=p.2002
p.ind[34]=p.2003
p.ind[35]=p.2004
p.ind[36]=p.2005
p.ind[37]=p.2006
p.ind[41]=p.2010
p.ind[43]=p.2012
p.ind[48]=p.2017
p.ind[49]=p.2018
p.ind[50]=p.2019


tau=0.026


inits=list(gamma=gamma,
           beta=beta,
           theta=theta,
           kappa=kappa,
           tau=tau,
           p.ind=p.ind,
           alpha=alpha
           # K=K,
           # K.g=K.g,
           # K.f=K.f
)

###
### Parameters to monitor
###

st.info$years.keep=c(1,6,11,16,21,26,31,36,41,46,51)

parameters=c("gamma",
             "beta",
             "kappa",
             "theta",
             "tau",
             "p",
             "n.tot",
             #n.tot.v",
             "lambda",
             #"lambda11",
             #"lambda.all",
             #"tuners",
             "alpha",
             "delta.bar",
             'pp',
            'N.t',
             'N')

###
### Output location
###

out.loc.l=list()
# for(i in 1:No.chains){
#   out.loc.l[[i]]=path.expand(paste0(getwd(),
#                                     "/MCMC.Output.seed.",
#                                     seed[i],".RData"))
# }

out.loc.l[[1]]=path.expand(paste0(getwd(),
                                  "/MCMC.Output.seed.",
                                  seed[i],".RData"))

Background=list(bath=bath,
                Boundary=Boundary,
                BoundaryNA=BoundaryNA,
                BoundaryInf=BoundaryInf,
                BoundaryInf.v=BoundaryInf.v,
                Boundary.us=Boundary.us)

rm(list=ls()[!ls() %in% c("data",
                          "priors",
                          "inits",
                          "parameters",
                          "st.info",
                          "Background",
                          "n.iter",
                          "checkpoint",
                          "savepoint",
                          "out.loc.l",
                          "run.chain.2pl.list")])
gc()

##
## save for stringing together chains
##

#save.image(file=paste0('otter_work_',data$seed,'.RData'))

########################################################################
### Run MCMC
########################################################################

source("MCMC_Logistic_reint_tuned.R")


##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sys.time=Sys.time()
run.chain.2pl.list(data=data,
                   priors=priors,
                   inits=inits,
                   parameters=parameters,
                   st.info=st.info,
                   Background=Background,
                   n.iter=n.iter,
                   checkpoint=checkpoint,
                   savepoint=savepoint,
                   out.loc=out.loc.l
)
Sys.time()-sys.time

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

