# Simulation function for MR-DOC2 paper

# Packages you will need


dep <- function(x) {
  ifelse(!all(x %in% installed.packages()),
    install.packages(x[x %in% installed.packages()[, "Package"] == F],
      repos = "http://cran.us.r-project.org"),
    lapply(x, require, character.only = T)
  )
}


dep(c("foreach", "ggplot2", "doParallel", "MASS", "umx"))


sim  <- function(
               b1 = sqrt(.05),
               b3 = sqrt(.05),
               g1 = sqrt(.05),
               g2 = sqrt(.05),
               abs = .25,
                ass = .25,
               cbs = .15,
               css = .15,
               ra = .3,
               rc = .25,
               re = .3,
               rf = .25,
               b2 = 0,
               b4 = 0,
               Nmz = 1000,
               Ndz = 1000,
               noC = FALSE,
               do_id = FALSE,
               Ve_pinned = FALSE,
               hp = FALSE,
               do_pow = "connor",
               group = NULL,
               est_out = FALSE){

ny=8
ne=20
alpha=0.05
#
# parameter
#
npar=18 #
partrue=rep(0,npar)
parfree=rep(F,npar)
parstart=rep(0,npar)
parnames=
c('x','y',
'ab','cb','eb','as','cs','es',
'ra','rc','re','rf',
'b1','b2','b3','b4','g1','g2')
parfree=rep(F,npar)
parstart=rep(0,npar)
#
partrue[1]=x=1 # allways 1
partrue[2]=y=1 # allways 1

ifelse(hp==T, partrue[14] <- b2, partrue[14] <- b2 <- 0)
ifelse(hp==T, partrue[16] <- b4, partrue[16] <- b4 <- 0)


cs=sqrt(css)
cb=sqrt(cbs)
as=sqrt(ass)
ab=sqrt(abs)

if (Ve_pinned) {
  es = eb = sqrt(0.1)
  } else {
  es=sqrt(1-ass-css)
  eb=sqrt(1-abs-cbs)
}

bgpars=c(b1,b3,g1,g2)^2
if (g1<0) bgpars[3]=-(g1^2)
if (g2<0) bgpars[4]=-(g2^2)


  keep  =  data.frame(b1 = bgpars[1], b3 = bgpars[2], g1 = bgpars[3], g2 = bgpars[4])


keep$ra = ra
keep$re = re
keep$rf = rf
keep$as = as
keep$es = es
keep$ab = ab
keep$eb = eb
keep$cs = cs
keep$cb = cb
keep$rc = rc


partrue[13]=b1  #=sqrt(.1)
partrue[15]=b3  #=sqrt(.1) #
partrue[17]=g1  #=.15 # -.25 # .20
partrue[18]=g2  #=-.102 # .25  #-.052
#
names(partrue)=
c('x','y','ab','cb','eb','as','cs','es',
'ra','rc','re','rf',
'b1','b2','b3','b4','g1','g2')

#
vnames_ne=c('B','S','PB','PS','B','S','PB','PS','AB','CB','EB','AS','CS','ES','AB','CB','EB','AS','CS','ES')
Ly=matrix(0,ny,ne)
diag(Ly[1:ny,1:ny])=c(1,1,x,y,1,1,x,y)
Be=matrix(0,ne,ne)
Iden=diag(ne)
Psmz= Psdz=matrix(0,ne,ne)
#
#                              B      S      B      S
#            ph  prs  ph  prs  a c e  a c e  a c e  a c e
diag(Psmz)=c(0,0,1,1, 0,0,1,1, 1,1,1, 1,1,1, 1,1,1, 1,1,1)
#
Psmz[12,9]=ra;
Psmz[15,12]=ra;
Psmz[18,9]=ra;
Psmz[18,15]=ra;
Psmz[13,10]=rc;
Psmz[19,16]=rc;
Psmz[16,13]=rc;
Psmz[19,10]=rc;
Psmz[20,17]=re;
Psmz[14,11]=re;
Psmz[15,9]=1;
Psmz[16,10]=1;
Psmz[17,11]=0;
Psmz[18,12]=1;
Psmz[19,13]=1;
Psmz[20,14]=0;
Psmz[7,3]=1;
Psmz[8,4]=1;
Psmz[4,3]=rf;
Psmz[8,7]=rf;
Psmz[8,3]=rf;
Psmz[7,4]=rf;
for (i in 1:ne) { for (j in 1:i) {
Psmz[j,i]=Psmz[i,j]}}
#
Psdz=Psmz
Psdz[15,12]=.5*ra;
Psdz[18,9]=.5*ra;
Psdz[15,9]=.5;
Psdz[16,10]=1;
Psdz[17,11]=0;
Psdz[18,12]=.5;
Psdz[19,13]=1;
Psdz[20,14]=0;
Psdz[7,3]=.5;
Psdz[8,4]=.5;
Psdz[4,3]=rf;
Psdz[8,7]=rf;
Psdz[8,3]=.5*rf;
Psdz[7,4]=.5*rf;
#
for (i in 1:ne) { for (j in 1:i) { Psdz[j,i]=Psdz[i,j]}}
#

Be[1,3]=b1;
Be[2,1]=g1;
Be[2,3]=b2;
Be[1,2]=g2;
Be[2,4]=b3;
Be[1,4]=b4;  #  added PRSS->B
#

Be[5,7]=b1;
Be[6,5]=g1;
Be[6,7]=b2;
Be[5,6]=g2;
Be[6,8]=b3;
Be[5,8]=b4; #
#
Be[1,9]=ab;
Be[1,10]=cb;
Be[1,11]=eb;
Be[2,12]=as;
Be[2,13]=cs;
Be[2,14]=es;
Be[5,15]=ab;
Be[5,16]=cb;
Be[5,17]=eb;
Be[6,18]=as;
Be[6,19]=cs;
Be[6,20]=es;

iBe=solve(Iden-Be);
Smz=Ly%*%iBe%*%Psmz%*%t(iBe)%*%t(Ly)
Sdz=Ly%*%iBe%*%Psdz%*%t(iBe)%*%t(Ly)
#
SiMZ=SiDZ=matrix(0,8,8)

#
# the exact sim data are convenient for the regression analyses
# evenif they are not used in modeling fitten (where we use Smz and Sdz)
#
datmz=mvrnorm(Nmz,rep(0,8),Sigma= Smz*(Nmz/(Nmz-1)),emp=T)  # or Sigma=Smz
datmz=as.data.frame(datmz)
datdz=mvrnorm(Ndz,rep(0,8),Sigma=Sdz*(Ndz/(Ndz-1)),emp=T) # or Sigma=Sdz
datdz=as.data.frame(datdz)
colnames(datmz)=colnames(datdz)=
 c('BMI1','SBP1','GVBMI1','GVSBP1', 'BMI2','SBP2','GVBMI2','GVSBP2')

# how the parameter values pan out interms of r^2
keep$PB_B=summary(lm(BMI1~GVBMI1, data=datdz))$r.squared
keep$PS_S=summary(lm(SBP1~GVSBP1, data=datdz))$r.squared
keep$PS_PB_S=summary(lm(SBP1~GVSBP1+GVBMI1, data=datdz))$r.squared
keep$PS_PB_B=summary(lm(BMI1~GVSBP1+GVBMI1, data=datdz))$r.squared
keep$S_B=summary(lm(BMI1~SBP1, data=datdz))$r.squared
keep$B_S=summary(lm(SBP1~BMI1, data=datdz))$r.squared
keep$all_B=summary(lm(BMI1~SBP1+GVBMI1+GVSBP1, data=datdz))$r.squared
keep$all_S=summary(lm(SBP1~BMI1+GVBMI1+GVSBP1, data=datdz))$r.squared


# prepare openMx input
#
mxfreeBE=matrix(F,ne,ne)
mxfreeBE[1,3]=T	#'b1';
mxfreeBE[2,1]=T	#'g1';

ifelse(hp, mxfreeBE[2,3]<-T, mxfreeBE[2,3]<-F)
mxfreeBE[1,2]=T    # 'g2'; zero!
mxfreeBE[2,4]=T    # 'b3';

ifelse(hp, mxfreeBE[1,4]<-T, mxfreeBE[1,4]<-F)
mxfreeBE[5,7]=T	#'b1';
mxfreeBE[6,5]=T	#'g1';

ifelse(hp, mxfreeBE[6,7]<-T, mxfreeBE[6,7]<-F)
mxfreeBE[5,6]=T    # g2;
mxfreeBE[6,8]=T    # g3;

ifelse(hp, mxfreeBE[5,8]<-T, mxfreeBE[5,8]<-F)

mxfreeBE[1,9]=T	#'ab';

ifelse(noC == T,
  mxfreeBE[1, 10] <- F, # T	#'cb';
  mxfreeBE[1, 10] <- T # T	#'cb';
)

mxfreeBE[1, 11] <- T #' eb';
mxfreeBE[2, 12] <- T #' as';


ifelse(noC == T,
  mxfreeBE[2, 13] <- F, # T	#'cs';
  mxfreeBE[2, 13] <- T # T	#'cs';
)

mxfreeBE[2, 14] <- T #' es';
mxfreeBE[5, 15] <- T #' ab';

ifelse(noC == T,
  mxfreeBE[5, 16] <- F, # T	#'cb';
  mxfreeBE[5, 16] <- T # T	#'cb';
)

mxfreeBE[5, 17] <- T #' eb';
mxfreeBE[6, 18] <- T #' as';

ifelse(noC == T,
  mxfreeBE[6, 19] <- F, # T	#'cs';
  mxfreeBE[6, 19] <- T # T	#'cs';
)


mxfreeBE[6,20]=T	#'es';

#
mxlabelsBE=matrix(NA,ne,ne)
mxlabelsBE[1,3]='b1';
mxlabelsBE[2,1]='g1';
mxlabelsBE[2,3]='b2';
mxlabelsBE[1,2]='g2';
mxlabelsBE[2,4]='b3';
mxlabelsBE[1,4]='b4';
#
mxlabelsBE[5,7]='b1';
mxlabelsBE[6,5]='g1';
mxlabelsBE[6,7]='b2';
mxlabelsBE[5,6]='g2';
mxlabelsBE[6,8]='b3';
mxlabelsBE[5,8]='b4';
#
mxlabelsBE[1,9]='ab';
mxlabelsBE[1,10]='cb';
mxlabelsBE[1,11]='eb';
mxlabelsBE[2,12]='as';
mxlabelsBE[2,13]='cs';
mxlabelsBE[2,14]='es';
mxlabelsBE[5,15]='ab';
mxlabelsBE[5,16]='cb';
mxlabelsBE[5,17]='eb';
mxlabelsBE[6,18]='as';
mxlabelsBE[6,19]='cs';
mxlabelsBE[6,20]='es';
#
mxvaluesBE=matrix(0,ne,ne)
mxvaluesBE[1,3]=b1  #'b1';
mxvaluesBE[2,1]=g1  #'g1';
mxvaluesBE[2,3]=b2  #'b2';
mxvaluesBE[1,2]=g2  #'g2';
mxvaluesBE[2,4]=b3 #  'b3';
mxvaluesBE[1,4]=b4 #  'b4';

#
mxvaluesBE[5,7]=b1  #'b1';
mxvaluesBE[6,5]=g1  #'g1';
mxvaluesBE[6,7]=b2  #'b2';
mxvaluesBE[5,6]=g2  #'g2';
mxvaluesBE[6,8]=b3  #'b3';
mxvaluesBE[5,8]=b4  #'b4';
#
#
mxvaluesBE[1,9]=ab  #'ab';
mxvaluesBE[1,10]=cb  #'cb';
mxvaluesBE[1,11]=eb  #'eb';
mxvaluesBE[2,12]=as  #'as';
mxvaluesBE[2,13]=cs  #'cs';
mxvaluesBE[2,14]=es; #'es';
mxvaluesBE[5,15]=ab; #'ab';
mxvaluesBE[5,16]=cb; #'cb';
mxvaluesBE[5,17]=eb; #'eb';
mxvaluesBE[6,18]=as; #'as';
mxvaluesBE[6,19]=cs; #'cs';
mxvaluesBE[6,20]=es; #'es';
#
mxfreePS=matrix(F,ne,ne)
mxfreePS[12,9]=T;	# 'ra';
mxfreePS[15,12]=T; # 'ra';
mxfreePS[18,9]=T; # 'ra';
mxfreePS[18,15]=T; # 'ra';

ifelse(noC == T,
mxfreePS[13,10] <-mxfreePS[19,16] <- mxfreePS[16,13] <- mxfreePS[19,10] <- F,
mxfreePS[13,10] <-mxfreePS[19,16] <- mxfreePS[16,13] <- mxfreePS[19,10] <- T
)


mxfreePS[20,17]=T; # re;
mxfreePS[14,11]=T; # re;
mxfreePS[4,3]=T;   # rf
mxfreePS[8,7]=T;   # rf
mxfreePS[7,4]=T;   # rf
mxfreePS[8,3]=T;   # rf
#
for (i in 1:ne) { for (j in 1:i) {
mxfreePS[j,i]=mxfreePS[i,j]}}
#
mxvaluesPS=matrix(0,ne,ne)
mxvaluesPS[15,9]=1;
mxvaluesPS[16,10]=1;
mxvaluesPS[17,11]=0;
mxvaluesPS[18,12]=1;
mxvaluesPS[19,13]=1;
mxvaluesPS[20,14]=0;
mxvaluesPS[7,3]=1;
mxvaluesPS[8,4]=1;
# free
mxvaluesPS[12,9]=ra  # 'ra';
mxvaluesPS[15,12]=ra  # 'ra';
mxvaluesPS [18,9]=ra  # 'ra';
mxvaluesPS[18,15]=ra   # 'ra';
mxvaluesPS[13,10]=rc # 'rc';
mxvaluesPS[19,16]=rc  # 'rc';
mxvaluesPS[16,13]=rc # 'rc';
mxvaluesPS[19,10]=rc  # 'rc';
mxvaluesPS[20,17]=re  # re
mxvaluesPS[14,11]=re  # re
#
mxvaluesPS[4,3]=rf
mxvaluesPS[8,7]=rf
mxvaluesPS[7,4]=rf
mxvaluesPS[8,3]=rf
#
diag(mxvaluesPS)=c(0,0,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
for (i in 1:ne) { for (j in 1:i) {
mxvaluesPS[j,i]=mxvaluesPS[i,j]}}
#
mxlabelsPS=matrix(NA,ne,ne)
mxlabelsPS[12,9]='ra';
mxlabelsPS[15,12]='ra';
mxlabelsPS[18,9]='ra';
mxlabelsPS[18,15]='ra';
mxlabelsPS[13,10]='rc';
mxlabelsPS[19,16]='rc';
mxlabelsPS[16,13]='rc';
mxlabelsPS[19,10]='rc';
mxlabelsPS[20,17]='re'; # labeled not free
mxlabelsPS[14,11]='re';  # labeled not free
mxlabelsPS[4,3]='rf'  # re
mxlabelsPS[8,7]='rf'  # re
mxlabelsPS[7,4]='rf'  # re
mxlabelsPS[8,3]='rf'  # re
#
for (i in 1:ne) { for (j in 1:i) {
mxlabelsPS[j,i]=mxlabelsPS[i,j]}}
#
mxvaluesPS2=matrix(1,ne,ne)
mxvaluesPS2[15,12]=.5
mxvaluesPS2[18,9]=.5
mxvaluesPS2[15,9]=.5;
mxvaluesPS2[18,12]=.5;
mxvaluesPS2[7,3]=.5;
mxvaluesPS2[8,4]=.5;
mxvaluesPS2[8,3]=.5;
mxvaluesPS2[7,4]=.5;
for (i in 1:ne) { for (j in 1:i) {
mxvaluesPS2[j,i]=mxvaluesPS2[i,j]}}
#
mxvaluesLY=matrix(0,ny,ne)
mxvaluesLY[1:ny,1:ny]=diag(ny)
mxfreeLY=matrix(F,ny,ne)
mxlabelsLY=matrix(NA,ny,ne)
mxfreeLY[3,3]=mxfreeLY[7,7]=T
mxvaluesLY[3,3]=mxvaluesLY[7,7]=x
mxlabelsLY[3,3]= mxlabelsLY[7,7]='x'
mxfreeLY[4,4]=mxfreeLY[8,8]=T
mxvaluesLY[4,4]=mxvaluesLY[8,8]=y
mxlabelsLY[4,4]= mxlabelsLY[8,8]='y'
ny1=ny-2
#
mxFilt=matrix(c(
1,0,0,0,0,0,0,0,
0,1,0,0,0,0,0,0,
0,0,1,0,0,0,0,0,
0,0,0,1,0,0,0,0,
0,0,0,0,1,0,0,0,
0,0,0,0,0,1,0,0),ny1,ny,byrow=T)
#
covmDZ=Sdz #
covmMZ=Smz[1:6,1:6]
covmDZ[lower.tri(covmDZ)] = t(covmDZ)[lower.tri(t(covmDZ))]
covmMZ[lower.tri(covmMZ)] = t(covmMZ)[lower.tri(t(covmMZ))]
meanmDZ=rep(0,8)
meanmMZ=rep(0,6)
vnames=
 c('BMI1','SBP1','GVBMI1','GVSBP1', 'BMI2','SBP2','GVBMI2','GVSBP2')
# ------------------------------ OpenMx proper starts here
dzselvars=vnames
mzselvars=vnames[1:6]
names(meanmDZ)=vnames
names(meanmMZ)=vnames[1:6]
rownames(covmMZ)=colnames(covmMZ)=vnames[1:6]
rownames(covmDZ)=colnames(covmDZ)=vnames[1:8]

Matmodel=mxModel('TwRM',
    mxMatrix(type='Full',nrow=ne,ncol=ne,
       free=mxfreeBE,value=mxvaluesBE,labels=mxlabelsBE,
       ubound=10,lbound=-10,dimnames=list(vnames_ne, vnames_ne), name='BE'),
    mxMatrix(type='Iden',nrow=ne,ncol=ne,name='I'),
    mxMatrix(type='Full',nrow=6,ncol=8,free=F,values=mxFilt, name='F'),
    mxMatrix(type='Full',nrow=ny,ncol=ne,
	 free=mxfreeLY,value=mxvaluesLY,labels=mxlabelsLY,
       ubound=10,lbound=-10,name='LY'),
    mxMatrix(type='Symm',nrow=ne,ncol=ne,
       free=mxfreePS,value=mxvaluesPS,labels=mxlabelsPS,
       ubound=10,lbound=-10,, dimnames=list(vnames_ne, vnames_ne),name='PS'),
    mxMatrix(type='Symm',nrow=ne,ncol=ne,
       free=F,value=mxvaluesPS2,labels=c(NA),name='PS2'),

    mxMatrix(type='Full',nrow=1,ncol=ny,
	 free=c(T,T,T,T,T,T,T,T),value=0,
       labels=c('meB','meS','meFB','meFS','meB','meS','meFB','meFS'),
       ubound=10,lbound=-10,name='dzmu'),

    mxAlgebra(expression=dzmu%*%t(F), name='mzmu'),
    mxAlgebra(expression=solve(I-BE),name='iBE'),
    mxAlgebra(expression=LY%*%iBE,  name='LYB'),
#
    mxAlgebra(expression=F%*%LYB%*%PS%*%t(LYB)%*%t(F),name='Smz'),
    mxAlgebra(expression=LYB%*%(PS*PS2)%*%t(LYB),name='Sdz')
        # PS2 is just ones and 1/2 in the appropriate positions.
 )

# a model the data, the fit function (MZ)
MZmodel=mxModel("MZ",
 #  mxData( observed=datmz[,1:6], type="raw"),   # the data
  mxData(observed=covmMZ, type="cov", means=meanmMZ,numObs = Nmz),
  mxExpectationNormal(covariance="TwRM.Smz",means="TwRM.mzmu",mzselvars),
  mxFitFunctionML()
                                       )
# a model the data, the fit function (DZ)
DZmodel=mxModel("DZ",
#   mxData(observed=datdz[,1:8], type="raw"),
     mxData(observed=covmDZ, type="cov", means= meanmDZ,numObs = Ndz),
     mxExpectationNormal(covariance="TwRM.Sdz",means="TwRM.dzmu",dzselvars),
     mxFitFunctionML()
                                    )
twinModel1 <-  mxModel("twinmodel1", Matmodel, MZmodel, DZmodel,
                       mxFitFunctionMultigroup( c("MZ","DZ") )
)
# fit the model
twinModel1out <- mxRun(twinModel1, silent=TRUE)

# Check identification?
if (do_id) {
  keep$id <- mxCheckIdentification(twinModel1,
                                        details = F)$status
}



if (do_pow == "connor") {
#
# alpha=0.05 # abovs 	# user specified: type I error prob.
#
  twinModel2=omxSetParameters(twinModel1,labels=c('g1'),free=F, values=0)
  twinModel2out <- mxRun(twinModel2, silent=TRUE)
  c1=mxCompare(twinModel1out,twinModel2out)
  twinModel2=omxSetParameters(twinModel1,labels=c('g2'),free=F, values=0)
  twinModel2out <- mxRun(twinModel2, silent=TRUE)
  c2=mxCompare(twinModel1out,twinModel2out)
  twinModel2=omxSetParameters(twinModel1,labels=c('g2','g1'),free=F, values=0)
  twinModel2out <- mxRun(twinModel2, silent=TRUE)
  c3=mxCompare(twinModel1out,twinModel2out)

# power g1
  ncp=c1[2,7]
  if (ncp<.00000001)  {ncp=0}   # ncp >= 0 but abs(ncp<.0000001) can happen
  dfs=c1[2,8]
  ca=qchisq(alpha,dfs,ncp=0,lower.tail=F)	 # critical value given alpha
  power=pchisq(ca,dfs,ncp=ncp,lower.tail=F)

#
  keep$ncpg1=ncp
  keep$dfg1=dfs
  keep$powg1=power

#
  ncp=c2[2,7]
  if (ncp<.000001)  {ncp=0}   # ncp < 0 but abs(ncp<.0000001) can happen
  dfs=c2[2,8]
  ca=qchisq(alpha,dfs,ncp=0,lower.tail=F)	 # critical value given alpha
  power=pchisq(ca,dfs,ncp=ncp,lower.tail=F)
#
  keep$ncpg2=ncp
  keep$dfg2=dfs
  keep$powg2=power

  ncp=c3[2,7]
  if (ncp<.000001)  {ncp=0}   # ncp < 0 but abs(ncp<.0000001) can happen
  dfs=c3[2,8]
  ca=qchisq(alpha,dfs,ncp=0,lower.tail=F)	 # critical value given alpha
  power=pchisq(ca,dfs,ncp=ncp,lower.tail=F)
#
  keep$ncpg12=ncp
  keep$dfg12=dfs
  keep$powg12=power


} else {
  keep$twins = umxPower(twinModel1out,  update = "g1", power = .8)[1]
  keep$group = group
}

# Twin correlations

keep$rmz1=diag(Smz[1:4,5:8])[1]
keep$rmz2=diag(Smz[1:4,5:8])[2]
keep$rdz1=diag(Sdz[1:4,5:8])[1]
keep$rmz2=diag(Sdz[1:4,5:8])[2]

parms_est <- parameters(twinModel1out) |>
  tidyr::gather(variable, value, -name) |>
  tidyr::spread(name, value)


keep$g1_diff <-  parms_est$g1 - g1


if (est_out) {
  keep$b1_hat = parms_est$b1
  keep$b3_hat = parms_est$b3
  keep$g1_hat = parms_est$g1
  keep$g2_hat = parms_est$g2
  keep$ra_hat = parms_est$ra
  keep$re_hat = parms_est$re
  keep$rf_hat = parms_est$rf
  keep$as_hat = parms_est$as
  keep$es_hat = parms_est$es
  keep$ab_hat = parms_est$ab
  keep$eb_hat = parms_est$eb
}
  return(keep)
}

## Simulation foreach loop, remember to save
## This one is the AE design
cl <- makeCluster(detectCores()/2)
registerDoParallel(cl)

ae <- foreach(b1=c(sqrt(.025),sqrt(.05),sqrt(.075)), .combine =rbind,
             .packages = c("umx", "MASS")) %:%
     foreach(b3=c(sqrt(.025),sqrt(.05),sqrt(.075)), .combine =rbind) %:%
      foreach(g1=c(sqrt(.020),sqrt(.040),sqrt(.060)), .combine =rbind) %:%
      foreach(g2=c(sqrt(.020),sqrt(.040),sqrt(.060)), .combine =rbind) %:%
      foreach(abs=c(.10,.25), .combine =rbind) %:%
      foreach(ass = c(.10,.25), .combine = rbind) %:%
      foreach(cbs=0, .combine =rbind) %:%
      foreach(css=0, .combine =rbind) %:%
      foreach(ra=c(.0, .25, .50), .combine =rbind) %:%
      foreach(rc=0, .combine =rbind) %:%
      foreach(re=c(.0, .25, .50), .combine =rbind) %:%
      foreach(rf=c(.0, .25, .50), .combine =rbind) %dopar% {
        sim(b1, b3, g1, g2, abs, ass, cbs, css, ra, rc, re, rf,,,,,noC = TRUE)
  }
stopCluster(cl)
