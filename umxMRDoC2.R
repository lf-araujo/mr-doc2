## Beta software, test before production use. Note rf requirement there,
## which is something we are trying to change
## This function should be added to umx() in the future
## It is a shorter version of the MR-DoC2 model

umxMRDoC2 <- function(pheno, prss,
                      mzData = NULL, dzData = NULL,
                      data = NULL, zyg = NULL,
                    sep = "_T", rf = NULL,
                    name = "MRDoC2", autoRun = getOption("umx_auto_run"),
                     tryHard = c("no", "yes", "ordinal", "search"),
                    optimizer = NULL) {

  tryHard <- match.arg(tryHard)

  # Managing data
  if (!is.null(data)) {
    if ("tbl" %in% class(data)) {
      data = as.data.frame(data)
    }
    mzData = data[data[, zyg] %in% ifelse(is.null(mzData), "MZ", mzData), ]
    dzData = data[data[, zyg] %in% ifelse(is.null(dzData), "DZ", dzData), ]
  } else {
    if ("tbl" %in% class(mzData)) {
      mzData = as.data.frame(mzData)
      dzData = as.data.frame(dzData)
    }
  }

  vnames = tvars(c(pheno, prss), sep = sep)

  xmu_twin_check(
    selDVs = c(pheno, prss),
    sep = sep, dzData = dzData, mzData = mzData, enforceSep = TRUE,
    nSib = 2, optimizer = optimizer
  )

  mzData = xmu_make_mxData(mzData, manifests = vnames)
  dzData = xmu_make_mxData(dzData, manifests = vnames)


  ny = 8
  ne_ = 4
  ny1 = ny-2

  x=1 # x
  y=1

  mxlabelsBE <- matrix(c(
    NA, "g2", "b1", "b4",
    "g1", NA, "b2", "b3",
    NA, NA, NA, NA,
    NA, NA, NA, NA
  ),4,4, byrow = T)

  mxfreeBE <- matrix(c(
    F, T, T, F,
    T, F, F, T,
    F, F, F, F,
    F, F, F,F
  ),4,4, byrow = T)

  mxfreeA <- matrix(c(
    T,T,F,F,
    T,T,F,F,
    F,F,T,T,
    F,F,T,T), 4, 4, byrow=T)

  mxvaluesA <- matrix(c(
    1, 1, 0,0,
    1, 1, 0,0,
    0,0, x^2,x*y*rf,
    0, 0, x*y*rf,y^2), 4,4, byrow = T)

  mxlabelsA <- matrix(c(
    "ab2","abraas",NA,NA,
    "abraas","as2",NA,NA,
    NA,NA,"x2","xrfy",
    NA,NA,"xrfy","y2"), 4, 4, byrow=T)

  mxfreeC <- matrix(c(
    T,T,F,F,
    T,T,F,F,
    F,F,F,F,
    F,F,F,F),4, 4, byrow=T)

  mxlabelsC <- matrix(c(
    "cb2","cbrccs",NA,NA,
    "cbrccs","cs2",NA,NA,
    NA,NA,NA,NA,
    NA,NA,NA,NA),4, 4, byrow=T)

  mxfreeE <- matrix(c(
    T,T,F,F,
    T,T,F,F,
    F,F,F,F,
    F,F,F,F),4, 4, byrow=T)

  mxlabelsE <- matrix(c(
    "eb2","ebrees",NA,NA,
    "ebrees","es2",NA,NA,
    NA,NA,NA,NA,
    NA,NA,NA,NA),4, 4, byrow=T)

  mxFilt <- matrix(c(
    1,0,0,0,0,0,0,0,
    0,1,0,0,0,0,0,0,
    0,0,1,0,0,0,0,0,
    0,0,0,1,0,0,0,0,
    0,0,0,0,1,0,0,0,
    0,0,0,0,0,1,0,0),ny1,ny,byrow=T)

top=mxModel('TwRM',
    mxMatrix(type='Full',nrow=ne_,ncol=ne_,
       free=mxfreeBE,labels=mxlabelsBE,
       ubound=10,lbound=-10, name='BE'),
    mxMatrix(type='Iden',nrow=ne_,ncol=ne_,name='I'),
    mxMatrix(type='Full',nrow=6,ncol=8,free=F,values=mxFilt, name='F'),

    mxMatrix(type='Full',nrow=ne_,ncol=ne_,
     free = F,value = diag(ne_),labels = NA,
       ubound=10,lbound=-10,name='LY'),

    mxMatrix(type='Symm',nrow=ne_,ncol=ne_,
       free=mxfreeA,labels=mxlabelsA,
       ubound=10,lbound=-10,name='A'),

    mxMatrix(type='Symm',nrow=ne_,ncol=ne_,
       free=mxfreeC,labels=mxlabelsC,
       ubound=10,lbound=-10,name='C'),

    mxMatrix(type='Symm',nrow=ne_,ncol=ne_,
       free=mxfreeE,labels=mxlabelsE,
       ubound=10,lbound=-10,name='E'),

    mxMatrix(type='Full',nrow=1,ncol=ny,
     free=T,value=0,
       labels=c('meB','meS','meFB','meFS','meB','meS','meFB','meFS'),
       ubound=10,lbound=-10,name='dzmu'),
    mxAlgebra(expression=dzmu%*%t(F), name='mzmu'),

    mxAlgebra(expression=solve(I-BE),name='iBE'),

    mxAlgebra(expression=iBE %*% A %*% t(iBE), name='A_'),
    mxAlgebra(expression=iBE %*% C %*%t(iBE), name='C_'),
    mxAlgebra(expression=iBE%*%E%*%t(iBE), name='E_'),
    mxAlgebra(expression=A_ + C_ + E_, name='SPh'),
    mxAlgebra(expression=rbind(
                         cbind(SPh, A_+C_),
                         cbind(A_+C_, SPh)),name='Smz_'),
    mxAlgebra(expression=rbind(
                         cbind(SPh, .5%x%A_+C_),
                         cbind(.5%x%A_+C_, SPh)),name='Sdz'),

    mxAlgebra(expression=F%*%Smz_%*%t(F),name='Smz')
 )

MZ=mxModel("MZ",mzData,
  mxExpectationNormal(covariance="TwRM.Smz",means="TwRM.mzmu",vnames[1:6]),
  mxFitFunctionML()
)

DZ=mxModel("DZ",dzData,
     mxExpectationNormal(covariance="TwRM.Sdz",means="TwRM.dzmu",vnames),
     mxFitFunctionML()
)


model =  mxModel("MRDoC2", top, MZ, DZ,
mxFitFunctionMultigroup( c("MZ","DZ") )
 )
  # model = as(model, "MxModel")
  model = mxAutoStart(model)
  # model = umxModify(model, update = c("x2", "y2"),  value = 1)

  model = xmu_safe_run_summary(model, autoRun = autoRun, tryHard = tryHard,
                               std = T)
}
