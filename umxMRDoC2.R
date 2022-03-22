## Beta software, test before production use. 
## This function is adapted to work with umx() 
## It is a shorter version, identical spec of the MR-DoC2 model

umxMRDoC2 <- function(pheno, prss,
                      mzData = NULL, dzData = NULL,
                      data = NULL, zyg = NULL,
                    sep = "_T", 
                    name = "MRDoC2", autoRun = getOption("umx_auto_run"),
                     tryHard = c("no", "yes", "ordinal", "search"),
                    optimizer = NULL) {

  tryHard <- match.arg(tryHard)
  options(mxByrow = TRUE)

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

  x=y=1
  rf=0.1

  top = mxModel("top",
    umxMatrix("BE", type = "Full",nrow=4,  ncol = 4,
              labels = c(NA, "g2", "b1", "b4",
                         "g1", NA, "b2", "b3",
                         NA, NA, NA, NA,
                         NA, NA, NA, NA),
              free = c(F, T, T, F,
                       T, F, F, T,
                       F, F, F, F,
                       F, F, F, F)),
    umxMatrix('I', type='Iden', nrow= 4,ncol= 4 ),
    umxMatrix('F', type='Full', nrow=6, ncol=8, free=F,
              values=c(1,0,0,0,0,0,0,0,
                       0,1,0,0,0,0,0,0,
                       0,0,1,0,0,0,0,0,
                       0,0,0,1,0,0,0,0,
                       0,0,0,0,1,0,0,0,
                       0,0,0,0,0,1,0,0)),
    umxMatrix('LY', type='Full',nrow=4, ncol = 4, free = F,values = diag(4),
              labels = NA),
    umxMatrix('A', type='Symm', nrow=4,ncol = 4,
              labels=c("ab2","abraas",NA,NA,
                       "abraas","as2",NA,NA,
                       NA,NA,"x2","xrfy", 
                       NA,NA,"xrfy","y2"),
              free=c(T,T,F,F,
                     T,T,F,F,
                     F,F,T,T,
                     F,F,T,T)),
    umxMatrix('C', type='Symm',nrow=4, ncol = 4,
              labels =c("cb2","cbrccs",NA,NA,
                       "cbrccs","cs2",NA,NA,
                       NA,NA,NA,NA,
                       NA,NA,NA,NA),
              free=c(T,T,F,F,
                     T,T,F,F,
                     F,F,F,F,
                     F,F,F,F)),
    umxMatrix('E', type='Symm', nrow=4, ncol = 4,
              labels =c("eb2","ebrees",NA,NA,
                       "ebrees","es2",NA,NA,
                       NA,NA,NA,NA,
                       NA,NA,NA,NA),
              free= c(T,T,F,F,
                      T,T,F,F,
                      F,F,F,F,
                      F,F,F,F)),
    umxMatrix('dzmu', type='Full', nrow=1, ncol=8, free=T, value=0,
       labels=c('meB','meS','meFB','meFS','meB','meS','meFB','meFS')),
    mxAlgebra('mzmu', expression = dzmu %*% t(F)),
    mxAlgebra('A_', expression = solve(I - BE) %&% A),
    mxAlgebra('C_', expression = solve(I - BE) %&% C),
    mxAlgebra('E_', expression = solve(I - BE) %&% E),
    mxAlgebra('SPh', expression= A_ + C_ + E_),
    mxAlgebra('Smz_', expression=rbind(
                         cbind(SPh,A_+C_),
                         cbind(A_+C_,SPh))),
    mxAlgebra('Sdz', expression=rbind(
                         cbind(SPh,.5%x%A_+C_),
                         cbind(.5%x%A_+C_,SPh))),
    mxAlgebra('Smz', expression= F%&%Smz_)
  )

  MZ = mxModel("MZ", mzData,
    mxExpectationNormal(covariance = "top.Smz",means = "top.mzmu", vnames[1:6]),
    mxFitFunctionML()
  )

  DZ = mxModel("DZ", dzData,
       mxExpectationNormal(covariance = "top.Sdz",means = "top.dzmu", vnames),
       mxFitFunctionML()
  )

  model = mxModel(name, top, MZ, DZ, 
                  mxFitFunctionMultigroup(c("MZ","DZ") ))

  model = mxAutoStart(model)
  model = mxRun(model)
  print(summary(model))
  return(model)
}


