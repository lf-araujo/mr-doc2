## Beta software, test before production use
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

 top = mxModel('top',
               umxMatrix('BE',  'Full', nrow = ne_, ncol = ne_,
                         free = mxfreeBE, labels = mxlabelsBE),
               umxMatrix('I', 'Iden', nrow = ne_, ncol = ne_),
               umxMatrix('F', 'Full', nrow = 6, ncol = 8, free = F,
                         values = mxFilt),
               umxMatrix('LY', 'Full', nrow = ne_, ncol = ne_,
                         free = F, value = diag(ne_)),
              umxMatrix('A', 'Symm', nrow = ne_, ncol = ne_, free = mxfreeA,
                        value = mxvaluesA, labels=mxlabelsA),
              umxMatrix('C', 'Symm', nrow=ne_, ncol=ne_, free = mxfreeC,
                        labels=mxlabelsC),
              umxMatrix('E', 'Symm',nrow=ne_,ncol=ne_, free = mxfreeE,
                        labels = mxlabelsE),
              umxMatrix('dzmu','Full',nrow = 1,ncol = ny, free = T,value = 0,
                 labels = c('meB', 'meS', 'meFB', 'meFS',
                            'meB', 'meS', 'meFB', 'meFS')),
              mxAlgebra(expression = dzmu%*%t(F), name = 'mzmu'),
              mxAlgebra(expression = solve(I-BE) %&% A, name = 'Av'),
              mxAlgebra(expression = solve(I-BE) %&% C, name = 'Cv'),
              mxAlgebra(expression = solve(I-BE) %&% E, name = 'Ev'),
              mxAlgebra(expression = A+C+E, name = 'SPh'),
              mxAlgebra(expression = rbind(
                                   cbind(SPh, Av+Cv),
                                   cbind(Av+Cv, SPh)),name='Smzv'),
              mxAlgebra(expression = rbind(
                                   cbind(SPh, .5 %x% Av + Cv),
                                   cbind(.5 %x% Av + Cv, SPh)),
                        name = 'expCovDZ'),
              mxAlgebra(expression = F %&% Smzv, name = 'expCovMZ')
   )

  MZ = mxModel("MZ", mzData,
               # mxExpectationRAM(M="mzmu", A = "BE",  F = "F",   vnames[1:6] ),
               mxExpectationNormal(covariance = "top.expCovMZ",
                                   means = "top.mzmu",vnames[1:6]),
               mxFitFunctionML())

  DZ = mxModel("DZ", dzData,
               # mxExpectationRAM(M="dzmu", A = "BE", F = "F", vnames ),
               mxExpectationNormal(covariance="top.expCovDZ",
                                   means="top.dzmu",vnames),
               mxFitFunctionML())

  model = mxModel("MRDoC2", top, MZ, DZ,
                  mxFitFunctionMultigroup(c("MZ","DZ")))
  # model = as(model, "MxModel")
  model = mxAutoStart(model)
  model = umxModify(model, update = c("x2", "y2"),  value = 1)

  model = xmu_safe_run_summary(model, autoRun = autoRun, tryHard = tryHard,
                               std = T)
}
