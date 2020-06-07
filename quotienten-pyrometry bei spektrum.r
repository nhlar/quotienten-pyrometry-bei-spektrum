# quotienten-pyrometry bei spektrumr
# 6.6.2020

rm(list = ls())

path <-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Parameters
real_curve = FALSE
real_curve = TRUE
data_file = "data.borchert/Metall_2s_1kW_nStelle_2_400.dat"
calib_messung_file = "data.borchert/B12A_Messung.txt"
calib_factor_file = "data.borchert/B12Aacalw.txt"
calib_emiss <- function(lambda,T) {
  eps_b0=0.5296; eps_b1=5.9531e-6; eps_b2=3.7677e-9
  eps_b0 - lambda * 10 * (eps_b1 + eps_b2*T)
} 
T_cal = 2448  # in K

d_lambda = 5  # Abstand in nm
T_range = seq(200,4000,10)
noise = 100
smooth_fac = 10000
lambda_restrict = c(100,1200)

#------------------------

source(file="quotienten-pyrometry bei spektrum - includes.r") # Definition of functions

#------------------------
# get data from measurement

  if (real_curve) {
    # get measurement data
    data.in <- read.table(data_file, header = FALSE, col.names=c("lambda","L"),	
      skipNul=TRUE, na.strings ="NaN",dec=".",sep=",")
    data.in$lambda = data.in$lambda/10
    # plot(data.in,type="l")
    lambda_range = data.in$lambda
    spec_data = data.in
    
    # get calibration measurement data
    calib_messung.in <- read.table(calib_messung_file, header = FALSE, col.names=c("lambda","L"), skipNul=TRUE,
                                   na.strings ="NaN",dec=",",sep="\t", comment.char=">")
    calib_messung = approx(calib_messung.in$lambda, calib_messung.in$L-0, xout=lambda_range)
    calib_messung = data.frame(lambda=calib_messung$x, L=calib_messung$y)
    # plot(calib_messung,type="l")
    spec_data = calib_messung
    
    # get emissivity for calibration
    calib_emiss.in = data.frame(lambda=lambda_range, emiss=calib_emiss(lambda_range,T_cal))
    # plot(calib_emiss.in,type="l")
    calib_emiss_func = approxfun(calib_emiss.in$lambda, calib_emiss.in$emiss)
    # plot(lambda_range, calib_emiss_func(lambda_range),type="l")
    
    # get calibration curve from file
    calib_func.in <- read.table(calib_factor_file, header = FALSE, col.names=c("lambda","factor"), skipNul=TRUE, 
                                na.strings ="NaN",dec=".",sep="", comment.char=">")
    calib_func.in$lambda = calib_func.in$lambda/10
    calib_func = approxfun(calib_func.in$lambda, calib_func.in$factor, ties=mean)
    # plot(lambda_range, calib_func(lambda_range),type="l")
  }
# windows()
# dev.list()["RStudioGD"]
# dev.set(dev.list()["RStudioGD"])

#------------------------
# make artifical data

  if (!real_curve) {
    lambda_range=seq(200,3000,1)
    calib_emiss_func = approxfun(c(300,600,10000),c(0.8,0.5,0.2), rule=2)
    calib_emiss_func = cobs(lambda_range,calib_emiss_func(lambda_range),lambda=10000)
    calib_emiss_func = approxfun(predict(calib_emiss_func,lambda_range))
    # plot(calib_emiss_func(lambda_range), type="l")

    calib_knots = list(x=c(100,300,500,1000,3500), y=c(1e-5,1e-5,1e-3,1e-5,1e-13))
    calib_sensibility_func = approxfun(calib_knots$x, calib_knots$y, rule=2)
    calib_sensibility_func = cobs(lambda_range,calib_sensibility_func(lambda_range),lambda=0, degree=2)
    calib_sensibility_func = approxfun(predict(calib_sensibility_func,lambda_range))
    plot(calib_sensibility_func(lambda_range), type="l")
    
    spec_data = spec_data_make(lambda_range, plot=TRUE, eps=calib_emiss_func, sensibility=calib_sensibility_func, T=T_cal, noise=noise)
    calib_messung = spec_data
  }

#------------------------
# get lambda-range
  lambda_range = range(spec_data$lambda)
  lambda_range = c(max(lambda_range[1],lambda_restrict[1]), min(lambda_range[2], lambda_restrict[2]))
  lambda_range = seq(lambda_range[1], lambda_range[2], d_lambda)

#------------------------
# calc reference black body matrix for T_range and lambda_range with radiation and ratio from adjacent wavelengths

  bb = black_body_ref(T_range, lambda_range, plot=TRUE)
  bb_L_group_T = group_by(bb$L,T)
  bb_L_group_lambda = group_by(bb$L,lambda)

#------------------------
# smooth calibrate data, adapt to lambda_range
  # spec_calib_messung_pre = smooth_curve(calib_messung, lambda=lambda_range, plot=TRUE, span=smooth_fac)
  title_y = eval(bquote(expression("L("*lambda*",T="*.(T_cal)*"K)") ))
  spec_calib_messung_pre = smooth_curve(calib_messung, lambda=lambda_range, title_y=title_y, plot=TRUE, fac=smooth_fac)
  colnames(spec_calib_messung_pre) = c("lambda", "L")
  spec_calib = spec_calib_messung_pre
  # plot(spec_calib, type = "l")

  # normalize to black body and calc normalisation factor
  spec_calib$L = spec_calib$L/calib_emiss_func(lambda_range)
  bb_ref = bb_calc(lambda_range, T_cal)
  # plot(bb_ref, type = "l")
  calib_norm = spec_calib$L/bb_ref$L
  # plot(calib_norm, type = "l")
  
#------------------------
# smooth data and adapt to lambda_range, calc quotient from adjacent wavelengths
  
  spec_data_pre = smooth_curve(spec_data, lambda=lambda_range, plot=TRUE, fac=smooth_fac)
  colnames(spec_data_pre) = c("lambda", "L")
  
  spec_data_pre2 = spec_data_pre
  spec_data_pre2$L = spec_data_pre$L/calib_norm
  # plot(spec_data_pre2, type = "l")
  
  spec_ratio = spec_data_ratio(spec_data_pre2,lambda=lambda_range, plot=TRUE)
  # plot(spec_ratio, type = "l")
  
#------------------------
# look for equivalent temperature

  # bb_L_group_lambda %$% plot(T~lambda) # operator aus magrittr
  
  T_find <- function(lambda_in, value, bb_) {
    bb_vals <- filter(bb_, near(lambda, lambda_in, tol = 2000*.Machine$double.eps))
    T = approx(bb_vals$bb_ratio,bb_vals$T, value)$y
  }

  T_res = list()
  for (xx in 1:length(spec_ratio$lambda)) {
    T_res=rbind(T_res, c(lambda=spec_ratio$lambda[xx], T_res=T_find(spec_ratio$lambda[xx], spec_ratio$ratio[xx], bb$ratio)))
  }
  T_res
  plot(T_res, type = "l")

# group_map(bb_group_lambda, ~min(.x$T))
# 
# # -----------------------------
# # 
# # library(purrr)
# 
# select(L_ref_ratio_group, lambda)
# L_ref=filter(L_ref_ratio_group, lambda == 0.000000325)
# filter(L_ref_ratio_group, T==200)
# approx(L_ref$value,L_ref$T,0.85)
# 
# filter(L_ref_ratio_group, map(L_ref_ratio_group, ~any(near(L_ref_ratio_group$lambda, 6.75e-07, tol = .Machine$double.eps^0.5)))[1] )
# L_ref_ratio_group[1,]
# 
# filter(map_lgl(x, ~ any(near(.x, c(0.5679, 5.6789), tol = 1e-4))))
# 
# 
# 
# # -----------------------------
# # Messung
# 
# T1 = 3463
# T2 = 3360
# sigma = 5.670400e-12 # Stefan-Boltzmann constant in W/(cm^2*K^4)
# q = sigma * (T1^4 - T2^4) # radiant emittance
# q
# 
# 
# # **********************************************************************
# 
# # Möglichkeiten zum Verbessern
# 
# L = expression(2*h*c^2/lambda^5*1/(exp(h*c/(lambda*k*T))-1));
# deriv(L_, "lambda");
# 
# library("Deriv")
# fx = Deriv("sin(x^2) * y", "x")
# fx_ = as.formula(paste("y ~ ", fx))
# eval(fx_)
# 
# deriv(~ x^2, "x")
# dx2x <- deriv(~ x^2, "x") ; dx2x
# mode(dx2x)
# eval(dx2x)
# dx2x.gradient
# x
# x <- -1:2
# 
# L = "2*h*c^2/lambda^5*1/(exp(h*c/(lambda*k*T))-1)"
# fx = Deriv(L, "lambda")
# 
# L1 = "2*h*c^2/lambda^5"
# Deriv(L1, "lambda")
# 
# deriv(~ 2*h*c^2/lambda^5*1/(exp(h*c/(lambda*k*T))-1), "lambda")
