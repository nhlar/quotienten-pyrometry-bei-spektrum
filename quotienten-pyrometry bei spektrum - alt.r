# quotienten-pyrometry bei spektrumr
# 31.5.2020

rm(list = ls())

# Parameter

T_fix = 1000
T_range = seq(293,3000,1)
T_range = seq(300,3000, 700)
d_lambda = 10  # Abstand in nm
lambda_range = seq(100,5000,d_lambda) # alles in nm

c = 299792458 # in m s
h = 6.62607015e-34   # in J s
k = 1.380649e-23  # in J s
lambda_range = lambda_range/1e9
  
# library("reshape2")
library(dplyr)

L <- function(lambda,T) 2*h*c^2/lambda^5*1/(exp(h*c/(lambda*k*T))-1)
L_ = deriv(~ 2*h*c^2/lambda^5*1/(exp(h*c/(lambda*k*T))-1), "lambda"); 

L_matrix = sapply(T_range, L, lambda=lambda_range)
L_set = melt(L_matrix, varnames=c("lambda","T"), value.name="value")
L_set$lambda = lambda_range[L_set$lambda]
L_set$T = T_range[L_set$T]

L_quote = function(lambda_range,T) L(lambda_range[-length(lambda_range)],T)/L(lambda_range[-1],T)
# L_quote = function(lambda_range,T) L(lambda_range[-length(lambda_range)],T)
L_quote_matrix = sapply(T_range, L_quote, lambda_range=lambda_range)
 
L_quote_set <- melt(L_quote_matrix, varnames=c("lambda","T"), value.name="value")
L_quote_set$lambda = lambda_range[L_quote_set$lambda]
L_quote_set$T = T_range[L_quote_set$T]

# -----------------------------
# Plots black body

library(ggplot2)
windows(height=4.5, width=7, xpos=10, ypos=100)
q = ggplot() + ggtitle("Black body radiation") + scale_x_continuous(labels = scales::label_number(scale = 1e9)) + theme_bw()
q = q + xlab(expression(lambda*" in nm")) + ylab(eval(bquote(expression("L("*lambda*",T="*.(T)*"K)") )))
q = q+geom_line(data=L_set, aes(x=lambda,y=value, group = T, colour = T))
q

windows(height=4.5, width=7, xpos=10, ypos=100)
q = ggplot() + ggtitle("Black body radiation") + scale_x_continuous(labels = scales::label_number(scale = 1e9)) + theme_bw()
q = q + xlab(expression(lambda*" in nm")) + ylab(eval(bquote(expression("Ratio") )))
q = q+geom_line(data=L_quote_set, aes(x=lambda,y=value, group = T, colour = T))
q

# -----------------------------
# 
select(L_quote_set, lambda)
aa=filter(L_quote_set, lambda == 1.01e-6)
approx(aa$value,aa$T,0.85)

# -----------------------------

# library(hrbrthemes)
# windows(height=4.5, width=7, xpos=10, ypos=100)
df <- data.frame(
  x = lambda_range,
  y = L(lambda_range,T)
)
q = ggplot() + ggtitle("Black body radiation") + scale_x_continuous(labels = scales::label_number(scale = 1e9))
q = q + xlab(expression(lambda*" in nm")) + ylab(eval(bquote(expression("L("*lambda*",T="*.(T)*"K)") )))

q = q+geom_line(data=df, aes(x,y), colour="blue")


#-----------------------------------------------------------------------------
# Need to define the second axis transformation to be the inverse of the data
# transformation to everything cancels out appropriately
#-----------------------------------------------------------------------------
lambda=lambda_range
ds <- data.frame(
  x = lambda_range,
  y = unname(attributes(eval(L_))[[1]])
)

# second axis needs treansformations. See more http://rstudio-pubs-static.s3.amazonaws.com/329613_f53e84d1a18840d5a1df55efb90739d9.html

scale.a      <- range(df$y)
scale.b      <- range(ds$y)
scale.factor <- diff(scale.a)/diff(scale.b)
trans <- ~ ((. - scale.a[1]) / scale.factor) + scale.b[1]

ds$y_        <- ((ds$y - scale.b[1]) * scale.factor) + scale.a[1]

rm(lambda)
q = q + scale_y_continuous(sec.axis = sec_axis(trans=trans, name = "dL"))
cols = hcl(c(15, 15+180), 100, 65)
q = q + theme_bw()
q= q+geom_line(data=ds, aes(x,y_), colour="red", show.legend=TRUE)  +
   theme( plot.title = element_text(hjust = 0.5,colour="red"),
        axis.text.y.right=element_text(colour="red"),
        axis.ticks.y.right=element_line(colour="red"),
        axis.title.y.right=element_text(colour="red"),
        axis.text.y=element_text(colour="blue"),
        axis.ticks.y=element_line(colour="blue"),
        axis.title.y=element_text(colour="blue")) +
        scale_colour_manual(values=cols) 
q

# fill matrix from 293 K to 3000 K for black body



#---------------------------------------------------------
# Messung





jitter(3)




f = expression(x^3)
dx2x <- D(f,"x")
x=1:10
eval(dx2x)

T1 = 3463
T2 = 3360
sigma = 5.670400e-12 # Stefan-Boltzmann constant in W/(cm^2*K^4)
q = sigma * (T1^4 - T2^4) # radiant emittance
q


# **********************************************************************

# Möglicvhkeiten zum Verbessern

L = expression(2*h*c^2/lambda^5*1/(exp(h*c/(lambda*k*T))-1));
deriv(L_, "lambda");

library("Deriv")
fx = Deriv("sin(x^2) * y", "x")
fx_ = as.formula(paste("y ~ ", fx))
eval(fx_)

deriv(~ x^2, "x")
dx2x <- deriv(~ x^2, "x") ; dx2x
mode(dx2x)
eval(dx2x)
dx2x.gradient
x
x <- -1:2

L = "2*h*c^2/lambda^5*1/(exp(h*c/(lambda*k*T))-1)"
fx = Deriv(L, "lambda")

L1 = "2*h*c^2/lambda^5"
Deriv(L1, "lambda")

deriv(~ 2*h*c^2/lambda^5*1/(exp(h*c/(lambda*k*T))-1), "lambda")
