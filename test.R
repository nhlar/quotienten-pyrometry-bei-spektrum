# test
# 7.6.2020

f1 <- function(x) print(x)
A <- 2
f1(A)
do.call("f1", list(A)) 
do.call(f1, list(A)) 

do.call("f1", list(A)) 
do.call(f1, list(A)) 

do.call("f1", as.list(spec_ratio$lambda))


A <- 2
f <- function(x) print(x^2)
env <- new.env()
assign("A", 10, envir = env)
assign("f", f, envir = env)
eval(print(A), envir = env)
ls(A, envir = env)
get("A", envir = env)

f <- function(x) print(x)
f(A)                                      # 2
do.call("f", list(A))                     # 2
do.call("f", list(A), envir = env)        # 4
do.call( f,  list(A), envir = env)        # 2
do.call("f", list(quote(A)), envir = env) # 100
do.call( f,  list(quote(A)), envir = env) # 10
do.call("f", list(as.name("A")), envir = env) # 100

eval(call("f", A))                      # 2
eval(call("f", quote(A)))               # 2
eval(call("f", A), envir = env)         # 4
eval(call("f", quote(A)), envir = env)  # 100

