a1 <- function(a,c) {
  print(a)
  pars <- as.list(match.call()[-1])
  print(pars$a)
  print(pars$c)
  # data[,as.character(pars$a)]
}

b = 2
a1(a=b,c=3)
a1()
formals(a1)$a

formals(a1)$a =3
formalArgs(a1)
formalArgs(args(a1))
