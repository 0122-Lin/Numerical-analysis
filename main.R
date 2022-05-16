#В18. Лабораторная работа
# 2.1. Вычисление собственных векторов и собственных значений матриц

e<-0.0001
er<-0.1^(6)
library("matlab")
library("matlib")
library("psych")
library("expm")
#проверка
Sys.setlocale(category = "LC_ALL", locale = "ru")
Checker<-function(A,l,rez,n){
  check <-(A-eye(n)*l[1])%*%rez[,1]
  correct<-1
  for (i in 1:n)
    if (abs(check[i,1])>e)
      correct<-0
  return (c(correct,check))
}
Checker1<-function(A,l,rez,n){
  check <-(A-eye(n)*l)%*%rez[,1]
  correct<-1
  for(i in 1:n)
    check[i,1] <- check[i,1]*0.0001
  if (abs(check[,1])>e)
    correct<-0
  return (c(correct,check))
}
#собственные многочлены
P<-function(a1,n){
  P <- c()
  for (i in 1:n)
    P <- c(P,-a1[i,n])
  P <- c(P, 1)
  return(P)
}
#метод Данилевского
D<-function(A,n){
  a1 <- A
  a2 <- 0
  pS <- eye(n)
  S <- matrix(0,n,n)
  fEye <- eye(n-1)
  for (i in 2:n)
    for (j in 1:n-1)
      S[i,j]<-fEye[i-1,j]
  check <- 1
  check1<-1
  while(check1<=(n-1)){
    if(a1[1,n]==0)
      check <- 0
    if(check)
    {
      for (i in 1:n)
        S[i,n]<-a1[i,n]
      pS<-pS%*%S
      So <- solve(S)
      a2<-a1
      a1<-So%*%a2%*%S
      P<-P(a1,n)}
    else{
      i <- 1
      a<-0
      while(i<=n-check1){
        if(a1[i,n]!=0)
        {
          a<-a1[1,]
          a1[1,]<-a1[i,]
          a1[i,]<-a
          i <- n+1
          check <- 1
        }
        i <- i+1
      }
      if(!check)
      {
        print(check1)
        F <- matrix(0,n-check1,n-check1)
        print(F)
        for (i in check1:n){
          for (j in check1:n){

            print(a1[i,j])
            F[i-check1+1,j-check1+1]<-a1[i,j]
          }
        }
        B <- matrix(0,check1-1,check1-1)
        for (i in 1:check1-1){
          for (j in 1:check1-1){
            B[i,j]<-a1[i,j]
          }
        }
        P <- P(B,n)%*%P(F,n)
      }
    }
    check1<-check1+1
  }
  #собственные числа матрицы
  l <- polyroot(P)
  X<-matrix(0,n,n)
  for (z in 1:n)
  {
    y<-matrix(0,n,1)
    t <- 0
    for (i in 1:n)
    {for(j in 1:(n-i+1)){
      y[i]<-y[i]+(l[z]^(j-i+t))*P[j+i]}
      t <- t+1}
    x <- pS%*%y
    for (i in 1:n)
      x[i,1] <- x[i,1]/x[n,1]
    for (i in 1:n)
      X[i,z]<-x[i,1]
  }

  return(c(l,X))
}
#метод Фаддеева
F<-function(A,n){
  I<-eye(n)
  tA<-A
  q <- matrix(0,n,1)#след матрицы
  B <- A
  allB <- list()
  allB[[1]]<-B
  for(i in 1:(n-1))
  {
    q[i,1]<-tr(tA)/i
    B<-tA-q[i,1]*I
    allB[[i]]<-B
    tA <- A%*%B
  }
  q[n,1]<-tr(tA)/n
  B<-tA-q[n,1]*I
  allB[[n]]<-B

  P<-c(1)
  for (i in 1:n+1)
    P[i]<-1
  for (i in 1:n+1)
    P[n+2-i]<-(-q[i-1,1])
  #собственные числа матрицы

  l <- polyroot(P)
  check <- 1
  for (i in 1:(n-1)){
    tempL<-l[i]
    for(j in (i+1):n){
        if(l[j] == tempL)
        {check <- 0}
    }}
  if(check){
    X<-zeros(n)
    for(j in 1:n){
      ans <- (l[j]^(n-1))*I
      for (i in 1:(n-1)){
        m <- allB[[i]]
        ans <- ans + (l[j]^(n-i-1))*m
      }
      for(i in 1:n){
        X[i,j]<-ans[i,1]/ans[n,1]
      }
    }
    return(c(l,X))
  }
  else{return(c(l,0))}

}
#метод Леверье
L<-function(A,n){
  S<-matrix(0,n,1)
  P<-matrix(0,(n),1)
  for(i in 1:n)
    S[i]<-tr(A%^%i)
  P[1]<-S[1]
  for(i in 2:n){
    P[i]<-S[i]
    for ( j in 1:(i-1)){
      P[i]<-P[i]-(P[i-j]*S[j])}
    P[i] <- (1/i)*P[i]
  }
  P_inv <-matrix(1,(n+1),1)
  for(i in 2:(n+1)){
    P_inv[i]<--P[i-1]
  }
  temp<-0
  for (i in 1:((n+1)%/%2)){
    temp <- P_inv[i]
    P_inv[i]<-P_inv[(n+1)-i+1]
    P_inv[(n+1)-i+1]<-temp
  }
  l <-polyroot(P_inv)
  X<-matrix(0,n,n)
  for(i in 1:n){
    m<-(A-l[i]*eye(n))
    x <- solve(m)
    X[,i]<-x[,1]
  }
  for(i in 1:n){
    for(j in 1:n){
      X[i,j] <- X[i,j]/X[n,j]
    }
  }
  return(c(l,X))
}

#метод Крылова
Kr<-function(A,n){
  #собственный многочлен
  C<-eye(n,1)
  m <- 0
  C0 <- matrix(0,n,0)
  while(m<rankMatrix(C)){
     m<-m+1
     C<- cbind(C,(A%^%m)%*%C[,1])
  }
  for (i in 1:m){
      C0 <- cbind(C0,C[,(m+1-i)])}
  Q <- solve(C0,C[,(m+1)])
  P<- c(1)
  for (i in 1:(n+1))
    P[i]<-1
  for (i in 2:(n+1))
    P[n+2-i]<-(-Q[i-1])
  #собственные числа матрицы
  l <- polyroot(P)
  P<-c(1)
  for(i in 1:m){
    P[i+1]<- -Q[i]
  }
  #собственные векторы матрицы
   X<-matrix(0,n,n)
  for(i in 1:n){
    m<-(A-l[i]*eye(n))
    x <- solve(m)
    X[,i]<-x[,1]
  }
  for(i in 1:n){
    for(j in 1:n){
      X[i,j] <- X[i,j]/X[n,j]
    }
  }
  return(c(l,X))
}
#метод вращений
Rot<-function(A,n,er){
  S<-eye(n)
  At<-A
  sum<-er+1
  while (sum>er){
    max <- At[1,2]
    imax<-1
    jmax<-2
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        if(abs(max) < abs(At[i,j])){
          max<-At[i,j]
          imax<-i
          jmax<-j
        }
        }
      }
    d<-(2*At[imax,jmax])/(At[imax,imax]-At[jmax,jmax])
    cosf<-sqrt((1/2)*(1+(1+d^2)^(-1/2)))
    sinf<-sqrt((1/2)*(1-(1+d^2)^(-1/2)))*sign(d)

    V<-eye(n)
    V[imax,imax]<-cosf
    V[imax,jmax]<--sinf
    V[jmax,imax]<-sinf
    V[jmax,jmax]<-cosf
    At<-t(V)%*%At%*%V
    S<-S%*%V
    sum<-0
    for(i in 1:(n-1)){
     for(j in (i+1):n){
        sum<-sum+(At[i,j]^2)
      }}
    sum <- 2*sum
    }
  l<-diag(At)
  X<-S
  for(i in 1:n){
    for(j in 1:n){
      X[i,j]<-X[i,j]/X[n,j]
    }
  }
  return(c(l,X))
  }
#Степенной метод
Pow<-function(A,n,k){
  y0<-matrix(0,n,1)
  y0[1]<-1
  yk1<-A%^%(k)%*%(y0)
  yk2<-A%^%(k+1)%*%(y0)
  l<-sum(yk2*yk1)/sum(yk1*yk1)
  X<-c()
  for(i in 1:n){
      X[i] <- yk1[i]/yk1[n]
    }
  return(c(round(l),X))
}
#метод лямбда-разности
LR<-function (A,n,l,k){
  y0<-matrix(1,n,1)
  yk1<-(A%^%(k-1))%*%(y0)
  yk2<-(A%^%(k))%*%(y0)
  yk3<-(A%^%(k+1))%*%(y0)
  l2<- (yk3[1]-l*(yk2[1]))/(yk2[1]-l*(yk1[1]))
  X<-c()
  for(i in 1:n){
      X[i] <- (yk3[i]-l*yk2[i])/(yk3[n]-l*yk2[n])
    }
  return(c((l2),X))
}
main<-function(){
  #num <- c(1.28,0.38,0.24,0.12,0.38,0.85,0.22,-0.11,0.24,0.22,1.35,0.32,0.12,-0.11,0.32,0.83)
  #n <- 4
  num <- c(1,2,1,2,-1,0,1,0,2)# регулярный
  n<-3
  num1 <- c(1,2,3,2,3,4,3,4,5)
  n1 <- 3
  free <- c(-4.86, 1.36,-3.19,0.38)

  A <- matrix(num,n,n)
  B <- matrix(num1,n1,n1)
  f <- matrix(free)
  ##############################
  print('Метод Данилевского: ')
  print('Рег. сл.: ')
  print("Собственные значения: ")
  D <- D(A,n)
  l <- matrix(D[1:n],n,1)
  print(l)
  print("Собственные векторы: ")
  rez <- matrix(D[n+1:(n+n*n)],n,n)
  print(rez)
  print("Проверка: ")
  check <- Checker(A,l,rez,n)
  tempCheck <-check[2:(n+1)]
  for (i in 1:n)
    print(tempCheck[i])
  if(check[1])
    print("Верно")
  else
    print("Неверно")
  print('Нерег. сл.: ')
  print("Собственные значения: ")
  D <- D(B,n1)
  l1 <- matrix(D[1:n1],n1,1)
  print(l1)
  print("Собственные векторы: ")
  rez1 <- matrix(D[n1+1:(n1+n1*n1)],n1,n1)
  print(rez1)
  print("Проверка: ")
  check <- Checker(B,l1,rez1,n1)
  tempCheck <-check[2:(n1+1)]
  for (i in 1:n1)
    print(tempCheck[i])
  if(check[1])
    print("Верно")
  else
    print("Неверно")
  ##############################
  print('Метод Фаддеева: ')
  print("Собственные значения: ")
  F <- F(A,n)
  l <- matrix(F[1:n],n,1)
  print(l)
  print("Собственные векторы: ")
  rez <- matrix(F[n+1:(n+n*n)],n,n)
  print(rez)
  print("Проверка: ")
  check <- Checker(A,l,rez,n)
  tempCheck <-check[2:(n+1)]
  for (i in 1:n)
    print(tempCheck[i])
  if(check[1])
    print("Верно")
  else
    print("Неверно")
  ##############################
  print('Метод Леверье: ')
  print("Собственные значения: ")
  L <- L(A,n)
  l <- matrix(L[1:n],n,1)
  print(l)
  print("Собственные векторы: ")
  rez <- matrix(L[(n+1):(n+n*n)],n,n)
  print(rez)
  print("Проверка: ")
  check <- Checker(A,l,rez,n)
  tempCheck <-check[2:(n+1)]
  for (i in 1:n)
    print(tempCheck[i])
  if(check[1])
    print("Верно")
  else
    print("Неверно")
 ##############################
 print('Метод Крылова: ')
 print("Собственные значения: ")
 K <- Kr(A,n)
  l <- matrix(K[1:n],n,1)
  print(l)
  print("Собственные векторы: ")
  rez <- matrix(K[(n+1):(n+n*n)],n,n)
  print(rez)
   print("Проверка: ")
  check <- Checker(A,l,rez,n)
  tempCheck <-check[2:(n+1)]
  for (i in 1:n)
    print(tempCheck[i])
  if(check[1])
    print("Верно")
  else
    print("Неверно")
   ##############################
 print('Метод вращ.: ')
 print("Собственные значения: ")
 R <- Rot(A,n,er)
  l <- matrix(R[1:n],n,1)
  print(l)
  print("Собственные векторы: ")
  rez <- matrix(R[(n+1):(n+n*n)],n,n)
  print(rez)
   print("Проверка: ")
  check <- Checker(A,l,rez,n)
  tempCheck <-check[2:(n+1)]
  for (i in 1:n)
    print(tempCheck[i])
  if(check[1])
    print("Верно")
  else
    print("Неверно")
 ##############################
 print('Степ. метод: ')
 print("Собственное число: ")
 P <- Pow(A,n,10)
  l <- P[1]
  print(l)
  print("Собств. вектор: ")
  rez1 <- P[2:(n+1)]
  rez<-matrix(0,n,1)
  for (i in 1:n)
    rez[i,1]<-rez1[i]
  print(rez)
  print("Проверка: ")
  check <- Checker1(A,l,rez,n)
  tempCheck <-check[2:(n+1)]
  for (i in 1:n)
    print(tempCheck[i])
  if(check[1])
    print("Верно")
  else
    print("Неверно")
   ##############################
 print('Метод лямбда-разности: ')
 print("Собственное число: ")
 LR <- LR(A,n,l,10)
   l <- LR[1]
  print(l)
  print("Собств. вектор: ")
  rez1 <- LR[2:(n+1)]
  rez<-matrix(0,n,1)
  for (i in 1:n)
    rez[i,1]<-rez1[i]
  print(rez)
  print("Проверка: ")
  check <- Checker1(A,l,rez,n)
  tempCheck <-check[2:(n+1)]
  for (i in 1:n)
    print(tempCheck[i])
  if(check[1])
    print("Верно")
  else
    print("Неверно")

}
main()
