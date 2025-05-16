program regularisation
  use fonctions
  implicit none
  integer :: N , i , k
  real*8 :: dx
  real*8,dimension(:),allocatable :: X , Y
  real*8 :: a , b
  print*,"Entrez le decoupage N de l'intervalle"
  print*,"  "
  read*,N
  print*,"  "
  print*,"Entrez k"
  print*,"  "
  read*,k

  print*,"Entrez a"
  print*,"  "
  read*,a
  print*,"  "
  print*,"Entrez b"
  print*,"  "
  read*,b

  open (unit = 100 , file = "data.txt", status = "new")
  
  
  dx = 2.0d0/(N - 1)
  allocate(X(N), Y(N))
  X(1) = -1.0d0
  Y(1) =  ((a + b)/2.0d0) + (tanh(k*X(1))/2.0d0) !(tanh(k*X(1)) + a + 1)*indicatrice_negatif(X(1)) + (tanh(k*X(1)) + b - 1)*indicatrice_positif(X(1))
  write(100,*) X(1),"  ",Y(1),"  ",discontinuite(X(1), a , b)
  do i = 2 , N
     X(i) = X(i - 1) + dx
     !print*,"X","   ","(", i ,")" , " = ","  ",X(i)
     Y(i) = ((a + b)/2.0d0) + (tanh(k*X(i))/2.0d0) !(tanh(k*X(i)) + a + 1)*indicatrice_negatif(X(i)) + (tanh(k*X(i)) + b - 1)*indicatrice_positif(X(i))
     write(100,*) X(i),"  ",Y(i),"  ",discontinuite(X(i), a , b)
  end do
  
  close(100)  
  deallocate(X)
  deallocate(Y)

end program regularisation

