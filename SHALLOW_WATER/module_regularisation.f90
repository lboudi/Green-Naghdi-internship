module fonctions
  implicit none

contains
  function indicatrice_negatif(x)
    real*8,intent(in) :: x
    real*8 :: indicatrice_negatif

    if ( x <= 0 ) then
       indicatrice_negatif = 1.0d0
    else
       indicatrice_negatif = 0.0d0
    end if

  end function indicatrice_negatif

  function indicatrice_positif(x)
    real*8,intent(in) :: x
    real*8 :: indicatrice_positif
    
    if ( x > 0 ) then
       indicatrice_positif = 1.0d0
    else
       indicatrice_positif = 0.0d0
    end if
  end function indicatrice_positif

  function discontinuite(x , a ,b)
    real*8,intent(in) :: x , a  , b

    real*8 :: discontinuite

    if( x <= 0 ) then
       discontinuite = a
    else
       discontinuite = b
    end if
  end function discontinuite
  

end module fonctions

