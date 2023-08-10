ft2m(feet::Real) = 0.3048 * feet
m2ft(meters::Real) = meters / 0.3048

knots2ms(knots::Real) = 0.5144444 * knots
ms2knots(mps::Real) = mps / 0.5144444

kel2cel(kelvin::Real) = kelvin - 273.15
cel2kel(celsius::Real) = celsius + 273.15
