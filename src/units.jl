"""
    ft2m(feet::Real)
Returns a copy of argument "feet" converted to meters
"""
ft2m(feet::Real) = 0.3048 * feet
"""
    m2ft(feet::Real)
Returns a copy of argument "meters" converted to feet
"""
m2ft(meters::Real) = meters / 0.3048

"""
    knots2ms(knots::Real)
Returns a copy of argument "knots" converted to meters per second
"""
knots2ms(knots::Real) = 0.5144444 * knots
"""
    ms2knots(mps::Real)
Returns a copy of argument "mps" (meters per second) converted to meters per second
"""
ms2knots(mps::Real) = mps / 0.5144444

"""
    kel2cel(kelvin::Real)
Returns a copy of argument "kelvin" converted to celsius
"""
kel2cel(kelvin::Real) = kelvin - 273.15
"""
    cel2kel(kelvin::Real)
Returns a copy of argument "celsius" converted to kelvin
"""
cel2kel(celsius::Real) = celsius + 273.15
