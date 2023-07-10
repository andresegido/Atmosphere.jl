meters(feet::Real) = 0.3048 * feet
feet(meters::Real) = meters / 0.3048

mps(knots::Real) = 0.5144444 * knots
knots(mps)::Real = mps / 0.5144444

celsius(kelvin::Real) = kelvin - 273.15
kelvin(celsius::Real) = celsius + 273.15
