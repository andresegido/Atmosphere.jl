using Roots: find_zero
RealOrNothing = Union{Real,Nothing}

"""
    FlightLevel(mach::Real, h::Real, eas::Real, cas::Real, tas::Real, a::Real, P::Real, T::Real, ρ::Real, 
                μ::Real, ν::Real, PDyn::Real, Re_x::Real, ΔT::Real, feet::Bool, knots::Bool, celsius::Bool)

Structure with several fields associated with a flight level and condition with atmospheric data
"""
struct FlightLevel
    "Mach number"
    mach::Real
    "Atmospheric altitude (in feet or meters; see Boolean `feet`)"
    h::Real
    "Equivalent air speed (in knots or meters per second; see Boolean `knots`)"
    eas::Real
    "Calibrated air speed (in knots or meters per second; see Boolean `knots`)"
    cas::Real
    "True air speed (in knots or meters per second; see Boolean `knots`)"
    tas::Real
    "Sound speed (in knots or meters per second; see Boolean `knots`)"
    a:: Real
    "Pressure (in Pascals)"
    P::Real
    "Temperature (in Celsius or Kelvin; see Boolean `celsius`)"
    T::Real
    "Density (in kilograms per cubic meter)"
    ρ::Real
    "Viscosity (in Pascals·second)"
    μ::Real
    "Kinematoc viscosity (in square meters·second)"
    ν::Real
    "Dynamic pressure (in Pascals)"
    PDyn::Real
    "Reynolds number per length (in 1/feet or 1/meters; see Boolean `feet`)"
    Re_x::Real
    "Temperature deviation from ISA"
    ΔT::Real
    "Show length magnitudes in feet (true) or meters (false)"
    feet::Bool
    "Show speed magnitudes in knots (true) or meters per second (false)"
    knots::Bool
    "Show temperature magnitudes in Celsius (true) or Kelvin (false)"
    celsius::Bool
end


function Base.show(io::IO, fl::FlightLevel)
    fl.feet ? distUnit = "ft" : distUnit = "m"
    fl.knots ? velUnit = "kts" : velUnit = "m/s"
    fl.celsius ? tempUnit = "C" : tempUnit = "K"

    pressUnit = "Pa"
    densUnit = "kg / m3"
    viscUnit = "Pa · s"
    kviscUnit = "m2 · s"

    print(io, 
    "Flight condition defined by:\n",
    "Mach number (.mach) = ", fl.mach, "\n",
    "Height (.h) = ", fl.h, " ", distUnit, "\n",
    "Equivalent Air Speed (.eas) = ", fl.eas, " ", velUnit, "\n",
    "Calibrated Air Speed (.cas) = ", fl.cas, " ", velUnit, "\n",
    "True Air Speed (.eas) = ", fl.tas, " ", velUnit, "\n",
    "Sound Speed (.a) = ", fl.a, " ", velUnit, "\n",
    "Pressure (.P) = ", fl.P, " ", pressUnit, "\n",
    "Temperature (.T) = ", fl.T, " ", tempUnit, "\n",
    "Density (.ρ) = ", fl.ρ, " ", densUnit, "\n",
    "Viscosity (.μ) = ", fl.μ, " ", viscUnit, "\n",
    "Kinematic viscosity (.ν) = ", fl.ν, " ", kviscUnit, "\n",
    "Dynamic pressure (.PDyn) = ", fl.PDyn, " ", pressUnit, "\n",
    "Reynolds / distance (.Re_x) = ", fl.Re_x, " ", distUnit, "^-1 \n",
    "Temperature deviation from ISA (.ΔT) = ", fl.ΔT, " ", tempUnit, "\n",
    )
end

"""
    FlightLevel(;mach::RealOrNothing=nothing, h::RealOrNothing=nothing, eas::RealOrNothing=nothing, 
    cas::RealOrNothing=nothing, tas::RealOrNothing=nothing, ΔT::Real=0, 
    feet::Bool=false, knots::Bool=false, celsius::Bool=false)

Create a FlightLevel structure based on ISA functions and transformstions between equivalent,
calibrated and true airspeed. 

From the following arguments, TWO AND ONLY TWO must be set:
- `mach::Real`: Mach number
- `h::Real`: Atmospheric altitude
- `eas::Real`: Equivalent air speed
- `cas::Real`: Calibrated air speed
- `tas::Real`: True air speed

The following arguments are optional:
- `ΔT::Real`: Temperature deviation from ISA
- `feet::Bool`: "Set length magnitudes in feet (true) or meters (false)"
- `knots::Bool`: "Set spped magnitudes in knots (true) or meters per second (false)"
- `celsius::Bool`: "Set temperature magnitudes in Celsius (true) or Kelvin (false)"
"""
function FlightLevel(;mach::RealOrNothing=nothing, h::RealOrNothing=nothing, eas::RealOrNothing=nothing, 
                     cas::RealOrNothing=nothing, tas::RealOrNothing=nothing, ΔT::Real=0, 
                     feet::Bool=false, knots::Bool=false, celsius::Bool=false)
    
    # Check which inputs are introduced and which ones are Nothing
    inputs = filter(p->p[1] in Set([:mach :h :eas :cas :tas]) , Base.@locals)
    filledInputs = [k for (k,v) in inputs if v !== nothing]

    # Exception if inputs defined are not 2
    i = length(filledInputs)
    i == 2 || error("""From the available input options for Mach number, height, EAS, CAS and TAS,
                    2 must be defined instead of $i.""")

    # Exception if pair of input values is not valid
    (
        (:mach in filledInputs && :h in filledInputs) || 
        (:mach in filledInputs && :eas in filledInputs) ||  
        (:mach in filledInputs && :cas in filledInputs) ||  
        (:h in filledInputs && :eas in filledInputs) ||
        (:h in filledInputs && :cas in filledInputs) ||
        (:h in filledInputs && :tas in filledInputs) ||
        error("Pair of values $filledInputs is not valid to define a flight level.")
    )

    # Convert inputs to International Units to work with ISA functions

    if :h in filledInputs && feet
        hFL = ft2m(h)
    elseif :h in filledInputs && !feet
        hFL = h
    end
    
    # If a speed was an input, get its FL field already
    if knots
        :eas in filledInputs && (easFL = knots2ms(eas))
        :cas in filledInputs && (casFL = knots2ms(cas))
        :tas in filledInputs && (tasFL = knots2ms(tas))
    else
        :eas in filledInputs && (easFL = eas)
        :cas in filledInputs && (casFL = cas)
        :tas in filledInputs && (tasFL = tas)
    end

    # Solve for height and ISA first

    if :h in filledInputs
        P_FL = P(hFL)
    else
        if :eas in filledInputs
            P_FL = (easFL / a0 / mach) ^ 2 * P0
        elseif :cas in filledInputs
            # Only for subsonic. Supersonic expression must be coded
            P_FL = P0 * ((1 + 0.2 * (casFL / a0) ^ 2) ^ 3.5 - 1) / ((1 + 0.2 * mach ^ 2) ^ 3.5 - 1)
        end
        hFL = find_zero(x -> P(x) - P_FL, [h_lim[1], h_lim[end]])
    end
    T_FL = T(hFL, ΔT)
    ρFL = ρ(hFL, ΔT)
    μFL = μ(hFL, ΔT)
    νFL = μFL / ρFL
    aFL = a(hFL, ΔT)

    # Solve for Mach number, speeds and other derived variables

    if :mach in filledInputs
        machFL = mach
    else
        if :tas in filledInputs
            machFL = tasFL / aFL
        elseif :eas in filledInputs
            machFL = eas / a0 * √(P0 / P_FL)
        elseif :cas in filledInputs
            machFL = √(5 * ((P0 / P_FL * ((1 + 0.2 * (casFL / a0) ^ 2) ^ 3.5 - 1) + 1) ^ (2/7) - 1))
        end
    end

    # If a speed was NOT an input, compute it
    :tas in filledInputs || (tasFL = aFL * machFL)
    :eas in filledInputs || (easFL = tasFL * √(ρFL/ρ0))
    :cas in filledInputs || (
        casFL = a0 * √(5 * ((P_FL / P0 * ((1 + 0.2 * machFL ^ 2) ^ 3.5 - 1) + 1) ^ (2/7) - 1))
        )
        # 5 * ((self.pressure / P0 * ((1 + 0.2 * self.mach ** 2) ** 3.5 - 1) + 1) ** (2/7) - 1)
    PDynFL = ρFL * tasFL ^ 2 / 2
    Re_xFL = tasFL / νFL

    # Convert units to requested. Careful with keeping original values.
    # Add new unit conversions in the future
    if feet
        :h in filledInputs ? hFL = h : hFL = m2ft(hFL)
        Re_xFL /= m2ft(1)
    end
    if knots
        :eas in filledInputs ? easFL = eas : easFL = ms2knots(easFL)
        :cas in filledInputs ? casFL = cas : casFL = ms2knots(casFL)
        :tas in filledInputs ? tasFL = tas : tasFL = ms2knots(tasFL)
        aFL = ms2knots(aFL)
    end
    if celsius
        T_FL = kel2cel(T_FL) 
    end

    # return FlightLevel struct
    return FlightLevel(
        machFL, hFL, easFL, casFL, tasFL, 
        aFL, P_FL, T_FL, ρFL, μFL, νFL, PDynFL, Re_xFL, 
        ΔT, 
        feet, knots, celsius
        )

end

"""
    FlightLevel(fl::FlightLevel; feet::Bool=fl.feet, knots::Bool=fl.knots, celsius::Bool=fl.celsius)

Return a copy of fl transforming the length, speed and temperature magnitudes to its corresponding units
"""
function FlightLevel(fl::FlightLevel; feet::Bool=fl.feet, knots::Bool=fl.knots, celsius::Bool=fl.celsius)
    mach = fl.mach
    h = fl.h
    eas = fl.eas
    cas = fl.cas
    tas = fl.tas
    a = fl.a
    P = fl.P
    T = fl.T
    ρ = fl.ρ
    μ = fl.μ
    ν = fl.ν
    PDyn = fl.PDyn
    Re_x = fl.Re_x
    ΔT = fl.ΔT

    if feet & !fl.feet
        h = m2ft(h)
        Re_x /= m2ft(1)
    elseif !feet & fl.feet
        h = ft2m(h)
        Re_x /= ft2m(1)
    end
    if knots & !fl.knots
        eas, cas, tas, a = ms2knots.([eas, cas, tas, a])
    elseif !knots & fl.knots
        eas, cas, tas, a = knots2ms.([eas, cas, tas, a])
    end
    if celsius & !fl.celsius
        T = kel2cel(T)
    elseif !celsius & fl.celsius
        T = cel2kel(T)
    end
    # TODO: Finish with conversions for pressure and (kinematic and dynamic) viscosity

    return FlightLevel(mach, h, eas, cas, tas, a, P, T, ρ, μ, ν, PDyn, Re_x, ΔT, feet, knots, celsius)
end