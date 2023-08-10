struct FlightLevel
    mach::Real
    h::Real
    eas::Real
    cas::Real
    tas::Real
    P::Real
    T::Real
    ρ::Real
    μ::Real
    ν::Real
    PDyn::Real
    Re_x::Real
    ΔT::Real
    feet::Bool
    knots::Bool
    celsius::Bool
end

function FlightLevel(fl::FlightLevel, feet::Bool, knots::Bool, celsius::Bool)
    mach = fl.mach
    h = fl.h
    eas = fl.eas
    cas = fl.cas
    tas = fl.tas
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
        eas = ms2knots(eas)
        cas = ms2knots(cas)
        tas = ms2knots(tas)
    elseif !knots & fl.knots
        eas = knots2ms(eas)
        cas = knots2ms(cas)
        tas = knots2ms(tas)
    end
    if celsius & !fl.celsius
        T = kel2cel(T)
    elseif !celsius & fl.celsius
        T = cel2kel(T)
    end
    # TODO: Finish with conversions for pressure and (kinematic and dynamic) viscosity

    return FlightLevel(mach, h, eas, cas, tas, P, T, ρ, μ, ν, PDyn, Re_x, ΔT, feet, knots, celsius)
end