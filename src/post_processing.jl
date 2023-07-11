"""
    calibrate_frequency(ω, slope, offset)

Calibrate the frequency axis by a linear fit.
"""
function calibrate_frequency(ω, slope, offset)
    ω_calibrated = zeros(length(ω))
    for i in eachindex(ω)
        ω_calibrated[i] = slope * ω[i] + offset
    end
    return ω_calibrated
end

"""
    cross_section(slice, ω1, ω3; integrate = :horizontal, ω1_limits=nothing, ω3_limits=nothing)

Integrate a cross section of the 2DIR spectrum,
either vertically or horizontally.
"""
function cross_section(slice, ω1, ω3; integrate = :horizontal, ω1_limits=nothing, ω3_limits=nothing)
    if ω3_limits === nothing
        ω3_limits = (ω3[end], ω3[1])
    end
    if ω1_limits === nothing
        ω1_limits = (ω1[1], ω1[end])
    end

    ω3_idx = findall(x -> ω3_limits[1] <= x <= ω3_limits[2], ω3)
    ω1_idx = findall(x -> ω1_limits[1] <= x <= ω1_limits[2], ω1)
    to_integrate = slice[ω3_idx, ω1_idx]

    if integrate == :horizontal
        slice = vec(sum(to_integrate, dims = 1)) ./ size(to_integrate, 1)
    elseif integrate == :vertical
        slice = vec(sum(to_integrate, dims = 2)) ./ size(to_integrate, 2)
    else
        error("integrate must be :horizontal or :vertical")
    end
    return ω1_idx, ω3_idx, slice
end

"""
    slopeoffset(x1, x2, x1_i, x2_i)

Scale and shift the frequency axes (ω1 and ω3) by a linear function.
"""
function slopeoffset(x1, x2, x1_i, x2_i)
    slope = (x2_i - x1_i) / (x2 - x1)
    offset = (-x2_i * x1 + x1_i * x2) / (x2 - x1)
    return slope, offset
end


"""
    pumpprobe(raw)

Get the pump-off, pump-on, and pump-probe spectra directly from raw data.
"""
function probe_spectra(raw)
    pump_off = raw[2, :] + raw[4, :]
    pump_on = raw[1, :] + raw[3, :]
    pump_probe = raw[1, :] + raw[3, :] - raw[2, :] - raw[4, :]
    return pump_off, pump_on, pump_probe
end