module Process2DIR

using DelimitedFiles
using StatsBase: mean
using FFTW

export Spectra,
       process_2dir,
       calibrate_frequency,
       PIXELS,
       cross_section

include("post_processing.jl")

struct Spectra
    frequency::Array{Float64}
    t2::Array{Int64}
    pump_off::Array{Float64, 2}
    pump_on::Array{Float64, 2}
    pump_probe::Array{Float64, 2}
    spectra::Array{Float64, 3}
end

const PIXELS::Int64 = 128  # number of pixels in a row of the MCT detector array

"""
    process_2dir(dir, prefix, f0=0.0; zeropad_multiple=8, extension=".2DIR",
    background=true, background_scale = 1, norm_scale=0)

Load all 2DIR data from a directory and process it,
returning the pump-off, pump-on, pump-probe spectra and 2DIR spectra (at each time step).

`f0` should be set appropriately with rotating frame.
It is zero without rotating frame (default).

`zeropad_multiple` determines the zero padding of the Fourier transform by
multiplying the number of time steps in the raw data.
"""
function process_2dir(dir, prefix, f0=0.0; zeropad_multiple=8, extension=".2DIR",
            background=true, background_scale = 1, norm_scale=0)
    
    files = readdir(dir)
    time = get_time(dir, files[1])
    t2 = get_t2(files)
    n_spectra = num_spectra(files)

    zero_padding = zeropad_multiple * length(time)  # simple Fourier transform zero padding
    frequency = calculate_frequency(time, zero_padding, f0)
    pump_off = zeros(length(t2), PIXELS)
    pump_on = zeros(length(t2), PIXELS)
    pump_probe = zeros(length(t2), PIXELS)
    spectra = zeros(zero_padding, PIXELS, length(t2))

    for i in eachindex(t2)
        pump_offs = zeros(n_spectra, PIXELS)
        pump_ons = zeros(n_spectra, PIXELS)
        pump_probes = zeros(n_spectra, PIXELS)
        spectras = zeros(length(time), PIXELS, n_spectra)

        for j in 1:n_spectra
            file = joinpath(dir, "$(prefix)_$(t2[i])fs_$(j - 1)$(extension)")

            if isfile(file)
                raw = readdlm(file)[:, 2:end]  # first column is the time axis, obtained above.
                pump_offs[j, :], pump_ons[j, :], pump_probes[j, :] = probe_spectra(raw)
                spectras[:, :, j] = fourframe(raw)
            end
        end

        # Retrieve pump-probe spectra and reshape Matrix to Vector
        pump_off[i, :] = vec(mean(pump_offs, dims=1))
        pump_on[i, :] = vec(mean(pump_ons, dims=1))
        pump_probe[i, :] = vec(mean(pump_probes, dims=1))

        spectras[1, :, :] ./= 2  # divide first row (time) elements of each spectrum by 2 for some reason
        transformed = fourier_transform(spectras, zero_padding)
        if i > 1
            if background == true
                spectra[:, :, i] = -(transformed - spectra[:, :, 1] .* background_scale)
            else
                spectra[:, :, i] = -transformed
            end
        elseif i == 1
            spectra[:, :, i] = -transformed
        end
    end

    # Normalize spectra by pump-off spectrum (with scale factor `norm_scale`)
    for (i, m) in enumerate(eachslice(spectra, dims=3))
        for (j, row) in enumerate(eachrow(m))
            if iszero(norm_scale)
                spectra[j, :, i] = spectra[j, :, i]
            else
                spectra[j, :, i] ./= (pump_off[1, :] ./ maximum(pump_off[1, :])) .* norm_scale
            end
        end
    end
    return Spectra(frequency, t2, pump_off, pump_on, pump_probe, spectra)
end

"""
    get_t2(files)

Return a reverse-sorted list of time values for t2 in femtoseconds by parsing the 
file names in the experiment directory.
"""
function get_t2(files; suffix="fs")
    t2s = Int[]
    for file in files
        if occursin(suffix, file)
            t2 = parse(Int, split(file, "_")[2][1:end-length(suffix)])
            if !(in(t2, t2s))
                push!(t2s, t2)
            end
        end
    end
    return sort(t2s, rev=true)  # The largest (positive) t2 is background 
end

"""
    num_spectra(files; extension=".2DIR")

Return the number of spectra at each time step by parsing the file names
in the experiment directory.
"""
function num_spectra(files; extension=".2DIR")
    numspectra = 0
    for file in files
        spectrum_number = parse(Int, chop(split(file, "_")[end], tail = length(extension)))
        if spectrum_number > numspectra
            numspectra = spectrum_number
        end
    end
    return numspectra + 1  # file names start at 0
end

"""
    get_time(file; extension=".2DIR")

Get the unique time steps (every 4) from the first column in a 2DIR file.
"""
function get_time(dir, file; extension=".2DIR")
    if occursin(extension, file)
        raw = readdlm(joinpath(dir, file))
        time = unique(raw[:, 1])
        return time
    else
        println("This is not a .2DIR file. Select a different file.")
    end
end

"""
    calculate_frequency(time, freq_length, f0 = 0.0)

Calculate the frequency domain axis along ω_1, convert to wavenumbers (cm⁻¹),
and shift by factor `f0`.
The length of the vector `freq_length` must be the same as
the number of rows in the Fourier transformed spectra.
"""
function calculate_frequency(time, freq_length, f0=0.0)
    frequency =  range(-1, 1, length = freq_length) ./ (2 * (time[2] - time[1]) * 1e-15)
    return f0 .+ frequency ./ 2.9997e10  # shift by f0 and convert to wavenumbers
end

"""
    fourframe(raw)

Phase-correct the 2DIR raw data.
"""
function fourframe(raw)
    return raw[1:4:end, :] - raw[2:4:end, :] + raw[3:4:end, :] - raw[4:4:end, :]
end

function fourier_transform(data, newlength)
    avg = dropdims(mean(data, dims=3), dims=3)
    zero_padded = zeros(newlength, PIXELS)
    zero_padded[1:size(avg)[1], 1:size(avg)[2]] = avg
    transformed = fftshift(fft(zero_padded, 1), 1)
    return real(transformed)
end

end # module Process2DIR
