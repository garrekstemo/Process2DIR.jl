# Process 2DIR

A simple package the loads and processes 2DIR data from a
128 x 128 pixel mercury cadmium telluride (MCT) infrared detector.

File naming convention is as follows:

`<file_prefix>_<#######>fs_##.2DIR`

where `file_prefix` is the prefix of the file name, `#######` is the delay time between the two pump pulses in femtoseconds, and `##` is the spectrum number at
that delay time. The extension is `.2DIR`, but
the files are simple text files and this extension
is an optional argument and can be modified by the user.

The pixel number is a global constant, `PIXELS`, and is set to 128.

The result is stored in the type `Spectra`, 
which includes
- frequecy in wavenumbers
- pump-probe spectra (on, off, and difference)
- 2DIR spectra at each delay time

The frequency will need to be further calibrated and
the second axis frequencies will need to be calculated to
produce a 2D plot.

Results can be retrieved in three lines if
the calibration constants are known:

```
using Process2DIR

result = process_2dir(dir, fileprefix, f0_shift)
ω1 = calibrate_frequency(result.frequency, cal1, cal2, shift1)
ω3 = calibrate_frequency(result.frequency, cal3, cal4, shift2)
```

Then plot using `result.spectra[:, :, i]`, where `i` is the index for each 2DIR spectrum.