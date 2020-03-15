# ECG Filter

An Octave/Matlab implementation of an electrocardiogram filter, tailored to filter a given electrocardiogram data, containing a 60 Hz powerline noise, as well as a few other unspecified noise. This code is written to fulfill a Biomedical Signal Processing (EB3102) assignment in Institut Teknologi Bandung. There is also a report written in Bahasa Indonesia detailing on my approach on solving this problem, along with the LaTeX code that was used to generate the report.

The code uses some built-in Octave functions, some functions from the signal package for Z-domain analysis, as well as some self-written functions used for:

* Designing a notch filter at a certain frequency
* Designing a comb filter at a certain frequency multiples
* Designing a nth order low-pass butterworth filter
* Designing a nth order high-pass butterworth filter
* Designing a nth order band-pass butterworth filter
