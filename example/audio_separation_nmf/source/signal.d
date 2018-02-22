/++
Signal processing package
 +/

module signal;

import mir.ndslice.topology : iota, map;

/++
Classic Blackman window slice

Params: n = length of window

Returns: window weight slice

See_Also:
    https://ccrma.stanford.edu/~jos/sasp/Blackman_Window_Family.html
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.blackman.html
 +/
pure nothrow @nogc @safe
auto blackmanWindow(size_t n, double a0=0.42, double a1=0.5, double a2=0.08)
{
    import mir.math.common : cos;
    import std.math : PI;
    immutable t = 2.0 * PI / (n - 1);
    auto ks = iota(n);
    return a0 - a1 * map!cos(t * ks) + a2 * map!cos(t * 2 * ks);
}

/++
Hanning window slice

Params: n = length of window

Returns: window weight slice

See_Also:
    https://ccrma.stanford.edu/~jos/sasp/Generalized_Hamming_Window_Family.html
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.hann.html
 +/
pure nothrow @nogc @safe
auto hanningWindow(size_t n, double a = 0.5, double b = 0.25)
{
    import mir.math.common : cos;
    import std.math : PI;
    immutable t = 2.0 * PI / (n - 1);
    return a - b / 2 * map!cos(t * iota(n));
}

/++
Split (window and stride) time frames for FFT or convolutions

Params:
     xs = input slice
     width = length of each segment
     stride = the number of skipped frames between the head of split slices
Returns: slice of split (windowed and strided) slices
 +/
pure nothrow @safe @nogc
auto splitFrames(Xs)(Xs xs, size_t width, size_t stride)
{
    import mir.ndslice.topology : windows; // , s = stride;
    import mir.ndslice.dynamic : strided;
    immutable nframes = (xs.length - width) / stride + 1;
    return xs.windows(width).strided!0(stride);
}

///
pure nothrow @safe @nogc
unittest
{
    static immutable ys = [[0,1,2], [2,3,4]];
    assert(iota(6).splitFrames(3, 2) == ys);
}

/++
Check whether the Constant OverLap Add (COLA) constraint (STFT inversible) is met
+/
bool checkCOLA(W)(W window, size_t nperseg, size_t noverlap, double tol=1e-10)
in
{
    import numir.core : Ndim;
    assert(nperseg > 1);
    assert(nperseg > noverlap);
    static assert(Ndim!window == 1);
    assert(window.length == nperseg);
} 
do
{
    auto stride = nperseg - noverlap;
    auto binsums = iota(nperseg / stride).map!(i => window[iota([stride], i*stride)]);
}


/++
Short time Fourie transform

Params:
    xs = input 1d sequence
    nperseg = (default 256) short-time frame width for each FFT segment
    noverlap = (default nperseg / 2) short-time frame overlapped length for each FFT segment
See_Also:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.stft.html
 +/
auto stft(alias windowFun=hanningWindow, Xs)(Xs xs, size_t nperseg=256, size_t noverlap=0)
in 
{
    assert(noverlap < nperseg);
}
do 
{
    import std.numeric : fft;
    import std.complex : Complex;
    import numir : empty;

    if (noverlap == 0) noverlap = nperseg / 2; // default value
    auto frames = splitFrames(xs, nperseg, nperseg - noverlap);
    immutable nfreq = nperseg; // for rfft: / 2 + 1;
    auto ret = empty!(Complex!double)(frames.length, nfreq);
    auto window = windowFun(nperseg);
    foreach (i; 0 .. frames.length) {
        fft(frames[i] * window, ret[i]);
    }
    return ret;
}
