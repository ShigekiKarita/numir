/++
Signal processing package
 +/

module signal;

import mir.ndslice.topology : iota, map;

/++
Classic Blackman window slice

Params: n = length of window

Returns: window weight slice

See_Also: https://ccrma.stanford.edu/~jos/sasp/Blackman_Window_Family.html
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

See_Also: https://ccrma.stanford.edu/~jos/sasp/Generalized_Hamming_Window_Family.html
 +/
pure nothrow @nogc @safe
auto hanningWindow(size_t n, double a = 0.5, double b = 0.25)
{
    import mir.math.common : cos;
    import std.math : PI;
    immutable t = 2.0 * PI / (n - 1);
    return a - b / 2 * map!cos(t * iota(n));
}

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
Short time Fourie transform

Params:
     xs = input 1d sequence
     width = short-time frame width for each FFT
     stride = (default framelen / 2) short-time frame stride length for each FFT
 +/
auto stft(alias windowFun=hanningWindow, Xs)(Xs xs, size_t width=1024, size_t stride=0)
{
    import std.numeric : fft;
    import std.complex : Complex;
    import numir : empty;

    if (stride == 0) stride = width / 2; // default value
    auto frames = splitFrames(xs, width, stride);
    immutable nfreq = width; // for rfft: / 2 + 1;
    auto ret = empty!(Complex!double)(frames.length, nfreq);
    auto window = windowFun(width);
    foreach (i; 0 .. frames.length) {
        fft(frames[i] * window, ret[i]);
    }
    return ret;
}
