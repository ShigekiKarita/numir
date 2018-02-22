/++
Audio separation using STFT and non-negative matrix factorization

See_Also: https://github.com/r9y9/julia-nmf-ss-toy
+/

import std.stdio;
import std.file : exists;
import std.net.curl : download;
import std.complex : abs;
import std.typecons : tuple;

import mir.math.common : log;
import mir.ndslice : sliced, ndarray, map, transposed, reversed, slice, isSlice, maxPos;
import lubeck : mtimes;
import ggplotd.ggplotd : GGPlotD;
import ggplotd.axes : xaxisLabel, yaxisLabel;
import ggplotd.colour : colourGradient;
import ggplotd.colourspace : XYZ;

import numir;
import plot : plot1d, plot2d;
import signal : blackmanWindow, hanningWindow, stft;
import dffmpeg : Audio;

/++
Non-negative matrix factorization
 
Params:
    y = input mixed matrix with the shape of (ntime x nfreq)
    nbasis = number of basis vectors
    maxiter = the maximum number of iterations
    eps = the minimum value of error

Returns:
    matrix tuple [h, u] where
    
    argmin_{h, u} ||y - h \times u||_2

    h = shape (ntime x nbasis)
    u = shape (nbasis x nfreq)
 +/
auto nmf(S)(S y, size_t nbasis, size_t maxiter=100, double eps=1e-21) if (isSlice!S)
{
    auto h = uniform(y.length!0, nbasis).slice;
    auto u = uniform(nbasis, y.length!1).slice;
    foreach (i; 0 .. maxiter)
    {
        h[] *= y.mtimes(u.transposed) / (h.mtimes(u).mtimes(u.transposed) + eps);
        u[] *= h.transposed.mtimes(y) / (h.transposed.mtimes(h).mtimes(u) + eps);
        u[] /= u.maxPos.first;
    }
    return tuple!("h", "u")(h, u);
}

void main()
{
    // prepare audio
    auto filename = "test10k.wav";
    if (!filename.exists)
        download("https://raw.githubusercontent.com/ShigekiKarita/torch-nmf-ss-toy/master/test10k.wav",
                 filename);
    auto wav = Audio!short().load(filename);
    writeln(wav.now);

    // plot waveform
    auto xs = wav.data.sliced;
    GGPlotD().plot1d(xs).save("wav.png");

    // STFT
    auto zs = xs.stft(512)[0..$, 0..257]; // take real part
    auto ys = zs.map!abs.slice;

    // plot STFT result
    auto logy = ys.map!log.reversed!(1).transposed;
    GGPlotD().plot2d(logy)
        .put("white-cornflowerBlue-crimson".colourGradient!XYZ)
        .put("time".xaxisLabel)
        .put("freq".yaxisLabel)
        .save("spectogram.png");

    // plot STFT window
    GGPlotD().plot1d(blackmanWindow(1024), 0.0)
        .plot1d(hanningWindow(1024), 1.0)
        .put("time".xaxisLabel)
        .put("gain".yaxisLabel)
        .put("cornflowerBlue-crimson".colourGradient!XYZ)
        .save("windows.png");


    // NMF
    auto nbasis = 4;
    auto factorized = nmf(ys, nbasis);

    // plot NMF time/freq basis
    GGPlotD hfig, ufig;
    hfig = hfig.plot2d(logy * 0.5);
    ufig = ufig.plot2d(logy.transposed.reversed!(1, 0) * 0.5);
    auto lmax = logy.maxPos.first;
    auto hmax = cast(double) ys.front.length / factorized.h.maxPos.first;
    auto umax = cast(double) ys.length / factorized.u.maxPos.first;
    foreach (i; 0..nbasis)
    {
        auto color = lmax * (i + 1) / nbasis + lmax;
        hfig = hfig.plot1d(hmax * factorized.h[0..$, i], color);
        ufig = ufig.plot1d(umax * factorized.u[i, 0..$], color);
    }
    auto cg = "white-orange-green-cornflowerBlue-crimson";
    hfig.put(cg.colourGradient!XYZ)
        .put("time".xaxisLabel)
        .put("gain".yaxisLabel)
        .save("time_basis.png");
    ufig.put(cg.colourGradient!XYZ)
        .put("freq".xaxisLabel)
        .put("gain".yaxisLabel)
        .save("freq_basis.png");
}
