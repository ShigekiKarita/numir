import numir;
import std.stdio;
import std.net.curl : get;

import numir;
import wavfile = numir.io.wavfile;

import plot : plot1d, plot2d;
import signal : blackmanWindow, hanningWindow, stft;
import numir.io.audio : loadAudio;


void main()
{
    // auto wav = get("https://raw.githubusercontent.com/ShigekiKarita/torch-nmf-ss-toy/master/test10k.wav");
    import std.array : array;
    import mir.ndslice : sliced, ndarray, map, transposed, reversed, slice;
    import ggplotd.ggplotd : GGPlotD, putIn;
    import ggplotd.axes : xaxisLabel, yaxisLabel;
    import ggplotd.colour : colourGradient;
    import ggplotd.colourspace : XYZ;

    plot1d(blackmanWindow(1024))
        .plot1d(hanningWindow(1024))
        .put("time".xaxisLabel)
        .put("gain".yaxisLabel)
        .save("windows.png");

    auto wav = loadAudio("test10k.wav", 10000);
    writeln(wav.data.length);
    wav.data.plot1d.save("wav.png");

    import std.complex : abs;
    import mir.math.common : log;

    // STFT
    auto xs = wav.data.sliced;
    auto zs = xs.stft(512)[0..$, 0..257]; // take real part
    auto ys = zs.map!abs.slice;
    ys.map!log.reversed!(1).transposed.plot2d
        .put("white-cornflowerBlue-crimson".colourGradient!XYZ)
        .save("spectogram.png");
}
