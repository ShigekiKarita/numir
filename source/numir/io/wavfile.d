module numir.io.wavfile;

struct RiffChunk
{
    uint id; // BE
    uint size;
    uint format; // BE
}

struct FormatChunk
{
    uint id; // BE
    uint size;
    ushort audioFormat;
    ushort numChannels;
    uint sampleRate;
    uint byteRate;
    ushort blockAlign;
    ushort bitsPerSample;
}

struct DataChunk
{
    uint id; // BE
    uint size;
}

/// wav information header
struct WavInfo
{
    RiffChunk riff;
    FormatChunk format;
    DataChunk data;
}

/// endian conversion little -> big (wav default)
auto littleToBigEndian(T)(T val)
{
    import std.bitmanip : nativeToLittleEndian, bigEndianToNative;
    return val.nativeToLittleEndian.bigEndianToNative!T;
}


/// binary wav data struct
struct WavData
{
    WavInfo info;
    ulong[] data;
    double bias = 0.0, scale = 1.0;
    double[] normalized;

    void fixEndian()
    {
        this.info.riff.id = this.info.riff.id.littleToBigEndian;
        this.info.riff.format = this.info.riff.format.littleToBigEndian;
        this.info.format.id = this.info.format.id.littleToBigEndian;
        this.info.data.id = this.info.data.id.littleToBigEndian;
    }

    this(string path)
    {
        import std.array : array;
        import std.algorithm : min, max, map;
        import std.bitmanip : bigEndianToNative;
        import std.range : enumerate;
        import std.stdio : File, fread, writeln;
        import std.format : format;


        scope f = File(path, "rb");

        info = f.rawRead(new WavInfo[1])[0];
        fixEndian();
        if (info.format.numChannels != 1) {
            throw new Exception(
                "Not implemented: Wavinfo.format.numChannels: %d != 1"
                .format(info.format.numChannels));
        }

        const numBytes = info.format.bitsPerSample / 8;
        data = new ulong[info.data.size / numBytes];
        ubyte[ulong.sizeof] bs;
        auto maxVal = ulong.min, minVal = ulong.max;
        foreach (ref d; data)
        {
            bs[$-numBytes .. $] = f.rawRead(new ubyte[numBytes]); // raw[n .. n + numBytes];
            d = bigEndianToNative!ulong(bs);
            minVal = min(d, minVal);
            maxVal = max(d, maxVal);
        }

        bias = data[0];
        scale = max(maxVal - bias, bias - minVal);
        normalized = data.map!(x => (x - bias) / scale).array;
    }
}

/// read wav file from path to data array
auto read(string path)
{
    return WavData(path);
}
