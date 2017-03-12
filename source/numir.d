/**
   License: $(LINK2 http://boost.org/LICENSE_1_0.txt, Boost License 1.0).

   Authors: $(LINK2 http://shigekikarita.github.io, karita)
*/

module numir;

/*
  this library is motivated by
  https://github.com/torch/torch7/wiki/Torch-for-Numpy-users
*/


import mir.ndslice;

import std.array : array;
import std.conv : to;
import std.format : format;
import std.meta : staticMap;
import std.range : ElementType, isInputRange;
import std.stdio : writeln;
import std.traits : CommonType, isArray, isFloatingPoint;


static if (__VERSION__ < 2073)
{
    import old : maxIndex; // not supported yet (2.071)
}
else
{
    import std.algorithm.searching: maxIndex;
}


///
unittest
{
    /*
      Types

      np.ndarray | mir
      -----------+--------
      np.float32 | float
      np.float64 | double
      np.int8    | byte
      np.uint8   | ubyte
      np.int16   | short
      np.int32   | int
      np.int64   | long
      np.astype  | as

      see also https://dlang.org/spec/type.html
    */
}


/++
 construct new uninitialized slice of an element type `E` and shape(`length ...`)

 Params:
 length = elements of shape
 Returns:
 new uninitialized slice
 +/
auto empty(E=double, size_t N)(size_t[N] length...)
{
    return uninitializedSlice!E(length);
}

/++
 construct new slice having the same element type and shape to given slice

 Params:
 initializer = template function(ElementType)(shape) that initializes slice
 slice = source slice to refer shape and element type
 Returns:
 new uninitialized slice
 +/
auto like(alias initializer, S)(S slice)
{
    return initializer!(DeepElementType!S)(slice.shape);
}

///
auto empty_like(S)(S slice)
{
    return slice.like!empty;
}

///
auto ones(E=double, size_t N)(size_t[N] length...)
{
    return slice!E(length, 1);
}

///
auto ones_like(S)(S slice)
{
    return slice.like!ones;
}

///
auto zeros(E=double, size_t N)(size_t[N] length...)
{
    return slice!E(length, 0);
}

///
auto zeros_like(S)(S slice)
{
    return slice.like!zeros;
}

///
auto eye(E=double)(size_t m, size_t n=0, long k=0)
{
    if (n == 0) n = m;
    auto z = zeros!E(m, n);
    z.diag(k)[] = 1;
    return z;
}

///
auto identity(E=double)(size_t n) {
    return eye!E(n);
}

///
unittest
{
    /* Constructors

       numpy            | numir
       -----------------+---------------------------------------
       np.empty([2,2])  | numir.empty(2, 2), numir.empty([2, 2])
       np.empty_like(x) | numir.empty_like(x)
       np.eye           | numir.eye
       np.identity      | numir.identity
       np.ones          | numir.ones
       np.ones_like     | numir.ones_like
       np.zeros         | numir.zeros
       np.zeros_like    | numir.zeros_like

       see also http://mir.dlang.io/mir_ndslice_allocation.html
    */

    // np.empty, empty_like
    assert(empty(2 ,3).shape == [2, 3]);
    assert(empty([2 ,3]).shape == [2, 3]);
    auto e0 = empty!int(2, 3);
    auto e1 = empty_like(e0);
    assert(e1.shape == e0.shape);
    e0[0, 0] += 1;
    assert(e1 != e0);
    alias E0 = DeepElementType!(typeof(e0));
    alias E1 = DeepElementType!(typeof(e1));
    static assert(is(E0 == E1));

    // np.ones, ones_like
    auto o = ones(2, 3);
    assert(o.all!(x => x == 1));
    assert(o.shape == [2, 3]);
    assert(o == ones([2, 3]));
    assert(o == o.ones_like);

    // np.zeros, np.zeros_like
    auto z = zeros(2, 3);
    assert(z.all!(x => x == 0));
    assert(z.shape == [2, 3]);
    assert(z == zeros([2, 3]));
    assert(z == z.zeros_like);

    // np.eye, identity
    assert(eye(2, 3, 1) == [[0.0, 1.0, 0.0],
                            [0.0, 0.0, 1.0]]);
    assert(identity(2) == [[1.0, 0.0],
                           [0.0, 1.0]]);
}


///
template rank(R)
{
    static if (isInputRange!R || isArray!R)
    {
        enum size_t rank = 1 + rank!(ElementType!R);
    }
    else
    {
        enum size_t rank = 0;
    }
}

///
template NestedElementType(T)
{
    static if (isArray!T)
    {
        alias NestedElementType = NestedElementType!(ElementType!T);
    }
    else
    {
        alias NestedElementType = T;
    }
}

///
size_t[rank!T] shapeNested(T)(T array)
{
    static if (rank!T == 0)
    {
        return [];
    }
    else
    {
        return to!(size_t[rank!T])(array.length ~ shapeNested(array[0]));
    }
}

unittest
{
    int[2][3] nested = [[1,2],[3,4],[5,6]];
    assert(nested.shapeNested == [3, 2]);
    assert([1].shapeNested == [1]);
    assert([1, 2].shapeNested == [2]);
    assert([[1,2],[3,4],[5,6]].shapeNested == [3, 2]);
    static assert(is(NestedElementType!(int[][]) == int));
}

///
auto nparray(E=void, T)(T a)
{
    static if (is(E == void))
    {
        alias E = NestedElementType!T;
    }
    auto m = slice!E(a.shapeNested);
    m[] = a;
    return m;
}

///
auto concatenate(int axis=0, Slices...)(Slices slices)
{
    enum int N = Ndim!(Slices[0]);
    static assert(-N <= axis, "out of bounds: axis(=%s) < %s".format(axis, -N));
    static assert(axis < N, "out of bounds: %s <= axis(=%s)".format(N, axis));
    static if (axis < 0) {
        enum axis = axis + N;
    }

    foreach (S; Slices) {
        static assert(Ndim!S == N,
                      "all the input arrays must have same number of dimensions: %s"
                      .format([staticMap!(Ndim, Slices)]));
    }

    import mir.ndslice.concatenation: concatenation;
    return concatenation!axis(slices).slice;
}

///
unittest
{
    /* From existing data

       numpy                     | numir
       --------------------------+------------------------
       np.array([ [1,2],[3,4] ]) | nparray([ [1,2],[3,4] ])
       np.ascontiguousarray(x)   | x.assumeContiguous
       np.copy(x)                | ????
       np.fromfile(file)         | ????
       np.concatenate            | concatenate
    */

    auto s = [[1,2],[3,4]].sliced; // mir's sliced
    // error: s[0, 0] = -1;

    auto m = nparray([[1,2],[3,4]]);
    m[0, 0] = -1;
    assert(m == [[-1,2], [3,4]]);
    static assert(is(DeepElementType!(typeof(m)) == int)); // maybe double?

    auto v = nparray([1, 2]);
    v[1] = -2;
    assert(v == [1, -2]);

    auto u = nparray([[5, 6]]);
    assert(concatenate(m, u) == [[-1, 2], [3, 4], [5, 6]]);
    assert(concatenate(u, m) == [[5, 6], [-1, 2], [3, 4]]);

    auto uT = u.universal.transposed;
    assert(concatenate!1(m, uT) == [[-1, 2, 5], [3, 4, 6]]);

    assert(concatenate!0([[0,1]].nparray, [[2,3]].nparray, [[4,5]].nparray) == iota(3, 2));
    // axis=-1 is the same to axis=$-1
    assert(concatenate!(-1)([[0,1]].nparray, [[2,3]].nparray, [[4,5]].nparray) == [iota(6)]);
    assert(concatenate!(-1)([[0,1]].nparray, [[2]].nparray) == [[0, 1, 2]]);
}



///
auto arange(size_t size)
{
    return size.iota;
}

///
auto linspace(E=double)(E start, E stop, size_t num=50)
{
    static if (!isFloatingPoint!E) {
        alias E = double;
    }
    return mir.ndslice.linspace([num].to!(size_t[1]), [[start, stop]].to!(E[2][1]));
}


version (DigitalMars) // FIXME: LDC fails
{
    ///
    auto steppedIota(E)(size_t num, E step, E start=0)
    {
        return iota(num).map!(i => E(i * step + start));
    }

    ///
    auto arange(E)(E start, E end, E step=1)
    {
        size_t num = to!size_t((end - start) / step) + 1;
        return num.steppedIota!E(step, start);
    }

    ///
    auto logspace(E=double)(E start, E stop, size_t num=50, E base=10)
    {
        return linspace(start, stop, num).map!(x => base ^^ x);
    }

    ///
    unittest
    {
        /* Numerical Ranges

           numpy                | numir
           ---------------------+--------
           np.arange(10)        | arange(10)
           np.arange(2, 3, 0.1) | arange(2, 3, 0.1)
           np.linspace(1, 4, 6) | linspace(1, 4, 6)
           np.logspace          | logspace

           see also: http://mir.dlang.io/mir_ndslice_topology.html#.iota
        */
        assert(arange(3) == [0, 1, 2]);
        assert(arange(2, 3, 0.3) == [2.0, 2.3, 2.6, 2.9]);
        assert(linspace(1, 2, 3) == [1.0, 1.5, 2.0]);
        assert(logspace(1, 2, 3, 10) == [10. ^^ 1.0, 10. ^^ 1.5, 10. ^^ 2.0]);
    }
}


/++ return diagonal slice +/
auto diag(S)(S s, long k=0)
{
    auto sk = k >= 0 ?  s[0 .. $, k .. $] : s[-k .. $, 0 .. $];
    return sk.diagonal;
}

///
unittest
{
    /* Building Matrices

       numpy    | numir
       ---------+---------
       np.diag  | diagonal
       np.tril  | <WIP>
       np.triu  | <WIP>
    */

    //  -------
    // | 0 1 2 |
    // | 3 4 5 |
    //  -------
    auto a = iota(2, 3);
    assert(a.diag == [0, 4]);
    assert(a.diag(1) == [1, 5]);
    assert(a.diag(-1) == [3]);
}

/// return
auto dtype(S)(S s)
{
    return typeid(DeepElementType!S);
}

///
template Ndim(S)
{
    enum Ndim = ndim(S());
}

///
size_t ndim(SliceKind kind, size_t[] packs, Iterator)(Slice!(kind, packs, Iterator) s)
{
    import mir.ndslice.internal: sum;
    return packs.sum;
}

/// return strides of byte size
size_t[] byteStrides(S)(S s)
{
    enum b = DeepElementType!S.sizeof;
    return s.strides.sliced.map!(n => n * b).array;
}

/// return size of raveled array
auto size(S)(S s)
{
    return s.elementsCount;
}

///
unittest
{
    /* Attributes

       numpy     | numir
       ----------+---------------
       x.shape   | x.shape
       x.strides | x.byteStrides (!) already <strides> function is defined in mir
       x.ndim    | x.ndim
       x.data    | ???
       x.size    | x.size
       len(x)    | x.length
       x.dtype   | x.dtype
    */

    auto e = empty!double(2, 3, 1, 3);
    assert(e.dtype == typeid(double));
    assert(e.dtype != typeid(float));
    assert(e.length == 2);
    assert(e.size == 2*3*1*3);
    assert(e.ndim == 4);
    assert(e.strides == [9, 3, 3, 1]);
    assert(e.byteStrides == [72, 24, 24, 8]);

    auto a = iota(3, 4, 5, 6);
    auto b = a.pack!2;
    assert(b.ndim == 4);
}


/// Shape Manipulation
unittest
{

}

/// Item selection and manipulation

/// Calculation

/// Arithmetic and comparison operations
