/**
   $(BIG $(LINK2 creation.html, numir.core.creation))
   $(BR)
   $(BIG $(LINK2 manipulation.html, numir.core.manipulation))
   $(BR)
   $(BIG $(LINK2 utility.html, numir.core.utility))
   $(BR)
   $(BIG $(LINK2 random.html, numir.random))
   $(BR)
   $(BIG $(LINK2 stats.html, numir.stats))
   $(BR)
   $(BIG $(LINK2 io.html, numir.io))
   $(BR)
   $(BIG $(LINK2 format.html, numir.format))
   $(BR)
   $(BIG $(LINK2 testing.html, numir.testing))
 */

/**
   License: $(LINK2 http://boost.org/LICENSE_1_0.txt, Boost License 1.0).

   Authors:
       $(LINK2 https://github.com/ShigekiKarita, Shigeki Karita),
       $(LINK2 https://github.com/jmh530, John Michael Hall),
       $(LINK2 https://github.com/9il, Ilya Yaroshenko)

   See_Also:
       https://github.com/libmir/numir
       $(BR)
       $(LINK2 http://docs.algorithm.dlang.io/latest/index.html, mir.algorithm)
       $(BR)
       $(LINK2 http://docs.random.dlang.io/latest/index.html, mir.random)
*/

module numir;

/**
  this library is motivated by
  https://github.com/torch/torch7/wiki/Torch-for-Numpy-users
*/


public import numir.core;
public import numir.random;
public import numir.io;
public import numir.testing;
public import numir.format;
public import numir.stats;
