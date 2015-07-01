# jpeg2png
## Silky smooth JPEG decoding - no more artifacts!
JPEG encoding loses information. But it is JPEG decoding that introduces artifacts by filling the missing information with noise.

jpeg2png is smarter and fills the missing information to create the smoothest possible picture.

## [Examples](/../../tree/images)

![Lena](/../images/lena_tiles.png?raw=true)

* Top left: [original](../images/lena.png) Lena image
* Top right: original, 64x64 detail
* Bottom left: [JPEG encoded](../images/lena.jpg) at 10% quality with 4:2:0 chroma subsampling using the GIMP, 64x64 detail
* Bottom right: [JPEG decoded](../images/lena_restored.jpg) with jpeg2png using the default settings, 64x64 detail

![Hige](/../images/deviantart_tiles.png?raw=true)

* Top left: [original](../images/deviantart.png) Hige image ([source](http://whitedovehemlock.deviantart.com/art/Hige-315700935))
* Top right: original, 64x64 detail
* Bottom left: [JPEG encoded](../images/deviantart.jpg) at 90% quality with 4:4:4 chroma subsampling using the GIMP, 64x64 detail
* Bottom right: [JPEG decoded](../images/deviantart_restored.png) with jpeg2png using the default settings, 64x64 detail

## Installation
A compiled Windows version is available on the ["Releases" page](../../releases).

jpeg2png is written in portable C, specifically C11. It relies on libjpeg and libpng.
You can compile the executable using GNU make and either GCC or Clang. Execute one of

    make
    CC=clang make

You may use ``./jpeg2png`` to execute it without installing, or install using

    sudo make install

jpeg2png is licensed GPLv3+.

## Usage
Just execute ``jpeg2png picture.jpg`` to create ``picture.png``. Execute ``jpeg2png --help`` to see all options.

Under Windows, you can also drag-and-drop JPEG files onto the program.

jpeg2png gives best results for pictures that should never be saved as JPEG.
Examples are charts, logo's, and cartoon-style digital drawings.

On the other hand, jpeg2png gives poor result for photographs or other finely textured pictures.

For pictures that have been reencoded to JPEG multiple times I recommend [shred](https://www.gnu.org/software/coreutils/manual/html_node/shred-invocation.html).

## What "smooth" means
jpeg2png finds the smoothest possible picture that encodes to the given JPEG file.
But what is "smooth"? A first approximation is [Total Variation](https://en.wikipedia.org/wiki/Total_variation_denoising).
It says that smoothness is the sum of the differences between neighboring pixels.
This works pretty well, but it can create very sharp transitions.
``0 1 2 3 4`` and ``0 0 0 4 4`` both have a Total Variation of 4.

What we need is to look not just at the differences, but also at the differences of the differences.
This called Total Generalized Variation[1]. The sum of differences of differences for
``0 1 2 3 4`` is 0, and for ``0 0 0 4 4`` it is 4.

To combine the first order sum of differences and the sum of second order differences we choose a weight w.
Then the total smoothness of ``0 1 2 3 4`` is ``4 + w * 0``, and of ``0 0 0 4 4`` is ``4 + w * 4``.

jpeg2png lets you choose the weight for the sum of second order differences with the ``-w`` parameter.

Optimizing purely for smoothness can sometimes go too far.
If we smooth ``0 3 3 .. 3 3 0`` with maximal deviation 0.5 the result is ``0.5 2.5 2.5 .. 2.5 2.5 0.5``.
The whole inner block gets changed to fit with the border in an improbable way.
In a real picture this can be seen as a slight change in brightness or color.
To prevent this we also optimize for the minimal sum of squared deviations of DCT coefficients, with a very small weight.

jpeg2png lets you choose the weight for the sum of squared deviations with the ``-p`` parameter.

## Finding the smoothest picture
Now we know what we are looking for. But how do we find it?
It turns out this is a non-linear convex optimization problem.
We use a method that gets to smoother decodings in steps.
The higher the number of steps, the more smooth the decoding.

jpeg2png lets you choose the number of steps with the ``-i`` parameter.

A low number of steps, like 10, will take only a few seconds.
The quality could be better, but compared to regular decoding it is already very good.

A high number of steps, like 1000, might take a few minutes.
The quality is very good, but such a high number is probably overkill.

## Nitty gritty
JPEG encoding goes something like
``convert colors to YCbCr -> chroma subsampling -> blockwise DCT -> quantization -> rounding -> entropy coding``

Standard JPEG decoding goes something like
``entropy decoding -> dequantization -> blockwise IDCT -> chroma upsampling -> convert colors to RGB``.

The crucial step that is missing in decoding is reversing the rounding.
Of course rounding is not one-to-one invertible, but we can unround x to the interval ``[x-0.5, x+0.5]``.
This gives us the set of possible pictures.

Our objective is to minimize ``sum i=1 to n (norm(gradient(u_i))) + w * sum i=1 to n (norm(gradient(gradient(u_i)))) + p * sum (DCT(u-original)/quant)^2)``.
To get the gradient for the TV term of the objective we use forward differences.
The norm is an Euclidean norm.
For the second order TGV term we use backward differences for the second gradient, giving us a 2x2 Hessian matrix.
We symmetrize the matrix, that is we average dxdy and dydx.
The norm here is a Frobenius norm.
We do not use any higher order TGV terms.
The deviations are normalized by the quantization factors.
We do not differentiate between deviations in the DC and AC coefficients.

Unfortunately if one of the norms is zero the gradient of our objective is not defined, so our objective is not smooth.
The subderivative chosen when a norm is zero is 0.

The objective is convex and our search space ``Q`` is convex too.
We can project onto our search space easily because DCT is orthogonal. So we can DCT, project the deviations onto the box ``[-0.5, +0.5]^n``, and IDCT back.

In conclusion, we use the subgradient method with projection and FISTA acceleration.
The step size chosen is ``radius(Q) / sqrt(1 + number of steps)``, where ``radius(Q)`` is ``sqrt(n) / 2``.

## Wishlist

* do more testing on different kinds of images
* make comparisons with known JPEG artifact reduction techniques
* ~~make it go faster~~
  * basically everything has SSE2 versions now, 2x speedup versus pure C
  * parallel (OpenMP)
    * almost linear speedup for multiple files
    * runs max 3x as fast with --separate-components
    * otherwise 1.1x speedup (2 cores)
    * not sure if it was worth the time in the end, but it made sense when --separate-components was the only mode
  * things that didn't work out
    * not boxing/unboxing: no performance difference (branch no_boxing)
    * using a z-order curve: no performance difference
    * turning TV inside out to make it parallel: 1.15x speedup in C version (2 cores) (branch parallel_tv)
      * makes the code very hard to read
      * too small of a difference to justify rewriting tgv and simd versions
  * ultimately, the progress bar was more important than improving time from 40 to 18 seconds
* ~~investigate optimizing all components together~~
  * implemented for TV, ~~don't know if it make sense for TGV~~
    * implemented for TGV
* ~~investigate better chroma upsampling~~
  * too late, too small to make a real difference
* ~~investigate automake / autoconf~~
  * too much work to learn, patches welcome
* ~~investigate smoothing methods~~
  * only accelerates the start, no improvement in the end
* ~~investigate other stop conditions than a fixed number of steps~~
  * no good criterion when using subgradient method
* ~~investigate dual methods, Bregman~~
  * too complicated and inflexible, primal-dual has a good stopping criterion but same complexity
* ~~support gray-scale, maybe other JPEG features~~
  * low interest, file an issue if you have a real-world use for this

## References with comments

[1] ["Total generalized variation" (2010) by Kristian Bredies, Karl Kunisch, Thomas Pock](http://gpu4vision.icg.tugraz.at/papers/2009/pock_tgv.pdf)

Maybe over-mathematical, but that's my impression of a lot of image papers.

[2] ["Introductory Lectures on Convex Programming, Volume I: Basic course" (1998) by Yurii Nestorov](http://enpub.fulton.asu.edu/cseml/Fall2008_ConvOpt/book/Intro-nl.pdf)

I believe this may be a draft of his 2004 book, but it's phenomenal nonetheless.
This is a solid foundation, well organized, with clear explanations and full proofs.
Strongly recommended.

[3] ["Adapted Total Variation for Artifact Free Decompression of JPEG Images." (2005) by François Alter, Sylvain Durand, Jacques Froment](http://www.mediafire.com/view/o9ya9gsdzyb0cwq/art10.1007s10851-005-6467-9.pdf)

TV objective, convex optimization with subgradient method with projection, in retrospect it's all quite obvious.
I found this after I figured it out but it's still a good read.

[4] ["Artifact-Free Decompression and Zooming of JPEG Compressed Images with Total Generalized Variation" (2013) by Kristian Bredies, Martin Holler](http://www.ma.tum.de/foswiki/pub/IGDK1754/ProceedingOther/BrediesHoller_2013.pdf)

More advanced than the above, considers subsampling, TGV, primal-dual algorithms. Promo material [[video](http://www.youtube.com/watch?v=GJG3B4X3eiQ)] [[presentation](http://www.uni-graz.at/~hollerm/presentations/presentation_tgv_jpeg.pdf)] [[TO](https://static.uni-graz.at/fileadmin/forschen/dokumente/technologietransfer/TO_JPEG_TGV.pdf)]

[5] ["DCT Quantization Noise in Compressed Images" by Mark A. Robertson and Robert L. Stevenson](https://www3.nd.edu/~lisa/mrobert2/csvt2001submit.pdf)

This is the source for the DCT deviations model. Uniform distribution for the errors is close enough, even if a generalized normal distribution with beta = 1/2 is more realistic. I ignore the HMRF stuff because I find the Huber function very ad hoc and inelegant.

## Links

[qjpegrest](http://viric.name/soft/qjpegrest/) is a tool that lets you try many different JPEG restoration methods (TV based, band-pass based and Huber MRF based). I learned from the code. See notes/qjpegrest.txt for installation help.

## License

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
