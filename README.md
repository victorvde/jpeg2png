# jpeg2png
## Silky smooth JPEG decoding - no more artifacts!
JPEG encoding loses information. But it is JPEG decoding that introduces artifacts by filling the missing information with noise.

jpeg2png is smarter and fills the missing information to create the smoothest possible picture. Just do

``$ jpeg2png input.jpg output.png``

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

## Availability

jpeg2png is written in portable C, specfically C11. It relies on libjpeg and libpng.
You can compile the executable using make and gcc.
jpeg2png is licensed GPLv3+.

## What "smooth" means
jpeg2png finds the smoothest possible picture that encodes to the given JPEG file.
But what is "smooth"? A first approximation is [Total Variation](https://en.wikipedia.org/wiki/Total_variation_denoising).
It says that smoothness is the sum of the differences between neighbouring pixels.
This works pretty well, but it can create very sharp transitions.
``0 1 2 3 4`` and ``0 0 0 4 4`` both have a Total Variation of 4.

What we need is to look not just at the differences, but also at the differences of the differences.
This called Total Generalized Variation[1]. The sum of differences of differences for
``0 1 2 3 4`` is 0, and for ``0 0 0 4 4`` it is 4.

To combine the first order sum of differences and the second order sum of differences we choose a weight w.
Then the total smoothness of ``0 1 2 3 4`` is ``4 + w * 0``, and of ``0 0 0 4 4`` is ``4 + w * 4``.

jpeg2png lets you choose the weight with the ``-w`` parameter.

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

Standand JPEG decoding goes something like
``entropy decoding -> unquantization -> blockwise IDCT -> chroma upsampling -> convert colors to RGB``.

The crucial step that is missing in decoding is reversing the rounding.
Of course rounding is not one-to-one invertible, but we can unround x to the interval ``[x-0.5, x+0.5]``.
This gives us the set of possible pictures.

Our objective is to minimize ``sum_i=1^n (norm(gradient(u_i))) + w * sum_i=1^n (norm(gradient(gradient(u_i))))``.
To get the gradient for the TV term of the objective we use forward differences.
The norm is an Euclidian norm.
For the second order TVG term we use backward differences for the second gradient, giving us a 2x2 Hessian matrix.
Currently we do not symmetrize the matrix.
The norm here is a Frobenius norm.
We do not use any higher order TVG terms.

The objective is convex and our search space ``Q`` is convex too.
Unfortunately if one of the norms is zero the gradient of our objective is not defined, so our objective is not smooth.
We can project onto our search space easily because DCT is orthogonal. So we can DCT, project onto the box ``[-0.5, +0.5]^n``, and IDCT back.

We use the subgradient method with projection and FISTA acceleration. The subderivative chosen when a norm is zero is 0.
The step size chosen is ``radius(Q) / sqrt(1 + number of steps)``, where ``radius(Q) is sqrt(n) / 2``.

## To Do

* do more testing on different kinds of images
* make comparisons with known JPEG artifact reduction techniques
* investigate automake / autoconf
* investigate smoothing methods
* investigate other stop conditions than a fixed number of steps
* investigate dual methods, Bregman
* optimize more
* investigate better chroma upsampling
* support gray-scale, maybe other JPEG features

## References with comments

[1] ["Total generalized variation" (2010) by Kristian Bredies, Karl Kunisch, Thomas Pock](http://gpu4vision.icg.tugraz.at/papers/2009/pock_tgv.pdf)

Maybe over-mathematical, but that's my impression of a lot of image papers.

[2] ["Introductory Lectures on Convex Programming, Volume I: Basic course" (1998) by Yurii Nestorov](http://enpub.fulton.asu.edu/cseml/Fall2008_ConvOpt/book/Intro-nl.pdf)

I believe this may be a draft of his 2004 book, but it's phenomonal nonetheless.
This is a solid foundation, well organized, with clear explanations and full proofs.
Strongly recommended.

[3] ["Adapted Total Variation for Artifact Free Decompression of JPEG Images." (2005) by François Alter, Sylvain Durand, Jacques Froment](http://www.mediafire.com/view/o9ya9gsdzyb0cwq/art10.1007s10851-005-6467-9.pdf)

TV objective, convex optimization with subgradient method with projection, in retrospect it's all quite obvious.
I found this after I figured it out but it's still a good read.

[4] ["Artifact-Free Decompression and Zooming of JPEG Compressed Images with Total Generalized Variation" (2013) by Kristian Bredies, Martin Holler](http://www.ma.tum.de/foswiki/pub/IGDK1754/ProceedingOther/BrediesHoller_2013.pdf)

More advanced than the above, considers subsampling, TGV, primal-dual algorithms. Promo material [[video](http://www.youtube.com/watch?v=GJG3B4X3eiQ)] [[presentation](http://www.uni-graz.at/~hollerm/presentations/presentation_tgv_jpeg.pdf)] [[TO](https://static.uni-graz.at/fileadmin/forschen/dokumente/technologietransfer/TO_JPEG_TGV.pdf)]

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
