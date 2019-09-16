# [A Fast Sinc Function Gridding Algorithm or Fourier Inversion in Computer Tomography (1985)](https://ieeexplore.ieee.org/abstract/document/4307723)

## Summary

### Gridding method

We consider a 2D function in real space $m(x,y)$ whose Fourier transform $M(u,v)$ is given by
\begin{equation}
M(u,v) = \int_{-\infty}^{+\infty} \int_{-\infty}^{+\infty} m(x,y) e^{2\pi i (ux+vy)} dx dy.
\end{equation}
In the tomography case, $M(u,v)$ is the Fourier transform of the image we want to reconstruct.

We also define a sampling function $S(u,v)$ in Fourier space, which is given by
\begin{equation}
S(u,v) = \sum_{j=i}^{k} \delta(u-u_j, v-v_j).
\end{equation}
The 2D delta function in the equation above \emph{samples} a function at points $(u_j, v_j)$, which form any grid (i.e.~not necessarily Cartesian) in Fourier space. If we multiply this function to $M(u,v)$, we get a sampled Fourier transform
\begin{equation}
M_S(u,v) = M(u,v) \cdot S(u,v).
\end{equation}
The fist step of the gridding method is to convolve the sampled Fourier transform with a suitable kernel $C(u,v)$:
\begin{equation}
M_{SC}(u,v) = [M(u,v) \cdot S(u,v)] \ast C(u,v),
\end{equation}
where the subscripts in $M_{SC}(u,v)$ denote the operations sampling and convolution, respectively. 

Next, the convolution shown above is sampled at Cartesian points
\begin{equation}
M_{SCS}(u,v) = \{[M(u,v) \cdot S(u,v)] \ast C(u,v)\} \cdot \Sha(u,v),
\end{equation}
where $\Sha(u,v)$ denotes the Shah or comb function
\begin{equation}
\Sha(u,v) = \sum_{k} \sum_{j} \delta(u-k,v-j),
\end{equation}
which samples points $(k,j)$ that are equally spaced.

The reconstructed image $m_{SCS}(x,y)$ is given by the 2D inverse Fourier transform of $M_{SCS}(u,v)$. Therefore,
\begin{equation}
m_{SCS}(x,y) = \{[m(x,y) \ast s(x,y)] \cdot c(x,y)\} \ast \Sha(x,y)
\end{equation}

If the Fourier transform $M(u,v)$ is not sufficiently sampled by $S(u,v)$, we cannot correct the aliasing of $m(x,y)$ via post-processing. However, we can add an additional correction for non-uniform sampling in Fourier space. For this, we consider an area density function given by
\begin{equation}
\rho(u,v) = S(u,v) \ast C(u,v).
\end{equation}
Areas that are oversampled have a large area density while those that are undersampled have a small area density. Including this contribution in (5), we get a sampled and weighted function
\begin{equation}
M_{SWCS}(u,v) = \Bigg\{\Bigg[M(u,v) \cdot \frac{S(u,v)}{S(u,v) \ast C(u,v)}\Bigg] \ast C(u,v)\Bigg\} \cdot \Sha(u,v)
\end{equation}
and its corresponding inverse Fourier transform
\begin{equation}
m_{SWCS}(x,y) =  \{[m(x,y) \ast [s(x,y) \ast^{-1} (s(x,y) \cdot c(x,y))]] \cdot c(x,y)\} \ast \Sha(x,y),
\end{equation}
where $\ast^{-1}$ denotes deconvolution.

The question that presents itself at this stage is how to choose a kernel $C(u,v)$. This is discussed in the following section.

### Ideal interpolation
We note that for an ideal interpolation, the inverse Fourier transform of our sampled, weighted and convolved function should be equal to the objective function at the sampled points. That is,
\begin{equation}
\Bigg(m_{SWC}(x,y)\Bigg)_S = \Bigg(m(x,y)\Bigg)_S.
\end{equation}
Thus, what we want is to be able to recover a \emph{continuous} function $m(x,y)$ from the \emph{discrete sequence} of values at the sampled points by doing a convolution.
Suppose that our Fourier grid sampling is done by a function of the form $\Sha(u/u_0,v/v_0)$, such that the grid spacing in Fourier domain is given by $(1/u_0, 1/v_0)$. The inverse Fourier transform of this function is also a comb function and has the form
\begin{equation}
\Sha(x,y) = u_0 v_0 \Sha(xu_0, yv_0).
\end{equation}
The sampled Fourier transform of our function is given by
\begin{equation}
M_{S}(u,v) = M(u,v) \cdot \Sha(u/u_0, v/v_0).
\end{equation}
Its inverse Fourier transform is given by
\begin{align}
m_{S}(x,y) &= m(x,y) \ast u_0v_0\Sha(xu_0, yv_0)\\
&= m(x,y) \ast u_0v_0\sum_j \sum_k \delta(xu_0-j, yv_0-k)\\
&=\sum_j \sum_k m\Big(x - \frac{k}{u_0}, y - \frac{j}{v_0}\Big).
\end{align}
Additional weights can be added as constants multiplied to the sampling function. Here, we consider all weights to be equal to 1.

The function we want to reconstruct, $m(x,y)$, has compact support. We define its support such that $m(x,y) = 0$ for $x>|x_0|$ or $y>|y_0|$. We make use of the Whittaker-Shannon interpolation formula, which is a method to reconstruct a continuous time, bandlimited function from a sequence of real numbers. Given $x[n] = x(nT)$, a sequence of real numbers with sampling period $T$, the formula states that 
\begin{equation}
x(t) = \sum_n x[n] \text{ sinc}\Bigg( \frac{t-nT}{T}\Bigg)
\end{equation}
is a perfect reconstruction of the original function $x(t)$ if the bandlimit of the function, $B$, is less than the Nyquist frequency, $\frac{1}{2T}$. The reconstruction is equivalently given by 
\begin{equation}
x(t) = \Bigg(\sum_n x[n] \cdot \delta(t-nT)\Bigg) \ast \text{ sinc}\Bigg( \frac{t}{T}\Bigg).
\end{equation}
The function we want to reconstruct, $m(x,y)$, has compact support in real space and its limits are given by $(x \leq x_0, y \leq y_0)$. Therefore, we can reconstruct it perfectly by multiplying its sampled Fourier transform with a sinc function, as shown in (17) and (18). To avoid aliasing, the sampling period in Fourier domain must obey the following
\begin{align}
\frac{1}{2u_0} > x_0 \\
\frac{1}{2v_0} > y_0
\end{align}
Following (18), we have
\begin{equation}
M_{SC}(u,v) =\Bigg( \frac{1}{u_0v_0} \sum_k \sum_j M(u,v) \cdot \delta(u-ku_0, v-jv_0)\Bigg) \ast \text{sinc}\Bigg(\frac{u}{u_0}\Bigg)\text{sinc}\Bigg(\frac{v}{v_0}\Bigg),
\end{equation}
where the term inside the first parentheses is $M_{S}(u,v)$. Thus, we get a form for the kernel $C(u,v)$ in Fourier space:
\begin{equation}
C(u,v) = \frac{1}{u_0v_0}\text{sinc}\Bigg(\frac{u}{u_0}\Bigg)\text{sinc}\Bigg(\frac{v}{v_0}\Bigg).
\end{equation}
Taking $x_1 = 1/2u_0$ and $y_1 = 1/2v_0$, we can rewrite this as
\begin{equation}
C(u,v) = 4x_1 y_1 \text{sinc}(2x_1 u) \text{sinc}(2y_1 v).
\end{equation}
$M_{SC}(u,v)$ is a perfect reconstruction of the continuous function $M(u,v)$ if
\begin{align}
x_0 < x_1 \\
y_0 < y_1.
\end{align}
In real space, the sinc function in (23) is an ideal brick-wall low pass filter $c(x,y)$, such that
\begin{align}
c(x,y) &= 1, |x| \leq x_1 \text{ \& } |y| \leq y_1 \\
&= 0, \text{otherwise}.
\end{align}
Then, 
\begin{equation}
m_{SC}(x,y) = m_{S}(x,y) \cdot c(x,y) = m(x,y)
\end{equation}
provided that
\begin{align}
x_0 \leq x_1 < \frac{1}{u_0} - x_0 \\
y_0 \leq y_1 < \frac{1}{v_0} - y_0
\end{align}
The interpolation is completed by resampling the continuous function onto a grid (in our case, Cartesian) with spacing $(u_1, v_1)$
\begin{equation}
m_{SCS}(x,y)=\sum_j \sum_k m_{SC}\Big(x - \frac{k}{u_1}, y - \frac{j}{v_1}\Big),
\end{equation}
provided \begin{align}
x_0 < \frac{1}{2u_1}\\
y_0 < \frac{1}{2v_1}.
\end{align}

Thus, \emph{exact} interpolation to the Cartesian grid is possible. Note that in the case of the sinc kernel, the contribution of the kernel does not have to be removed after the inverse Fourier transform (a process known as deapodisation) because this kernel is chosen such that the sampled Fourier transform is \emph{exactly equal} to the Fourier transform of the objective function on the Cartesian grid.

\subsection{Finite kernel}
The sinc kernel has infinite extent and therefore convolution with the sinc kernel is computationally infeasible. Now, we consider the case that the kernel $C(u,v)$ has finite extent, which implies that its inverse Fourier transform $c(x,y)$ no longer has compact support. We also relax the condition that $m(x,y)$ has compact support (because of the finite sampling of $M_S(u,v)$ due to measurement). Thus, the function
\begin{equation}
m_{SC}(x,y) = m_S(x,y) \cdot c(x,y)
\end{equation}
is of infinite extent and we take into account only a truncated version of this function for $|x|<x_1, |y|<y_1$. In general, the sampling of this function
\begin{equation}
m_{SCS}(x,y) = \sum_j \sum_k m_{SC}\Bigg(x - \frac{k}{u_1}, y - \frac{j}{v_1} \Bigg)
\end{equation}
is not an accurate representation of $m_{SC}(x,y)$ and is affected by aliasing errors. In general, the gridding method for functions of non-compact support does not correspond to an exact interpolation and the values of the sampled function do not match the values of the objective function exactly, even at the measurement points. However, the objective function can be estimated by deapodisation:
\begin{equation}
m_e(x,y) = \frac{m_{SCS}(x,y)}{c(x,y)}, |x|<x_1 \text{ \& } |y|<y_1,
\end{equation}
where $m_e(x,y)$ is the estimate of the objective function $m(x,y)$. The central point ($k=0, j=0$) is computed exactly by deapodisation. At the other points, the error in the estimated function is given by
\begin{align}
m_e(x,y) - m_{S}(x,y) &= \frac{m_{SCS}(x,y)}{c(x,y)} - m_{S}(x,y) \\
&= \frac{1}{c(x,y)}\sum_j \sum_k m_{SC}\Bigg(x - \frac{k}{u_1}, y - \frac{j}{v_1} \Bigg) - m_S(x,y)\\
&= \frac{1}{c(x,y)}\sum_j \sum_k m_{S}\Bigg(x - \frac{k}{u_1}, y - \frac{j}{v_1} \Bigg) \cdot c\Bigg( x - \frac{k}{u_1}, y - \frac{j}{v_1}\Bigg) - m_S(x,y) \\
&= \frac{1}{c(x,y)}\sum_{j\neq0} \sum_{k\neq0} m_{S}\Bigg(x - \frac{k}{u_1}, y - \frac{j}{v_1} \Bigg) \cdot c\Bigg( x - \frac{k}{u_1}, y - \frac{j}{v_1}\Bigg)
\end{align}

The accuracy of the gridding method depends on the rate of decay of $c(x,y)$ outside the region of interest ($|x| \leq x_1, |y| \leq y_1$). An optimality criterion for the selection of a suitable kernel was proposed by Schwab in a 1983 paper:
\begin{equation}
R  = \ddfrac{\iint_A |c(x,y)|^2 w(x,y) dx dy}{\int_{-\infty}^{\infty} \int_{-\infty}^{\infty} |c(x,y)|^2 w(x,y) dx dy},
\end{equation}
where $A$ is the region in which $c(x,y)$ should be concentrated and $w(x,y)$ are weight functions. Schwab also showed that $R$ is maximised by zero-order prolate spheroidal wavefunctions, which can be approximated by Kaiser-Bessel window functions.

