# Gaussian Wave Packet Analytical Solution Derivation


# Background Introduction

The Gaussian wave packet is a solution of the free-particle Schrodinger equation,

$$
  i \hbar \frac{ \partial }{ \partial t } \psi( x, t ) = - \frac{ \hbar^2 }{ 2 m } \frac{ \partial^2 }{ \partial x^2 } \psi( x, t )
$$

where $\psi( x, t )$ is the wave function in x-space, $m$ is the particle mass, and $\hbar$ is the reduced Planck constant.

Its wave function can be seen as a superposition of the plane waves, $A( k ) e^{ i ( k x - \omega t ) }$, where the angular frequency $\omega$ and the wavenumber $k$ follow the dispersion relation

$$
  \omega = \frac{ \hbar k^2 }{ 2 m }
$$

and the amplitudes of different plane wave wavelengths $A( k )$ follow a Gaussian distribution.


---

# k-space Distribution

## Center of the Gaussian distribution in k-space

$$
  k_0 = \mathrm{ center\ of\ the\ initial\ density\ distribution\ in\ k\ space\ } \tilde{ \rho }( k, t=0 ) = | \tilde{ \psi }( k, t=0 ) |^2
$$

- The velocity of the center of the wave packet in x-space = the group velocity of the central k-mode =

$$
  v_0 = \frac{ \partial \omega }{ \partial k } {\Large \mid}_{ k = k_0 } = \frac{ \hbar }{ m } k_0
$$


## Width of the Gaussian distribution in k-space

$$
  \sigma_{ k, 0 } = \mathrm{ standard\ deviation\ of\ the\ initial\ density\ distribution\ in\ k\ space\ } \tilde{ \rho }( k, t=0 ) = | \tilde{ \psi }( k, t=0 ) |^2
$$

## Wave function in k-space

At $t=0$, the density in k-space is a Gaussian distribution with the center $k_0$ and width $\sigma_{k, 0}$.

The corresponding wave function is

$$
  \tilde{ \psi }( k, t=0 ) = \frac{ 1 }{ \sqrt{ \sigma_{ k, 0 } \sqrt{ 2 \pi } } } \exp \Big[ - \frac{ ( k - k_0 )^2 }{ 4 \sigma_{ k, 0 }^2 } \Big]
$$

- The normalization of the wave function in k-space: $\int_{-\infty}^{\infty} |\tilde{\psi}(k, t=0)|^2 dk = 1$

At time $t$, we have

$$
  \begin{aligned}
  \tilde{ \psi }( k, t )
  &=
  \tilde{ \psi }( k, t=0 ) e^{ -i \omega t } \\
  &=
  \frac{ 1 }{ \sqrt{ \sigma_{ k, 0 } \sqrt{ 2 \pi } } } \exp \Big[ - \frac{ ( k - k_0 )^2 }{ 4 \sigma_{ k, 0 }^2 } - i \frac{ \hbar k^2 }{ 2 m } t \Big]
  \end{aligned}
$$


---

# x-space Distribution

## Center of the Gaussian distribution in x-space

$$
  x_0 = \mathrm{ center\ of\ the\ initial\ density\ distribution\ } \rho( x, t=0 ) = | \psi( x, t=0 ) |^2
$$

## Width of the Gaussian distribution in x-space

$$
  \sigma_{ x, 0 } = \mathrm{ standard\ deviation\ of\ the\ initial\ density\ distribution\ } \rho( x, t=0 ) = | \psi( x, t=0 ) |^2
$$

- The width of the Gaussian distribution in x-space and in k-space are inversely proportional to each other, which is the uncertainty principle

$$
  \sigma_{ x, 0 } \times \sigma_{ k, 0 } = \frac{ 1 }{ 2 }
$$


## Wave function in x-space

The wave function in the x-space is the inverse Fourier transform of the wave function in the k-space

$$
  \psi( x^\prime, t ) = \frac{ 1 }{ \sqrt{ 2 \pi } } \int_{ -\infty }^{ \infty } \tilde{ \psi }( k, t ) e^{ i k x^\prime } dk
$$

where $x^\prime = x - x_0$.

With the Gaussian k-space wave function, we get

$$
  \begin{aligned}
  \psi( x, t )
  &=
  \sqrt{ \frac{ 1 }{ 2 \pi \sigma_{ k, 0 } \sqrt{ 2 \pi } } } \int_{ -\infty }^{ \infty } \exp[ -\sigma_{ x, 0 }^2 ( k - k_0 )^2 - i \frac{ \hbar t }{ 2 m } k^2 + i k ( x - x_0 ) ] dk \\
  &=
  \sqrt{ \frac{ \sigma_{ x, 0 } }{ \pi \sqrt{ 2 \pi } } } \int_{ -\infty }^{ \infty } \exp{ [ -( ( \sigma_{ x, 0 }^2 + i \frac{ \hbar t }{ 2 m } ) k^2 - ( 2 \sigma_{ x, 0 }^2 k_0 + i ( x - x_0 ) ) k + \sigma_{ x, 0 }^2 k_0^2 ) ] } dk \\
  \end{aligned}
$$

With the [Gaussian integral formula](https://en.wikipedia.org/wiki/Gaussian_integral#The_integral_of_a_Gaussian_function), $\int_{ -\infty }^{ \infty } e^{ -( a x^2 + b x + c ) } dx = \sqrt{ \frac{ \pi }{ a } } e^{ \frac{ b^2 }{ 4 a } - c }$, we can integrate it to get

$$
  \begin{aligned}
  \psi( x, t )
  &=
  \sqrt{ \frac{ \sigma_{ x, 0 } }{ \pi \sqrt{ 2 \pi } } } \sqrt{ \frac{ \pi }{ ( \sigma_{ x, 0 }^2 + i \frac{ \hbar t }{ 2 m } ) } } \exp{ \Big[ \frac{ ( -i ( -i 2 \sigma_{ x, 0 }^2 k_0 + ( x - x_0 ) ) )^2 }{ 4 ( \sigma_{ x, 0 }^2 + i \frac{ \hbar t }{ 2 m } ) } - \sigma_{ x, 0 }^2 k_0^2 \Big] } \\
  &=
  \frac{ 1 }{ \sqrt{ \sigma_{ x, 0 } \sqrt{ 2 \pi } } } \sqrt{ \frac{ \sigma_{ x, 0 }^2 }{ ( \sigma_{ x, 0 }^2 + i \frac{ \hbar t }{ 2 m } ) } } \exp{ \Big[ - \frac{ ( ( x - x_0 ) - i 2 \sigma_{ x, 0 }^2 k_0 )^2 }{ 4 ( \sigma_{ x, 0 }^2 + i \frac{ \hbar t }{ 2 m } ) } - \sigma_{ x, 0 }^2 k_0^2 \Big] } \\
  &=
  \frac{ 1 }{ \sqrt{ \sqrt{ 2 \pi \sigma_{ x, 0 }^2 } } } \sqrt{ \frac{ 2 \sigma_{ x, 0 }^2 }{ ( 2 \sigma_{ x, 0 }^2 + i \frac{ \hbar t }{ m } ) } } \exp{ \Big[ - \frac{ ( ( x - x_0 ) - i 2 \sigma_{ x , 0 }^2 k_0 )^2 }{ 2 ( 2 \sigma_{ x, 0 }^2 + i \frac{ \hbar t }{ m } ) } - \frac{ 2 \sigma_{ x, 0 }^2 k_0^2 }{ 2 } \Big] }
  \end{aligned}
$$

Here, we introduce the new parameter, $\alpha = w_0^2 = 2 \sigma_{ x, 0 }^2$, to simplify the equation

$$
  \psi( x, t ) =
  \frac{ 1 }{ \sqrt{ \sqrt{ \pi \alpha } } }
  \underbrace{ \sqrt{ \frac{ \alpha }{ \alpha + i \frac{ \hbar }{ m } t } } } _\mathrm{ Part.\ I }
  \underbrace{ \exp \Big[ - \frac{ ( x - x_0 - i k_0 \alpha )^2 }{ 2 ( \alpha + i \frac{ \hbar }{ m } t ) } \Big] } _\mathrm{ Part.\ II }
  \underbrace{ \exp \Big[ - \frac{ \alpha k_0^2 }{ 2 } \Big] } _\mathrm{ Part.\ III }
$$

Let’s break it down and simplify each part

### Part. I

$$
  \begin{aligned}
  \sqrt{ \frac{ \alpha }{ \alpha + i \frac{ \hbar }{ m } t } }
  &=
  \sqrt{ \frac{ 1 }{ 1 + i \frac{ \hbar t }{ m \alpha } } } \\
  &=
  \sqrt{ \frac{ ( 1 - i \frac{ \hbar t}{ m \alpha } ) }{ ( 1 + i \frac{ \hbar t }{ m \alpha } )( 1 - i \frac{ \hbar t }{ m \alpha } ) } } \\
  &=
  \sqrt{ \frac{ 1 }{ 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } } } \sqrt{ 1 - i \frac{ \hbar t }{ m \alpha } } \\
  &=
  \frac{ 1 }{ ( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } )^{ 1 / 2 } } \Big[ \sqrt{ ( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } ) } \exp \Big( -i { \cos^{ -1 } \big ( \frac{ 1 }{ \sqrt{ 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } } } \big) } \Big) \Big]^{ 1 / 2 } \\
  &=
  \frac{ 1 }{ ( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } ) ^{ 1 / 4 } } \exp \Big[ -i \frac{ 1 }{ 2 }{ \cos^{ -1 } \big( \frac{ 1 }{ \sqrt{ 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } } } \big) }\Big]
  \end{aligned}
$$

### Part. II

$$
  \begin{aligned}
  \exp \Big[ - \frac{ ( x - x_0 - i k_0 \alpha )^2 }{ 2 ( \alpha + i \frac{ \hbar }{ m } t ) } \Big]
  &=
  \exp \Big[ - \frac{ ( x - x_0 )^2 - k_0^2 \alpha^2 - 2 i ( x - x_0 ) k_0 \alpha }{ 2 \alpha ( 1 + i \frac{ \hbar t }{ m \alpha } ) } \Big] \\
  &=
  \exp \Big[ - \frac{ \big( ( x - x_0 )^2 - k_0^2 \alpha^2 - 2 i ( x - x_0 ) k_0 \alpha \big) \big( 1 - i \frac{ \hbar t }{ m \alpha } \big) }{ 2 \alpha ( 1 + i \frac{ \hbar t }{ m \alpha } )( 1 - i \frac{ \hbar t }{ m \alpha } ) } \Big] \\
  &=
  \exp \Big[ - \frac{ \big( ( x - x_0 )^2 -k_0^2 \alpha^2 - 2 ( x - x_0 ) k_0 \alpha \frac{ \hbar t }{ m \alpha } \big) - i \big( ( x - x_0 )^2 \frac{ \hbar t }{ m \alpha } - k_0^2 \alpha^2 \frac{ \hbar t }{ m \alpha } + 2 ( x - x_0 ) k_0 \alpha \big) }{ 2 \alpha ( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } ) } \Big] \\
  &=
  \underbrace{ \exp \Big[ - \frac{ ( x - x_0 )^2 - k_0^2 \alpha^2 - 2 ( x - x_0 ) k_0 \frac{ \hbar t }{ m } }{ 2 \alpha ( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } ) } \Big] } _\mathrm{ Part.\ II-1 }
  \times
  \underbrace{ \exp \Big[ i \frac{ ( x - x_0 )^2 \frac{ \hbar t }{ m \alpha } - k_0^2 \alpha^2 \frac{ \hbar t }{ m \alpha } + 2 ( x - x_0 ) k_0 \alpha }{ 2 \alpha ( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } ) } \Big] } _\mathrm{ Part.\ II-2 }
  \end{aligned}
$$

- Part. II-1

$$
  \begin{aligned}
  \exp \Big[ - \frac{ ( x - x_0 )^2 - k_0^2 \alpha^2 - 2 ( x - x_0 ) k_0 \frac{ \hbar t }{ m } }{ 2 \alpha ( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } ) } \Big]
  &=
  \exp \Big[ - \frac{ ( x - x_0 - k_0 \frac{ \hbar t }{ m } )^2 - k_0^2 \frac{ \hbar^2 t^2 }{ m^2 } - k_0^2 \alpha^2 }{ 2 \alpha ( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } ) } \Big]
  \end{aligned}
$$

- Part. II-2

$$
  \begin{aligned}
  \exp \Big[ i \frac{ ( x - x_0 )^2 \frac{ \hbar t }{ m \alpha } - k_0^2 \alpha^2 \frac{ \hbar t }{ m \alpha } + 2 ( x - x_0 ) k_0 \alpha }{ 2 \alpha( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } ) } \Big]
  &=
  \exp \Big[ i \frac{ \big( ( x - x_0 )^2 \frac{ \hbar t }{ m \alpha } - ( 2 ( x - x_0 ) k_0 \alpha \frac{ \hbar t }{ m \alpha } ) \frac{ \hbar t }{ m \alpha } + ( k_0^2 \alpha^2 \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } ) \frac{ \hbar t }{ m \alpha } \big) + ( 2 ( x - x_0 ) k_0 \alpha \frac{ \hbar t }{ m \alpha } ) \frac{ \hbar t }{ m \alpha } - ( k_0^2 \alpha^2 \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } ) \frac{ \hbar t }{ m \alpha } -k_0^2 \alpha^2 \frac{ \hbar t }{ m \alpha } + 2 ( x - x_0 ) k \alpha }{ 2 \alpha ( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } ) } \Big] \\
  &=
  \exp \Big[ i \frac{ \big( ( x - x_0 -k_0 \alpha \frac{ \hbar t }{ m \alpha } )^2 \frac{ \hbar t }{ m \alpha } \big) + 2 ( x - x_0 ) k_0 \alpha \big( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } \big) - k_0^2 \alpha^2 \frac{ \hbar t }{ m \alpha } \big( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } \big) }{ 2 \alpha ( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } ) } \Big] \\
  &=
  \exp \Big[ i \frac{ ( x - x_0 - k_0 \alpha \frac{ \hbar t }{ m \alpha } )^2 \frac{ \hbar t }{ m \alpha } }{ 2 \alpha ( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } ) } \Big] \times \exp \Big[ i \Big( k_0 ( x - x_0 ) - \frac{ 1 }{ 2 } k_0^2 \alpha \frac{ \hbar t }{ m \alpha } \Big) \Big] \\
  &=
  \exp \Big[ i \frac{ ( x - x_0 - k_0 \frac{ \hbar t }{ m } )^2 \frac{ \hbar t }{ m \alpha } }{ 2 \alpha ( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } ) } \Big] \times \exp \Big[ i \Big( k_0 ( x - x_0 - \frac{ 1 }{ 2 } k_0 \frac{ \hbar t }{ m } ) \Big) \Big]
  \end{aligned}
$$


### Part. III

$$
  \exp \Big[ - \frac{ \alpha k_0^2 }{ 2 } \Big] =
  \exp \Big[ - \frac{ \alpha^2 k_0^2 ( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } ) }{ 2 \alpha ( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } ) } \Big] =
  \exp \Big[ - \frac{ k_0^2 \alpha^2 + k_0^2 \frac{ \hbar^2 t^2 }{ m^2 } }{ 2 \alpha ( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } ) } \Big]
$$

Finally, we can combine all the parts together to get

$$
  \begin{aligned}
  \psi( x, t )
  &=
  \frac{ 1 }{ ( \pi \alpha )^{ 1 / 4 } } \frac{ 1 }{ ( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } )^{ 1 / 4 } } \exp \Big[ -i \frac{ 1 }{ 2 } { \cos^{ -1 } \big( \frac{ 1 }{ \sqrt{ 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } } } \big) } \Big] \exp \Big[ - \frac{ ( x - x_0 -k_0 \frac{ \hbar t }{ m } )^2 - k_0^2 \frac{ \hbar^2 t^2 }{ m^2 } - k_0^2 \alpha^2 }{ 2 \alpha ( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } ) } \Big] \exp \Big[ i \frac{ ( x - x_0 - k_0 \frac{ \hbar t }{ m } )^2 \frac{ \hbar t }{ m \alpha } }{ 2 \alpha ( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } ) } \Big] \exp \Big[ i \Big( k_0 ( x - x_0 - \frac{ 1 }{ 2 } k_0 \frac{ \hbar t }{ m } ) \Big) \Big] \exp \Big[ - \frac{ k_0^2 \alpha^2 + k_0^2 \frac{ \hbar^2 t^2 }{ m^2 } }{ 2 \alpha ( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } ) } \Big] \\
  &=
  \frac{ 1 }{ ( \pi \alpha )^{ 1 / 4 } } \frac{ 1 }{ ( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } )^{ 1 / 4 } } \exp \Big[ -i \frac{ 1 }{ 2 } { \cos^{ -1 } \big( \frac{ 1 }{ \sqrt{ 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } } } \big) } \Big] \exp \Big[ - \frac{ ( x - x_0 - k_0 \frac{ \hbar }{ m } t )^2 }{ 2 \alpha ( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } ) } \Big] \exp \Big[ i \frac{ ( x - x_0 - k_0 \frac{ \hbar t }{ m } )^2 \frac{ \hbar t }{ m \alpha } }{ 2 \alpha ( 1 + \frac{ \hbar^2 t^2 }{ m^2 \alpha^2 } ) } \Big] \exp \Big[ i \Big( k_0 ( x - x_0 - \frac{ 1 }{ 2 } k_0 \frac{ \hbar t }{ m } ) \Big) \Big] \\
  \end{aligned}
$$

$$
  \begin{aligned}
  \psi( x, t )
  &=
  \frac{ 1 }{ ( \pi w_0^2 )^{ 1 / 4 } } \frac{ 1 }{ [ 1 + \frac{ \hbar^2 t^2 }{ m^2 ( w_0^2 )^2 } ]^{ 1 / 4 } } \exp \Big[ -i \frac{ 1 }{ 2 } { \cos^{ -1 } \big( \frac{ 1 }{ \sqrt{ 1 + \frac{ \hbar^2 t^2 }{ m^2 ( w_0^2 )^2 } } } \big) } \Big] \exp \Big[ - \frac{ ( x - x_0 - v_0 t )^2 }{ 2 ( w_0^2 ) ( 1 + \frac{ \hbar^2 t^2 }{ m^2 ( w_0^2 )^2 } ) } \Big] \exp \Big[ i \frac{ ( x - x_0 - v_0 t )^2 \frac{ \hbar t }{ m ( w_0^2 ) } }{ 2 ( w_0 )^2 ( 1 + \frac{ \hbar^2 t^2 }{ m^2 ( w_0^2 )^2 } ) } \Big] \exp \Big[ i \Big( v_0 \frac{ m }{ \hbar } ( x - x_0 - \frac{ 1 }{ 2 } v_0 t \Big) \Big] \\
  &=
  \frac{ 1 }{ \{ \pi w_0^2 [ 1 + ( \frac{ \hbar t }{ m w_0^2 } )^2 ] \}^{ 1 / 4 } } \exp \Big[ { - \frac{ 1 }{ 2 } \frac{ ( x - x_0 - v_0 t )^2 }{ w_0^2 [ 1 + ( \frac{ \hbar t }{ m w_0^2 } )^2 ] } } \Big] \exp \Big[ i \Big( - \frac{ 1 }{ 2 } \cos^{ -1 } \big( \frac{ 1 }{ \sqrt{ 1 + ( \frac{ \hbar t }{ m w_0^2 } )^2 } } \big ) + \frac{ 1 }{ 2 } ( x - x_0 - v_0 t )^2 \frac{ ( \frac{ \hbar t }{ m w_0^2 } ) }{ w_0^2 [ 1 + ( \frac{ \hbar t }{ m w_0^2 } )^2 ] } + v_0 \frac{ m }{ \hbar }( x - x_0 - \frac{ 1 }{ 2 } v_0 t ) \Big) \Big] \\
  \end{aligned}
$$

$$
  \begin{aligned}
  \psi( x, t )
  &=
  \frac{ 1 }{ \{ 2 \pi \sigma_{ x, 0 }^2 [ 1 + ( \frac{ \hbar t }{ 2 m \sigma_{ x, 0 }^2 } )^2 ] \}^{ 1 / 4 } } \exp \Big[ { - \frac{ 1 }{ 2 } \frac{ ( x - x_0 - v_0 t )^2 }{ 2 \sigma_{ x, 0 }^2 [ 1 + ( \frac{ \hbar t }{ 2 m \sigma_{ x , 0 }^2 } )^2 ] } } \Big] \exp \Big[ i \Big( - \frac{ 1 }{ 2 } \cos^{ -1 } \big( \frac{ 1 }{ \sqrt{ 1 + ( \frac{ \hbar t }{ 2 m \sigma_{ x, 0 }^2 } )^2 } } \big) + \frac{ 1 }{ 2 }( x - x_0 - v_0 t )^2 \frac{ ( \frac{ \hbar t }{ 2 m \sigma_{ x, 0 }^2 } ) }{ 2 \sigma_{ x, 0 }^2 [ 1 + ( \frac{ \hbar t }{ 2 m \sigma_{ x, 0 }^2 } )^2 ] } + v_0 \frac{ m }{ \hbar }( x - x_0 - \frac{ 1 }{ 2 } v_0 t ) \Big) \Big]
  \end{aligned}
$$

## Density Distribution in x-space

The density distribution is the square of the wave function

$$
  \begin{aligned}
  \rho( x, t ) = | \psi( x, t ) |^2
  &=
  \frac{ 1 }{ \{ 2 \pi \sigma_{ x, 0 }^2 [ 1 + ( \frac{ \hbar t }{ 2 m \sigma_{ x, 0 }^2 } )^2 ] \}^{ 1 / 2 } } \exp \Big[ { - \frac{ ( x - x_0 - v_0 t )^2 }{ 2 \sigma_{ x, 0 }^2 [ 1 + ( \frac{ \hbar t }{ 2 m \sigma_{ x, 0 }^2 } )^2 ] } } \Big] \\
  &=
  \frac{ 1 }{ \sigma_{ x, t } \sqrt{ 2 \pi } } \exp \Big[ { - \frac{ ( x - x_0 - v_0 t )^2 }{ 2 \sigma_{ x, t }^2 } } \Big]
  \end{aligned}
$$

- The normalization of the density distribution: $\int_{ -\infty }^{ \infty } \rho( x, t ) dx = 1$
- The center of the Gaussian distribution will be at $x_0 + v_0 t$ and move with velocity $v_0$.
- The standard deviation of the Gaussian distribution $\sigma_{ x, t } = \sigma_{ x, 0 } \sqrt{ 1 + ( \frac{ \hbar t }{ 2 m \sigma_{ x, 0 }^2 } )^2 }$ will increase with time.


---

# In GAMER ELBDM `GaussianWavePacket` Test Problem

## Parameters

- `Gau_v0`

$$
  v_0
$$

- `Gau_Width`

$$
  w_0 = \sqrt{ \alpha } = \sqrt{ 2 } \times \sigma_{ x, 0 }
$$

- `Gau_Center`

$$
  x_0
$$

- `ELBDM_ETA`

$$
  \frac{ m }{ \hbar }
$$

- `Gau_Const1`

$$
  \begin{aligned}
  C_1
  &=
  1 + \Big( \frac{ t }{ \frac{ m }{ \hbar } w_0^2 } \Big)^2 \\
  &=
  1 + \Big( \frac{ \hbar t }{ m w_0^2 } \Big)^2
  \end{aligned}
$$

- `Gau_Theta1`

$$
  \begin{aligned}
  \theta_1
  &=
  -\frac{ 1 }{ 2 } \cos^{ -1 } ( C_1^{ - \frac{ 1 }{ 2 } } ) \\
  &=
  -\frac{ 1 }{ 2 } \cos^{ -1 } \Big( \frac{ 1 }{ \sqrt{ 1 + ( \frac{ \hbar t }{ m w_0^2 } )^2 } } \Big)
  \end{aligned}
$$

- `Gau_Const2`

$$
  \begin{aligned}
  C_2
  &=
  ( w_0^2 \pi C_1 )^{-1/4} \exp \Big[ - \frac{ 1 }{ 2 } \Big( \frac{ ( x - v_0 t - x_0 )^2 }{ w_0^2 C_1 } \Big) \Big] \\
  &=
  \frac{ 1 }{ \{ \pi w_0^2 [ 1 + ( \frac{ \hbar t }{ m w_0^2 } )^2 ] \}^{ 1 / 4 } } \exp \Big[ { - \frac{ 1 }{ 2 } \frac{ ( x - x_0 - v_0 t )^2 }{ w_0^2 [ 1 + ( \frac{ \hbar t }{ m w_0^2 } )^2 ] } } \Big]
  \end{aligned}
$$

- `Gau_Theta2`

$$
  \begin{aligned}
  \theta_2
  &=
  \frac{ 1 }{ 2 }( x - v_0 t - x_0 )^2 \frac{ \frac{ m }{ \hbar} t }{ ( \frac{ m }{ \hbar } w_0^2 )^2 + t^2 } + v_0 \frac{ m }{ \hbar }( x - \frac{ 1 }{ 2 } v_0 t - x_0 ) \\
  &=
  \frac{ 1 }{ 2 } ( x - x_0 - v_0 t )^2 \frac{ ( \frac{ \hbar t }{ m w_0^2 } ) }{ w_0^2 [ 1 + ( \frac{ \hbar t }{ m w_0^2 } )^2 ] } + v_0 \frac{ m }{ \hbar }( x - x_0 - \frac{ 1 }{ 2 } v_0 t )
  \end{aligned}
$$

- Real part of the wave function

$$
  \mathrm{ Re }[ \psi ] = C_2 \cos( \theta_1 + \theta_2 )
$$

- Imaginary part of the wave function

$$
  \mathrm{ Im }[ \psi ] = C_2 \sin( \theta_1 + \theta_2 )
$$

## Wave function

$$
  \begin{aligned}
  \psi
  &=
  \mathrm{ Re }[ \psi ] + i \mathrm{ Im }[ \psi ] \\
  &=
  C_2 e^{ i ( \theta_1 + \theta_2 ) }\\
  &=
  \frac{ 1 }{ \{ \pi w_0^2 [ 1 + ( \frac{ \hbar t }{ m w_0^2 } )^2 ] \}^{ 1 / 4 } } \exp \Big[ { - \frac{ 1 }{ 2 } \frac{ ( x - x_0 - v_0 t )^2 }{ w_0^2 [ 1 + ( \frac{ \hbar t }{ m w_0^2 } )^2 ] } } \Big] \exp \Big[ i \Big( - \frac{ 1 }{ 2 } \cos^{ -1 } \big( \frac{ 1 }{ \sqrt{ 1 + ( \frac{ \hbar t }{ m w_0^2 } )^2 } } \big ) + \frac{ 1 }{ 2 }( x - x_0 - v_0 t )^2 \frac{ ( \frac{ \hbar t }{ m w_0^2 } ) }{ w_0^2 [ 1 + ( \frac{ \hbar t }{ m w_0^2 } )^2 ] } + v_0 \frac{ m }{ \hbar }( x - x_0 - \frac{ 1 }{ 2 } v_0 t ) \Big) \Big] \\
  &=
  \frac{ 1 }{ \{ 2 \pi \sigma_{ x, 0 }^2 [ 1 + ( \frac{ \hbar t }{ 2 m \sigma_{ x, 0 }^2 } )^2 ] \}^{ 1 / 4 } } \exp \Big[ { - \frac{ 1 }{ 2 } \frac{ ( x - x_0 - v_0 t )^2 }{ 2 \sigma_{ x, 0 }^2 [ 1 + ( \frac{ \hbar t }{ 2 m \sigma_{ x, 0 }^2 } )^2 ] } } \Big] \exp \Big[ i \Big( - \frac{ 1 }{ 2 } \cos^{ -1 } \big( \frac{ 1 }{ \sqrt{ 1 + ( \frac{ \hbar t }{ 2 m \sigma_{ x , 0 }^2 } )^2 } } \big) + \frac{ 1 }{ 2 } ( x - x_0 - v_0 t )^2 \frac{ ( \frac{ \hbar t }{ 2 m \sigma_{ x , 0 }^2 } ) }{ 2 \sigma_{ x, 0 }^2 [ 1 + ( \frac{ \hbar t }{ 2 m \sigma_{ x, 0 }^2 } )^2 ] } + v_0 \frac{ m }{ \hbar } ( x - x_0 - \frac{ 1 }{ 2 } v_0 t ) \Big) \Big]
  \end{aligned}
$$


---
