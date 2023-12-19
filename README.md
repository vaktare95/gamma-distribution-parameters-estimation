This is my implementation of estimation of two parameters of Gamma distribution: $\alpha$ (shape) and $\theta$ (scale). The estimations of parameters $\hat{\alpha}$ and $\hat{\theta}$ are calculated based on equations derived from maximum likelihood estimation (MLE):<br>
$\hat{\theta} = \frac{\overline{X}_n}{\hat{\alpha}}$,<br>
$\frac{\Gamma\'(\hat{\alpha})}{\Gamma(\hat{\alpha})} - \ln(\hat{\alpha}) - C_n = 0$,<br>
where $\overline{X}_n$ is a mean based on the $n$-element probe, $\Gamma(\alpha)$ is the Gamma function:<br>
$\Gamma(\alpha) =  \int_0^{\infty} u^{\alpha-1} e^{-u} du  $<br>
and constant $C_n = \overline{(\ln{X})}_n - \ln \overline{X}_n$.<br>

The second equation for $\hat{\alpha}$ calculation is solved numerically using Newton-Raphson method:<br>
$f(\alpha) = \frac{\Gamma\'(\alpha)}{\Gamma(\alpha)} - \ln(\alpha) - C_n$,<br>
$\hat{\alpha} _{i+1} = \hat{\alpha} _i - \frac{f(\hat{\alpha} _i)}{f\'(\hat{\alpha} _i)}$.

Derivative $\Gamma\'(\hat{\alpha})$ is calculated using central difference method and integral  $\int_0^{\infty} u^{\alpha-1} e^{-u} du$  is calculated using trapezoid method.

The code is adjusted to Python3.8.

