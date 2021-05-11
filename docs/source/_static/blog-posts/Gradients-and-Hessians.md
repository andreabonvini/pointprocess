[TOC]

## Inverse Gaussian

$$
\mu(H_{u_{i}},\theta_t)=\theta_0+\sum_{j=1}^{p}\theta_ju_{i-j+1}
$$

$$
p(u_i|k,\theta_t) = \left[\frac{k}{2\pi\cdot{u_i}^3}\right]^{1/2}e^{-\frac{1}{2}\frac{k\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)^2}{\mu(H_{u_{i-1}},\theta_t)^2\cdot u_i}}
$$

$$
\mathcal{L}(u_{t-n:t}|k,\mathbf{\eta}_t,\theta_t) = -\sum_{i=1}^{n}\eta_i\log p({u}_i|k,\theta_t)\\
$$

### Gradient `Vector`

$$
\nabla\mathcal{L}(u_{t-m:t}|k,\mathbf{\eta}_t,\theta_t)=
\begin{bmatrix}
\frac{\part}{\part k}\mathcal{L}(u_{t-m:t}|k,\mathbf{\eta}_t,\theta_t)\\
\frac{\part}{\part \theta_0}\mathcal{L}(u_{t-m:t}|k,\mathbf{\eta}_t,\theta_t)\\
\vdots\\
\frac{\part}{\part \theta_{p}}\mathcal{L}(u_{t-m:t}|k,\mathbf{\eta}_t,\theta_t)
\end{bmatrix}
$$

For $k$ we have: 
$$
\frac{\part}{\part k}\mathcal{L}(u_{t-m:t}|k,\mathbf{\eta}_t,\theta_t) =\\

\frac{\part}{\part k}\left(-\sum_{i=0}^{m-1}\eta_i\log p(u_i|k,\theta_t)\right) = \\

\frac{\part}{\part k}
\left(
-
\sum_{i=0}^{m-1}
\eta_i
\log
\left(
\left[\frac{k}{2\pi\cdot{u_i}^3}\right]^{1/2}e^{-\frac{1}{2}\frac{k\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)^2}{\mu(H_{u_{i-1}},\theta_t)^2\cdot u_i}}
\right)
\right)=\\

\frac{\part}{\part k}
\left(
\sum_{i=0}^{m-1}
-
\frac{1}{2}
\eta_i
\log
\left(
\left[\frac{k}{2\pi\cdot{u_i}^3}\right]
\right)
+\frac{1}{2}
\eta_i
\frac{k\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)^2}{\mu(H_{u_{i-1}},\theta_t)^2\cdot u_i}
\right)=\\

\frac{\part}{\part k}
\left(
\sum_{i=0}^{m-1}
-
\frac{1}{2}
\eta_i
\log
\left(
k
\right)
+\frac{1}{2}\eta_i\log
\left(2\pi\cdot{u_i}^3\right)
+\frac{1}{2}
\eta_i
\frac{k\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)^2}{\mu(H_{u_{i-1}},\theta_t)^2\cdot u_i}
\right)=\\

\color{blue}
{
\sum_{i=0}^{m-1}
\frac{\eta_i}{2}
\left(-\frac{1}{k}+\frac{\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)^2}{\mu(H_{u_{i-1}},\theta_t)^2\cdot u_i}\right)
}
$$
For $\theta_j$ we have:
$$
\frac{\part}{\part \theta_j}\mathcal{L}(u_{t-m:t}|k,\mathbf{\eta}_t,\theta_t) =\\

\frac{\part}{\part \theta_j}\left(-\sum_{i=0}^{m-1}\eta_i\log p(u_i|k,\theta_t)\right) = \\

\frac{\part}{\part \theta_j}
\left(
-
\sum_{i=0}^{m-1}
\eta_i
\log
\left(
\left[\frac{k}{2\pi\cdot{u_i}^3}\right]^{1/2}e^{-\frac{1}{2}\frac{k\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)^2}{\mu(H_{u_{i-1}},\theta_t)^2\cdot u_i}}
\right)
\right)=\\

\frac{\part}{\part \theta_j}
\left(
\sum_{i=0}^{m-1}
-
\frac{1}{2}
\eta_i
\log
\left(
\left[\frac{k}{2\pi\cdot{u_i}^3}\right]
\right)
+\frac{1}{2}
\eta_i
\frac{k\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)^2}{\mu(H_{u_{i-1}},\theta_t)^2\cdot u_i}
\right)=\\

\frac{\part}{\part \theta_j}
\left(
\sum_{i=0}^{m-1}
-
\frac{1}{2}
\eta_i
\log
\left(
k
\right)
+\frac{1}{2}\eta_i\log
\left(2\pi\cdot{u_i}^3\right)
+\frac{1}{2}
\eta_i
\frac{k\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)^2}{\mu(H_{u_{i-1}},\theta_t)^2\cdot u_i}
\right)=\\


\frac{\part}{\part \theta_j}
\left(
\sum_{i=0}^{m-1}
\frac{\eta_ik}{2u_i}
\frac
{
\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)^2
}
{
\mu(H_{u_{i-1}},\theta_t)^2
}
\right)=\\


\sum_{i=0}^{m-1}
\frac{\eta_ik}{2u_i}
\cdot
\frac{
2\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)(-u_{i-j})\cdot\mu(H_{u_{i-1}},\theta_t)^2
-
2\mu(H_{u_{i-1}},\theta_t)(u_{i-j})\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)^2
}
{\mu(H_{u_{i-1}},\theta_t)^4}=\\


\sum_{i=0}^{m-1}
\frac{\eta_ik}{u_i}
\cdot
\frac{
-\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)(u_{i-j})\cdot\mu(H_{u_{i-1}},\theta_t)
-
(u_{i-j})\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)^2
}
{\mu(H_{u_{i-1}},\theta_t)^3}=\\


\sum_{i=0}^{m-1}
\frac{\eta_ik}{u_i}
\cdot
\frac{
\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)(u_{i-j})
\left(
-\mu(H_{u_{i-1}},\theta_t)
-u_i
+\mu(H_{u_{i-1}},\theta_t)
\right)
}
{\mu(H_{u_{i-1}},\theta_t)^3}=\\


-\sum_{i=0}^{m-1}
\frac{\eta_ik}{u_i}
\cdot
\frac{
\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)(u_{i-j})
\left(
u_i
\right)
}
{\mu(H_{u_{i-1}},\theta_t)^3}=\\

\color{blue}
{
-\sum_{i=0}^{m-1}
{\eta_ik}
\cdot
\frac{
\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)(u_{i-j})
}
{\mu(H_{u_{i-1}},\theta_t)^3}
}
$$
Since $\mu(H_{u_{i-1}},\theta_t)=\theta_0+\sum_{j=1}^{n-1}\theta_ju_{i-j}$  we note that in case $\color{black}{j=0}$, we must use $1$ instead of $u_{i-j} = u_i$

(`Note: In code this case is handled by the fact the u[i] = 1`)

### Hessian `Matrix`

$$
\nabla^2\mathcal{L}(u_{t-n:t}|k,\mathbf{\eta}_t,\theta_t)=

\begin{bmatrix}
\frac{\part^2}{\part^2 k}\mathcal{L}(\cdot)
&
\frac{\part}{\part k\part\theta_0}\mathcal{L}(\cdot)
&
\cdots
&
\frac{\part}{\part k\part\theta_{p}}\mathcal{L}(\cdot)\\


\\

\cdot &
\frac{\part^2}{\part^2 \theta_0}\mathcal{L}(\cdot)
&
\cdots
&
\frac{\part}{\part \theta_0\theta_{p}}\mathcal{L}(\cdot)
\\


\cdot & \cdot &  \ddots & \vdots

\\

\cdot & \cdot & \cdots &
\frac{\part^2}{\part^2 \theta_{p}}\mathcal{L}(\cdot)
\end{bmatrix}
$$

$$
\frac{\part^2}{\part^2 k}\mathcal{L}(u_{t-n:t}|k,\mathbf{\eta}_t,\theta_t) =\\

\frac{\part}{\part k}
\left(
{
\sum_{i=1}^{n}
\frac{\eta_i}{2}
\left(-\frac{1}{k}+\frac{\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)^2}{\mu(H_{u_{i-1}},\theta_t)^2\cdot u_i}\right)
}
\right)=\\


\frac{\part}{\part k}
\left(
-
\sum_{i=1}^{n}
\frac{\eta_i}{2}k^{-1}
\right)=\\


\color{blue}{
\sum_{i=1}^{n}
\frac{\eta_i}{2}k^{-2}
}
$$


$$
\frac{\part}{\part k\part\theta_{j}}\mathcal{L}(u_{t-m:t}|k,\mathbf{\eta}_t,\theta_t) =\\


\frac{\part}{\part \theta_j}
\left(
{
\sum_{i=1}^{n}
\frac{\eta_i}{2}
\left(-\frac{1}{k}+\frac{\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)^2}{\mu(H_{u_{i-1}},\theta_t)^2\cdot u_i}\right)
}
\right)=\\


\frac{\part}{\part \theta_j}
\left(
{
\sum_{i=1}^{n}
\frac{\eta_i}{2u_i}
\frac{\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)^2}{\mu(H_{u_{i-1}},\theta_t)^2}
}
\right)=\\
\cdots\\
\color{blue}
{
-\sum_{i=1}^{n}
{\eta_i}
\cdot
\frac{
\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)(u_{i-j})
}
{\mu(H_{u_{i-1}},\theta_t)^3}
}
$$

$$
\frac{\part}{\part k\part\theta_{j}}\mathcal{L}(u_{t-n:t}|k,\mathbf{\eta}_t,\theta_t) =
\color{blue}{
\begin{cases}
-\sum_{i=1}^{n}
{\eta_i}
\cdot
\frac{
\left(u_i-\mu_i\right)(u_{i-j})
}
{\mu_i^3}
\ \ \ \ \ \ \ \ j\neq0\\
-\sum_{i=1}^{n}
{\eta_i}
\cdot
\frac{
\left(u_i-\mu_i\right)(1)
}
{\mu_i^3}
\ \ \ \ \ \ \ \ \ \ \  j=0
\end{cases}}
$$

$$
\frac{\part}{\part \theta_j\part k}\mathcal{L}(u_{t-m:t}|k,\mathbf{\eta}_t,\theta_t) =\\

\frac{\part}{\part k}
\left(
-\sum_{i=1}^{n}
{\eta_ik}
\cdot
\frac{
\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)(u_{i-j})
}
{\mu(H_{u_{i-1}},\theta_t)^3}
\right)
=\\


\color{purple}
{
-\sum_{i=1}^{n}
{\eta_i}
\cdot
\frac{
\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)(u_{i-j})
}
{\mu(H_{u_{i-1}},\theta_t)^3}
}
$$

We have that $\frac{\part}{\part \theta_j\part k}\mathcal{L}(\cdot) = \frac{\part}{\part k\part \theta_j}\mathcal{L}(\cdot)$ as *expected*.
$$
\frac{\part}{\part \theta_j \part \theta_q}\mathcal{L}(u_{t-n:t}|k,\mathbf{\eta}_t,\theta_t) =\\


\frac{\part}{\part \theta_q}
\left(
-\sum_{i=1}^{n}\eta_i{k}\cdot
\left(
\frac{
\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)\cdot(u_{i-j})
}
{\mu(H_{u_{i-1}},\theta_t)^3}
\right)
\right)=\\


-\sum_{i=1}^{n}\eta_i{k}
(u_{i-j})\cdot
\frac{\part}{\part \theta_q}
\left(
\frac{
\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)\
}
{\mu(H_{u_{i-1}},\theta_t)^3}
\right)=\\

-\sum_{i=1}^{n}\eta_i{k}
(u_{i-j})\cdot
\left(
\frac{
(-u_{i-q})\mu(H_{u_{i-1}},\theta_t)^3
-
3\mu(H_{u_{i-1}},\theta_t)^2\left(u_{i-q}\right)
\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)\
}
{\mu(H_{u_{i-1}},\theta_t)^6}
\right)=\\

-\sum_{i=1}^{n}\eta_i{k}
(u_{i-j})\cdot
\left(
\frac{
(-u_{i-q})\mu(H_{u_{i-1}},\theta_t)
-
3\left(u_{i-q}\right)
\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)\
}
{\mu(H_{u_{i-1}},\theta_t)^4}
\right)=\\

\sum_{i=1}^{n}\eta_i{k}
(u_{i-j})\cdot
\left(
\frac{
(u_{i-q})\mu(H_{u_{i-1}},\theta_t)
+
3\left(u_{i-q}\right)
\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)\
}
{\mu(H_{u_{i-1}},\theta_t)^4}
\right)=\\

\sum_{i=1}^{n}\eta_i{k}
(u_{i-j})
\left(u_{i-q}\right)
\cdot
\left(
\frac{
\mu(H_{u_{i-1}},\theta_t)
+
3
\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)\
}
{\mu(H_{u_{i-1}},\theta_t)^4}
\right)=\\

\color{blue}
{
\sum_{i=1}^{n}\eta_i{k}(u_{i-j})(u_{i-q})\cdot
\left(
\frac{
3\cdot u_i
-
2 \cdot \mu(H_{u_{i-1}},\theta_t)
}
{\mu(H_{u_{i-1}},\theta_t)^4}
\right)
}
$$

$$
\frac{\part}{\part \theta_j \part \theta_q}\mathcal{L}(u_{t-n:t}|k,\mathbf{\eta}_t,\theta_t) =
\color{blue}{
\begin{cases}
\sum_{i=1}^{n}\eta_i{k}(u_{i-j})(u_{i-q})\cdot
\left(
\frac{
3\cdot u_i
-
2 \cdot \mu_i
}
{\mu_i^4}
\right)
\ \ \ \ \ \ j\neq0,q\neq0
\\
\sum_{i=1}^{n}\eta_i{k}(1)(u_{i-q})\cdot
\left(
\frac{
3\cdot u_i
-
2 \cdot \mu_i
}
{\mu_i^4}
\right)
\ \ \ \ \ \ \ \ \ \ \ j=0,q\neq0
\\
\sum_{i=1}^{n}\eta_i{k}(u_{i-j})(1)\cdot
\left(
\frac{
3\cdot u_i
-
2 \cdot \mu_i
}
{\mu_i^4}
\right)
\ \ \ \ \ \ \ \ \ \ \ j\neq0,q=0
\\
\sum_{i=1}^{n}\eta_i{k}(1)(1)\cdot
\left(
\frac{
3\cdot u_i
-
2 \cdot \mu_i
}
{\mu_i^4}
\right)
\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ j=0,q=0
\end{cases}
}
$$



### Right-Censoring

Some mathematical trick to treat the *pdf* and *cdf* of the Inverse Gaussian:
$$
p(u_i|k,\theta_t) = \\

\left[\frac{k}{2\pi\cdot{u_i}^3}\right]^{1/2}e^{-\frac{1}{2}\frac{k\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)^2}{\mu(H_{u_{i-1}},\theta_t)^2\cdot u_i}}=\\

\text{exp}\left(\log\left(\left[\frac{k}{2\pi\cdot{u_i}^3}\right]^{1/2}e^{-\frac{1}{2}\frac{k\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)^2}{\mu(H_{u_{i-1}},\theta_t)^2\cdot u_i}}\right)\right) = \\

\text{exp}\left(\log\left(\left[\frac{k}{2\pi\cdot{u_i}^3}\right]^{1/2}\right)+\log\left(e^{-\frac{1}{2}\frac{k\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)^2}{\mu(H_{u_{i-1}},\theta_t)^2\cdot u_i}}\right)\right) = \\

\text{exp}\left(
\frac{1}{2}\log\left(\frac{k}{2\pi\cdot{u_i}^3}\right)
-
\frac{1}{2}\frac{k\left(u_i-\mu(H_{u_{i-1}},\theta_t)\right)^2}{\mu(H_{u_{i-1}},\theta_t)^2\cdot u_i})\right) = \\
$$

$$
P(u_i|k,\theta_t) )= \\

\Phi\left(\sqrt{\frac{k}{u_i}}\left(\frac{u_i}{\mu(H_{u_{i-1}},\theta_t)}-1\right)\right)+e^{\frac{2k}{\mu(H_{u_{i-1}},\theta_t)}}\Phi\left(-\sqrt{\frac{k}{u_i}}\left(\frac{u_i}{\mu(H_{u_{i-1}},\theta_t)}+1\right)\right) = \\

\Phi\left(\sqrt{\frac{k}{u_i}}\left(\frac{u_i}{\mu(H_{u_{i-1}},\theta_t)}-1\right)\right)
+
e^{\frac{2k}{\mu(H_{u_{i-1}},\theta_t)}+\log\left[\Phi\left(-\sqrt{\frac{k}{u_i}}\left(\frac{u_i}{\mu(H_{u_{i-1}},\theta_t)}+1\right)\right)\right]}
$$

Where $\Phi$ is the standard normal distribution *cdf*.

-----

The additional term we add to our cost function when we apply right censoring is:
$$
- \color{green}{\eta_n\log\left(1-P(t-u_{n})\right)}
$$


$$
t_n = t-u_n\\
P(t) = \int_{0}^{t}p(v)dv\\
\mathcal{L}(t)=-\eta_n\log\left(1-P(t)\right)\\
\frac{\part}{\part k}
\left(-\eta_n\log\left(1-P(t_n)\right)\right) = \frac{\part \left(-\eta_n\log\left(1-P(t_n)\right)\right)}{\part \left(1-P(t_n)\right)}
\cdot
\frac{\part \left(1-P(t_n)\right)}{\part P(t_n)}
\cdot
\frac{\part P(t_n)}{\part k}=\\
-\frac{\eta_n}{1-P(t_n)}\cdot(-1)\cdot\frac{\part P(t_n)}{\part k}=\\
\frac{\eta_n}{1-P(t_n)}\cdot\frac{\part P(t_n)}{\part k}
$$

$$
P(t) =Pr(X< t)= \Phi\left(\sqrt{\frac{k}{t}}\left(\frac{t}{\mu}-1\right)\right)+e^{\frac{2k}{\mu}}\Phi\left(-\sqrt{\frac{k}{t}}\left(\frac{t}{\mu}+1\right)\right)
$$

Where $\Phi$ is the standard normal distribution *cdf*.
$$
\frac{\part}{\part k}P(t)=\\

 \Phi^{'}\left(\sqrt{\frac{k}{t}}\left(\frac{t}{\mu}-1\right)\right)
 \color{blue}{\frac{\part}{\part k}\left(\sqrt{\frac{k}{t}}\left(\frac{t}{\mu}-1\right)\right)}
 
 + \\
 
 \color{green}{\frac{\part}{\part k}\left(\frac{2k}{\mu}\right)}e^{\frac{2k}{\mu}}\Phi\left(-\sqrt{\frac{k}{t}}\left(\frac{t}{\mu}+1\right)\right) 
 
 + \\
 
 \Phi^{'}\left(-\sqrt{\frac{k}{t}}\left(\frac{t}{\mu}+1\right)\right)
 \color{purple}{\frac{\part}{\part k}\left(-\sqrt{\frac{k}{t}}\left(\frac{t}{\mu}+1\right)\right)}
 e^{\frac{2k}{\mu}}
 
 
 \\= \\
 
 \Phi^{'}\left(\sqrt{\frac{k}{t}}\left(\frac{t}{\mu}-1\right)\right)

 \color{blue}{\left(\frac{1}{2t\sqrt{\frac{k}{t}}}\left(\frac{t}{\mu}-1\right)\right)}
 
 +\\
 
 \color{green}{\frac{2}{\mu}}e^{\frac{2k}{\mu}}\Phi\left(-\sqrt{\frac{k}{t}}\left(\frac{t}{\mu}+1\right)\right)
 
 + \\
 
 \Phi^{'}\left(-\sqrt{\frac{k}{t}}\left(\frac{t}{\mu}+1\right)\right)
 \color{purple}{
 \left(-\frac{1}{2t
 \sqrt{\frac{k}{t}}}\left(\frac{t}{\mu}+1\right)\right)}
 e^{\frac{2k}{\mu}}
$$

Which can be rewritten as
$$
\frac{\part}{\part k}P(t)= \\
 
 \Phi^{'}\left(\sqrt{\frac{k}{t}}\left(\frac{t}{\mu}-1\right)\right)

 \color{blue}{\left(\frac{1}{2t\sqrt{\frac{k}{t}}}\left(\frac{t}{\mu}-1\right)\right)}
 
 +\\
 
 \color{green}{\frac{2}{\mu}}e^{\frac{2k}{\mu}+\log\left(\Phi\left(-\sqrt{\frac{k}{t}}\left(\frac{t}{\mu}+1\right)\right)\right)}
 
 + \\
 

 \color{purple}{
 \left(-\frac{1}{2t
 \sqrt{\frac{k}{t}}}\left(\frac{t}{\mu}+1\right)\right)}
 e^{\frac{2k}{\mu} + \log\left( \Phi^{'}\left(-\sqrt{\frac{k}{t}}\left(\frac{t}{\mu}+1\right)\right)\right)}
$$
We prefer this formulation since computing `exp(2k/mu)` alone would be impracticable since the values for $k$ may be greater than $10^3$, moreover there exists valid approximations of exact formulations for the `logcdf` and `logpdf` of the normal distribution.

$$
\frac{\part}{\part \mu}P(t)=\\

\Phi^{'}\left(\sqrt{\frac{k}{t}}\left(\frac{t}{\mu}-1\right)\right)
 \color{blue}{\frac{\part}{\part \mu}\left(\sqrt{\frac{k}{t}}\left(\frac{t}{\mu}-1\right)\right)}
 
 + \\
 
 \color{green}{\frac{\part}{\part \mu}\left(\frac{2k}{\mu}\right)}e^{\frac{2k}{\mu}}\Phi\left(-\sqrt{\frac{k}{t}}\left(\frac{t}{\mu}+1\right)\right) 
 
 + \\
 
 \Phi^{'}\left(-\sqrt{\frac{k}{t}}\left(\frac{t}{\mu}+1\right)\right)
 \color{purple}{\frac{\part}{\part \mu}\left(-\sqrt{\frac{k}{t}}\left(\frac{t}{\mu}+1\right)\right)}
 e^{\frac{2k}{\mu}}

\\ = \\
 
 \Phi^{'}\left(\sqrt{\frac{k}{t}}\left(\frac{t}{\mu}-1\right)\right)
 \color{blue}{\left(-\sqrt{\frac{k}{t}}\cdot \frac{t}{\mu^2}\right)}
 
 +\\
 
 \color{green}{\left(-\frac{2k}{\mu^2}\right)}e^{\frac{2k}{\mu}}\Phi\left(-\sqrt{\frac{k}{t}}\left(\frac{t}{\mu}+1\right)\right)
 +\\
 
 \Phi^{'}\left(-\sqrt{\frac{k}{t}}\left(\frac{t}{\mu}+1\right)\right)
 \color{purple}{\left(\sqrt{\frac{k}{t}}\cdot\frac{t}{\mu^2}\right)}
 e^{\frac{2k}{\mu}}
$$

Which can be rewritten as
$$
\frac{\part}{\part \mu}P(t)= \\
 
 \Phi^{'}\left(\sqrt{\frac{k}{t}}\left(\frac{t}{\mu}-1\right)\right)
 \color{blue}{\left(-\sqrt{\frac{k}{t}}\cdot \frac{t}{\mu^2}\right)}
 
 +\\
 
 \color{green}{\left(-\frac{2k}{\mu^2}\right)}e^{\frac{2k}{\mu}+\log\left(\Phi\left(-\sqrt{\frac{k}{t}}\left(\frac{t}{\mu}+1\right)\right)\right)}
 +\\
 
 
 \color{purple}{\left(\sqrt{\frac{k}{t}}\frac{t}{\mu^2}\right)}
 e^{\frac{2k}{\mu}+\log\left(\Phi^{'}\left(-\sqrt{\frac{k}{t}}\left(\frac{t}{\mu}+1\right)\right)\right)}
$$
Then, we can add to the gradient the following terms
$$
\frac{\part}{\part k}
\left(-\eta_n\log\left(1-P(t_n)\right)\right) = 
\frac{\eta_n}{1-P(t_n)}\cdot\frac{\part P(t_n)}{\part k}
\\
\frac{\part}{\part \theta_j}
\left(-\eta_n\log\left(1-P(t_n)\right)\right) = 
\frac{\eta_n}{1-P(t_n)}\cdot\frac{\part P(t_n)}{\part \mu}\cdot\frac{\part \mu}{\part \theta_j}
$$

(Remember that $\mu(H_{u_{i}},\theta_t)=\theta_0+\sum_{j=1}^{p}\theta_ju_{i-j+1}$)

## Gaussian

$$
\mu(H_{u_{i}},\theta_t)=\theta_0+\sum_{j=1}^{p}\theta_ju_{i-j+1}
$$

$$
p(u_i|\sigma,\theta_t) = \frac{1}{\sigma\sqrt{2\pi}}e^{-\frac{1}{2}\frac{\left(u_i-\mu_i\right)^2}{\sigma}}
$$

$$
\mathcal{L}(u_{t-n:t}|\sigma,\mathbf{\eta}_t,\theta_t) = -\sum_{i=1}^{n}\eta_i\log p({u}_i|\sigma,\theta_t)\\
$$

### Gradient `Vector`

$$
\nabla\mathcal{L}(u_{t-m:t}|\sigma,\mathbf{\eta}_t,\theta_t)=
\begin{bmatrix}
\frac{\part}{\part \sigma}\mathcal{L}(u_{t-m:t}|\sigma,\mathbf{\eta}_t,\theta_t)\\
\frac{\part}{\part \theta_0}\mathcal{L}(u_{t-m:t}|\sigma,\mathbf{\eta}_t,\theta_t)\\
\vdots\\
\frac{\part}{\part \theta_{p}}\mathcal{L}(u_{t-m:t}|\sigma,\mathbf{\eta}_t,\theta_t)
\end{bmatrix}
$$

For $\sigma$ we have: 
$$
\frac{\part}{\part \sigma}\mathcal{L}(u_{t-m:t}|\sigma,\mathbf{\eta}_t,\theta_t) =\\

\frac{\part}{\part \sigma}\left(-\sum_{i=0}^{m-1}\eta_i\log p(u_i|\sigma,\theta_t)\right) = \\

\frac{\part}{\part \sigma}\left(-\sum_{i=0}^{m-1}\eta_i\log \left(\frac{1}{\sigma\sqrt{2\pi}}e^{-\frac{1}{2}\frac{\left(u_i-\mu_i\right)^2}{\sigma}}\right)\right) = \\

-\sum_{i=0}^{m-1}\eta_i\left(\frac{\part}{\part\sigma}\log \left(\frac{1}{\sigma\sqrt{2\pi}}\right)+\frac{\part}{\part\sigma}\log\left(e^{-\frac{1}{2}\frac{\left(u_i-\mu_i\right)^2}{\sigma}}\right)\right) = \\

-\sum_{i=0}^{m-1}\eta_i\left(\frac{\part}{\part\sigma}\log \left(\frac{1}{\sigma\sqrt{2\pi}}\right)+\frac{\part}{\part\sigma}\left(-\frac{1}{2}\frac{\left(u_i-\mu_i\right)^2}{\sigma}\right)\right) = \\

-\sum_{i=0}^{m-1}\eta_i
\left(-\sigma\sqrt{2\pi}\frac{1}{2\pi\sigma^2}
+
\frac{1}{2}\frac{\left(u_i-\mu_i\right)^2}{\sigma^2}\right) = \\

\color{blue}{
-\sum_{i=0}^{m-1}\eta_i
\left(
-\frac{1}{\sigma}
+
\frac{1}{2}\frac{\left(u_i-\mu_i\right)^2}{\sigma^2}
\right)}
$$
For $\theta_j$ we have:
$$
\frac{\part}{\part \theta_j}\mathcal{L}(u_{t-m:t}|\sigma,\mathbf{\eta}_t,\theta_t) =\\

\frac{\part}{\part \theta_j}\left(-\sum_{i=0}^{m-1}\eta_i\log p(u_i|\sigma,\theta_t)\right) = \\

\frac{\part}{\part \theta_j}\left(-\sum_{i=0}^{m-1}\eta_i\log \left(\frac{1}{\sigma\sqrt{2\pi}}e^{-\frac{1}{2}\frac{\left(u_i-\mu_i\right)^2}{\sigma}}\right)\right) = \\

-\sum_{i=0}^{m-1}\eta_i\left(\frac{\part}{\part\theta_j}\log\left(\frac{1}{\sigma\sqrt{2\pi}}\right)+\frac{\part}{\part\theta_j}\log\left(e^{-\frac{1}{2}\frac{\left(u_i-\mu_i\right)^2}{\sigma}}\right)\right) = \\

-\sum_{i=0}^{m-1}\eta_i\left(\frac{\part}{\part\theta_j}\log\left(\frac{1}{\sigma\sqrt{2\pi}}\right)+\frac{\part}{\part\theta_j}\left(-\frac{1}{2}\frac{\left(u_i-\mu_i\right)^2}{\sigma}\right)\right) = \\

-\sum_{i=0}^{m-1}\eta_i\left(\frac{\part}{\part\theta_j}\left(-\frac{1}{2}\frac{\left(u_i-\mu_i\right)^2}{\sigma}\right)\right) = \\

-\sum_{i=0}^{m-1}\eta_i\left(
-\frac{1}{2\sigma}\frac{\part}{\part\theta_j}\left(u_i-\mu_i\right)^2\right) = \\

\color{blue}{
-\sum_{i=0}^{m-1}\eta_i\left(
\frac{1}{\sigma}\left(u_i-\mu_i\right)u_{i-j}\right)}
$$
Since $\mu(H_{u_{i-1}},\theta_t)=\theta_0+\sum_{j=1}^{n-1}\theta_ju_{i-j}$ we note that in case $\color{black}{j=0}$, instead of $u_{i-j} = u_i$ we have to use $ 1$

### Hessian `Matrix`

$$
\nabla^2\mathcal{L}(u_{t-n:t}|\sigma,\mathbf{\eta}_t,\theta_t)=

\begin{bmatrix}
\frac{\part^2}{\part^2 \sigma}\mathcal{L}(\cdot)
&
\frac{\part}{\part k\part\theta_0}\mathcal{L}(\cdot)
&
\cdots
&
\frac{\part}{\part \sigma\part\theta_{p}}\mathcal{L}(\cdot)\\


\\

\cdot &
\frac{\part^2}{\part^2 \theta_0}\mathcal{L}(\cdot)
&
\cdots
&
\frac{\part}{\part \theta_0\theta_{p}}\mathcal{L}(\cdot)
\\


\cdot & \cdot &  \ddots & \vdots

\\

\cdot & \cdot & \cdots &
\frac{\part^2}{\part^2 \theta_{p}}\mathcal{L}(\cdot)
\end{bmatrix}
$$

$$
\frac{\part^2}{\part^2 \sigma}\mathcal{L}(u_{t-n:t}|\sigma,\mathbf{\eta}_t,\theta_t) =\\

\frac{\part}{\part\sigma}
\left(
-\sum_{i=0}^{m-1}\eta_i
\left(-\frac{1}{\sigma}
+
\frac{1}{2}\frac{\left(u_i-\mu_i\right)^2}{\sigma^2}\right)
\right)=\\

-\sum_{i=0}^{m-1}\eta_i\left(
\frac{\part}{\part\sigma}
\left(-\frac{1}{\sigma}\right)
+
\frac{1}{2}\left(u_i-\mu_i\right)^2
\frac{\part}{\part\sigma}\frac{1}{\sigma^2}
\right)

=\\

\color{purple}{
-\sum_{i=0}^{m-1}\eta_i
\left(\frac{1}{\sigma^2}
-
\left(u_i-\mu_i\right)^2
\frac{1}{\sigma^3}
\right)
}
$$

$$
\frac{\part^2}{\part\sigma\part\theta_j}\mathcal{L}(u_{t-n:t}|\sigma,\mathbf{\eta}_t,\theta_t) =\\

\frac{\part}{\part\theta_j}
\left(
-\sum_{i=0}^{m-1}\eta_i
\left(-\frac{1}{\sigma}
+
\frac{1}{2}\frac{\left(u_i-\mu_i\right)^2}{\sigma^2}\right)
\right)=\\

-\sum_{i=0}^{m-1}\eta_i\left(
\frac{1}{2\sigma^2}\frac{\part}{\part\theta_j}\left(u_i-\mu_i\right)^2
\right)

=\\

-\sum_{i=0}^{m-1}\eta_i\left(
-\frac{1}{\sigma^2}\left(u_i-\mu_i\right)u_{i-j}
\right)

=\\

\color{purple}{
\sum_{i=0}^{m-1}\eta_i\left(
\frac{1}{\sigma^2}\left(u_i-\mu_i\right)u_{i-j}
\right)
}
$$

Consistency check:
$$
\frac{\part^2}{\part\theta_j\part\sigma}\mathcal{L}(u_{t-n:t}|\sigma,\mathbf{\eta}_t,\theta_t) =\\

\frac{\part}{\part\sigma}
\left(
-\sum_{i=0}^{m-1}\eta_i\left(
\frac{1}{\sigma}\left(u_i-\mu_i\right)u_{i-j}\right)
\right)
= \\
\color{purple}{
\sum_{i=0}^{m-1}\eta_i\left(
\frac{1}{\sigma^2}\left(u_i-\mu_i\right)u_{i-j}
\right)
}
$$
We have that $\frac{\part}{\part \theta_j\part \sigma}\mathcal{L}(\cdot) = \frac{\part}{\part \sigma\part \theta_j}\mathcal{L}(\cdot)$ as *expected*.
$$
\frac{\part^2}{\part\theta_j\theta_k}\mathcal{L}(u_{t-n:t}|\sigma,\mathbf{\eta}_t,\theta_t) =\\

\frac{\part}{\part\theta_k}
\left(
-\sum_{i=0}^{m-1}\eta_i\left(
\frac{1}{\sigma}\left(u_i-\mu_i\right)u_{i-j}\right)
\right)= \\

\frac{\part}{\part\theta_k}
\left(
\sum_{i=0}^{m-1}\eta_i\left(
\frac{1}{\sigma}u_{i-j}\mu_i\right)
\right)= \\

\color{purple}{
\sum_{i=0}^{m-1}\eta_i\left(
\frac{1}{\sigma}u_{i-j}u_{i-k}
\right)
}
$$



### Right Censoring

The additional term we add to our cost function when we apply right censoring is:
$$
- \color{green}{\eta_n\log\left(1-P(t-u_{n})\right)}
$$


$$
t_n = t-u_n\\
P(t) = \int_{0}^{t}p(v)dv\\
\mathcal{L}(t)=-\eta_n\log\left(1-P(t)\right)
$$

$$
\frac{\part}{\part \sigma}
\left(-\eta_n\log\left(1-P(t_n)\right)\right) = \frac{\part \left(-\eta_n\log\left(1-P(t_n)\right)\right)}{\part \left(1-P(t_n)\right)}
\cdot
\frac{\part \left(1-P(t_n)\right)}{\part P(t_n)}
\cdot
\frac{\part P(t_n)}{\part\sigma}=\\
-\frac{\eta_n}{1-P(t_n)}\cdot(-1)\cdot\frac{\part P(t_n)}{\part \sigma}=\\

\frac{\eta_n}{1-P(t_n)}\cdot\frac{\part P(t_n)}{\part \sigma}
$$

$$
\frac{\part}{\part \theta_j}
\left(-\eta_n\log\left(1-P(t_n)\right)\right)= \frac{\eta_n}{1-P(t_n)}\cdot\frac{\part P(t_n)}{\part \theta_j}
$$

Given that $\phi$ is the standard normal distribution pdf,  we know that:
$$
\frac{\part P(t_n)}{\part \sigma} = -\left(\frac{t_n-\mu}{\sigma^2}\right)\phi\left(\frac{t_n-\mu}{\sigma}\right)
$$

$$
\frac{\part P(t_n)}{\part \mu} = -\frac{1}{\sigma}\phi\left(\frac{t_n-\mu}{\sigma}\right)
$$

Which yields:
$$
\frac{\part}{\part \sigma}
\left(-\eta_n\log\left(1-P(t_n)\right)\right) = \color{green}{\frac{-\eta_n}{1-P(t_n)}\left(\frac{t_n-\mu}{\sigma^2}\right)\phi\left(\frac{t_n-\mu}{\sigma}\right)}
$$

$$
\frac{\part}{\part \mu}
\left(-\eta_n\log\left(1-P(t_n)\right)\right) = \color{green}{\frac{-\eta_n}{1-P(t_n)}\frac{1}{\sigma}\phi\left(\frac{t_n-\mu}{\sigma}\right)}
$$

So
$$
\frac{\part}{\part \theta_j}
\left(-\eta_n\log\left(1-P(t_n)\right)\right) = \color{green}{\frac{-\eta_n}{1-P(t_n)}\frac{1}{\sigma}\phi\left(\frac{t_n-\mu}{\sigma}\right)u_{i-j}}
$$
The hessian is approximated.

We can use the gaussian distribution parameters as a starting point for the *Inverse Gaussian* parameters.

We can use $\theta_{\textit{InverseGaussian}} = \theta_{\textit{Gaussian}}$ and $k_{\textit{InverseGaussian}} = \frac{\mu_{\text{mean}}^3}{\sigma^2}$

## LogNormal
$$
\mu(H_{u_{i}},\theta_t)= \theta_0+\sum_{j=1}^{p}\theta_ju_{i-j+1}
$$

$$
p(u_i|\sigma,\theta_t) = \frac{1}{u_i\sigma\sqrt{2\pi}}e^{-\frac{1}{2}\frac{\left(\log(u_i)-\log(\mu_i)\right)^2}{\sigma^2}}
$$

$$
\mathcal{L}(u_{t-n:t}|\sigma,\mathbf{\eta}_t,\theta_t) = -\sum_{i=1}^{n}\eta_i\log p({u}_i|\sigma,\theta_t)\\
$$

### Gradient `Vector`

$$
\nabla\mathcal{L}(u_{t-m:t}|\sigma,\mathbf{\eta}_t,\theta_t)=
\begin{bmatrix}
\frac{\part}{\part \sigma}\mathcal{L}(u_{t-m:t}|\sigma,\mathbf{\eta}_t,\theta_t)\\
\frac{\part}{\part \theta_0}\mathcal{L}(u_{t-m:t}|\sigma,\mathbf{\eta}_t,\theta_t)\\
\vdots\\
\frac{\part}{\part \theta_{p}}\mathcal{L}(u_{t-m:t}|\sigma,\mathbf{\eta}_t,\theta_t)
\end{bmatrix}
$$

For $\sigma$ we have: 
$$
\frac{\part}{\part \sigma}\mathcal{L}(u_{t-m:t}|\sigma,\mathbf{\eta}_t,\theta_t) =\\

\frac{\part}{\part \sigma}\left(-\sum_{i=0}^{m-1}\eta_i\log p(u_i|\sigma,\theta_t)\right) = \\

\frac{\part}{\part \sigma}
\left(-\sum_{i=0}^{m-1}\eta_i\log
\left(

\frac{1}{u_i\sigma\sqrt{2\pi}}e^{-\frac{1}{2}\frac{\left(\log(u_i)-\log(\mu_i)\right)^2}{\sigma^2}}

\right)
\right) = \\

-\sum_{i=0}^{m-1}\eta_i
\left(
\frac{\part}{\part\sigma}
\log
\left(\frac{1}{u_i\sigma\sqrt{2\pi}}\right)
+\frac{\part}{\part\sigma}\log\left(e^{-\frac{1}{2}\frac{\left(\log\left(u_i\right)-\log(\mu_i)\right)^2}{\sigma^2}}\right)\right) = \\

-\sum_{i=0}^{m-1}\eta_i
\left(
\frac{\part}{\part\sigma}
\log
\left(\frac{1}{u_i\sigma\sqrt{2\pi}}\right)
-\frac{\part}{\part\sigma}\left({\frac{1}{2}\frac{\left(\log\left(u_i\right)-\log(\mu_i)\right)^2}{\sigma^2}}\right)\right) = \\

-\sum_{i=0}^{m-1}\eta_i
\left(
-\frac{u_i\sigma\sqrt{2\pi}}{u_i\sigma^2\sqrt{2\pi}}
-\frac{\left(\log\left(u_i\right)-\log(\mu_i)\right)^2}{2}\frac{\part}{\part\sigma}\left(\frac{1}{\sigma^2}\right)\right) = \\

\color{blue}{
-\sum_{i=0}^{m-1}\eta_i
\left(
-\frac{1}{\sigma}
+\left(\log\left(u_i\right) - \log(\mu_i)\right)^2\left(\frac{1}{\sigma^3}\right)\right)
}
$$
For $\theta_j$ we have:
$$
\frac{\part}{\part \theta_j}\mathcal{L}(u_{t-m:t}|\sigma,\mathbf{\eta}_t,\theta_t) =\\

\frac{\part}{\part \theta_j}\left(-\sum_{i=0}^{m-1}\eta_i\log p(u_i|\sigma,\theta_t)\right) = \\

\frac{\part}{\part \theta_j}
\left(-\sum_{i=0}^{m-1}\eta_i\log
\left(

\frac{1}{u_i\sigma\sqrt{2\pi}}e^{-\frac{1}{2}\frac{\left(\log(u_i)-\log(\mu_i)\right)^2}{\sigma^2}}

\right)
\right) = \\

-\sum_{i=0}^{m-1}\eta_i
\left(
\frac{\part}{\part\theta_j}
\log
\left(\frac{1}{u_i\sigma\sqrt{2\pi}}\right)
+\frac{\part}{\part\theta_j}\log\left(e^{-\frac{1}{2}\frac{\left(\log\left(u_i\right)-\log(\mu_i)\right)^2}{\sigma^2}}\right)\right) = \\

-\sum_{i=0}^{m-1}\eta_i
\left(
\frac{\part}{\part\theta_j}
\log
\left(\frac{1}{u_i\sigma\sqrt{2\pi}}\right)
-\frac{\part}{\part\theta_j}\left({\frac{1}{2}\frac{\left(\log\left(u_i\right)-\log(\mu_i)\right)^2}{\sigma^2}}\right)\right) = \\

-\sum_{i=0}^{m-1}\eta_i
\left(
-\frac{1}{2\sigma^2}\frac{\part}{\part\theta_j}\left(\left(\log\left(u_i\right)-\log(\mu_i)\right)^2\right)\right) = \\

\color{blue}{
-\sum_{i=0}^{m-1}\eta_i
\left(
\frac{1}{\sigma^2}\left(\log\left(u_i\right)-\log(\mu_i)\right)\frac{1}{\mu_i}u_{i-j}\right)
}
$$
Since $\mu(H_{u_{i-1}},\theta_t)=\theta_0+\sum_{j=1}^{n-1}\theta_ju_{i-j}$ we note that in case $\color{black}{j=0}$, instead of $u_{i-j} = u_i$ we have to use $ 1$

### Hessian `Matrix`

$$
\nabla^2\mathcal{L}(u_{t-n:t}|\sigma,\mathbf{\eta}_t,\theta_t)=

\begin{bmatrix}
\frac{\part^2}{\part^2 \sigma}\mathcal{L}(\cdot)
&
\frac{\part}{\part k\part\theta_0}\mathcal{L}(\cdot)
&
\cdots
&
\frac{\part}{\part \sigma\part\theta_{p}}\mathcal{L}(\cdot)\\


\\

\cdot &
\frac{\part^2}{\part^2 \theta_0}\mathcal{L}(\cdot)
&
\cdots
&
\frac{\part}{\part \theta_0\theta_{p}}\mathcal{L}(\cdot)
\\


\cdot & \cdot &  \ddots & \vdots

\\

\cdot & \cdot & \cdots &
\frac{\part^2}{\part^2 \theta_{p}}\mathcal{L}(\cdot)
\end{bmatrix}
$$

$$
\frac{\part^2}{\part^2 \sigma}\mathcal{L}(u_{t-n:t}|\sigma,\mathbf{\eta}_t,\theta_t) =\\

\frac{\part}{\part\sigma}
\left(
-\sum_{i=0}^{m-1}\eta_i
\left(
-\frac{1}{\sigma}
+\left(\log\left(u_i\right) - \log(\mu_i)\right)^2\left(\frac{1}{\sigma^3}\right)\right)
\right)=\\

\color{purple}{
-\sum_{i=0}^{m-1}\eta_i
\left(
\frac{1}{\sigma^2}
-\frac{3\left(\log\left(u_i\right) - \log(\mu_i)\right)^2}{\sigma^4}\right)
}
$$

$$
\frac{\part^2}{\part\sigma\part\theta_j}\mathcal{L}(u_{t-n:t}|\sigma,\mathbf{\eta}_t,\theta_t) =\\

\frac{\part}{\part\theta_j}
\left(
-\sum_{i=0}^{m-1}\eta_i
\left(
-\frac{1}{\sigma}
+\left(\log\left(u_i\right) - \log(\mu_i)\right)^2\left(\frac{1}{\sigma^3}\right)\right)
\right)=\\

\color{purple}{
\sum_{i=0}^{m-1}\eta_i
\left(
\frac{2\left(\log\left(u_i\right) - \log(\mu_i)\right)}{\sigma^3}\frac{1}{\mu_i}u_{i-j}\right)}
$$

Consistency check:
$$
\frac{\part^2}{\part\theta_j\part\sigma}\mathcal{L}(u_{t-n:t}|\sigma,\mathbf{\eta}_t,\theta_t) =\\

\frac{\part}{\part\sigma}
\left(
-\sum_{i=0}^{m-1}\eta_i
\left(
\frac{1}{\sigma^2}\left(\log\left(u_i\right)-\log(\mu_i)\right)\frac{\part\mu_i}{\part\theta_j}\right)
\right)
= \\

\color{purple}{
\sum_{i=0}^{m-1}\eta_i
\left(
\frac{2\left(\log\left(u_i\right) - \log(\mu_i)\right)}{\sigma^3}\frac{1}{\mu_i}u_{i-j}\right)}
$$
We have that $\frac{\part}{\part \theta_j\part \sigma}\mathcal{L}(\cdot) = \frac{\part}{\part \sigma\part \theta_j}\mathcal{L}(\cdot)$ as *expected*.
$$
\frac{\part}{\part\theta_j\theta_k}\mathcal{L}(u_{t-n:t}|\sigma,\mathbf{\eta}_t,\theta_t) =\\

\frac{\part}{\part\mu}
\left(
-\sum_{i=0}^{m-1}\eta_i
\left(
\frac{1}{\sigma^2}\left(\log\left(u_i\right)-\log(\mu_i)\right)\frac{1}{\mu_i}u_{i-j}\right)
\right)
\frac{\part \mu}{\part \theta_k}
= \\



-\sum_{i=0}^{m-1}
\eta_i
\frac{u_{i-j}}{\sigma^2}
\frac{\part }{\part \mu}
\left(
\frac{\log(u_i)-\log(\mu_i)}{\mu_i}
\right)
\frac{\part \mu}{\part \theta_k}
= \\

-\sum_{i=0}^{m-1}
\eta_i
\frac{u_{i-j}}{\sigma^2}
\left(
\frac{
- \frac{1}{\mu_i}\mu_i
- 
(\log(u_i) -\log(\mu_i))
}{\mu_i^2}
\right)

\frac{\part \mu}{\part \theta_k}
= \\

\color{purple}{
\sum_{i=0}^{m-1}
\eta_i
\left(
\frac{
1
+
(\log(u_i) -\log(\mu_i))
}{\mu_i^2\sigma^2}
\right)
u_{i-j}
u_{i-k}
}
$$





### Right Censoring

The additional term we add to our cost function when we apply right censoring is:
$$
- \color{green}{\eta_n\log\left(1-P(t-u_{n})\right)}
$$


$$
t_n = t-u_n\\
P(t) = \int_{0}^{t}p(v)dv\\
\mathcal{L}(t)=-\eta_n\log\left(1-P(t)\right)
$$

$$
\frac{\part}{\part \sigma}
\left(-\eta_n\log\left(1-P(t_n)\right)\right) = \frac{\part \left(-\eta_n\log\left(1-P(t_n)\right)\right)}{\part \left(1-P(t_n)\right)}
\cdot
\frac{\part \left(1-P(t_n)\right)}{\part P(t_n)}
\cdot
\frac{\part P(t_n)}{\part\sigma}=\\
-\frac{\eta_n}{1-P(t_n)}\cdot(-1)\cdot\frac{\part P(t_n)}{\part \sigma}=\\

\frac{\eta_n}{1-P(t_n)}\cdot\frac{\part P(t_n)}{\part \sigma}
$$

$$
\frac{\part}{\part \theta_j}
\left(-\eta_n\log\left(1-P(t_n)\right)\right)= \frac{\eta_n}{1-P(t_n)}\cdot\frac{\part P(t_n)}{\part \theta_j}
$$

Given that $\phi$ is the standard normal distribution pdf,  we know that:
$$
\frac{\part P(t_n)}{\part \sigma} = \phi\left(\frac{\log t_n-\log(\mu)}{\sigma}\right)\frac{\part}{\part \sigma}\left(\frac{\log t_n-\log(\mu)}{\sigma}\right) = \\

-\phi\left(\frac{\log t_n-\log(\mu)}{\sigma}\right)\left(\frac{\log t_n-\log(\mu)}{\sigma^2}\right)
$$

$$
\frac{\part P(t_n)}{\part \mu} = \phi\left(\frac{\log t_n-\log(\mu)}{\sigma}\right)\frac{\part}{\part \mu}\left(\frac{\log t_n-\log(\mu)}{\sigma}\right) = \\

-\phi\left(\frac{\log t_n-\log(\mu)}{\sigma}\right)\left(\frac{1}{\mu\sigma}\right)
$$

Which yields:
$$
\frac{\part}{\part \sigma}
\left(-\eta_n\log\left(1-P(t_n)\right)\right) = 
\color{green}{
-\frac{\eta_n}{1-P(t_n)}
\phi\left(\frac{\log t_n-\log(\mu)}{\sigma}\right)\left(\frac{\log t_n-\log(\mu)}{\sigma^2}\right)}
$$

$$
\frac{\part}{\part \mu}
\left(-\eta_n\log\left(1-P(t_n)\right)\right) =
\color{green}{
-\frac{\eta_n}{1-P(t_n)}
\phi\left(\frac{\log t_n-\log(\mu)}{\sigma}\right)\left(\frac{1}{\mu\sigma}\right)
}
$$

The hessian is approximated.