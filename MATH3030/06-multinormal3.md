# Part III:  Inference using the MVN {-}


Part III of these lecture notes covers statistical inference based on the multivariate normal (MVN) distribution, which is of relevance when there are several measurements per experimental unit and observations consist of random vectors.

Chapter \@ref(#multinormal) focuses on classical distribution theory relating to the MVN distribution, including the Wishart distribution, which is defined on the set of symmetric positive definite matrices and is a natural generalisation of the $\chi^2$ distribution.  Another important distribution related to the MVN distribution is the Hotelling $T^2$ distribution, which is a multivariate analogue of the Student $t$-distribution.

Chapter 7 is concerned with testing hypotheses concerning vector means in $1$-sample and $2$-sample settings.  There is a close connection with the classical $1$-sample and $2$-sample $t$-tests in univariate statistics, but here we are dealing with random vectors rather than random variables.

Chapter 8 is concerned with the multivariate linear model, in which the responses consist of random vectors rather than single random variables.  Errors in this setting take the form of random vectors.

The results in Part III turn out to be natural but non-trivial generalisations of the results in the univariate case.

# Multivariate Normal Distribution Theory {#multinormal}



The multivariate normal distribution (MVN) is important for a number of reasons:

1. It is a straightforward generalisation of the univariate normal distribution.
2.  It is entirely defined by its mean vector $\bmu$ and its covariance matrix $\bSigma$.
3. Zero correlation implies independence.
4. Linear functions of multivariate normal vectors are also multivariate normal vectors.
5. There is a multivariate version of the Central Limit Theorem.
6. It has simple geometric properties.


## Definition and Properties of the MVN

```{definition, mvn} \quad A random vector $\bx=(x_1, \ldots , x_p)^\top$ has a $p$-dimensional MVN distribution if and only if $\ba^\top \bx$ is univariate normal for all fixed $p \times 1$ vectors $\ba$.
```

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:six1"><strong>(\#prp:six1) </strong></span>If $\stackrel{p \times 1}{\bx}$  is MVN then for each constant matrix $\bA$ ($q \times p$) and constant vector $\bc$ ($q \times 1$), $\by = \bA \bx + \bc$ has a $q$-dimensional MVN.</div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}Let $\bb$ ($q \times 1)$ be a fixed vector. Then
$$ \bb^\top \by = \bb^\top \bA \bx + \bb^\top \bc = \ba^\top \bx + \bb^\top \bc $$
where $\ba^\top = \bb^\top \bA$.  Now $\ba^\top \bx$ is univariate normal for all $\ba$ since $\bx$ is MVN.  Therefore $\bb^\top \by$ is univariate normal for all $\bb$, so $\by$ is MVN.</div>\EndKnitrBlock{proof}


\BeginKnitrBlock{corollary}<div class="corollary"><span class="corollary" id="cor:csix1"><strong>(\#cor:csix1) </strong></span>Any subset of the components of a MVN vector $\bx$ is also MVN.</div>\EndKnitrBlock{corollary}

**Notation** \quad If $\bx$ ($p \times 1$) is MVN with mean $\bmu$ and covariance matrix $\bSigma$ then we can write
$$ \bx \sim N_p (\bmu, \bSigma).$$

This notation is a direct extension of the univariate notation where $p=1$.









If the population covariance matrix $\bSigma$ ($p \times p$) is positive definite, so that $\bSigma^{-1}$ exists,
then the  **probability density function** (pdf) of the MVN distribution is given by
$$ f(\bx) = \frac{1}{| 2 \pi \bSigma |^{1/2}} \exp \lb -\frac{1}{2}(\bx - \bmu)^\top \bSigma^{-1} (\bx - \bmu) \rb.$$

If $p=1$, so that $\bx = x$, $\bmu = \mu$ and $\bSigma = \sigma^2$, say, then the pdf simplifies to
\begin{eqnarray*}
f(\bx) &=& \frac{1}{|2 \pi \sigma^2|^{1/2}} \exp \lb -\frac{1}{2}(x - \mu) (\sigma^2)^{-1} (x - \mu) \rb \\
&=& \frac{1}{(2 \pi \sigma^2)^{1/2}} \exp \lb -\frac{1}{2 \sigma^2}(x - \mu)^2 \rb
\end{eqnarray*}
which is the familiar pdf of the univariate normal distribution $N(\mu,\sigma^2)$.

If $p>1$ and $\bSigma = \bI_p$ then
\begin{eqnarray*}
f(\bx) &=& \frac{1}{(2 \pi)^{p/2}} \exp \lb -\frac{1}{2}(\bx - \bmu)^\top (\bx - \bmu) \rb \\
&=& \frac{1}{(2 \pi)^{p/2}} \exp \lb -\frac{1}{2} \sum_{i=1}^p (x_i - \mu_i)^2 \rb \\
&=& \lb \frac{1}{\sqrt{2 \pi}} \exp \lb -\frac{1}{2} (x_1 - \mu_1)^2 \rb \rb \\
 && \qquad \qquad \times \ldots \lb \frac{1}{\sqrt{2 \pi}} \exp \lb -\frac{1}{2} (x_p - \mu_p)^2 \rb \rb
\end{eqnarray*}
This means that, by the factorisation theorem for probability densities, the components of $\bx$ have independent univariate normal distributions.

If $p=2$ we can plot $f(\bx)$ using contour plots, where each contour shows values of $\bx$ for which $f(\bx) = c$ for some constant, $c > 0$.  Four examples are shown below for $c = 0.02, 0.04$.

FIX
<!--
%\begin{center}
%$\bmu = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \quad \bSigma = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} \qquad
%\bmu = \begin{pmatrix} 1 \\ 1 \end{pmatrix} \quad \bSigma = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$
%\includegraphics[width=6cm,angle=270]{mvncontours1.ps}
%\includegraphics[width=6cm,angle=270]{mvncontours2.ps} \\
%\includegraphics[width=6cm,angle=270]{mvncontours3.ps}
%\includegraphics[width=6cm,angle=270]{mvncontours4.ps}
%$\bmu = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \quad \bSigma = \begin{pmatrix} 1 & -1 \\ -1 & 2 \end{pmatrix} \qquad
%\bmu = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \quad \bSigma = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$
%\end{center}

\begin{center}
\includegraphics[width=6cm,angle=0]{mvncontours1.pdf}\hfill
\includegraphics[width=6cm,angle=0]{mvncontours2.pdf}
$\bmu = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \quad \bSigma = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} \qquad
\bmu = \begin{pmatrix} 1 \\ 1 \end{pmatrix} \quad \bSigma = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$
\end{center}

%\newpage
\begin{center}
\includegraphics[width=6cm,angle=0]{mvncontours3.pdf}\hfill
\includegraphics[width=6cm,angle=0]{mvncontours4.pdf}
$\bmu = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \quad \bSigma = \begin{pmatrix} 1 & -1 \\ -1 & 2 \end{pmatrix} \qquad
\bmu = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \quad \bSigma = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$
\end{center}
-->

## Transformations

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:six2"><strong>(\#prp:six2) </strong></span> If $\bx \sim N_p(\bmu,\bSigma)$ and $\by = \bA \bx + \bc$, where $\bA$ ($q \times p$) and $\bc$ ($q \times 1$) are constant, then
$$\by \sim N_q(\bA \bmu + \bc, \bA \bSigma \bA^\top).$$</div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}We know $\by$ is MVN by Proposition \@ref(prp:six1).  We also know $E(\by)$ and $\text{Var}(\by)$ from results (1) and (3) of \S 1.4.</div>\EndKnitrBlock{proof}

The above result implies that a linear transformation of a MVN random variable is also MVN.  We can use this result to prove two important corollaries.  The first corollary is useful for simulating data from a general MVN distribution.

\BeginKnitrBlock{corollary}<div class="corollary"><span class="corollary" id="cor:csix2"><strong>(\#cor:csix2) </strong></span> If $\bx \sim N_p(\bzero,\bI_p)$ and $\by = \bSigma^{1/2} \bx + \bmu$ then $\by \sim N_p(\bmu,\bSigma)$.</div>\EndKnitrBlock{corollary}

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}We apply \@ref(prp:six2)  with $\bA = \bSigma^{1/2}$ and $\bc = \bmu$.  Therefore $E(\by) = \bSigma^{1/2} \bzero_p + \bmu = \bmu$ and $\text{Var}(\by) = \bSigma^{1/2} \bI_p \bSigma^{1/2} = \bSigma$.</div>\EndKnitrBlock{proof}


The second corollary says that any MVN random variable can be transformed into standard form.

\BeginKnitrBlock{corollary}<div class="corollary"><span class="corollary" id="cor:csix3"><strong>(\#cor:csix3) </strong></span>If $\bx \sim N_p(\bmu,\bSigma)$, $\bSigma$ has full rank and we define $\by = \bSigma^{-1/2}(\bx - \bmu)$ then $\by \sim N_p(\bzero,\bI_p)$.</div>\EndKnitrBlock{corollary}

```{proof} Apply Proposition \@ref(prp:six2) with $\bA = \bSigma^{-1/2}$ and $\bc = - \bSigma^{-1/2} \bmu$.  Then $E(\by) = \bSigma^{-1/2} \bmu - \bSigma^{-1/2} \bmu = \bzero_p$ and $\text{Var}(\by) = \bSigma^{-1/2} \bSigma \bSigma^{-1/2} = \bI_p$. 
```

The moment generating function of a random vector $\bx$ ($p \times 1$) is given by
$$
M({\mathbf t})=E[e^{{\mathbf t}^\top \bx}],
$$
and is defined for all $t \in \mathbb{R}^p$ for which $M({\mathbf t})$ is finite.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:six3"><strong>(\#prp:six3) </strong></span> The moment generating function of $\bx \sim N_p(\bmu , \bSigma)$ is given by
\begin{equation}
M({\mathbf t})=\exp \left (\bmu^\top  {\mathbf t} + \frac{1}{2} {\mathbf t}^\top \bSigma {\mathbf t} \right).
(\#Mt)
\end{equation}</div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}For fixed $\mathbf t$, define the random variable $Y=\bx^\top {\mathbf t}$.  From Proposition \@ref(prp:six2), $Y \sim N(\mu_{\mathbf t}, \sigma_{\mathbf t}^2)$, where $\mu_{\mathbf t}=\bmu^\top {\mathbf t}$ and $\sigma_{\mathbf t}^2={\mathbf t}^\top \bSigma {\mathbf t}$.  If $\sigma_{\mathbf t}\equiv {\mathbf t}^\top \bSigma {\mathbf t}=0$ then $Y=\bmu^\top {\mathbf t}$ with probability one, and then $M({\mathbf t})=e^{\bmu^\top {\mathbf t}}$ which agrees with \@ref(eq:Mt).  So from now on we assume $\sigma_{\mathbf t}>0$.   Then
\begin{align*}
M({\mathbf t})&=E[e^{\bx^\top {\mathbf t}}]\\
&=E[e^{Y}]=\int_{-\infty}^\infty \exp(y) \frac{1}{\sqrt{2\pi \sigma_{\mathbf t}^2}}
\exp\left (-\frac{1}{2}\frac{(y-\mu_{\mathbf t})^2}{\sigma_{\mathbf t}^2} \right )dy.
\end{align*}
The integral above can be evaluated by completing the square in the exponent, using the identity
$$
y-\frac{1}{2}\frac{(y-\mu_{\mathbf t})^2}{\sigma_{\mathbf t}^2}=\mu_{\mathbf t}
+\frac{1}{2}\sigma_{\mathbf t}^2-\frac{1}{2}\frac{(y-\mu_{\mathbf t}-\sigma_{\mathbf t}^2)^2}{\sigma_{\mathbf t}^2}.
$$
Consequently
\begin{align*}
M({\mathbf t})&=\int_{-\infty}^\infty \exp \left \{\mu_{\mathbf t} +\frac{1}{2}\sigma_{\mathbf t}^2 \right \}
\frac{1}{\sqrt{2 \pi \sigma_{\mathbf t}^2}}\exp \left \{ -\frac{1}{2} \frac{(y-\mu_{\mathbf t}-\sigma_{\mathbf t}^2)^2}
{\sigma_{\mathbf t}^2}\right \}dy\\
&=\exp\left ( \mu_{\mathbf t} + \frac{1}{2}\sigma_{\mathbf t}^2  \right )\\
&=\exp\left (
\bmu^\top {\mathbf t} + \frac{1}{2}{\mathbf t}^\top \bSigma {\mathbf t}\right ),
\end{align*}
as required.</div>\EndKnitrBlock{proof}

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:six4"><strong>(\#prp:six4) </strong></span>Two vectors $\bx$ ($p \times 1$) and $\by$ ($q \times 1$) which are jointly multivariate normal are independent if and only if they are uncorrelated (i.e. $\text{Cov}(\bx,\by) = \bzero_{p,q}$).</div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}We prove this result using the factorisation theorem for moment generating functions (MGFs), which is now stated.
Let ${\mathbf t}=({\mathbf t}_1^\top , {\mathbf t}_2^\top)^\top$ where ${\mathbf t}_1 \in \mathbb{R}^p$, ${\mathbf t}_2 \in \mathbb{R}^q$ and ${\mathbf t} \in \mathbb{R}^{p+q}$.  The joint MGF of two arbitrary random vectors $\stackrel{p \times 1}{\bx}$ and $\stackrel{q \times 1}{\by}$ is defined by
  $$
 M({\mathbf t}_1, {\mathbf t}_2)=E[e^{{\mathbf t}_1^\top \bx + {\mathbf t}_2^\top \by}],
  $$
  for all ${\mathbf t}=({\mathbf t}_1^\top , {\mathbf t}_2^\top )^\top$ at which $M({\mathbf t}_1, {\mathbf t}_2)$ is finite.
The factorisation theorem for MGFs states that
$\bx$ and $\by$ are independent if an only if $M({\mathbf t}_1 , {\mathbf t}_2)$ factorises, i.e.
$$
M({\mathbf t}_1 , {\mathbf t}_2)=M_1({\mathbf t}_1)M_2({\mathbf t}_2)
$$
for some functions $M_1$ and $M_2$, in which case $M_1$ and $M_2$ are the marginal MGFs of $\bx$ and $\by$.  Now we focus on the MVN case.  Suppose
\begin{equation}
E[\bx]=\bmu_{\bx}, \qquad \qquad E[\by]=\bmu_{\by}, \quad  \text{Var}(\bx)=\bSigma_{\bx \bx},
\quad  \text{Var}(\by)=\bSigma_{\by \by},
(\#def1)
\end{equation}
 and
\begin{equation}
\text{Cov}(\bx,\by)=\bSigma_{\bx \by}=\bSigma_{\by \bx}^\top = \text{Cov}(\by, \bx)^\top.
(\#def2)
\end{equation}
Using Proposition \@ref(prp:six3) and definitions \@ref(eq:def1) and \@ref(eq:def2),
\begin{align*}
M({\mathbf t}_1, {\mathbf t}_2)&=\exp\left ( \bmu^\top {\mathbf t} + \frac{1}{2}{\mathbf t}^\top \bSigma {\mathbf t} \right )\\
&=\exp\bigg (\bmu_{\bx}^\top {\mathbf t}_1 +\bmu_{\by}^\top {\mathbf t}_2+\frac{1}{2}{\mathbf t}_1^\top \bSigma_{\bx \bx}{\mathbf t_1}\\
& \qquad \qquad +\frac{1}{2}{\mathbf t}_2^\top  \bSigma_{\by \by}{\mathbf t}_2+\frac{1}{2} 2{\mathbf t}_1^\top \bSigma_{\bx \by}{\mathbf t}_2 \bigg)\\
&=M_1({\mathbf t}_1)M_2({\mathbf t}_2)M_3({\mathbf t}_1, {\mathbf t}_2),
\end{align*}
where $M_1({\mathbf t}_1)$ and $M_2({\mathbf t}_2)$ are the marginal MGFs of $\bx$ and $\by$ respectively, and
$$
M_3({\mathbf t}_1, {\mathbf t}_2)=\exp\left ({\mathbf t}_1^\top \bSigma_{\bx \by}{\mathbf t}_2 \right ).
$$
The factorisation theorem holds if an only if $M_3({\mathbf t}_1, {\mathbf t}_2)$ is constant with respect to
${\mathbf t}_1$ and ${\mathbf t}_2$, which is the case if and only if $\bSigma_{\bx \by}={\mathbf 0}_{p,q}$. </div>\EndKnitrBlock{proof}
 
 In words: Proposition \@ref(prp:six4) means that zero correlation implies independence for the MVN distribution.  This is not generally true for other distributions.

**Note** that  Propositions \@ref(prp:six1)  -   \@ref(prp:six4) each holds regardless of whether the covariance matrix $\bSigma$ is positive definite or not.

The term $(\bx-\bmu)^\top \bSigma^{-1} (\bx-\bmu)$ appears in the exponent of the pdf and we derive its distribution in  Proposition \@ref(prp:six5).

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:six5"><strong>(\#prp:six5) </strong></span>If $\bx \sim N_p(\bmu, \bSigma)$ and $\bSigma$ is positive definite then
$$(\bx-\bmu)^\top \bSigma^{-1} (\bx-\bmu) \sim \chi_p^2.$$</div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}Define $\by = \bSigma^{-1/2} (\bx-\bmu)$ so
\begin{align*}
(\bx-\bmu)^\top \bSigma^{-1} (\bx-\bmu) &= \lb \bSigma^{-1/2} (\bx-\bmu) \rb^\top \lb \bSigma^{-1/2} (\bx-\bmu) \rb \\
&= \by^\top \by = \sum_{i=1}^p y_i^2
\end{align*}
By Corollary \@ref(cor:csix3), $\by \sim N_p (\bzero, \bI_p)$, and so the components of $\by$ have independent univariate normal distributions with mean 0 and variance 1.  Recall from univariate statistics that if $z \sim N(0,1)$ then $z^2 \sim \chi^2_1$ and if $z_1, \ldots, z_n$ are iid $N(0,1)$ then $\sum_{i=1}^n z_i^2 \sim \chi_n^2$.  It therefore follows that $\sum_{i=1}^p y_i^2 \sim \chi^2_p$. </div>\EndKnitrBlock{proof}

We saw earlier in this chapter that the MVN distribution in $p$ dimensions has constant density on ellipses or ellipsoids given by $f(\bx) = c$ for some constant $c > 0$.  We can rearrange this equation to be of the form
$$U(\bx) = (\bx-\bmu)^\top \bSigma^{-1} (\bx-\bmu) = k$$
where $k = - 2 \log(c) - \log |2 \pi \bSigma| > 0$ is a combination of the constant, $c$,  and the normalising constant in the pdf.   Proposition \@ref(prp:six5) means we can calculate the probability, $P(U(\bx)<k)$, which is the probability of $\bx$ lying within a particular ellipsoid.

## Two important results for the MVN

In this section we present two important results which are natural
generalisations of what happens in the univariate case.


\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:six6"><strong>(\#prp:six6) </strong></span>If $\bx_1, \ldots, \bx_n$ is an IID random sample from $N_p(\bmu, \bSigma)$, then the sample mean $\bar{\bx} = \frac{1}{n} \sum_{i=1}^n \bx_i$ and the sample variance matrix $\bS = \frac{1}{n} \sum_{i=1}^n (\bx_i - \bar{\bx})(\bx_i-\bar{\bx})^\top$ are independent. </div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}Define $\bx =\bar{\bx}$ and $\by_i=\bx_i-\bar{\bx}$, $i=1, \ldots , n$. From  Proposition \@ref(prp:six1) and  Proposition \@ref(prp:six2) we can see that if $\bx_1, \ldots, \bx_n$ is a random sample from $N_p(\bmu, \bSigma)$ then $\bar{\bx} \sim N_p (\bmu, n^{-1}\bSigma)$.  Then
\begin{align*}
\text{Cov}(\bar{\bx},\by_i)&=\text{Cov}(\bar{\bx}, \bx_i -\bar{\bx})\\
&=\text{Cov}(\bar{\bx}, \bx_i) - \text{Cov}(\bar{\bx}, \bar{\bx})\\
&=n^{-1}\sum_{j=1}^n \left \{E[(\bx_j -\bmu)(\bx_i-\bmu)^\top]\right \}\\
& \qquad \qquad -E[(\bar{\bx}-\bmu)(\bar{\bx}-\bmu)^\top]\\
&=n^{-1}\bSigma - n^{-1}\bSigma= {\mathbf 0}_{p,p}.
\end{align*}
Now define the $(np) \times 1$ vector $\by=(\by_1^\top, \ldots , \by_n^\top)^\top$.  Then $\text{Cov}(\bx,\by)={\mathbf 0}_{p, np}$,
so we may apply  Proposition \@ref(prp:six4) to conclude that $\bar{\bx}$ and $\by$ are independent.  Therefore $\bar{\bx}$ and
$$
\bS=\frac{1}{n}\sum_{i=1}^n \by_i \by_i^\top =n^{-1}\sum_{i=1}^n (\bx_i -\bar{\bx})(\bx_i -\bar{\bx})^\top
$$
are independent, becasue $\bS$ is a function of $\by$ alone.  </div>\EndKnitrBlock{proof}

 Recall from above that if $\bx_1, \ldots, \bx_n$ is a random sample from $N_p(\bmu, \bSigma)$ then $\bar{\bx} \sim N_p (\bmu, \frac{1}{n}\bSigma)$.  This result is also approximately true for a large random sample from a non-normal population, as is now stated in the multivariate central limit theorem.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:clt"><strong>(\#prp:clt) </strong></span>**Central limit theorem**  Let $\bx_1, \bx_2, \ldots$ be a sample of independent and identically distributed random vectors from a distribution with mean vector $\bmu$ and finite variance matrix $\bSigma$. Then asymptotically as $n \rightarrow \infty$, $\sqrt{n}(\bar{\bx}-\bmu)$ converges in distribution to $N_p ({\mathbf 0}_p, \bSigma )$.</div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}Beyond the scope of this module. </div>\EndKnitrBlock{proof}



## The Wishart distribution

The Wishart distribution is a multivariate generalisation of the univariate $\chi^2$ distribution.  In univariate statistics the $\chi^2$ distribution plays an important role in inference related to the univariate normal, e.g. in the definition of Student's $t$ distribution.  An analogous role is played by the Wishart distribution in multivariate statistics.

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:wishart"><strong>(\#def:wishart) </strong></span> Let $\bx_1, \ldots, \bx_n$ be an IID random sample from $N_p (\bzero, \bSigma)$. Then $\bM = \sum_{i=1}^n \bx_i \bx_i^\top$
is said to have a Wishart distribution with $n$ degrees of freedom and scale matrix $\bSigma$.  We write this as $\bM \sim W_p(\bSigma, n)$.</div>\EndKnitrBlock{definition}

Note that $W_p(\bSigma,n)$ is a probability distribution on the set of $p \times p$ symmetric non-negative definite random matrices.

We call $W_p(\bI_p,n)$ a standard Wishart distribution.

When $p=1$, $W_1(1,n)$ is the $\chi_n^2$ distribution and $W_1(\sigma^2,n)$ is the $\sigma^2 \chi_n^2$ distribution.  This claim follows from \@ref(prp:six9) below.



We now use the definition of $W_p(\bSigma, n)$ to prove some important\\
 results.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:six8"><strong>(\#prp:six8) </strong></span>If $\bM \sim W_p(\bSigma,n)$ and $\bA$ is a fixed $q \times p$ matrix, then
$$ \bA \bM \bA^\top \sim W_q \lb \bA \bSigma \bA^\top, n \rb.$$</div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}From the definition, let $\bM = \sum_{i=1}^n \bx_i \bx_i^\top$, where the $\bx_i$ are
IID $N_p(\bzero,\bSigma)$.   Then
\begin{align*}
\bA \bM \bA^\top &= \bA \lb \sum_{i=1}^n \bx_i \bx_i^\top \rb \bA^\top\\
 &= \sum_{i=1}^n (\bA \bx_i)(\bA \bx_i)^\top = \sum_{i=1}^n \by_i \by_i^\top
\end{align*}
where $\by_i = \bA \bx_i \sim N_q(\bzero,\bA \bSigma \bA^\top)$, by  Proposition \@ref(prp:six2).  Now we apply the definition of the Wishart distribution to $\by_1,\ldots,\by_n$ and, hence, $\sum_{i=1}^n \by_i \by_i^\top \sim W_q\lb \bA \bSigma \bA^\top, n \rb$. </div>\EndKnitrBlock{proof}

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:six9"><strong>(\#prp:six9) </strong></span>If $\bM \sim W_p(\bSigma,n)$ and $\ba$ is a fixed $p \times 1$ vector then
$$ \ba^\top \bM \ba \sim \lb \ba^\top \bSigma \ba \rb \chi_n^2.$$</div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}Apply Proposition \@ref(prp:six8) with $\bA = \ba^\top$ then $\ba^\top \bM \ba \sim W_1( \ba^\top \bSigma \ba, n)$. But $W_1(1,n)$ is equal in distribution to $\sum_{i=1}^n z_i^2$ where the $z_i$ are IID $N(0,1)$, and so has $\chi_n^2$ distribution.  Moreover,  using Proposition \@ref(prp:six8)  with $p=q=1$ and $A=\sigma$, it is seen that $W_1(\sigma^2,n)$ is equal in distribution to $\sigma^2 \chi_n^2$. </div>\EndKnitrBlock{proof}

Note that an alternative form of the above result is
$$\frac{ \ba^\top \bM \ba }{ \ba^\top \bSigma \ba } \sim \chi_n^2.$$

\BeginKnitrBlock{corollary}<div class="corollary"><span class="corollary" id="cor:csix4"><strong>(\#cor:csix4) </strong></span>Let $m_{ii}$ and $\sigma_{ii}$ be the $i$th diagonal entry for $\bM$ and $\bSigma$ respectively, then $m_{ii} \sim \sigma_{ii} \chi^2_n$ for $i=1,\ldots,p$.</div>\EndKnitrBlock{corollary}

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}Let $\ba = (a_1,\ldots,a_p)^\top$ where $a_j = 1$ if $j=i$ and $a_j = 0$ otherwise. Then $\ba^\top \bM \ba = m_{ii}$ and $\ba^\top \bSigma \ba = \sigma_{ii}$.  Now apply Proposition \@ref(prp:six9)  </div>\EndKnitrBlock{proof}

Note, however, that the $m_{ii}$, $i=1,\ldots,p$, are not, in general, independent.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:six10"><strong>(\#prp:six10) </strong></span>If $\bM_1 \sim W_p(\bSigma,n_1)$ and $\bM_2 \sim W_p(\bSigma,n_2)$ are independent
then
$$\bM_1 + \bM_2 \sim W_p(\bSigma,n_1 + n_2).$$</div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}From the definition, let $\bM_1 = \sum_{i=1}^{n_1} \bx_i \bx_i^\top$ and let $\bM_2 = \sum_{i=n_1+1}^{n_1+n_2} \bx_i \bx_i^\top$, where $\bx_i \sim N_p(\bzero,\bSigma)$, then $\bM_1+\bM_2 = \sum_{i=1}^{n_1+n_2} \bx_i \bx_i^\top \sim W_p(\bSigma,n_1 + n_2)$ by the definition of the Wishart distribution. </div>\EndKnitrBlock{proof}

Our next result is known as Cochran's theorem.  Recall the definition of projection matrices at the end of \S 2.4.

\BeginKnitrBlock{theorem}<div class="theorem"><span class="theorem" id="thm:six11"><strong>(\#thm:six11) </strong></span>**(Cochran's Theorem)**  Suppose $\stackrel{n \times n}{\mathbf P}$ is a projection matrix of rank $r$.  Assume that $\bX$ is an $n \times p$ data matrix with IID rows that have a common $N_p({\mathbf 0}_p, \bSigma)$ distribution, where $\Sigma$ has full rank $p$, and note the identity
\begin{equation}
\bX^\top \bX = \bX^\top {\mathbf P} \bX + \bX^\top ({\mathbf I}_n -{\mathbf P})\bX.
(\#Cochran1)
\end{equation}
Then
\begin{equation}
\bX^\top {\mathbf P} \bX \sim W_p(\bSigma, r), \qquad  \bX^\top ({\mathbf I}_n -{\mathbf P})\bX \sim W_p(\bSigma, n-r),
(\#Cochran2)
\end{equation}
and $\bX^\top {\mathbf P} \bX$ and $\bX^\top ({\mathbf I}_n -{\mathbf P})\bX$
are independent.</div>\EndKnitrBlock{theorem}


\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}We first of all prove the result in the particular case $\bSigma = {\mathbf I}_p$ and
then  consider the general case.  Using the Spectral Decomposition Theorem \@ref(prp:spectraldecomp) and noting Proposition \@ref(prp:two1), we may write
$$
{\mathbf P}=\sum_{j=1}^r \bq_j \bq_j^\top \qquad  \hbox{and} \qquad 
(\bI_n-{\mathbf P})=\sum_{j=r+1}^n \bq_j \bq_j^\top
$$
where $\bq_1, \ldots , \bq_n$ are mutually orthogonal unit vectors.  Then
\begin{align}
\bX^\top \bP \bX &=  \bX^\top \left (\sum_{j=1}^r \bq_j \bq_j^\top \right) \bX \nonumber \\
& =\sum_{j=1}^r \bX^\top \bq_j \bq_j^\top \bX =\sum_{j=1}^r \by _j \by_j^\top,
(\#Prep)
\end{align}
and similarly,
\begin{align}
\bX^\top (\bI_n -\bP) \bX &=  \bX^\top \left (\sum_{j=r+1}^n \bq_j \bq_j^\top \right) \bX \nonumber \\
& =\sum_{j=r+1}^n \bX^\top \bq_j \bq_j^\top \bX =\sum_{j=r+1}^n \by _j \by_j^\top,
(\#I-Prep)
\end{align}
where, for $j=1, \ldots , n$,  $\by_j=\bX^\top \bq_j$  is a $p \times 1$ vector.  We shall now prove that the $\by_j$ are
IID $N_p({\mathbf 0}_p, \bI_p)$. Write $\bX=[\bx_{[1]}, \ldots , \bx_{[p]}]$, where $\bx_{[u]}$ is column $u$ of
$\bX$.  Then $\bx_{[1]}, \ldots , \bx_{[p]}$ are IID $N_n({\mathbf 0}_n ,\bI_n)$.  Moreover,
$$
\by_j=\bX^\top \bq_j = \left [  \begin{array}{c}
\bq_j^\top \bx_{[1]}\\
\bq_j^\top \bx_{[2]}\\
..\\
..\\
..\\
\bq_j^\top \bx_{[p]}
\end{array} \right ].
$$
But
\begin{align}
E[\bq_j^\top \bx_{[u]} \bq_k^\top \bx_{[v]}] &=E[\bq_j^\top \bx_{[u]} \bx_{[v]}^\top \bq_k]\nonumber \\
&=\bq_j^\top E[\bx_{[u]} \bx_{[v]}^\top ]\bq_k \nonumber\\
&=\bq_j^\top \left (\delta_{uv}\bI_n \right)\bq_k \nonumber\\
&=\bq_j^\top \bq_k \delta_{uv}\nonumber\\
&=\delta_{jk}\delta_{uv},
(\#jkuv)
\end{align}
where $\delta$ is the Kronecker $\delta$ defined by
$$
\delta_{ab}=\begin{cases} 0 &\text{if} \quad  a \neq b\\
1 &\text{if} \quad  a=b \end{cases}.
$$
If follows immediately from \@ref(eq:jkuv) that
$$
\text{Var}(\by_j)=\bI_p \qquad  \text{Cov}(\by_j , \by_k)={\mathbf 0}_{pp} \quad 
\text{if} \quad  j \neq k.
$$
By Proposition \@ref(prp:six4), the $\by_j$, $j=1,\ldots , n$, are IID $N_p({\mathbf 0}_p , \bI_p)$, and therefore it follows from the
definition of the Wishart distribution that, when $\bSigma=\bI_p$,  \@ref(eq:Prep) has a Wishart $W_p(\bI_p,r)$ distribiton, \@ref(eq:I-Prep) has a Wishart $W_p(\bI_p, n-r)$ distrubtion.  Moreover, these random Wishart matrices are independent becasue the $\by_j$ are all independent.

Finally, we consider the case of a general covariance matrix $\bSigma$.  We have proved that \@ref(eq:Cochran1) holds
when $\bSigma=\bI_p$, so pre-multiply both sides by the matrix square root $\bSigma^{1/2}$, and post-multiply both sides by $\bSigma^{1/2}$.  This corresponds to the case where the $\bx_i$ are IID $N_p({\mathbf 0}_p, \bSigma)$.  Then, using Proposition \@ref(prp:six8),
$$
\bSigma^{1/2} W_p(\bI_p, t)\bSigma^{1/2} \stackrel{d}{=} W_p(\bSigma^{1/2} \bSigma^{1/2}, t)
\stackrel{d}{=}W_p(\bSigma,t),
$$
when $t=r$ and $t=n-r$.  Moreover, since $\bSigma^{1/2}$ is a non-random matrix, independence is preserved
when we pre- and post-multiply by $\bSigma^{1/2}$, and the result is proved. </div>\EndKnitrBlock{proof}

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:six12"><strong>(\#prp:six12) </strong></span>If $\bx_1,\ldots,\bx_n$ is an IID sample from $N_p(\bmu,\bSigma)$, then
$$ n \bS = \sum_{i=1}^n (\bx_i - \bar{\bx})(\bx_i - \bar{\bx})^\top \sim W_p(\bSigma,n-1).$$</div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{} Define $\bP= {\mathbf H}\equiv \bI_n - n^{-1}{\mathbf 1}_n {\mathbf 1}_n^\top$ where ${\mathbf 1}_n$ is the $n \times 1$ vector of ones.  Note that $\bH$ is the $n \times n$ centering matrix and, from Property (i) of \S 2.7, $\bH$ is a projection matrix.     Clearly, $\bI_n - \bP=n^{-1} {\mathbf 1}_n {\mathbf 1}_n^\top$ has rank $1$, so $\bH$ has rank $n-1$.   Therefore, using Theorem \@ref(thm:six11) ,
$$
\bX^\top \bH \bX \sim W_p(\bSigma, n-1).
$$
But from Property (vi) in \S 2.7, $\bX^\top \bH \bX =n\bS$,
and consequently, $n\bS \sim   W_p(\bSigma, n-1)$, as required. </div>\EndKnitrBlock{proof}

## Hotelling's $T^2$ distribution

Hotelling's $T^2$ distribution is a multivariate analogue of the Student $t$ distribution.  It plays an important role in multivariate hypothesis testing and confidence region construction, just as the Student $t$ distribution does in the univariate setting.

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:Hotelling"><strong>(\#def:Hotelling) </strong></span>Suppose $\bx \sim N_p(\bzero,\bI_p)$ and $\bM \sim W_p(\bI_p,n)$ are independent, then the quantity $\tau ^2 = n \bx^\top \bM^{-1} \bx$ is said to have Hotelling's $T^2$ distribution with parameters $p$ and $n$.  We write this as $\tau^2 \sim T^2(p,n)$.</div>\EndKnitrBlock{definition}

We can generalise the definition with the following result.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:six13"><strong>(\#prp:six13) </strong></span>Suppose $\bx \sim N_p(\bmu,\bSigma)$ and $\bM \sim W_p(\bSigma,n)$ are independent and
$\bSigma$ has full rank $p$.  Then
$$ n (\bx-\bmu)^\top \bM^{-1} (\bx-\bmu) \sim T^2(p,n). $$</div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}Define $\by = \bSigma^{-1/2}(\bx-\bmu)$.   Then, by Corollary \@ref(cor:csix3), $\by \sim N_p(\bzero,\bI_p)$.   Further, let $\bZ = \bSigma^{-1/2} \bM \bSigma^{-1/2}$ then $\bZ \sim W_p(\bI_p,n)$ by applying \@ref(prp:six8)  with $\bA = \bSigma^{-1/2}$.  From the definition, $n \by^\top \bZ^{-1} \by \sim T^2(p,n)$ and
\begin{eqnarray*}
n \by^\top \bZ^{-1} \by &=& n (\bx-\bmu)^\top \bSigma^{-1/2} \bSigma^{1/2} \bM^{-1} \bSigma^{1/2} \bSigma^{-1/2} (\bx-\bmu) \\
&=& n(\bx-\bmu)^\top \bM^{-1}(\bx-\bmu)
\end{eqnarray*}
so the result is proved. </div>\EndKnitrBlock{proof}

This result gives rise to an important corollary used in hypothesis testing when $\bSigma$ is unknown.


\BeginKnitrBlock{corollary}<div class="corollary"><span class="corollary" id="cor:csix5"><strong>(\#cor:csix5) </strong></span>If $\bar{\bx}$ and $\bS$ are the mean and covariance matrix based on a sample of size $n$ from $N_p(\bmu,\bSigma)$ then
$$ (n-1)(\bar{\bx}-\bmu)^\top \bS^{-1} (\bar{\bx}-\bmu) \sim T^2(p,n-1).$$
</div>\EndKnitrBlock{corollary}

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}We have seen earlier that $\bar{\bx} \sim N_p(\bmu,\frac{1}{n}\bSigma)$. Let $\bx^\ast = n^{1/2} \bar{\bx}$ and let $\bmu^\ast = n^{1/2} \bmu$.  Then $\bx^\ast=n^{1/2} \bar{\bx} \sim N_p(\bmu^\ast, \bSigma)$.

From Proposition \@ref(prp:six12)  we know $n\bS \sim W_p(\bSigma,n-1)$, and from Theorem \@ref(prp:six6)  we know $\bar{\bx}$ and $\bS$ are independent.  Applying Proposition \@ref(prp:six13)  with $\bx = \bx^\ast$ and $\bM = n\bS$ we obtain
$$ (n-1)(\bx^\ast - \bmu^\ast)^\top (n\bS)^{-1} (\bx^\ast - \bmu^\ast) \sim T^2(p,n-1),$$
and given $\bx^\ast - \bmu^\ast = n^{1/2} (\bx-\bmu)$ then
\begin{align*}
&(n-1)(\bx^\ast - \bmu^\ast)^\top (n\bS)^{-1} (\bx^\ast - \bmu^\ast)\\
& \qquad \qquad = (n-1)n^{1/2}(\bar{\bx}-\bmu)^\top n^{-1} \bS^{-1} n^{1/2}(\bar{\bx}-\bmu) \\
&\qquad \qquad = (n-1)(\bar{\bx}-\bmu)^\top \bS^{-1} (\bar{\bx}-\bmu).
\end{align*}</div>\EndKnitrBlock{proof}

Hotelling's $T^2$ distribution is not often included in statistical tables but the next result  tells us that Hotelling's $T^2$ is a scale transformation of an $F$ distribution.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:six14"><strong>(\#prp:six14) </strong></span>If $\tau^2 \sim T^2(p,n)$ then
$$\gamma^2 = \frac{n-p+1}{np} \tau^2 \sim F_{p,n-p+1}.$$</div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{} Beyond the scope of the module.</div>\EndKnitrBlock{proof}

We can apply this result to the previous corollary.

\BeginKnitrBlock{corollary}<div class="corollary"><span class="corollary" id="cor:csix6"><strong>(\#cor:csix6) </strong></span> If $\tau^2 = (n-1)(\bar{\bx}-\bmu)^\top \bS^{-1} (\bar{\bx}-\bmu)$ then
$$ \gamma^2 = \frac{n-p}{p} (\bar{\bx}-\bmu)^\top \bS^{-1} (\bar{\bx}-\bmu) \sim F_{p,n-p}. $$</div>\EndKnitrBlock{corollary}

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}From Corollary \@ref(cor:csix6) we know $\tau^2 \sim T^2(p,n-1)$.  Applying Proposition \@ref(prp:six14) we get
\begin{align*}
&\gamma^2\\
& = \frac{(n-1)-p+1}{(n-1)p}(n-1)(\bar{\bx}-\bmu)^\top \bS^{-1} (\bar{\bx}-\bmu) \sim F_{p,(n-1)-p+1} \\
&= \frac{n-p}{p}(\bar{\bx}-\bmu)^\top \bS^{-1} (\bar{\bx}-\bmu) \sim F_{p,n-p}
\end{align*}</div>\EndKnitrBlock{proof}

