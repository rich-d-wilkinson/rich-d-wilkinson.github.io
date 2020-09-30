# Part IV:  Classification and Clustering {-}


In Part IV, we focus on different methods of classification,
i.e. allocating the observations in a sample to different subsets (or
groups). In Chapter \@ref(lda), we focus on an approach called discriminant
analysis, in which we have a training sample available, and we use
this training sample to set up a suitable classification rule.   An important type of situation where
discriminant analysis is used is in screening tests.  Here, several variables may be measured on each of a number of
individuals, and we want to decide whether each individual is ''negative'', in which case no further investigations are required,
or ''positive'', in which case further tests are required.


In Chapter \@ref(cluster), we consider an alternative approach, known as cluster analysis,
in which we allocate the observations into clusters (or similar subsets).
Here,  a training sample is not available and typically the number of clusters will
not be known in advance.  The idea is to form clusters in such a way that experimental units within clusters are as similar as possible, in a suitable sense,
and experimental units in different clusters are as dissimilar as possible.


# Discriminant analysis - FIX THE FIGURES {#lda}


Consider $g$ populations $\Pi_1, \ldots , \Pi_g$. Each population is described by a pdf $f_j(\bx)$, $\bx \in \mR^p$, $j=1,\ldots , g$.  Let $\bz \in \mR^p$ be a `new' observation assumed to come from one of $\Pi_1, \ldots , \Pi_g$.  The aim of discriminant analysis is to allocate $\bz$ to one of $\Pi_1, \ldots , \Pi_g$ with 'as small a \\
probability of error as possible'.

For example, $\bz$ might contain numerical measures of a person's \\
financial history.  A credit rating agency might then want to \\
classify this customer as ''safe'' or ''risky'' based on knowledge of \\
previous customers.

A **discriminant rule**, $d$, corresponds to a division of $\mR^p$ into disjoint regions $\cR_1, \ldots, \cR_g$, where
$$\bigcup_{j=1}^g \cR_j = \mR^p, \qquad \cR_j \cap \cR_k = \emptyset, j \neq k.$$
The rule $d$ is then defined by
\begin{center}
$d$: allocate $\bz$ to $\Pi_j$ if and only if  $\bz \in \cR_j$.
\end{center}

Three possible approaches to this problem are considered.  We examine the simplest case, where the $f_j(\bx)$ are known exactly, in Section 9.2.  If we have to estimate the parameters of $f_j(\bx)$ then we use the sample version in Section 9.3.  Finally, if we do not know the distribution of the populations, $\Pi_j$, then an alternative approach, described in Section 9.4, may be used.

## Maximum likelihood discriminant rule

Suppose initially that $f_1(\bx), \ldots, f_g(\bx)$ are **known** pdf's.  This is the simplest case but it is an unrealistic assumption in practice unless the sample sizes are *very* large.

\BeginKnitrBlock{definition}
<span class="definition" id="def:maxlik"><strong>(\#def:maxlik) </strong></span>The **maximum likelihood** (ML) discriminant rule
allocates $\bz$ to the population with the largest likelihood at $\bz$, i.e. it allocates $\bz$ to $\Pi_j$ where
$$ f_j(\bz)= \max\limits_{1 \leq k \leq g} f_k(\bz) $$
\EndKnitrBlock{definition}

\BeginKnitrBlock{example}
<span class="example" id="exm:exnine1"><strong>(\#exm:exnine1) </strong></span>Consider the univariate case with $g=2$ where $\Pi_1$ is the $N(\mu_1,\sigma_1^2)$ distribution and $\Pi_2$ is the $N(\mu_2,\sigma_2^2)$ distribution.  The ML discriminant rule allocates $z$ to $\Pi_1$ if and only if
$$
f_1(z) > f_2(z) ,
$$
which is equivalent to
%\begin{align*}
$$
\frac{1}{(2\pi\sigma_1^2)^{1/2}} \exp \lb -\frac{1}{2\sigma_1^2} (z-\mu_1)^2 \rb \\
> \frac{1}{(2\pi\sigma_2^2)^{1/2}} \exp \lb -\frac{1}{2\sigma_2^2} (z-\mu_2)^2 \rb.
$$
%\end{align*}
Collecting terms together on the left hand side (LHS) gives
\begin{eqnarray*}
&&  \qquad \frac{\sigma_2}{\sigma_1} \exp \lb -\frac{1}{2\sigma_1^2} (z - \mu_1)^2 +\frac{1}{2\sigma_2^2} (z - \mu_2)^2 \rb > 1 \\
&\iff& \qquad \log \lb \frac{\sigma_2}{\sigma_1} \rb -\frac{1}{2\sigma_1^2} (z - \mu_1)^2 + \frac{1}{2\sigma_2^2} (z - \mu_2)^2 > 0 \\
& \iff & \qquad z^2 \lb \frac{1}{\sigma_2^2} - \frac{1}{\sigma_1^2} \rb
+ z \lb \frac{2 \mu_1}{\sigma_1^2} - \frac{2 \mu_2}{\sigma_2^2} \rb + \frac{\mu_2^2}{\sigma_2^2} - \frac{\mu_1^2}{\sigma_1^2} + 2 \log \frac{\sigma_2}{\sigma_1} > 0.
\end{eqnarray*}
Suppose, for example, that $\mu_1 = \sigma_1 = 1$ and $\mu_2 = \sigma_2 = 2$, then this reduces to the quadratic expression
$$ -\frac{3}{4}z^2 + z + 2 \log 2 > 0.$$
Suppose that our new observation is $z=0$, say.  Then the LHS is $2 \log 2$ which is greater than zero and so we would allocate $z$ to population 1.

For more general values of $z$ we can solve the quadratic equation to find $z$ such that $f_1(z)=f_2(z)$.  Using the quadratic equation formula we find that
$$z = \frac{-1 \pm \sqrt{1+6 \log 2}}{-3/2} = \frac{2}{3} \pm \frac{2}{3} \sqrt{1 + 6 \log 2}.$$  Hence the solutions are $z = -0.85$ and $z = 2.18$. Our discriminant rule would then be to allocate $z$ to $\Pi_1$ if $-0.85 < z < 2.18$ and allocate it to $\Pi_2$ otherwise.  The situation is illustrated below.

FIX
<!--\begin{center}
\includegraphics[width=12cm,angle=0]{mldr_univariateB.pdf}
\end{center}-->
\EndKnitrBlock{example}
  
Now we consider the case of $g$ multivariate normal populations which, for simplicity, have the same covariance matrix.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:nine1"><strong>(\#prp:nine1) </strong></span>If $\Pi_k$ is the $N_p(\bmu_k,\bSigma)$ population, $k=1,\ldots,g$, then the ML discriminant rule allocates $\bz$ to $\Pi_j$ where $j$ is the value of $k$ which minimises
$$(\bz-\bmu_k)^\top \bSigma^{-1} (\bz-\bmu_k).$$
\EndKnitrBlock{proposition}

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{} The $k$th likelihood is
$$f_k(\bz) = | 2 \pi \bSigma |^{-1/2} \exp \lb -\frac{1}{2} (\bz - \bmu_k)^\top \bSigma^{-1} (\bz - \bmu_k) \rb.$$
This is maximised when the exponent is minimised, due to the minus sign in the exponent and the fact that $\bSigma$ is positive definite. 
\EndKnitrBlock{proof}


\BeginKnitrBlock{corollary}
<span class="corollary" id="cor:nine2c"><strong>(\#cor:nine2c) </strong></span>When $g=2$, the rule allocates $\bz$ to $\Pi_1$ if and only if
$$\ba^\top (\bz - \bh) > 0, $$
where $\ba = \bSigma^{-1} (\bmu_1 - \bmu_2)$ and $\bh = \frac{1}{2} (\bmu_1 + \bmu_2)$.
\EndKnitrBlock{corollary}

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}First, note that
\begin{eqnarray*}
(\bz-\bmu_k)^\top \bSigma^{-1} (\bz-\bmu_k) &=& \bz^\top \bSigma^{-1} \bz + \bmu_k^\top \bSigma^{-1} \bmu_k - \bmu_k^\top \bSigma^{-1} \bz - \bz^\top \bSigma^{-1} \bmu_k \\
&=& \bz^\top \bSigma^{-1} \bz + \bmu_k^\top \bSigma^{-1} \bmu_k - 2\bmu_k^\top \bSigma^{-1} \bz.
\end{eqnarray*}

From Proposition \@ref{prp:nine1} we know that $f_1(\bz) > f_2(\bz)$ if and only if
$$
(\bz-\bmu_1)^\top \bSigma^{-1} (\bz-\bmu_1) < (\bz-\bmu_2)^\top \bSigma^{-1} (\bz-\bmu_2).
$$
Expanding both sides we find find that
\begin{eqnarray*}
&&\bz^\top \bSigma^{-1} \bz + \bmu_1^\top \bSigma^{-1} \bmu_1 - 2\bmu_1^\top \bSigma^{-1} \bz \\
 && \qquad \qquad < \bz^\top \bSigma^{-1} \bz + \bmu_2^\top \bSigma^{-1} \bmu_2 - 2\bmu_2^\top \bSigma^{-1} \bz \\
& \iff& \qquad - 2\bmu_1^\top \bSigma^{-1} \bz  + 2\bmu_2^\top \bSigma^{-1} \bz < \bmu_2^\top \bSigma^{-1} \bmu_2 - \bmu_1^\top \bSigma^{-1} \bmu_1 \\
& \iff & \qquad - 2 \lb \bmu_1^\top \bSigma^{-1} \bz - \bmu_2^\top \bSigma^{-1} \bz \rb < (\bmu_2-\bmu_1)^\top \bSigma^{-1} (\bmu_1+\bmu_2) \\
& \iff & \qquad (\bmu_1 - \bmu_2)^\top \bSigma^{-1} \bz > \frac{1}{2} \lsb (\bmu_1-\bmu_2)^\top \bSigma^{-1} (\bmu_1+\bmu_2) \rsb \\
& \iff & \qquad \ba^\top \bz > \frac{1}{2} \ba^\top (\bmu_1+\bmu_2) \\
& \iff & \qquad \ba^\top \lb \bz - \frac{1}{2}  (\bmu_1+\bmu_2) \rb > 0 \\
& \iff & \qquad \ba^\top (\bz-\bh) > 0,
\end{eqnarray*}
where $\ba$ and $\bh$ are defined above.
\EndKnitrBlock{proof}

Note that the discriminant rule is *linear* in $\bz$.

\BeginKnitrBlock{example}
<span class="example" id="exm:exnine2"><strong>(\#exm:exnine2) </strong></span>Consider the bivariate case ($p=2$) with $g=2$ groups, where $\Pi_1$ is the $N_2(\bmu_1,\bI_2)$ distribution and $\Pi_2$ is the $N_2(\bmu_2,\bI_2)$ distribution.  Suppose $\bmu_1 = \begin{pmatrix} c \\ 0 \end{pmatrix}$ and $\bmu_2 =  \begin{pmatrix} -c \\ 0 \end{pmatrix}$ for some constant $c>0$.  Here, $\ba = \bSigma^{-1} (\bmu_1 - \bmu_2) = \begin{pmatrix} 2c \\ 0 \end{pmatrix}$ and $\bh = \frac{1}{2}( \bmu_1 + \bmu_2 ) = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$.

The ML discriminant rule allocates $z$ to $\Pi_1$ if $\ba^\top (\bz - \bh) = \ba^\top \bz > 0$.  If we write $\bz = \begin{pmatrix} z_1 \\ z_2 \end{pmatrix}$ then $\ba^\top \bz = 2cz_1$, which is greater than zero if $z_1 > 0$.  Hence we allocate $\bz$ to $\Pi_1$ if $z_1 > 0$ and allocate $\bz$ to $\Pi_2$ if $z_1 \leq 0$.

FIX
<!--\begin{center}
\includegraphics[width=12cm,angle=0]{mldr_bivariate.pdf}
\end{center}
-->
\EndKnitrBlock{example}
  
\BeginKnitrBlock{example}
<span class="example" id="exm:exnine3"><strong>(\#exm:exnine3) </strong></span>A slightly more complicated version of the previous example: we still assume $\bmu_1=-\bmu_2$ but make no assumption about $\bSigma$.  Write $\ba = \begin{pmatrix} a_1 \\ a_2 \end{pmatrix}$ and $\bh = \frac{1}{2}( \bmu_1 + \bmu_2 ) = \bzero$. The ML discriminant rule allocates $z$ to $\Pi_1$ if $\ba^\top (\bz - \bh) = \ba^\top \bz > 0$.  If we write $\bz = \begin{pmatrix} z_1 \\ z_2 \end{pmatrix}$ then the boundary separating $\cR_1$ and $\cR_2$ is given by $\ba^\top \bz =  \begin{pmatrix} a_1 & a_2 \end{pmatrix} \begin{pmatrix} z_1 \\ z_2 \end{pmatrix} = a_1 z_1 + a_2 z_2 = 0$, i.e. $z_2 = -\frac{a_1}{a_2} z_1$.  This is a straight line through the origin with gradient $-a_1/a_2$.
\EndKnitrBlock{example}

Note that when $g>2$, the boundaries for the ML rule will be piece-wise linear rather than linear.

## The sample ML discriminant rule

To use the ML discriminant rule, above, we need to know the model parameters for each group.  In reality, we often do not know these parameters but we can estimate them from ``training'' data.  Training data typically consists of samples $\bx_{1,j}, \ldots, \bx_{n_j,j}$ known to be from population $\Pi_j$ ($j=1,\ldots ,g$).  Note that there are $n_j$ observations from population $\Pi_j$.

For simplicity, we shall assume that the populations have multivariate normal distributions with different means $\bmu_j$, $j=1,\ldots,g$, and the same covariance matrix, $\bSigma$.  Let $\bar{\bx}_j$ and $\bS_j$ be the sample mean and sample covariance matrix for the $j$th group.  Then $\bar{\bx}_j$ is an unbiased estimate of $\bmu_j$ and
$$\bS_u = \frac{1}{n-g} \sum_{k=1}^g n_k \bS_k$$
is an unbiased estimate of $\bSigma$ where $n = n_1 + n_2 + \ldots + n_g$.  The sample ML discriminant rule is then defined by substituting these estimates into \@ref(prp:nine1).

\BeginKnitrBlock{definition}
<span class="definition" id="def:sampleML"><strong>(\#def:sampleML) </strong></span>If $\Pi_k$ is the $N_p(\bmu_k,\bSigma)$ population, $k=1,\ldots,g$, then the sample ML discriminant rule allocates $\bz$ to $\Pi_j$ where $j $ is the value of $k$ which minimises
$$(\bz-\bar{\bx}_k)^\top \bS_u^{-1} (\bz-\bar{\bx}_k).$$
\EndKnitrBlock{definition}

In the case $g=2$, the rule allocates $\bz$ to $\Pi_1$ if and only if
$$\hat{\ba}^\top (\bz - \hat{\bh}) > 0$$
where $\hat{\ba} = \bS_u^{-1} (\bar{\bx}_1 - \bar{\bx}_2)$, $\hat{\bh} = \frac{1}{2} (\bar{\bx}_1 + \bar{\bx}_2)$ and $\bS_u$,  the pooled estimate of $\bSigma$, is given by
$$ \bS_u =  \frac{1}{n_1 + n_2 -2} (n_1 \bS_1 + n_2 \bS_2 ).$$

This is analogous to the Corollary following Proposition \@ref(prp:nine1).

\BeginKnitrBlock{example}
<span class="example" id="exm:exnine4"><strong>(\#exm:exnine4) </strong></span>Consider the G11PRB and G11STA module marks for $n_1 = 98$ students on G100 and $n_2 = 46$ students on G103.  The sample means and variances for each group are given by
\begin{eqnarray*}
\bar{\bx}_1 = \begin{pmatrix} 60.582 \\ 62.786 \end{pmatrix} &\qquad& \bar{\bx}_2 = \begin{pmatrix} 64.761 \\ 60.457 \end{pmatrix} \\
\bS_1 = \begin{pmatrix} 201.04 & 129.56 \\ 129.56 & 316.21 \end{pmatrix} &\qquad& \bS_2 = \begin{pmatrix} 229.88 & 177.02 \\ 177.02 & 354.16 \end{pmatrix}
\end{eqnarray*}
Hence,
\begin{eqnarray*}
\bS_u &=& \frac{1}{98+46-2} \lb 98 \bS_1 + 46 \bS_2 \rb = \begin{pmatrix} 213.21 & 146.76 \\ 146.76 & 332.96 \end{pmatrix}, \\
\bar{\bx}_1 - \bar{\bx}_2 &=& \begin{pmatrix} -4.179 \\ 2.329 \end{pmatrix}, \\
\hat{\bh} &=& \frac{1}{2} (\bar{\bx}_1 + \bar{\bx}_2) = \begin{pmatrix} 62.671 \\ 61.621 \end{pmatrix},
\end{eqnarray*}
and
$$\hat{\ba} = \bS_u^{-1} (\bar{\bx}_1 - \bar{\bx}_2) = \begin{pmatrix} 0.0067 & -0.0030 \\ -0.0030 & 0.0043 \end{pmatrix} \begin{pmatrix} -4.179 \\ 2.329 \end{pmatrix} = \begin{pmatrix} -0.035 \\ 0.022 \end{pmatrix}.$$

The sample ML discriminant rule allocates a new observation\\
 $\bz = (z_1, z_2)^\top$ to $\Pi_1$ if and only if
$$ \hat{\ba}^\top (\bz - \hat{\bh}) = \begin{pmatrix} -0.035 & 0.022 \end{pmatrix} \begin{pmatrix} z_1 - 62.671 \\ z_2 - 61.621 \end{pmatrix} > 0.$$

For example, if a student on this year's course scores 80 on G11PRB and 60 on G11STA then
$$ \hat{\ba}^\top (\bz - \hat{\bh}) = \begin{pmatrix} -0.035 & 0.022 \end{pmatrix} \begin{pmatrix} 80 - 62.671 \\ 60 - 61.621 \end{pmatrix} = -0.644 < 0,$$
and so we would allocate this student to G103.  The boundary, where $\hat{\ba}^\top (\bz - \hat{\bh}) = 0$, is shown below.

FIX
<!--
\begin{center}
\includegraphics[width=12cm,angle=0]{sample_mldr1.pdf}
\end{center}
-->

Note that the boundary line passes half-way between the two sample means.  In this example it is difficult to discriminate accurately between G100 and G103 because there is a large overlap between the two populations.

We could extend the example to include, say, students on GL11.  Here the boundary between the three populations is piece-wise linear and they meet at a common point.

FIX
<!--
\begin{center}
\includegraphics[width=12cm,angle=0]{sample_mldr2.pdf}
\end{center}
-->
\EndKnitrBlock{example}

## Fisher's linear discriminant rule

Recall that when the $\Pi_j$ are $N_p ( \bs{\mu}_j, {\mathbf \Sigma})$ populations, the ML
discriminant rule is linear if $g=2$ but not when $g>2$.  An alternative approach due to Fisher is to look for a linear discriminant function without assuming that the $\Pi_j$ are multivariate normal.

Suppose we have a training sample $\bx_{1,j}, \ldots, \bx_{n_j,j}$ from $\Pi_j$ ($j=1,\ldots,g$).  Calculate the `within' sum of squares matrix:
$$ \bW = \sum_{j=1}^g \sum_{i=1}^{n_j} (\bx_{ij} - \bar{\bx}_j) (\bx_{ij} - \bar{\bx}_j)^\top  = \sum_{j=1}^g n_j \bS_j $$
where $\bar{\bx}_j= \frac{1}{n_j} \sum_{i=1}^{n_j} \bx_{ij}$ is the sample mean of the $j$th group.  Also, calculate the 'between' sum of squares matrix
$$ \bB = \sum_{j=1}^g n_j (\bar{\bx}_j - \bar{\bx}) (\bar{\bx}_j - \bar{\bx})^\top ,$$
where $\bar{\bx} = \frac{1}{n} \sum_{j=1}^g \sum_{i=1}^{n_j}
\bx_{ij}$ is the overall mean, and $n=\sum_{j=1}^g n_j$.

Fisher's Criterion is to choose a unit vector, $\blambda$, to maximise
$$ \frac{\blambda^\top \bB \blambda}{\blambda^\top \bW \blambda}, $$
the ratio of the `between' sum of squares to the `within' sum of squares along $\blambda$.

The function $L(\bz)=\blambda^\top \bz$ is called Fisher's linear discriminant function.  Once $L(\bz)$ has been obtained, we allocate $\bz$ to the $\Pi_j$ whose discriminant score $L(\bar{\bx}_j)$ is closest to $L(\bz)$, that is, allocate $\bz$ to $\Pi_j$ iff
$$ | \blambda^\top \bz - \blambda^\top \bar{\bx}_j | = \min\limits_{1 \leq k \leq g} | \blambda^\top \bz - \blambda^\top \bar{\bx}_k |. $$

How do we find $\blambda$?

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:nine2"><strong>(\#prp:nine2) </strong></span>A vector $\blambda$ that maximises $$\frac{\blambda^\top \bB \blambda}{\blambda^\top \bW \blambda}$$ is an eigenvector of $\bW^{-1}\bB$ corresponding to the largest eigenvalue.
\EndKnitrBlock{proposition}

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}Assume $\bW$ is positive definite and note that $\bW$ is symmetric, so we can use the spectral decomposition theorem and write $\bW = \bQ \bLambda \bQ^\top$.

Define $\bgamma = \bW^{1/2} \blambda$.  Then $\blambda = \bW^{-1/2} \bgamma$ where
$\bW^{-1/2}=\bQ \bLambda^{-1/2} \bQ^\top$ and
\begin{eqnarray*}
\max_{\blambda \colon \blambda^\top \blambda=1} \left\{ \frac{\blambda^\top \bB \blambda}{\blambda^\top \bW \blambda} \right\}
&=& \max_{\bgamma \colon \bgamma \neq \bzero} \left\{ \frac{\bgamma^\top \bW^{-1/2} \bB \bW^{-1/2} \bgamma} {\bgamma^\top \bW^{-1/2} \bW \bW^{-1/2} \bgamma} \right\} \\
&=& \max_{\bgamma \colon \bgamma \neq \bzero} \left\{ \frac{ \bgamma^\top \bW^{-1/2} \bB \bW^{-1/2} \bgamma}{\bgamma^\top \bI_p \bgamma} \right\} \\
&=& \max_{\bgamma \colon \bgamma^\top \bgamma =1} \left\{ \bgamma^\top \bW^{-1/2} \bB \bW^{-1/2} \bgamma \right\}
\end{eqnarray*}

This is similar to the PCA situation in \S 3.2 where we chose $\bu$ to be the eigenvector corresponding to the largest eigenvalue of $\bS$ to maximise $\bu^\top \bS \bu$.  Hence, we choose $\bgamma$ to be the eigenvector corresponding to the largest eigenvalue of $\bW^{-1/2} \bB \bW^{-1/2}$.

If $\bgamma$ is an eigenvector of $\bW^{-1/2} \bB \bW^{-1/2}$ then, by definition,
$$\bW^{-1/2} \bB \bW^{-1/2} \bgamma = \rho \bgamma$$
 where $\rho$ is the corresponding eigenvalue.  Pre-multiplying both sides by $\bW^{-1/2}$ gives
\begin{eqnarray*}
\bW^{-1} \bB (\bW^{-1/2} \bgamma) &=& \rho \bW^{-1/2} \bgamma \\
\bW^{-1} \bB \blambda &=& \rho \blambda.
\end{eqnarray*}

So, the $\blambda$ we require is the unit eigenvector corresponding to the largest
eigenvalue of $\bW^{-1} \bB$. 
\EndKnitrBlock{proof}


When $g=2$, Fisher's rule and the sample ML rule with $\bSigma_1=\bSigma_2=\bSigma$ turn out to be the same.  Note that in
the sample ML rule we assumed that the two groups are from $N_p(\bmu_i, \bSigma)$
populations, but Fisher's rule makes no such assumption.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:nine3"><strong>(\#prp:nine3) </strong></span>If $g=2$ then Fisher's rule and the sample ML rule described in \S 9.3 are equivalent.
\EndKnitrBlock{proposition}

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}First, note that
\begin{eqnarray*}
\bar{\bx}_1 - \bar{\bx} &=& \bar{\bx}_1 - \lb \frac{n_1 \bar{\bx}_1 + n_2 \bar{\bx}_2}{n_1+n_2} \rb
= \frac{ (n_1+n_2) \bar{\bx}_1 - n_1 \bar{\bx}_1 - n_2 \bar{\bx}_2 }{n_1+n_2} \\
&=& \frac{ n_2 (\bar{\bx}_1 - \bar{\bx}_2) }{n_1 + n_2} = \frac{n_2 \bd}{n_1+n_2}
\end{eqnarray*}
where $\bd = \bar{\bx}_1 - \bar{\bx}_2$.  By analogy $\bar{\bx}_2 - \bar{\bx} = \frac{n_1 (-\bd)}{n_1+n_2}$.  Therefore,
\begin{eqnarray*}
\bB &=& n_1 (\bar{\bx}_1 - \bar{\bx})(\bar{\bx}_1 - \bar{\bx})^\top + n_2 (\bar{\bx}_2 - \bar{\bx})(\bar{\bx}_2 - \bar{\bx})^\top \\
&=& \frac{n_1 n_2^2}{(n_1+n_2)^2} \bd\bd^\top + \frac{n_2 n_1^2}{(n_1+n_2)^2} (-\bd)(-\bd)^\top \\
&=& \frac{n_1 n_2 (n_1 + n_2)}{(n_1+n_2)^2} \bd\bd^\top = \frac{n_1 n_2}{n_1+n_2} \bd\bd^\top.
\end{eqnarray*}
Let $c = \frac{n_1 n_2}{n_1+n_2}$.  Now $\blambda$ is an eigenvector of $\bW^{-1} \bB = c \bW^{-1} \bd \bd^\top$.  Also, the non-zero eigenvalues of $c \bW^{-1} \bd \bd^\top$ are the same as the non-zero eigenvalues of $c \bd^\top \bW^{-1} \bd$, which is scalar and so itself is the only non-zero eigenvalue.  The eigenvector, $\blambda$, must then satisfy
$$ c \bW^{-1} \bd \bd^\top \blambda = c \bd^\top \bW^{-1} \bd \blambda. $$
If we choose $\blambda = \bW^{-1} \bd$ then the equation is satisfied.  Hence $\blambda = \bW^{-1} (\bar{\bx}_1 - \bar{\bx}_2)$.

Let $r = \blambda^\top \bz$, $s = \blambda^\top \bar{\bx}_1$ and $t = \blambda^\top \bar{\bx}_2$, then Fisher's rule allocates $\bz$ to $\Pi_1$ if and only if
\begin{eqnarray*}
&&| r-s | < | r-t | \\
&\iff & (r-s)^2 < (r-t)^2 \\
&\iff & r^2 - 2rs + s^2 < r^2 - 2rt + t^2 \\
&\iff & 0 < 2r(s-t) + t^2 - s^2 \\
& \iff & 0 < 2r(s-t) + (t-s)(t+s) \\
& \iff & 0 < (s-t)(2r-t-s)
\end{eqnarray*}
Now $s-t = \blambda^\top (\bar{\bx}_1 - \bar{\bx}_2) = \bd^\top \bW^{-1} \bd$ which is a quadratic form and must therefore be positive, because $\bW$ is assumed to be positive definite.  Hence Fisher's rule allocates $\bz$ to $\Pi_1$ if
\begin{eqnarray*}
&& (2r-s-t) > 0\\
&\iff & r - \frac{1}{2}(s+t) > 0 \\
&\iff & \blambda^\top \lb \bz - \frac{1}{2}(\bar{\bx}_1 + \bar{\bx}_2) \rb > 0 \\
&\iff & (\bar{\bx}_1 - \bar{\bx}_2)^\top \bW^{-1} \lb \bz - \frac{1}{2}(\bar{\bx}_1 + \bar{\bx}_2) \rb > 0 \\
&\iff & (\bar{\bx}_1 - \bar{\bx}_2)^\top \bS_u^{-1} \lb \bz - \frac{1}{2}(\bar{\bx}_1 + \bar{\bx}_2) \rb > 0
\end{eqnarray*}
where the last line follows since $\bW = (n_1 + n_2 - 2)\bS_u$.  This is equivalent to the sample ML rule for $g=2$.
\EndKnitrBlock{proof}

For $g > 2$, the sample ML rule and Fisher's linear rule will not, in general, be the same.  Fisher's rule is linear when $g>2$ and is easier to implement than ML rules when there are several populations.  It is often reasonable to use Fisher's rule for non-normal populations.  In particular, Fisher's rule requires fewer assumptions than ML rules.  However, the ML rule is `optimal' in some sense when its assumptions are valid.

## Probability of misclassification

Let $p_{jk}$ denote the probability of allocating an observation to population $\Pi_j$, when in fact it comes from $\Pi_k$.  Therefore $p_{kk}$ is the probability of correctly classifying this observation and $1-p_{kk}$ is the probability of misclassification.

One way of estimating $p_{jk}$ is to consider the number of observations from the training data that are misclassified.  For example, if $n_k$ observations come from population $k$ and $n_{jk}$ is the number of observations from population $k$ classified as from population $j$, then
$$ \hat{p}_{jk} = \frac{n_{jk}}{n_k} $$
is an estimate of $p_{jk}$.

When $g=2$, $\Pi_j$ is $N_p(\bmu_j, \bSigma)$ and we use the ML rule, we obtain an explicit expression for $p_{12}$ and $p_{21}$ as follows.

Recall that we allocate $\bz$ to $\Pi_1$ if and only if  $U = \ba^\top (\bz-\bh)>0$ where $\ba = \bSigma^{-1} (\bmu_1 - \bmu_2)$ and $\bh = \frac{1}{2} (\bmu_1 + \bmu_2)$.

Suppose $\bz$ is from $\Pi_2$. Then $\bz \sim N_p(\bmu_2,\bSigma)$, so
\begin{eqnarray*}
E[U] &=& E[\ba^\top (\bz-\bh)] = \ba^\top (E[\bz] -\bh) = \ba^\top (\bmu_2-\bh); \\
\text{Var}(U) &=& \text{Var}(\ba^\top \bz - \ba^\top \bh) = \text{Var}(\ba^\top \bz) = \ba^\top \bSigma \ba.
\end{eqnarray*}
Hence, when $\bz$ is from $\Pi_2$, $U \sim N(\ba^\top (\bmu_2 - \bh), \ba^\top \bSigma \ba)$.

Define $\Delta^2 = (\bmu_1-\bmu_2)^\top \bSigma^{-1} (\bmu_1-\bmu_2)$, Then
\begin{eqnarray*}
\ba^\top (\bmu_2-\bh) &=& \ba^\top \lb \bmu_2 - \frac{1}{2}\bmu_1 - \frac{1}{2}\bmu_2 \rb = \frac{1}{2} \ba^\top (\bmu_2 - \bmu_1) \\
&=& \frac{1}{2} (\bmu_1 - \bmu_2)^\top \bSigma^{-1} (\bmu_2 - \bmu_1) \\
&=& -\frac{1}{2} (\bmu_1 - \bmu_2)^\top \bSigma^{-1} (\bmu_1 - \bmu_2) = -\frac{1}{2}\Delta^2.
\end{eqnarray*}
and
$$ \ba^\top \bSigma \ba = (\bmu_1 - \bmu_2)^\top \bSigma^{-1} \bSigma \bSigma^{-1} (\bmu_1 - \bmu_2) = \Delta ^2. $$
Hence  $U \sim N( -\frac{1}{2}\Delta^2, \Delta^2)$ when $\bz$ is from $\Pi_2$.


Now $\bz$ is allocated to $\Pi_1$ if $U>0$.  The probability of this event when $\bz$ is in fact from $\Pi_2$ is
\begin{eqnarray*}
p_{12} = P(U>0) &=& P \lb \frac{ U - (-\Delta^2/2) }{\Delta} > \frac{0 - (-\Delta^2/2) }{\Delta} \rb \\
&=& P \lb Z > \frac{\Delta^2}{2\Delta} = \frac{\Delta}{2}  \bigg| Z \sim N(0,1) \rb \\
&=& P \lb Z < -\frac{\Delta}{2} \rb
\end{eqnarray*}
which can be found from statistical tables.  A similar argument shows that $p_{21}=p_{12}$.

If we are using the sample ML rule then we can replace $\Delta^2$ with
$$D^2 = (\bar{\bx}_1 - \bar{\bx}_2)^\top \bS_u^{-1} (\bar{\bx}_1 - \bar{\bx}_2)$$
and $\hat{p}_{12} = P(Z<-D/2)$.
