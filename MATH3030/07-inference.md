# Inference in $1$ and $2$ samples based on MVN

In the univariate case the Student $t$ distribution plays a key role when we are dealing with normal random samples (i) in hypothesis testing for means (especially the 'paired' $t$-test and the $t$-test for comparing means in two independent samples), and (ii) the construction of confidence intervals for means.  In this chapter we develop the analogous results in the multivariate case, where observations are now assumed to be random vectors.  The role of the Student  $t$ distribution will be played by Hotelling's $T^2$, and the role of the $\chi^2$ is played by the Wishart distribution.

Despite there being important technical differences, the results in the multivariate case are seen to be  natural generalisations of the univariate results.

## Hypothesis testing: $\bSigma$ known

Let $\bx_1,\ldots,\bx_n$ be a random sample from $N_p(\bmu,\bSigma)$ where $\bSigma$ is assumed known and $\bmu = (\mu_1,\ldots,\mu_p)^\top$.  We wish to test the null hypothesis $\bmu = \ba$, where $\ba$ is fixed and pre-specified, against the alternative $\bmu \neq \ba$.  We could conduct $p$ separate univariate tests with null hypotheses $$H_0: \mu_i = a_i \qquad \text{vs.} \qquad H_1: \mu_i \neq a_i,$$ $i=1,\ldots,p.$  However, this ignores possible correlations between variables involved in different tests.  An alternative approach is to conduct a single hypothesis test using Proposition \@ref(prp:six5).

Recall that $\bar{\bx} \sim N_p(\bmu, \frac{1}{n} \bSigma)$.   Therefore $n^{1/2} \bar{\bx} \sim N_p(n^{1/2} \bmu, \bSigma)$.  Applying Proposition  we see that
$$
(n^{1/2}\bar{\bx}-n^{1/2}\bmu)^\top \bSigma^{-1} (n^{1/2}\bar{\bx}-n^{1/2}\bmu)
=n(\bar{\bx}-\bmu)^\top \bSigma^{-1} (\bar{\bx}-\bmu) \sim \chi_p^2
$$

We use this as our test statistic in a test at the $\alpha$\% significance level of the following hypotheses:
$$
H_0: \bmu = \ba \qquad \text{vs} \qquad H_1: \bmu \neq \ba,
$$
where $\ba$ is a fixed, pre-specified vector.  The relevant test statistic is
$$
\zeta^2 = n(\bar{\bx}-\ba)^\top \bSigma^{-1} (\bar{\bx}-\ba),
$$
and when $H_0$ is true, then $\zeta^2 \sim \chi_p^2$.

Critical value: We reject $H_0$ if $\zeta^2 > \chi^2_{p,\alpha}$.

Alternatively, we can state the result as a $p$-value where $p = P(\chi^2_p > \zeta_{\text{obs}}^2)$, and
$\zeta_{\text{obs}}^2$ is the observed value of the statistic $\zeta^2$.

The multivariate equivalent of a confidence interval is a confidence region and the $100(1-\alpha)$\% confidence region for $\bmu$ is $\{ \ba: \zeta^2 \leq \chi^2_{p,\alpha} \}$.

This confidence region will be the interior of an ellipse or ellipsoid.

\BeginKnitrBlock{example}
<span class="example" id="exm:unnamed-chunk-1"><strong>(\#exm:unnamed-chunk-1) </strong></span>The scatterplot shows the module marks for $n=209$ students on probability (PRB, $x_1$) and statistics (STA, $x_2$).


FIX FIGURE
<!--
\begin{center}
\includegraphics[width=12cm,angle=0]{prbsta1.pdf}
\end{center}
-->

 The observations $\bx_1,\ldots,\bx_{209}$ are assumed to be a random sample from $N_2(\bmu,\bSigma)$ where
$$
\bSigma = \begin{pmatrix} 200 & 150 \\ 150 & 300 \end{pmatrix},
 $$
and the sample mean vector is $\bar{\bx} = ( 61.957 , 62.632)^\top$.
 The target for the module mean for a large population of students should be exactly 60 for both modules.  We now conduct a hypothesis test of $H_0$ versus $H_1$ at the 5\% level to see if the lecturers have missed their target,
 where
$$
H_0: \bmu = \begin{pmatrix} 60 \\ 60 \end{pmatrix} \qquad \text{and} \qquad H_1: \bmu \neq \begin{pmatrix} 60 \\ 60 \end{pmatrix}.
$$

The test statistic is
$$
\zeta^2 = 209 \begin{pmatrix} 61.957 - 60 \\ 62.632 - 60 \end{pmatrix}^\top \begin{pmatrix} 200 & 150 \\ 150 & 300 \end{pmatrix}^{-1} \begin{pmatrix} 61.957 - 60 \\ 62.632 - 60 \end{pmatrix}.
$$
Now $|\bSigma| = 200 \times 300 - 150^2 = 37500$, so
$$\bSigma^{-1} = \frac{1}{37500} \begin{pmatrix} 300 & -150 \\ -150 & 200 \end{pmatrix} =  \begin{pmatrix} 0.008 & -0.004 \\ -0.004 & 0.016/3 \end{pmatrix}$$
and
$$\zeta^2 = 209 \begin{pmatrix} 1.957 \\ 2.632 \end{pmatrix}^\top \begin{pmatrix} 0.008 & -0.004 \\ -0.004 & 0.016/3 \end{pmatrix} \begin{pmatrix} 1.957 \\ 2.632 \end{pmatrix} = 5.512.
$$
The critical value is $\chi^2_{2,0.05} = 5.991$ so $\zeta^2 < \chi^2_{p,0.05}$ and we do not reject the null hypothesis at the 5\% level.

Note that if we had conducted separate univariate hypothesis tests of $H_0: \mu_1 = 60$ and $H_0: \mu_2 = 60$ then the test statistics would have been:
\begin{eqnarray*}
z_1 &=& \frac{\bar{x}_1 - \mu_1}{\sqrt{\sigma_1^2/n}} = \frac{61.957-60}{\sqrt{200/209}} = 2.000 \\
z_2 &=& \frac{\bar{x}_2 - \mu_2}{\sqrt{\sigma_2^2/n}} = \frac{62.632-60}{\sqrt{300/209}} = 2.196.
\end{eqnarray*}
The critical value would have been $Z_{0.025} = 1.960$ and both null hypotheses would have been rejected.  Therefore a multivariate hypothesis can be accepted when each of its univariate components is rejected and vice-versa.

Returning to the multivariate test, the $p$-value is $0.064$ and the $95$\% confidence region is the interior of an ellipse, centred on $\bar{\bx}$, with the angle of the major-axis governed by $\bSigma$.  We can see from the plot below that $(60,60)^\top$, marked with a cross, lies just inside the confidence region.


FIX FIGURE
<!--
\begin{center}
\includegraphics[width=12cm,angle=0]{prbsta2.pdf}
\end{center}
-->
\EndKnitrBlock{example}
  
  
  
## Hypothesis testing - 1 sample case

In \S 7.2 we considered a hypothesis test of $H_0: \bmu = \ba$ vs. $H_1: \bmu \neq \ba$ based on an IID sample from $N_p(\bmu,\bSigma)$ when $\bSigma$ was known.  In reality, we rarely know $\bSigma$, so we replace it with the sample covariance matrix, $\bS$, and the $\chi^2_p$ distribution is replaced with the $F_{p,n-p}$ distribution using Corollary 6.6.  The details are as follows.

Hypotheses: $H_0: \bmu = \ba$ vs $H_1: \bmu \neq \ba$.

Test statistic: $\gamma^2 = \frac{n-p}{p} (\bar{\bx}-\ba)^\top \bS^{-1} (\bar{\bx}-\ba)$, where $\gamma^2 \sim F_{p,n-p}$ under the the assumption that $H_0$ is true.

Critical value: We reject $H_0$ if $\gamma^2 > F_{p,n-p,\alpha}$, where $\alpha$ is the significance level.

Alternatively, we can state the result as a $p$-value where $p = P(F_{p,n-p} > \gamma^2)$.

The $100(1-\alpha)$\% confidence region for $\bmu$ is $\{ \ba: \gamma^2 \leq F_{p,n-p,\alpha} \}$, which will again be the interior of an ellipse or ellipsoid, but the confidence region is now determined by $\bS$ rather than $\bSigma$.

\BeginKnitrBlock{example}
<span class="example" id="exm:unnamed-chunk-2"><strong>(\#exm:unnamed-chunk-2) </strong></span>We return to the example of \S 7.2, with the module marks for $n=209$ students on probability (PRB, $x_1$) and statistics (STA, $x_2$).

The observations $\bx_1,\ldots,\bx_{209}$ are assumed to be a random sample from $N_2(\bmu,\bSigma)$, but now we assume that $\bSigma$ is unknown.  We calculate the sample mean and sample covariance matrix as,
$$\bar{\bx} = \begin{pmatrix} 61.957 \\ 62.632 \end{pmatrix} \qquad \bS = \begin{pmatrix} 215.29 & 157.19 \\ 157.19 & 333.56 \end{pmatrix}.$$
We conduct a hypothesis test at the 5\% level of:
$$H_0: \bmu = \begin{pmatrix} 60 \\ 60 \end{pmatrix} \qquad \text{vs.} \qquad H_1: \bmu \neq \begin{pmatrix} 60 \\ 60 \end{pmatrix}.$$

The test statistic is
\begin{align*}
\gamma^2 &= \frac{209-2}{2} \begin{pmatrix} 61.957 - 60 \\ 62.632 - 60 \end{pmatrix}^\top \begin{pmatrix} 215.29 & 157.19 \\ 157.19 & 333.56 \end{pmatrix}^{-1} \begin{pmatrix} 61.957 - 60 \\ 62.632 - 60 \end{pmatrix} \\
&= \frac{207}{2} \begin{pmatrix} 1.957 \\ 2.632 \end{pmatrix}^\top \begin{pmatrix} 0.0071 & -0.0033 \\ -0.0033 & 0.0046 \end{pmatrix} \begin{pmatrix} 1.957 \\ 2.632 \end{pmatrix} = 2.525.
\end{align*}
The critical value is $F_{2,207,0.05} = 3.040$ so $\gamma^2 < F_{p,n-p,0.05}$ and we do not reject the null hypothesis at the 5\% level.

The $p$-value is $0.082$ and the $95$\% confidence region is the interior of an ellipse, centred on $\bar{\bx}$, with the angle of the major-axis governed by $\bS$.  The confidence region is slightly larger than when $\bSigma$ was known.


FIX FIGURE
<!--\begin{center}
\includegraphics[width=12cm,angle=270]{prbsta3.pdf}
\end{center}
-->
\EndKnitrBlock{example}

## Hypothesis testing - 2 sample case

As with the univariate case, we may wish to test the difference between two population means.    As with univariate statistics, there are two cases to consider:

**Paired case** If $m=n$ and there exists some experimental link between $\bx_i$ and $\by_i$ then we can look at the differences $\bz_i = \by_i - \bx_i$ for $i=1,\ldots,n$.  For example, $\bx_i$ and $\by_i$ could be vectors of pre-treatment and post-treatment measurements, respectively, of the same variables. The crucial assumption is that the differences $\bz_i$ are IID $N_p(\bmu, \bSigma)$. To examine the null hypothesis of no difference between the means we would test $H_0: \bmu ={\mathbf 0}_p$ against $H_1: \bmu \neq {\mathbf 0}_p$.

We then base our inference on $\bar{\bz} = \frac{1}{n} \sum_{i=1}^n \bz_i = \bar{\by} - \bar{\bx}$, and proceed exactly as in the 1 sample case, using the test in \S 7.2 if $\bSigma$ is known, or the test in \S 7.3 if $\bSigma$ if unknown with $\bS = \frac{1}{n} \sum_{i=1}^n (\bz_i - \bar{\bz})(\bz_i - \bar{\bz})^\top$.

**Unpaired case** The unpaired case is where $\bx_i$ and $\by_i$ are independent and not connected to each other.  For example, in a clinical trial we may have two separate groups of patients, where one group receives a placebo and the other group receives an active treatment. Let $\bx_1,\ldots,\bx_n$ be an IID sample from $N_p(\bmu_1,\bSigma)$ and let $\by_1,\ldots,\by_m$ be an IID sample from $N_p(\bmu_2,\bSigma)$. In this case, we can base our inference on the following result.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:seven1"><strong>(\#prp:seven1) </strong></span>Let $\bx_1,\ldots,\bx_n$ be a random sample  from $N_p(\bmu_1,\bSigma_1)$ and let $\by_1,\ldots,\by_m$ be a random sample  from $N_p(\bmu_2,\bSigma_2)$. Then when $\bmu_1 = \bmu_2$ and $\bSigma_1 = \bSigma_2$,
$$\frac{nm}{n+m} (\bar{\by} - \bar{\bx})^\top \bS_u^{-1} (\bar{\by} - \bar{\bx}) \sim T^2(p,n+m-2),$$
where
$$\bS_u = \frac{n\bS_1 + m\bS_2}{n+m-2}$$
is the pooled unbiased variance matrix estimator
and $\bS_j$ is the sample covariance matrix for group $j$, $j=1,2$.
\EndKnitrBlock{proposition}

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}From Result 1.1 FIX and Proposition \@ref(prp:six1) we know that $\bar{\bx} \sim N_p \lb \bmu_1,n^{-1}\bSigma_1 \rb$ and $\bar{\by} \sim N_p \lb \bmu_2,m^{-1}\bSigma_2 \rb$, and $\bar{\bx}$ and $\bar{\by}$ are independent, so
$$\bar{\by} - \bar{\bx} \sim N_p \lb \bmu_2 - \bmu_1, \frac{1}{n}\bSigma_1 + \frac{1}{m} \bSigma_2 \rb.$$
 If $\bmu_1 = \bmu_2$ and $\bSigma_1 = \bSigma_2 = \bSigma$, then $\bar{\by} - \bar{\bx} \sim N_p \lb \bzero_p, \lb \frac{1}{n} + \frac{1}{m} \rb \bSigma \rb$ and
$$\bz = \lb \frac{1}{n} + \frac{1}{m} \rb^{-1/2} (\bar{\by} - \bar{\bx}) \sim N_p(\bzero_p,\bSigma).$$

From Proposition \@ref(prp:six12) we know that $n\bS_1 \sim W_p(\Sigma_1,n-1)$ and $m\bS_2 \sim W_p(\Sigma_2,m-1)$.  Therefore when $\bSigma_1 = \bSigma_2 = \bSigma$,
\begin{eqnarray*}
\bM = (n+m-2)\bS_u &=& (n+m-2)\frac{n\bS_1 + m\bS_2}{n+m-2} \\
&=& n\bS_1 + m\bS_2 \sim W_p(\bSigma,n+m-2)
\end{eqnarray*}
by Proposition \@ref(prp:six10), using the fact that $\bS_1$ and $\bS_2$ are independent.

Now $\bz$ is independent of $\bM$, since $\bar{\bx}$ and $\bar{\by}$ are independent of $\bS_1$ and $\bS_2$, respectively, by Proposition \@ref(prp:six6).  Therefore, applying Proposition \@ref(prp:six13) with $\bx = \bz$ and $\bM = (n+m-2)\bS_u$, we have
$$(n+m-2) \bz^\top ((n+m-2)\bS_u)^{-1} \bz = \bz^\top \bS_u^{-1} \bz \sim T^2(p,n+m-2)$$ and
\begin{eqnarray*}
\bz^\top \bS_u^{-1} \bz &=& \lb \frac{1}{n} + \frac{1}{m} \rb^{-1/2} (\bar{\by} - \bar{\bx})^\top \bS_u^{-1} \lb \frac{1}{n} + \frac{1}{m} \rb^{-1/2} (\bar{\by} - \bar{\bx}) \\
&=& \lb \frac{1}{n} + \frac{1}{m} \rb^{-1} (\bar{\by} - \bar{\bx})^\top \bS_u^{-1} (\bar{\by} - \bar{\bx}).
\end{eqnarray*}
Finally,
$$\lb \frac{1}{n} + \frac{1}{m} \rb^{-1} = \lb \frac{m}{nm} + \frac{n}{nm} \rb^{-1} = \lb \frac{n+m}{nm} \rb^{-1} = \frac{nm}{n+m},$$
so Proposition \@ref(prp:seven1) is proved.  
\EndKnitrBlock{proof}
  
  
  
As in the one sample case, we can convert Hotelling's two-sample $T^2$ statistic to the $F$ distribution using Proposition \@ref(prp:six14).

\BeginKnitrBlock{corollary}
<span class="corollary" id="cor:cseven1"><strong>(\#cor:cseven1) </strong></span>Using the notation of Proposition \@ref(prp:seven1), it follows that
$$\delta^2 = \frac{(n+m-p-1)}{(n+m-2)p} \frac{nm}{(n+m)} (\bar{\by} - \bar{\bx})^\top \bS_u^{-1} (\bar{\by} - \bar{\bx}) \sim F_{p,n+m-p-1}.$$
\EndKnitrBlock{corollary}

```Simply apply Proposition \@ref(prp:six14) to the statistic in Proposition \@ref(prp:seven1) (replace $n$ with $n+m-2$). 
```

\BeginKnitrBlock{example}
<span class="example" id="exm:unnamed-chunk-4"><strong>(\#exm:unnamed-chunk-4) </strong></span>For the probability and statistics marks in Example 3.5, is there a significant difference between students registered on G100 and G103?  The data is shown below, together with the sample means.


FIX FIGURE
<!--
\begin{center}
\includegraphics[width=12cm,angle=0]{prbsta4.pdf}
\end{center}
-->

Let $\bmu_1$ and $\bmu_2$ be the population means for G100 and G103 respectively.  Our hypotheses are
$$H_0: \bmu_1 = \bmu_2 \qquad \text{vs.} \qquad H_1: \bmu_1 \neq \bmu_2.$$

Let $\bx_1,\ldots,\bx_{n}$ be the marks for G100 students, which we assume are a random sample from $N_2(\bmu_1,\bSigma_1)$.  Similarly, let $\by_1,\ldots,\by_m$ be the marks for G103 students, which we assume are a random sample from $N_2(\bmu_2,\bSigma_2)$.  The sample summary statistics are:
\begin{eqnarray*}
n = 98 &\qquad& m = 46 \\
\bar{\bx} = \begin{pmatrix} 60.582 \\ 62.786 \end{pmatrix} &\qquad& \bar{\by} = \begin{pmatrix} 64.761 \\ 60.457 \end{pmatrix} \\
\bS_1 = \begin{pmatrix} 201.04 & 129.56 \\ 129.56 & 316.21 \end{pmatrix} &\qquad& \bS_2 = \begin{pmatrix} 229.88 & 177.02 \\ 177.02 & 354.16 \end{pmatrix}
\end{eqnarray*}

The assumption $\bSigma = \bSigma_1 = \bSigma_2$ does not look unreasonable given the sample covariance matrices, so we compute
\begin{eqnarray*}
\bS_u &=& \frac{1}{98+46-2} \lb 98 \begin{pmatrix} 201.04 & 129.56 \\ 129.56 & 316.21 \end{pmatrix} + 46 \begin{pmatrix} 229.88 & 177.02 \\ 177.02 & 354.16 \end{pmatrix} \rb \\
&=& \begin{pmatrix} 213.21 & 146.76 \\ 146.76 & 332.96 \end{pmatrix}
\end{eqnarray*}
and, therefore, ${\ds \bS_u^{-1} = \begin{pmatrix} 0.0067 & -0.0030 \\ -0.0030 & 0.0043 \end{pmatrix}}.$

The test statistic is
$$\delta^2 = \frac{141}{284} \times \frac{4508}{144} \begin{pmatrix} 4.179 \\ -2.329 \end{pmatrix}^\top \begin{pmatrix} 0.0067 & -0.0030 \\ -0.0030 & 0.0043 \end{pmatrix} \begin{pmatrix} 4.179 \\ -2.329 \end{pmatrix}
= 3.089$$

The critical value for $\alpha=0.05$ is
$$F_{2,98+46-2-1,\alpha} = F_{2,141,0.05} = 3.060.$$

Therefore $\delta^2 > F_{p,n+m-p-1}$, so we reject the null hypothesis at the 5\% level.  The $p$-value is 0.049.  So there is moderate evidence against $H_0$, the null hypothesis that $\bmu_1=\bmu_2$ and $\bSigma_1=\bSigma_2$.
\EndKnitrBlock{example}
