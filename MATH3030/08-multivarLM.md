
# The Multivariate Linear Model

In the standard linear model, the response variable, $y$, is univariate and the mean of $y$, $\mu=E[y]$, is modelled as a linear function of the elements of a covariate vector $\bx=(x_1, \ldots , x_q)^\top \in \mathbb{R}^q$, i.e. it is assumed that
$$
\mu = {\pmb \beta}^\top \bx,
$$
where ${\pmb \beta} \in \mathbb{R}^q$ is an unknown parameter vector to be estimated.  In this chapter we consider a generalisation of the standard linear model in which the response, $\by$, is now a $p \times 1$ vector.  In this setting, the linear model takes the form
$$
{\pmb \mu}=\bB^\top \bx,
$$
where the mean vector $\bmu$ is $p \times 1$, the covariate vector $\bx$ is $q \times 1$ and the parameter vector $\bB$ is $q \times p$.  The reason for having $\bB^\top$ above rather than simply $\bB$ will become clear later in the chapter.

## The standard univariate linear model 

In this section we give a brief review of the standard linear model in which the response is univariate.  Then, in subsequent sections, we discuss the multivariate linear model, i.e. multiple regression with a vector response which has a general covariance matrix.

Consider the univariate linear model in which
\begin{equation}
y_i = \bx_i^\top {\pmb \beta}+\epsilon_i, \qquad i=1, \ldots , n,
(\#eq:slm1)
\end{equation}
where $\stackrel{q \times 1}{\pmb \beta}$ is a parameter vector and $\bx_i$ is a $q \times 1$  covariate vector for experimental unit $i$.  We can also write \@ref(eq:slm1) in equivalent vector-matrix form
\begin{equation}
\by=\bX {\pmb \beta} +{\pmb \epsilon},
(\#eq:slm2)
\end{equation}
where $\stackrel{n \times q}{\bX}=[\bx_1 , \ldots , \bx_n]^\top$ is the matrix of covariates,
 $\stackrel{n \times 1}{\by}$ is the vector of univariate responses and $\stackrel{n \times 1}{\pmb \epsilon}$ is the vector of univariate `error' terms.

 It is assumed throughout this chapter that $\bX$ has full rank $p$, so that $(\bX^\top \bX)^{-1}$ exists.

 Most or all of the following assumptions are usually made in the standard linear model:


1. For each $i=1, \ldots , n$, $E[\epsilon_i]=0$.

2.  The $\epsilon_i$ are uncorrelated, i.e. for $i \neq j$, $\text{Cov}(\epsilon_i, \epsilon_j)=0$.

3. The $\epsilon_i$ have constant variance, i.e. $\text{Var}(\epsilon_i)=\sigma^2$, i.e. $\sigma^2$ does not depend on $i$.
4. The $\epsilon_i$ are IID $N(0, \sigma^2)$.
 
 It is clear that assumption 4. implies each of assumptions 1-3.  However, note that the least squares approach to be discussed below makes sense under assumptions 1-3 alone.  The attraction of assumption 4. is that it enables us to perform exact inference since the relevant distributions are known exactly, and the estimators of ${\pmb \beta}$ and $\sigma^2$ are maximum likelihood estimators (MLEs).  We shall assume 4. and its multivariate analogue, see \@ref(eq:MVNassumption) below, throughout this chapter.

The log-likelihood for models \@ref(eq:slm1) and \@ref(eq:slm2) under the Gaussian assumption 4.  is given by
\begin{align*}
\ell({\pmb \beta}, \sigma^2)&=-\frac{n}{2}\log (\sigma^2)--\frac{n}{2}\log(2\pi)-\frac{1}{2\sigma^2} \sum_{i=1}^n (y_i-\bx_i^\top {\pmb \beta})^2\\
& \qquad = -\frac{n}{2}\log (\sigma^2)-\frac{n}{2}\log(2\pi)-\frac{1}{2\sigma^2} (\by - \bX {\pmb \beta})^\top (\by - \bX {\pmb \beta}).
\end{align*}
Applying the results in \S 2.10 to the second expression above,
$$
\frac{\partial \ell}{\partial {\pmb \beta}}({\pmb \beta}, \sigma^2)=\frac{1}{\sigma^2}\bX^\top (\by - \bX {\pmb \beta})
$$
and
$$
\frac{\partial \ell}{\partial \sigma^2}({\pmb \beta}, \sigma^2)=-\frac{n}{2}\frac{1}{\sigma^2}+\frac{1}{2\sigma^4}
 (\by - \bX {\pmb \beta})^\top (\by - \bX {\pmb \beta}).
$$
Setting $\partial \ell(\hat{\pmb \beta}, \hat{\sigma}^2))/\partial {\pmb \beta}={\mathbf 0}_q$, the zero vector, and assuming that $\bX^\top \bX$ is invertible, implies that
\begin{equation}
\hat{\pmb \beta}=\left (\bX^\top \bX \right )^{-1}\bX^\top \by.
(\#eq:uni1)
\end{equation}
Also,  setting $\partial \ell(\hat{\pmb \beta}, \hat{\sigma}^2))/\partial \sigma^2=0$ gives
\begin{equation}
\hat{\sigma}^2 = \frac{1}{n}\by^\top \bP \by,
(\#eq:uni2)
\end{equation}
where
\begin{equation}
\bP=\bI_n - \bX \left ( \bX^\top \bX\right)^{-1}\bX^\top
(\#eq:defP)
\end{equation}
is a projection matrix.  The maximised log-likelihood is given by
\begin{align*}
\ell(\hat{\pmb \beta}, \hat{\sigma}^2)&= -\frac{n}{2}\log(\hat{\sigma}^2)-\frac{n}{2}\log(2\pi) -\frac{1}{2\hat{\sigma}^2}(\by - \bX \hat{\pmb \beta})^\top (\by - \bX \hat{\pmb \beta})\\
&= -\frac{n}{2}\log(\hat{\sigma}^2)-\frac{n}{2}\log(2\pi)-\frac{n}{2\by^\top \bP \by}\by ^\top \bP \by\\
&=-\frac{n}{2}\log(\hat{\sigma}^2)-\frac{n}{2}\log(2\pi)-\frac{n}{2}.
\end{align*}

Under the IID Gaussian assumption for the $\epsilon_i$,
$$
\hat{\pmb \beta} \sim N_q\left \{{\pmb \beta}, \sigma^2 (\bX^\top \bX)^{-1}\right\},
$$
and
$$
n\hat{\sigma}^2/\sigma^2 \sim \chi_{n-q}^2.
$$

Before moving on, we briefly discuss an important case of the standard linear model - the one-way analysis of variance.
Here, there are $t$ groups, say, with $n_j$ observations in group $j$, and the model takes the form
\begin{equation}
y_{ij}=\mu_j +\epsilon_{ij}, \qquad i=1, \ldots , n_j; \quad j=1, \ldots t,
(\#eq:ANOVA1)
\end{equation}
where the $\epsilon_{ij}$ are IID $N(0,\sigma^2)$ random variables, and $\mu_j$ is the mean of population $j$.  Note that we are assuming that $\text{Var}(y_{ij})=\sigma^2$
is constant across populations $j=1, \ldots , t$.

Define
$$
\bar{y}_{+j}=n_j^{-1}\sum_{i=1}^{n_j}y_{ij}, \qquad j=1, \ldots , t; \qquad \bar{y}_{++}=n^{-1}\sum_{j=1}^t \sum_{i=1}^{n_j}y_{ij},
$$
where $n=\sum_{j=1}^t n_j$.  It is easily checked that, under model \@ref(eq:ANOVA1), the MLE of $\mu_j$ is $\bar{y}_{+j}$.
Consider the null hypothesis
$$
H_0:\, \mu_1= \cdots = \mu_t.
$$
The total sum of squares, $T$, defined below,
has the following decomposition:
\begin{align}
T&\equiv \sum_{j=1}^t \sum_{i=1}^{n_j} (y_{ij}-\bar{y}_{++})^2 \nonumber\\
&=\sum_{j=1}^t \sum_{i=1}^{n_j} (y_{ij}-\bar{y}_{+j} + \bar{y}_{+j} -\bar{y}_{++})^2 \nonumber\\
&\sum_{j=1}^t \sum_{i=1}^{n_j} \big \{(y_{ij}-\bar{y}_{+j})^2 + (\bar{y}_{+j} -\bar{y}_{++})^2 + 2 (y_{ij}-\bar{y}_{+j}) (\bar{y}_{+j} -\bar{y}_{++})\big \} \nonumber\\
&=\sum_{j=1}^t \sum_{i=1}^{n_j} (y_{ij}-\bar{y}_{+j})^2 +\sum_{j=1}^t n_j (\bar{y}_{+j}-\bar{y}_{++})^2 \nonumber\\
&=W + B,
(\#eq:T=B+W)
\end{align}
using the fact that sum over $i$ of the product term above is $0$.  In the above,
$W$ stands for `within' sum of squares, also known as the residual sum of squares, and $B$ stands for `between' sum of squares.  The MLE of $\sigma^2$ is given by $W/n$.   Standard theory tells us that under model \@ref(eq:ANOVA1), $W \sim \sigma^2 \chi_{n-t}^2$.  Moreover, under $H_0$, $B \sim \sigma^2 \chi_{t-1}^2$, $W$ and $B$ are independent, and therefore $T \sim \sigma^2 \chi_{n-1}^2$.
To test $H_0$ we use the fact that, under $H_0$,
\begin{equation}
f=\frac{B/(t-1)}{W/(n-t)}
(\#eq:Fstat)
\end{equation}
has an $F_{t-1, n-t}$ distribution. We reject $H_0$ when $f$ is `large' compared with an $F_{t-1,n-t}$ random variable.

In Example Sheet 3 you are asked to show that the `twice the log-likelihood ratio statistic' for testing $H_0$ against the general alternative with $\mu_1, \ldots , \mu_t$ is given by
\begin{equation}
n_+ \log \left ( 1+ \frac{B}{W}  \right )=n_+ \log (1+W^{-1}B),
(\#eq:ANOVA2)
\end{equation}
which is an increasing function of the statistic \@ref(eq:Fstat).  Consequently, we can think of the classical $F$ test as being equivalent to a likelihood ratio test.



## Multivariate Linear Model

In the standard linear model the responses $y_i$ are univariate.  In the multivariate linear model, the responses are $p \times 1$ vectors $\by_i$.  Here, the linear model takes the form
\begin{equation}
\by_i=\bB^\top \bx_i +{\pmb \epsilon}_i, \qquad i=1, \ldots , n,
(\#eq:mlm1)
\end{equation}
where $\stackrel{q \times p}{\bB}$ is a parameter matrix, the $\bx_i$ are $q \times 1$ covariate vectors as in \S 8.2, and the ${\pmb \epsilon}_i$ are $p \times 1$ error vectors.  The model \@ref(eq:mlm1) may be written in matrix form as
\begin{equation}
\bY = \bX \bB +\bE,
(\#eq:mlm2)
\end{equation}
where $\stackrel{n \times p}{\bY}=[\by_1, \ldots , \by_n]^\top$ is the data matrix for the $y$-variables, $\bX=[\bx_1, \ldots , \bx_n]^\top$ is the $n \times q$ data matrix for the $x$-variables, defined as in \S 8.2, and $\stackrel{n \times p}{\bE}=[{\pmb \epsilon}_1, \ldots , {\pmb \epsilon}_n]^\top$.

By analogy with assumption 4. in \S 8.2, it is assumed that
\begin{equation}
{\pmb \epsilon_1}, \ldots , {\pmb \epsilon}_n \quad \text{are
IID}\quad N_p({\mathbf 0}_p, \bSigma).
(\#eq:MVNassumption)
\end{equation}

Under this MVN assumption, the  log-likelihood function $\ell(\bB, \bSigma)$  for the parameter matrices $\bB$ and $\bSigma$ is given by
\begin{align}
\ell(\bB, \bSigma)&=-\frac{np}{2}\log(2\pi) -\frac{n}{2}\log(\vert \bSigma \vert) \nonumber \\
& \qquad \qquad -\frac{1}{2}\text{tr}\left \{
(\bY-\bX\bB) \bSigma^{-1} (\bY - \bX \bB)^\top\right \}.
(\#eq:MVNlik)
\end{align}


\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:eight1"><strong>(\#prp:eight1) </strong></span>The maximum likelihood estimators of $\bB$ and $\bSigma$ in \@ref(eq:MVNlik)
are given by
\begin{equation}
\hat{\bB}= (\bX^\top \bX)^{-1}\bX^\top \bY
(\#eq:MVbeta)
\end{equation}
and
\begin{equation}
\hat{\bSigma}=\frac{1}{n}\bY^\top \bP \bY,
(\#eq:RSS)
\end{equation}
where $\bP$ is the matrix defined in \@ref(eq:defP).  The maximised log-likelihood is given by
\begin{equation}
\ell(\hat{\bB}, \hat{\bSigma})=- \frac{n}{2} \log(\vert \hat{\bSigma} \vert )-\frac{np}{2}\log(2\pi) -\frac{np}{2}.
(\#eq:maxlik)
\end{equation}
\EndKnitrBlock{proposition}

**Remark:** note how similar \@ref(eq:MVbeta) and \@ref(eq:RSS) are to their univariate counterparts \@ref(eq:uni1) and \@ref(eq:uni2), respectively; the only thing that is different is that $\bY$ is now an $n \times p$ matrix rather than an $n \times 1$ vector.

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}**of \@ref(eq:MVbeta)}.** Recall that $\bP$ defined in \@ref(eq:defP) is a projection matrix, i.e. $\bP^\top =\bP$ and $\bP^2 = \bP$.
Moreover,
\begin{equation}
\bP \bX=\bX-\bX (\bX^\top \bX)^{-1}\bX^\top \bX=\bX-\bX={\mathbf 0}_{n,q}
(\#eq:orthog27a)
\end{equation}
and
\begin{equation}
\bX^\top \bP=\bX^\top - \bX^\top \bX (\bX^\top \bX)^{-1}\bX^\top =\bX^\top -\bX^\top ={\mathbf 0}_{q,n}.
(\#eq:orthog27b)
\end{equation}
Now write
\begin{align*}
\bY-\bX \bB&=\bY-\bX \hat{\bB} +\bX \hat{\bB}-\bX \bB\\
&=\bY-\bX (\bX^\top \bX)^{-1}\bX^\top \bY+\bX(\hat{\bB}-\bB)\\
&=\bP \bY +\bX (\hat{\bB}-\bB).
\end{align*}
Then, using \@ref(eq:orthog27a) and \@ref(eq:orthog27b),
\begin{align}
(\bY &- \bX \bB)^\top (\bY - \bX \bB) \nonumber \\
&=\{\bP \bY + \bX (\hat{\bB}-\bB)\}^\top \{\bP \bY +\bX(\hat{\bB}-\bB)\}\nonumber \\
&=\bY^\top \bP \bY + \bY^\top \bP\bX (\hat{\bB}-\bB)\nonumber \\
& \qquad \qquad + (\hat{\bB}-\bB)^\top \bX^\top \bP \bY +(\hat{\bB}-\bB)^\top \bX^\top \bX (\hat{\bB}-\bB)
\nonumber\\
&=\bY^\top \bP \bY + (\hat{\bB}-\bB)^\top \bX^\top \bX (\hat{\bB}-\bB).
(\#eq:decomp87)
\end{align}
 As noted in Chapter 2, for any compatible matrices $\bW$ and $\bZ$, we have $\text{tr}(\bW \bZ)=\text{tr}(\bZ \bW)$.   So,  using \@ref(eq:decomp87), it follows that the trace term in \@ref(eq:MVNlik) may be written
\begin{align}
&\text{tr}\left \{ (\bY-\bX \bB)\bSigma^{-1} (\bY - \bX \bB)^\top \right \}\nonumber \\
&=\text{tr}\left \{\bSigma^{-1} (\bY-\bX \bB)^\top (\bY - \bX \bB) \right \}\nonumber \\
&=\text{tr}\left \{  \bSigma^{-1} ( \bY -\bX \hat{\bB}+\bX\hat{\bB} -\bX \bB  )^\top
( \bY -\bX \hat{\bB}+\bX\hat{\bB} -\bX \bB  )       \right \}\nonumber \\
&= \text{tr}\left [\bSigma^{-1} \{\bY^\top \bP \bY +(\hat{\bB}-\bB)^\top \bX^\top \bX (\hat{\bB}-\bB)  \}\right ]
\nonumber\\
&=\text{tr} \left (\bSigma^{-1}\bY^\top \bP \bY\right ) + \text{tr}\left \{
\bX (\hat{\bB}-\bB) \bSigma^{-1} (\hat{\bB}-\bB)^\top  \bX^\top  \right \}
(\#eq:twobits)
\end{align}
The first term on the RHS of \@ref(eq:twobits) does not depend on $\bB$.  In the second term, the matrix inside the trace is non-negative definite for all $\bB$, and therefore the second term in \@ref(eq:twobits) is non-negative, and has a minimum value of $0$, achieved uniquely when $\bB=\hat{\bB}$.  Therefore, for fixed $\bSigma$, \@ref(eq:MVNlik) is maximised when $\bB=\hat{\bB}$, as defined in \@ref(eq:MVbeta).  \hfill $\square$

To prove \@ref(eq:RSS), we require the following result, whose proof is omitted because it is a bit too technical.

\EndKnitrBlock{proof}
\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}Suppose that $\bA$ is a fixed symmetric positive definite  $p \times p$ matrix and define the function
$$
G(\bSigma) \equiv -\frac{n}{2} \log \left (\vert \bSigma \vert\right)  -\frac{1}{2}\text{tr} \left (\bSigma^{-1} \bA \right).
$$
Then $G(\bSigma)$ is maximised over symmetric positive definite $p \times p$ matrices $\bSigma$ at
$\hat{\bSigma}=n^{-1}\bA$, and the maximum value of $G$ is given by
$$
G(\hat{\bSigma})= -\frac{n}{2} \log \left (\vert n^{-1} \bA \vert \right ) -np/2.
$$
\EndKnitrBlock{proof}
  
\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}**of \@ref(eq:RSS) and \@ref(eq:maxlik)}.**  Again using the result $\text{tr}(\bW \bZ)=\text{tr}(\bZ \bW)$ for compatible matrices $\bW$ and $\bZ$, we have
$$
\text{tr}\{(\bY -\bX \hat{\bB}) \bSigma^{-1}(\bY -\bX \hat{\bB})^\top\}=\text{tr}(\bSigma^{-1}(\bP \bY)^\top \bP \bY)=\text{tr}(\bSigma^{-1} \bY^\top \bP \bY),
$$
because $\bP$ is a projection matrix.    Consequently,
$$
\ell(\hat{\bB},\bSigma)=-\frac{np}{2}\log(2\pi)-\frac{n}{2}\log \vert \bSigma \vert -\frac{1}{2} \text{tr}
\left \{ \bSigma^{-1}\bY^\top \bP \bY \right \},
$$
and we want to maximise $\ell(\hat{\bB}, \bSigma)$ over symmetric positive definite matrices $\bSigma$.  So take $\bA=\bY^\top \bP \bY$.  Then, using Proposition \@ref(prp:eight2), $\hat{\bSigma} = n^{-1}\bY^\top \bP \bY$, which agrees with \@ref(eq:RSS).

To prove \@ref(eq:maxlik), note that
\begin{align*}
\text{tr}\{(\bY-\bX \hat{\bB}) \hat{\bSigma}^{-1}(\bY -\bX \hat{\bB})^\top \}
&=\text{tr}\{\hat{\bSigma}^{-1} (\bP \bY)^\top \bP \bY\}\\
&= n \text{tr}\{(\bY^\top \bP \bY)^{-1} \bY^\top \bP \bY\}\\
&=n\text{tr}(\bI_p)=np,
\end{align*}
so \@ref(eq:maxlik) follows after substitution of $\bB=\hat{\bB}$ and $\bSigma=\hat{\bSigma}$ into \@ref(eq:MVNlik).
\EndKnitrBlock{proof}

We now determine the joint distribution of $\hat{\bB}$ and $\hat{\bSigma}$.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:eight3"><strong>(\#prp:eight3) </strong></span>Assume that \@ref(eq:mlm2) and \@ref(eq:MVNassumption) both hold and assume that $\bX$ has full rank $q$, where $q < n$.  Then the following results hold.


1. The MLE of $\bB$, $\hat{\bB}$ defined in \@ref(eq:MVbeta), and the MLE of $\bSigma$, $\hat{\bSigma}$ defined in \@ref(eq:RSS),
are independent.  Moreover, the elements of $\hat{\bB}$ are jointly multivariate normal, while $n\hat{\bSigma}
\sim W_p(\bSigma, n-q)$.

2. The MLE $\hat{\bB}$ satisfies $E[\hat{\bB}]=\bB$, i.e. $\hat{\bB}$ is unbiased for $\bB$.

3.Write $\hat{\bb}_{[j]}$ for column $j$ of $\hat{\bB}$, $j=1, \ldots , p$, so that $\hat{\bB}=[\hat{\bb}_{[1]}, \ldots , \hat{\bb}_{[p]}]$.
Then for $j,k=1, \ldots , p$,
$$
\text{Cov}(\hat{\bb}_{[j]}, \hat{\bb}_{[k]})=\sigma_{jk} (\bX^\top \bX)^{-1}.
$$
\EndKnitrBlock{proposition}
  
\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}*of Proposition \@ref(prp:eight3)*
First, observe that $\hat{\bSigma}$ is a function of the elements of $\bP \bY$, namely
$$
\hat{\bSigma}=n^{-1} \bY^\top \bP^\top \bP \bY=n^{-1}\bY^\top \bP \bY,
$$
because $\bP$ is a projection matrix.  Moreover, under the MVN assumption \@ref(eq:MVNassumption), $\bP \bY$ and
$\hat{\bB}$ are jointly MVN, so to prove part 1. it will be sufficient to prove that each column of $\hat{\bB}$ is uncorrelated with
each column of $\bP \bY$.  Write $\by_{[j]}$ for column $j$ of $\bY$.  Then
$$
\bY=[\by_{[1]}, \ldots , \by_{[p]}], \qquad \bP \bY=[\bP\by_{[1]}, \ldots , \bP \by_{[p]}],
$$
and
$$
\hat{\bB}\equiv [\hat{\bb}_{[1]}, \ldots , \hat{\bb}_{[p]}] =[(\bX^\top \bX)^{-1}\bX^\top \by_{[1]}, \ldots , (\bX^\top \bX)^{-1}\bX^\top \by_{[p]}],
$$
where, as before,  $\hat{\bb}_{[j]}$ is column $j$ of $\hat{\bB}$.  The covariance between column $j$ of $\hat{\bB}$ and column $k$ of $\bP \bY$ is given by
\begin{align*}
\text{Cov}(\hat{\bb}_{[j]}, \bP \by_{[k]})&= \text{Cov}\{(\bX^\top \bX)^{-1}\bX^\top \by_{[j]}, \bP \by_{[k]}\}\\
&=(\bX ^\top \bX)^{-1} \bX^\top \text{Cov}(\by_{[j]}, \by_{[k]})\bP\\
&=(\bX ^\top \bX)^{-1} \bX^\top \sigma_{jk}\bI_n\bP\\
&=\sigma_{jk}(\bX^\top \bX)^{-1} \bX^\top (\bI_n - \bX (\bX^\top \bX)^{-1}\bX^\top)={\mathbf 0}_{q,n},
\end{align*}
where $\bSigma=(\sigma_{jk})_{j,k=1}^p$. In the above we have used the definition of $\bP$ in \@ref(eq:defP) and the fact that
\begin{equation}
\text{Cov}(\by_{[j]}, \by_{[k]})=\sigma_{jk}\bI_n
(\#eq:energise)
\end{equation}
due to the independence of the rows of $\bY$.  Therefore $\hat{\bB}$ and $\bP\bY$ are independent by Proposition \@ref(prp:six4), and
 by Proposition \@ref(thm:six11) (Cochran's theorem), $n\hat{\bSigma} \sim W_p(\bSigma, n-q)$ which concludes the proof of part 1.

Part 2. follows because, taking the expectation of \@ref(eq:MVbeta) we find that
\begin{align*}
E[\hat{\bB}]&=E[(\bX^\top \bX)^{-1}\bX^\top \bY]\\
&=(\bX^\top \bX)^{-1}\bX^\top E[\bY]\\
&=(\bX^\top \bX)^{-1}\bX^\top \bX \bB =\bB.
\end{align*}

For part 3., using \@ref(eq:energise) again,
\begin{align*}
\text{Cov}(\hat{\bb}_{[j]}, \hat{\bb}_{[k]})&=\text{Cov}\{(\bX^\top \bX)^{-1}\bX^\top \by_{[j]}, (\bX^\top \bX)^{-1} \bX^\top \by_{[k]}\}\\
&=(\bX^\top \bX)^{-1}\bX^\top \text{Cov}(\by_{[j]}, \by_{[k]}) \bX (\bX^\top \bX)^{-1}\\
&=(\bX^\top \bX)^{-1} \bX^\top \sigma_{jk}\bI_n\bX (\bX^\top \bX)^{-1}=\sigma_{jk}(\bX^\top \bX)^{-1},
\end{align*}
as required.  
\EndKnitrBlock{proof}

We now investigate the situation where we wish to test \@ref(eq:mlm2) against a sub-model.  Specifically, we consider the null hypothesis
$$
H_0: \, \bB=\begin{pmatrix}\bB^\ast \\{\mathbf 0}_{q-r,p} \end{pmatrix}, \qquad \bSigma \quad \text{unrestricted},
$$
where $\bB^\ast$ is $r \times p$ with $1 \leq r <q$. The alternative hypothesis is
$$
H_1: \, \stackrel{q \times p}{\bB} \quad \text{unrestricted}, \qquad \bSigma \quad
\text{unrestricted},
$$
i.e. model \@ref(eq:mlm2) under the MVN assumption \@ref(eq:MVNassumption).  Note that $H_0$ is nested within $H_1$.
Let us write $\hat{\bB}_j$ and $\hat{\bSigma}_j$ for the MLEs of $\bB$ and $\bSigma$ under hypothesis $H_j$, $j=0,1$.   Then, using \@ref(eq:maxlik) in Proposition \@ref(prp:eight1) to calculate the maximised likelihood under $H_0$ and $H_1$, we obtain the Wilks statistic, $\omega_{0,1}$ (equals ``twice the different in maximised log likelihoods''):
\begin{align*}
\omega_{0,1}&=2\{\ell(\hat{\bB}_1, \hat{\bSigma}_1)-\ell(\hat{\bB}_0, \hat{\Sigma}_0)\}\\
&= - n\log(\vert \hat{\bSigma}_1\vert) + n\log (\vert \hat{\bSigma}_0\vert)
=n\log \left ( \frac{\vert \hat{\bSigma}_0 \vert}{\vert \hat{\bSigma}_1\vert}  \right).
\end{align*}

In the univariate case there is a simple, explicit transformation from the Wilks statistic $\omega_{0,1}$ to the classical $F_{q-r,n-q}$ distribution for testing $H_0$ against $H_1$.  In the multivariate case, the exact distribution of $\omega_{0,1}$ under $H_0$ does not have a simple relationship with an $F$-distribution, and the situation is a lot more complicated.

When $n$ is large, however, we can use the large sample log-likelihood ratio test, which implies that, under $H_0$,
$\omega_{0,1}$ is approximately $\chi^2$, and we should reject $H_0$ when $\omega_{0,1}$ is sufficiently large.

The relevant degrees of freedom of the $\chi^2$ are now calculated.

Under $H_0$, the number of free parameters is
$$
rp + \frac{1}{2}p(p+1),
$$
while under $H_1$ there are
$$
qp +\frac{1}{2}p(p+1)
$$
free parameters.  So the difference is $(q-r)p$, and so we should refer $\omega_{0,1}$ to $\chi_{(q-r)p}^2$.

## One-way MANOVA

 We now consider the multivariate version of the one-way ANOVA considered at the end of \S 8.2, known as one-way MANOVA.  The  model is defined by
 \begin{equation}
 \by_{ij} = \bmu_j + {\pmb \epsilon}_{ij}, \qquad  i=1, \ldots , n_j; \quad j=1, \ldots, t,
 (\#eq:MANOVA1)
 \end{equation}
 where the ${\pmb \epsilon}_{ij}$ are IID $N_p({\mathbf 0}_p, \bSigma)$, and the $\by_{ij}$ and $\bmu_j$ are also
 $p \times 1$ vectors.  It is assumed that $\text{Var}({\pmb \epsilon}_i)=\bSigma$ is constant across the $t$ populations.

 One possibility is to use the linear model framework developed in \S 8.3.    In this case, $\bX$ is a matrix with each element equal to zero or one.  However, it is also feasible to do the calculations directly.  This is what we shall do.

 The log-likelihood of model \@ref(eq:MANOVA1) is given by
 \begin{align}
 \ell(\bmu_1, \ldots , \bmu_t)&=-\frac{np}{2}\log(2\pi)-\frac{n}{2}\log(\vert \bSigma\vert)\nonumber\\
 & \qquad \qquad -\frac{1}{2}\sum_{j=1}^t \sum_{i=1}^{n_j}
 (\by_{ij}-\bmu_j)^\top \bSigma^{-1} (\by_{ij}-\bmu_j),
 (\#eq:MANOVAlik)
 \end{align}
 where $n=n_1 + n_2 +\cdots + n_t$.

 Using results in \S 2.10, the partial derivative, or gradient,  of \@ref(eq:MANOVAlik) with respect to the vector $\bmu_k$ is given by
 \begin{align*}
 \frac{\partial \ell}{\partial \bmu_k}(\bmu_1, \ldots , \bmu_t, \bSigma)&=\sum_{i=1}^{n_j} \bSigma^{-1}(\by_{ik}-\bmu_k)\\
 &=n_k \bSigma^{-1} (\bar{\by}_{+k}-\bmu_k).
 \end{align*}
where $\bar{\by}_{+k}$ is the sample mean of group $k$, i.e.
$$
\bar{\by}_{+k}=n_k^{-1}\sum_{i=1}^{n_k} \by_{ik}.
$$
So, setting $\partial \ell/\partial \bmu_k={\mathbf 0}_p$ implies $\hat{\bmu}_k=\bar{\by}_{+k}$.
Therefore
\begin{align*}
\ell(\hat{\bmu}_1, \ldots , \hat{\bmu}_t, \bSigma)&= -\frac{np}{2}\log(2\pi)-\frac{n}{2}\log(\vert \bSigma \vert)\\
&\qquad \qquad - \frac{1}{2} \sum_{j=1}^t \sum_{i=1}^{n_j} (\by_{ij}-\bar{\by}_{+j})^\top \bSigma^{-1}
(\by_{ij}-\bar{\by}_{+j})\\
&=-\frac{n p}{2}\log(2\pi) -\frac{n}{2}\log(\vert \bSigma\vert)-\frac{1}{2}\text{tr}(\bSigma^{-1}\bW),
\end{align*}
where $\bW$, defined by
$$
\bW=\sum_{j=1}^t \sum_{i=1}^{n_j} (\by_{ij}-\bar{\by}_{+j})(\by_{ij}-\bar{\by}_{+j})^\top,
$$
is the matrix version of the `within' sum of squares considered in \S 8.2.  Using Proposition \@ref(prp:eight2), we deduce that
$$
\hat{\bSigma}=n^{-1}\bW.
$$
Consequently,
$$
\ell(\hat{\bmu}_1, \ldots , \hat{\bmu}_t,\hat{\bSigma})=-\frac{n}{2}\log(\vert \hat{\bSigma}\vert)-\frac{np}{2}- \frac{np}{2}\log(2\pi).
$$

Now consider the key null hypothesis for a one-way MANOVA:
$$
H_0: \, \bmu_1=\cdots =\bmu_t, \qquad \bSigma \quad \text{unrestricted}.
$$
Under $H_0$, the $\by_{ij}$ are IID $N_p(\bmu, \bSigma)$, so the log-likelihood under $H_0$ is
\begin{align*}
\ell_0(\bmu, \bSigma) & \equiv \ell (\bmu, \ldots ,\bmu, \bSigma)\\
&=-\frac{np}{2}\log(2\pi) - \frac{n}{2}\log(\vert \bSigma\vert)-\frac{1}{2}\sum_{j=1}^t \sum_{i=1}^{n_j} (\by_{ij}-\bmu)^\top \bSigma^{-1}(\by_{ij}-\bmu),
\end{align*}
and so the MLE of $\bmu$ under $H_0$ is given by
$$
\hat{\bmu}_0=\frac{1}{n} \sum_{j=1}^t \sum_{i=1}^{n_j} \by_{ij} =\bar{\by}_{++},
$$
and, using Proposition \@ref(prp:eight2) again, it is seen that the MLE of $\bSigma$ under $H_0$ is
$$
\hat{\bSigma}_0=n^{-1}\bT,
$$
where $\bT$ is the matrix analogue of the total sum of squares, i.e.
$$
\bT=\sum_{j=1}^t \sum_{i=1}^{n_j} (\by_{ij}-\bar{\by}_{++})(\by_{ij}-\bar{\by}_{++})^\top.
$$

The Wilks statistic, $\omega_0$, for testing $H_0$ against the general alternative \@ref(eq:MANOVA1) is then
\begin{align}
\omega_0=2\{\ell(\hat{\bmu}_1, \ldots , \hat{\bmu}_t, \hat{\bSigma})-\ell_0(\hat{\bmu}_0, \hat{\bSigma}_0)\}&=n\log \left (\vert \hat{\bSigma}_0 \vert / \vert \hat{\bSigma}\vert  \right ) \nonumber \\
&=n \log\left ( \vert \bT\vert/\vert \bW \vert \right ).
(\#eq:MANOVA2)
\end{align}
The degrees of freedom under $H_0$ are $p+p(p+1)/2$ and the degrees of freedom under \@ref(eq:MANOVA1) are $pt+p(p+1)/2$, so the difference is $p(t-1)$. Consequently, when the $n_j$ are all large, we should refer $\omega_0$ to
$\chi_{p(t-1)}^2$ and reject $H_0$ when $\omega_0$ is sufficiently large.

It is not immediately obvious that \@ref(eq:MANOVA2) is a natural generalisation of \@ref(eq:ANOVA2).  However, in Example Sheet 4 you are asked to prove that
$$
\bT=\bW+\bB,
$$
where $\bB$ is the matrix analogue of the `between' sum of squares $B$ in \@ref(eq:T=B+W), i.e.
$$
\bB = \sum_{j=1}^t  n_j (\bar{\by}_{+j}-\bar{\by}_{++})(\bar{\by}_{+j}-\bar{\by}_{++})^\top.
$$
Consequently,
\begin{align*}
\omega_0&=n \log (\vert \bT \vert /\vert \bW \vert)
=n \log (\vert \bW^{-1} \vert \vert \bW + \bB\vert)\\
&=n \log (\vert \bW^{-1}(\bW+\bB)\vert)
=n \log (\vert \bI_p +\bW^{-1} \bB \vert),
\end{align*}
which \textit{is} a natural generalisation of \@ref(eq:ANOVA2).

