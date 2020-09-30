# Canonical Correlation Analysis {#cca}


Suppose we observe a random sample of $n$ bivariate observations
$$
\bz_1=(x_1,y_1)^\top , \ldots , \bz_n=(x_n,y_n)^\top.
$$
If we are interested in exploring possible dependence between the $x_i$'s and $y_i$'s then among the first things we would do would be to obtain a scatterplot of the $x_i$'s against the $y_i$'s and calculate the correlation coefficient.  Recall that the sample  correlation coefficient is defined by
\begin{equation}
r=r[x,y]=\frac{n^{-1}\sum_{i=1}^n (x_i-\bar{x})(y_i-\bar{y})}{\left ( n^{-1}\sum_{i=1}^n (x_i-\bar{x})^2  \right )^{1/2}  \left ( n^{-1}\sum_{i=1}^n (y_i-\bar{y})^2 \right )^{1/2}}
(\#eq:scr)
\end{equation}
where $\bar{x}=n^{-1}\sum_{i=1}^n x_i$ and $\bar{y}=n^{-1}\sum_{i=1}^n y_i$ are the sample means.  Note that the sample correlation is a **scale-free measure** of the strength of **linear dependence** between the $x_i$'s and the $y_i$'s.

In this chapter we investigate the multivariate analogue of this question.  Suppose
$$
\bz_i =(\bx_i^\top,\by_i^\top)^\top, \qquad i=1,\ldots, n,
$$
 is a random sample of vectors. What is a sensible way to assess and describe the strength of the linear dependence between the $\bx_i$ vectors and the $\by_i$ vectors?  That is what this chapter is about.  A key role is played by the singular valued decomposition (SVD) introduced in Result 2.13 in Chapter 2.

\BeginKnitrBlock{example}
<span class="example" id="exm:prem"><strong>(\#exm:prem) </strong></span>From time to time we will return to the Premier League example in this chapter.  We shall treat $W$ and $D$, the number of wins and draws, respectively, as the $x$-variables; and $F$ and $A$, the number of goals for and against, will be treated as the $y$-variables.  The number of losses, $L$, is omitted as it provides no additional information when we know $W$ and $D$.  A question we shall consider is: how strongly associated are the match outcome variables, $W$ and $D$,  with the goals for and against variables, $F$ and $A$?
\EndKnitrBlock{example}
  
  
  
## Canonical Correlation Analysis

Assume we are given a random sample of vectors
$$
\bz_i=(\bx_i^\top , \by_i^\top )^\top: \, i=1,\ldots, n,
$$
where
the $\bx_i$ are $p \times 1$, the $\by_i$ are $q \times 1$ and, consequently, the $\bz_i$ are $(p+q)\times 1$. We are interested in determining the strength of linear association between the $\bx_i$ vectors and the $\by_i$ vectors.

Write
$$
\bar{\bz}=n^{-1}\sum_{i=1}^n \bz_i, \qquad \bar{\bx}=n^{-1} \sum_{i=1}^n \bx_i \qquad \text{and} \qquad \bar{\by}=n^{-1}\sum_{i=1}^n \by_i
$$
for the sample mean vectors of the $\bz_i$, $\bx_i$ and $\by_i$ respectively.

We formulate this task as an optimisation problem (cf. PCA). First, we introduce some notation. Let $\bS_{\bz\bz}$ denote the sample covariance matrix of the $\bz_i$, $i=1,\ldots, n$.  Then $\bS_{\bz\bz}$ can be written in block matrix form
$$
\bS_{\bz\bz}=\left [\begin{array}{cc}
\bS_{\bx \bx} & \bS_{\bx\by}\\
\bS_{\by \bx} & \bS_{\by \by} \end{array} \right ],
$$
where $\bS_{\bx \bx}$ ($p \times p$) is the sample covariance matrix of the $\bx_i$, $\bS_{\by \by}$ ($q \times q$) is the sample covariance of the $\by_i$, and the cross-covariance matrices are given by
$$
\stackrel{p \times q}{\bS}_{\bx \by}=n^{-1} \sum_{i=1}^n (\bx_i -\bar{\bx})(\by_i-\bar{\by})^\top
\qquad \text{and} \qquad \stackrel{q \times p}{\bS}_{\by \bx}=\bS_{\bx \by}^\top.
$$

**Example \@ref(exm:prem) (continued)**.  The relevant covariance matrix here is given in \@ref(eq:PLES),
but we need to delete the middle row and middle column because this relates to the variable $L$, the number of losses,
which we are omitting.  So we are left with
\begin{equation}
\bS_{\bx \bx}=\begin{pmatrix} 39.4 & -8.3\\ -8.3 & 8.1   \end{pmatrix} , \qquad
\bS_{\by \by}=\begin{pmatrix} 392.2 & -208.7\\ -208.7 & 230.9   \end{pmatrix}
(\#eq:bSxy1)
\end{equation}
and
\begin{equation}
\bS_{\bx \by}=\bS_{\by \bx}^\top =
\begin{pmatrix} 115.7  & -81.9\\ -29.4 & 6.0   \end{pmatrix}.
(\#eq:bSxy2)
\end{equation}
We shall return to this example in a little while.

We want to find the linear combination of the $x$-variables and the linear combination of the $y$-variables which is most highly correlated.

 One version of the optimisation problem we want to solve is: find non-zero vectors $\stackrel{p \times 1}{\ba}$ and $\stackrel{q \times 1}{\bb}$ which maximise the correlation coefficient
$$
r[\ba^\top \bx,\bb^\top \by]=\frac{\ba^\top \bS_{\bx \by}\bb}{(\ba^\top \bS_{\bx \bx}\ba)^{1/2}(\bb^\top \bS_{\by \by}\bb)^{1/2}}.
$$
In other words:
\begin{align}
  &\mbox{Find non-zero vectors }\quad  \ba \;\; (p \times 1)\mbox{ and  } \bb \;\; (q \times 1) \nonumber\\
  &\mbox{to maximise} \qquad  r[\ba^\top \bx,\bb^\top \by],
(\#eq:opt26)
\end{align}

where $r[.,.]$ is defined in \@ref(eq:scr).
  Intuitively, this makes sense, because we want to find the linear combination of the $x$-variables and the linear combination of the $y$-variables which are most highly correlated.

  However, note that for any $\gamma>0$ and $\delta>0$,
  \begin{equation}
  r[\gamma\ba^\top \bx, \delta \bb^\top \by]= \frac{\gamma \delta}{\sqrt{\gamma^2 \delta^2}}r[\ba^\top \bx,\bb^\top \by]=r[\ba^\top \bx,\bb^\top \by],
  (\#eq:invar)
  \end{equation}
  i.e. $r[\ba^\top \bx,\bb^\top \by]$ is invariant with respect to positive scalar multiplication of $\ba$ and $\bb$.  Consequently there will be an infinite number of solutions to this optimisation problem, because if $\ba$ and $\bb$ are solutions to optimization problem \@ref(eq:opt26), then so are $\gamma \ba$ and $\delta \bb$, for any $\gamma>0$ and $\delta>0$.

A more  useful way to formulate this optimisation problem is the following: find
\begin{equation}
\max_{\ba, \bb} \ba^\top \bS_{\bx \by}\bb
(\#eq:opt27a)
\end{equation}
subject to the constraints
\begin{equation}
\ba^\top \bS_{\bx \bx}\ba=1 \qquad \text{and} \qquad \bb^\top \bS_{\by \by}\bb=1.
(\#eq:opt27b)
\end{equation}

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-1"><strong>(\#prp:unnamed-chunk-1) </strong></span>Assume that $\bS_{\bx \bx}$ and $\bS_{\by \by}$ both are non-singular.  Then the following holds.

1. If $\ba=\hat{\ba}$ and $\bb=\hat{\bb}$ maximise \@ref(eq:opt26), then
$$
\ba=\check{\ba}\equiv\hat{\ba}/(\hat{\ba}^\top \bS_{\bx \bx}\hat{\ba})^{1/2} \qquad \text{and} \qquad
\bb=\check{\bb}\equiv \hat{\bb}/(\hat{\bb}^\top \bS_{\by \by}\hat{\bb})^{1/2}
$$
maximise  \@ref(eq:opt27a) subject to the constraints \@ref(eq:opt27b).  Moreover, if $\ba=\check{\ba}$ and $\bb=\check{\bb}$ maximise \@ref(eq:opt27a) subject to constraints
\@ref(eq:opt27b) then, for any $\gamma>0$ and $\delta>0$, $\ba=\gamma \check{\ba}$ and $\bb=\delta \check{\bb}$ maximise \@ref(eq:opt26).

2. The optimum solution to \@ref(eq:opt27a) and \@ref(eq:opt27b) is obtained when $\ba =\bS_{\bx \bx}^{-1/2}{\mathbf q}_1$ and $\bb=\bS_{\by \by}^{-1/2}
{\mathbf r}_1$, where $\bS_{\bx \bx}^{-1/2} \bS_{\bx \by}\bS_{\by \by}^{-1/2}$ has SVD
\begin{equation}
\bA\equiv \bS_{\bx \bx}^{-1/2}\bS_{\bx \by}\bS_{\by \by}^{-1/2}= \sum_{j=1}^t \xi_j {\mathbf q}_j {\mathbf r}_j^\top \equiv {\mathbf Q}{\pmb \Xi} {\mathbf R}^\top,
(\#eq:svdcca)
\end{equation}
where $\bA$ has rank $t$ and $\xi_1 \geq \cdots \geq \xi_t >0$.

3. The maximum value of the correlation coefficient is given by the largest singular value $\xi_1$.

\EndKnitrBlock{proposition}

Note: the matrix square roots $\bS_{\bx \bx}^{-1/2}$  and $\bS_{\by \by}^{-1/2}$ of $\bS_{\bx \bx}^{-1}$ and $\bS_{\by \by}^{-1}$, respectively, are defined using the definition of matrix square roots of symmetric non-negative definite matrices given in Chapter 2.

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}(i)  In \@ref(eq:invar) it was noted that,   for $\ba \neq {\mathbf 0}_p$ and $\bb \neq {\mathbf 0}_q$, the expression for $r[\ba^\top \bx, \bb^\top \by]$ is invariant when we change $\ba$ to $\gamma \ba$ and change $\bb$ to $\delta \bb$, where $\gamma>0$ and $\delta>0$ are scalars, so the second statement in Result 4.1(i) follows imnmdeiately.  Suppose now a solution to problem \@ref(eq:opt26) is achieved when $\ba = \hat{\ba}$ and $\bb=\hat{\bb}$.  Then, due to the invariance with respect to rescaling, the optimum is also achieved when $\ba =\check{\ba}\equiv\hat{\ba}/(\hat{\ba}^\top \bS_{\bx \bx} \hat{\ba})^{1/2}$ and $\bb =\check{\bb}\equiv \hat{\bb}/(\hat{\bb}^\top \bS_{\by \by} \hat{\bb})^{1/2}$. But by definition of $\check{\ba}$ and $\check{\bb}$, they satisfy the constraints \@ref(eq:opt27b) because
$$
\check{\ba}^\top \bS_{\bx \bx} \check{\ba}=\frac{\hat{\ba}^\top \bS_{\bx \bx}\hat{\ba}}{\left \{ \left (\hat{\ba}^\top \bS_{\bx \bx}\hat{\ba}\right )^{1/2}\right \}^2}
=\frac{\hat{\ba}^\top \bS_{\bx \bx}\hat{\ba}}{\hat{\ba}^\top \bS_{\bx \bx}\hat{\ba}}=1
$$
and, similarly,
$$
\check{\bb}^\top \bS_{\by \by} \check{\bb}=\frac{\hat{\bb}^\top \bS_{\by \by}\hat{\bb}}{\hat{\bb}^\top \bS_{\by \by}\hat{\bb}}=1.
$$
So $\ba=\check{\ba}$ and $\bb=\check{\bb}$ maximises \@ref(eq:opt27a) subject to the constraints \@ref(eq:opt27b).

\noindent (ii) \& (iii)  We may write the constraints \@ref(eq:opt27b) as
$$
\tilde{\ba}^\top \tilde{\ba}=1 \qquad \text{and} \qquad \tilde{\bb}^\top \tilde{\bb}=1
$$
where
$$
\tilde{\ba}=\bS_{\bx \bx}^{1/2} \ba \qquad \text{and} \qquad \tilde{\bb}=\bS_{\by \by}^{1/2}\bb.
$$

Recall that $\bS_{\bx \bx}$ and $\bS_{\by \by}$ are assumed to be non-singular.  Then, using results from Chapter 2, $\bS_{\bx \bx}^{1/2}$ and $\bS_{\by \by}^{1/2}$ will also be non-singular, and so
$$
(\bS_{\bx \bx}^{1/2})^{-1}=\bS_{\bx \bx}^{-1/2} \qquad \text{and} \qquad (\bS_{\by \by}^{1/2})^{-1}=\bS_{\by \by}^{-1/2}
$$
both exist and so we may write
$$
\ba=\bS_{\bx \bx}^{-1/2}\tilde{\ba} \qquad \text{and} \qquad \bb=\bS_{\by \by}^{-1/2} \tilde{\bb},
$$
and optimisation problem \@ref(eq:opt27a) subject to \@ref(eq:opt27b) becomes
$$
\max_{\tilde{\ba}, \tilde{\bb}}
\tilde{\ba}^\top \bS_{\bx \bx}^{-1/2}\bS_{\bx \by}\bS_{\by \by}^{-1/2} \tilde{\bb}
$$
subject to
$$
\vert \vert \tilde{\ba} \vert \vert =1 \qquad \text{and} \qquad \vert \vert \tilde{\bb}\vert \vert=1.
$$
From the properties of the SVD, and in particular Result 2.15  in Chapter 2, we know that the maximum correlation is $\xi_1$. Moreover, using the SVD again, this is achieved when $\tilde{\ba}={\mathbf q}_1$ and $\tilde{\bb}={\mathbf r}_1$ or, equivalently, $\ba=\bS_{\bx \bx}^{-1/2}{\mathbf q}_1$ and $\bb=\bS_{\by \by}^{-1/2}{\mathbf r}_1$.  
\EndKnitrBlock{proof}


**Example \@ref(exm:prem) (continued)** 
We now want to calculate the matrix $\bA$ in \@ref(eq:svdcca) and then find its singular valued decomposition.  We first need to find $\bS_{\bx \bx}^{-1/2}$ and $\bS_{\by \by}^{-1/2}$.
Using R to do the calculations, we obtain the following:
\begin{align*}
\bS_{\bx \bx}&=\bQ_{\bx} \bLambda_{\bx} \bQ_{\bx}^\top\\
&= \begin{pmatrix}  -0.970 & -0.241\\ 0.241 & -0.970 \end{pmatrix} \begin{pmatrix} 41.46 & 0 \\
 0 & 6.04\end{pmatrix} \begin{pmatrix} -0.970 & -0.241\\ 0.241 & -0.970  \end{pmatrix}^\top,
\end{align*}
and so
\begin{align*}
\bS_{\bx \bx}^{-1/2}&=\bQ_{\bx} \bLambda_{\bx}^{-1/2} \bQ_{\bx}^\top\\
&= \begin{pmatrix}  -0.970 & -0.241\\ 0.241 & -0.970 \end{pmatrix} \begin{pmatrix} 41.46^{-1/2} & 0 \\
 0 & 6.04^{-1/2}\end{pmatrix} \begin{pmatrix} -0.970 & -0.241\\ 0.241 & -0.970  \end{pmatrix}^\top\\
 &=\begin{pmatrix} 0.170 & 0.059\\0.059  &  0.392 \end{pmatrix};
\end{align*}
and, omitting details of the calculations this time,
$$
\bS_{\by \by}^{-1/2}=\bQ_{\by} \bLambda_{\by}^{-1/2} \bQ_{\by}^\top=\begin{pmatrix} 0.064 & 0.030\\0.030 &  0.086 \end{pmatrix}.
$$
Consequently,
\begin{align*}
\bA&=\bS_{\bx \bx}^{-1/2}\bS_{\bx \by}\bS_{\by \by}^{-1/2}\\
&=\begin{pmatrix} 0.170 & 0.059\\0.059  &  0.392 \end{pmatrix}
\begin{pmatrix} 115.7  & -81.9\\ -29.4 & 6.0   \end{pmatrix}
\begin{pmatrix} 0.064 & 0.030\\0.030 &  0.086 \end{pmatrix}\\
&=\begin{pmatrix} 0.741 & -0.628\\-0.374 &  -0.351\end{pmatrix}.
\end{align*}
The SVD of $\bA$ is given by
\begin{align}
\bA&=\bQ {\pmb \Xi} \bR^\top \nonumber \\
&=\begin{pmatrix} -0.997 & 0.082\\ 0.082 & 0.997   \end{pmatrix}
\begin{pmatrix}  0.974 & 0 \\0 & 0.508 \end{pmatrix}
\begin{pmatrix} -0.790 & -0.613\\ 0.613 & -0.790 \end{pmatrix}^\top.
(\#eq:SVDanalysis)
\end{align}
So the 1st CC coefficient is $0.974$, which is close to its maximum value of $1$.  The 1st CC weight vectors are
given by
\begin{align*}
\ba_1&=\bS_{\bx \bx}^{-1/2}\bq_1\\
&=\begin{pmatrix} 0.170 & 0.059\\0.059  &  0.392 \end{pmatrix} \begin{pmatrix} -0.997 \\ 0.082\end{pmatrix}\\
&=\begin{pmatrix}  -0.165, \\ - 0.027 \end{pmatrix}.
\end{align*}
Similar calculations show that
$$
\bb_1=\bS_{\by \by}^{-1/2}\br_1=\begin{pmatrix}  -0.032 \\ 0.029\end{pmatrix}.
$$

In order to make interpretation easier:


- We change $\ba_1$  to $-\ba_1$ and $\bb_1$ to $-\bb_1$.  [This entails changing $\bq_1$ to $-\bq_1$ and $\br_1$ to $-\br_1$; note that, provided we change the sign of **both** $\bq_1$ and $\br_1$, we do not change the matrix $\bA$.]
- We rescale $\ba_1$ and $\bb_1$ so that they are unit vectors.

This leads to the standardised 1st CC weight vectors
$$
\ba_1=\begin{pmatrix} 0.987\\0.160  \end{pmatrix} \qquad \text{and} \qquad
\begin{pmatrix} 0.743\\ -0.670\end{pmatrix}
$$
and the 1st CC variables, obtained by using these weights, are
$$
\eta_1 =0.987*(W-\bar{W}) +0.160*(D -\bar{D})
$$
and
$$
 \psi_1 = 0.743*(F-\bar{F}) - 0.670*(A-\bar{A}),
$$
where the bars are used to denote sample means.

We can see that $\psi_1$ is measuring something similar to goal difference $F-A$, as usually defined, but it gives slightly higher weight to goals scored than goals conceded ($0.743$ versus $0.670$).

It is also seen that $\eta_1$ is measuring something similar to number of points $3*W+D$, as usually defined, but the ratio of points for a win to points for a draw is somewhat higher, at around 6:1, as opposed to the usual ratio 3:1.




## The full set of canonical correlations

Let us first recap what we did in the previous section: we found the choices linear combinations of the $x$-variables and linear combinations of $y$-variables which
maximise the correlation, and expressed the answer in terms of quantities which arise in the SVD of $\bA$, where
$$
\bA\equiv \bS_{\bx \bx}^{-1/2} \bS_{\bx \by}\bS_{\by \by}^{-1/2}=\bQ {\pmb \Xi} \bR^\top=\sum_{j=1}^t \xi_j \bq_j \br_j^\top,
$$
with $t$ the rank of $\bA$, which in most examples is given by $t=\min(p,q)$, and singular values $\xi_1 \geq \xi_2 \geq \cdots \geq \xi_t>0$.
Specifically, the maximum value of the correlation is $\xi_1$, the optimal weights for the $x$-variables are given by $\ba=\bS_{\bx \bx}^{-1/2}\bq_1=\ba_1$, say, and
the optimal weights for the $y$-varables are given by $\bb=\bS_{\by \by}^{-1/2}\br_1 = \bb_1$, say.

Can we repeat this process, as we did with PCA?  Yes, we can.  To obtain the second canonical correlation coefficient, plus the associated sets of weights, we need to solve the following optimisation problem:
\begin{equation}
\max_{\ba,\, \bb} \ba^\top \bS_{\bx \by}\bb
(\#eq:cc2)
\end{equation}
subject to the constraints
\begin{equation}
\ba^\top \bS_{\bx \bx}\ba = 1, \qquad \bb^\top \bS_{\by \by}\bb=1,
(\#eq:conny21)
\end{equation}
\begin{equation}
\ba_1^\top \bS_{\bx \bx} \ba=0 \qquad \text{and} \qquad \bb_1^\top \bS_{\by \by}\bb=0.
(\#eq:conny22)
\end{equation}
Note that maximising \@ref(eq:cc2) subject to \@ref(eq:conny21) is very similar to the optimisation problem \@ref(eq:opt27a) and \@ref(eq:opt27b) considered in the previous section.  What is
new are the constraints \@ref(eq:conny22), which take into account that we have already found the first canonical correlation.  If for $j=1,2$ we write $\tilde{\ba}_j =\bS_{\bx \bx}^{1/2} \ba_j$ and $\tilde{\bb}_j=\bS_{\by \by}^{1/2} \bb_j$, then it is seen from \@ref(eq:conny22) that
$$
\tilde{\ba}_1^\top \tilde{\ba}_2=0 \qquad \text{and} \qquad \tilde{\bb}_1^\top \tilde{\bb}_2=0.
$$
Consequently, we may view constraints \@ref(eq:conny22) as corresponding to orthogonality constraints (cf. PCA) in  modified coordinate systems.

We now discuss the optimisation of \@ref(eq:cc2), \@ref(eq:conny21) and \@ref(eq:conny22).  At first glance it looks complex.  However, using arguments very similar to those used to prove
Result 2.15 in Chapter 2, we may deduce the following:

- The maximum of \@ref(eq:cc2) subject to constraints \@ref(eq:conny21) and \@ref(eq:conny22) is equal to $\xi_2$, the second largest singular value of $\bA$.
- The optimal weights for the $x$-variables for the second canonical correlation are given by $\ba_2=\bS_{\bx \bx}^{-1/2} \bq_2$.
- The optimal weights for the $y$-variables for the second canonical correlation are given by $\bb_2=\bS_{\by \by}^{-1/2}\br_2$.

Consider now the general case of the $k$th canonical correlation where $2 \leq k \leq t$.  In this case we replace \@ref(eq:conny21) and \@ref(eq:conny22) by, respectively,
\@ref(eq:connyk1) and \@ref(eq:connyk2) below, where
\begin{equation}
\ba^\top \bS_{\bx \bx}\ba= 1, \qquad \bb^\top \bS_{\by \by}\bb=1,
(\#eq:connyk1)
\end{equation}
\begin{equation}
\ba_j^\top \bS_{\bx \bx} \ba=0 \qquad \text{and} \qquad \bb_j^\top \bS_{\by \by}\bb=0, \qquad j=1, \ldots , k-1.
(\#eq:connyk2)
\end{equation}
Then the optimisation problem is
\begin{equation}
\max_{\ba, \, \bb} \ba^\top \bS_{\bx \by}\bb
(\#eq:cck)
\end{equation}
subject to constraints \@ref(eq:connyk1) and \@ref(eq:connyk2).  The solution in the general case is as follows.

- The maximum of \@ref(eq:cck) subject to constraints \@ref(eq:connyk1) and \@ref(eq:connyk2) is equal to $\xi_k$, the $k$th largest singular value of $\bA$.
- The optimal weights for the $x$-variables for the $k$th canonical correlation are given by $\ba_k=\bS_{\bx \bx}^{-1/2} \bq_k$.
- The optimal weights for the $y$-variables for the $k$th canonical correlation are given by $\bb_k=\bS_{\by \by}^{-1/2}\br_k$.


Terminology: we call $\ba_k$ and $\bb_k$ the $k$th cc (weight) vectors for the $x$-variables and $y$ variables, respectively.

We call $\eta_{ik}=\ba_k^\top (\bx_i - \bar{\bx})$ and $\psi_{ik}=\bb_k^\top (\by_i -\bar{\by})$, $i=1, \ldots , n$, the $k$th cc scores for the $x$-variables and the $y$-variables, respectively.

Define the CC score vectors ${\pmb \eta}_k=(\eta_{1k}, \ldots , \eta_{nk})^\top$ and ${\pmb \psi}_{k}=(\psi_{1k}, \ldots , \psi_{nk})^\top$.  Then we have the following result.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-3"><strong>(\#prp:unnamed-chunk-3) </strong></span>Assume that $\bS_{\bx \bx}$ and $\bS_{\by \by}$ both have full rank.  Then for $1 \leq k,\ell \leq t$,
$$
r[\eta_k,  \psi_{\ell}]=\begin{cases} \xi_k &\text{if} \quad k=\ell\\
0 & \text{if} \quad k \neq \ell, \end{cases}
$$
where $t$ is the rank of $\bA=\bS_{\bx \bx}^{-1/2}\bS_{\bx \by} \bS_{\by \by}^{-1/2}$ and $\xi_1 \geq \xi_2 \geq \cdots \xi_t > 0$ are the strictly positive singular values of $\bA$.
\EndKnitrBlock{proposition}


**Example \@ref(exm:prem) (continued)**  From \@ref(eq:SVDanalysis), it is seen that the 2nd CC coefficient is given by $\xi_2=0.508$.  So the correlation between the second pair of CC variables is a lot smaller than the 1st CC coefficient, though still appreciably different from $0$.  We now calculate the 2nd CC weight vectors:
$$
\ba_2=\bS_{\bx \bx}^{-1/2} \bq_2 = \begin{pmatrix} 0.073 \\ 0.396 \end{pmatrix}
\qquad \text{and} \qquad
\bb_2=\bS_{\by \by}^{-1/2}\br_2=-\begin{pmatrix}0.062\\ 0.086  \end{pmatrix},
$$
with standardised version (without the sign changes this time)
$$
\ba_2=\begin{pmatrix}0.181 \\ 0.984  \end{pmatrix}
\qquad \text{and} \qquad
\bb_2=-\begin{pmatrix}0.589 \\ 0.808  \end{pmatrix},
$$
and new variables
$$
\eta_2=0.181*(W-\bar{W}) +0.984*(D -\bar{D})
$$
and
$$
\psi_2=-\{0.589*(F-\bar{F})+0.808*(A-\bar{A})\}.
$$
Note that, to a good approximation, $\eta_2$ is measuring something similar to the number of draws and, approximately, $\psi_2$ is something related to the negative of total number of goals in a team's games.  So large $\psi_2$ means relatively few goals in a team's games, and small (i.e. large negative) $\psi_2$ means a relatively large number of goals in a team's games.

Interpretation of the 2nd CC: teams that have a lot of draws tend to be in low-scoring games and/or teams that have few draws tend to be in high-scoring games.

## Connection  with linear regression when $q=1$

Although CCA analysis is clearly a different technique to linear regression, it turns out that when either $p=1$ or $q=1$, there is a close  connection between the two approaches.

Without loss of generality we assume that $q=1$ and  $p>1$.  Hence there is only a single $y$-variable but we still have $p>1$ $x$-variables.

We also make the following assumptions:

1. The $\bx_i$ have been centred so that $\bar{\bx}={\mathbf 0}_p$, the zero vector.
2. The covariance matrix for the $x$-variables, $\bS_{\bx \bx}$, has full rank $p$.

Both of these are weak assumptions in the multiple linear regression context.

Since $q=1$,
$$
\bA=\bS_{\bx \bx}^{-1/2} \bS_{\bx y}\bS_{yy}^{-1/2}
$$
is a $p \times 1$ vector.  Consequently, in this rather special case,
the SVD tells us that
$$
\bA=\xi_1 \bq_1,
$$
where
$$
\xi_1=\vert \vert \bA \vert \vert \qquad \text{and} \qquad \bq_1=\bA /\vert \vert \bA \vert \vert=\tilde{\ba},
$$
and $\tilde{\ba}=\bS_{\bx \bx}^{1/2} \ba$.

Consequently,
\begin{align*}
\ba&=\bS_{\bx \bx}^{-1/2}\bq_1\\
&=\bS_{\bx \bx}^{-1/2} \frac{1}{\vert \vert \bS_{\bx \bx}^{-1/2}\bS_{\bx y}S_{yy}^{-1/2}\vert \vert}\bS_{\bx \bx}^{-1/2}\bS_{\bx \by}S_{yy}^{-1/2}\\
&=\frac{1}{\vert \vert \bS_{\bx \bx}^{-1/2}\bS_{\bx y}\vert \vert}\bS_{\bx \bx}^{-1/2}\bS_{\bx \bx}^{-1/2}\bS_{\bx y}\\
&=\frac{1}{\vert \vert \bS_{\bx \bx}^{-1/2}\bS_{\bx y}\vert \vert}\bS_{\bx \bx}^{-1}\bS_{\bx y}.
\end{align*}
But since $\bar{\bx}={\mathbf 0}_p$  and $\bS$ has full rank by the assumptions above, it follows that
$$
n\bS_{\bx \bx} =\sum_{i=1}^n \bx_i \bx_i^\top =\bX^\top \bX
$$
and
$$
 n\bS_{\bx \by}=\sum_{i=1}^n y_i \bx_i=\bX^\top \by,
$$
where $\by=(y_1, \ldots ,y_n)^\top$ is the $n \times 1$ data matrix for the $y$-variable and $\bX=[\bx_1, \ldots , \bx_n]^\top$ is the data matrix for the $x$-variables.
Consequently, the optimal $\ba$ is a scalar multiple of
$$
\bS_{\bx \bx}^{-1}\bS_{\bx y}=\left ( \bX^\top \bX \right )^{-1} \bX^\top \by=\hat{\pmb \beta},
$$
say, which is the classical expression for least squares estimator.  Therefore the least squares estimator $\hat{\pmb \beta}$ solves \@ref(eq:opt26).  However, it does not usually solve the optimisation problem defined by problems \@ref(eq:opt27a) and \@ref(eq:opt27b) because typically it will not be the case that $\hat{\pmb \beta}^\top \bS_{\bx \bx}\hat{\pmb \beta}=1$, so that \@ref(eq:opt27b) will  not be satisfied.


## Population CCA

So far in this chapter we have based CCA on the sample covariance matrix
$$
\bS_{\bz\bz}=\left [\begin{array}{cc}
\bS_{\bx \bx} & \bS_{\bx\by}\\
\bS_{\by \bx} & \bS_{\by \by} \end{array} \right ],
$$
However, just as there is a population analogue of PCA, so there is a population analogue of CCA.

Given random vectors $\stackrel{p \times 1}{\bx}$ and $\stackrel{q \times 1}{\by}$, define the random vector $\bz=(\bx^\top, \by^\top)^\top$ with population covariance matrix
$$
\text{Var}(\bz)=\bSigma_{\bz\bz}=\left [\begin{array}{cc}
\bSigma_{\bx \bx} & \bSigma_{\bx\by}\\
\bSigma_{\by \bx} & \bSigma_{\by \by} \end{array} \right ].
$$
Then, by analogy with what we have seen in the sample CCA, the population CCA is based on the
matrix
$$
\check{\bA}=\bSigma_{\bx \bx}^{-1/2}\bSigma_{\bx \by}\bSigma_{\by \by}^{-1/2},
$$
where, as in \S 3.4,  the check symbol has been used above and below to indicate population quantities.
If $\check{\bA}$ has SVD
$$
\check{\bA}=\sum_{j=1}^t \check{\xi}_j\check {\mathbf q}_j \check{\mathbf r}_j^\top \equiv \check{\mathbf Q}\check{\pmb \Xi} \check{\mathbf R}^\top,
$$
where $\check{\xi}_1 \geq \cdots \geq \check{\xi}_t \geq 0$ and $t=\min(p,q)$, and the $\check{\bq}_j$ and $\check{\mathbf r}_j$ are unit vectors, then the first population CC coefficient is given by $\check{\xi}_1$,
and the associated weights are given by
$$
\check{\ba}=\bSigma_{\bx \bx}^{-1/2}\check{\bq}_1=\check{\ba}_1 \qquad \text{and} \qquad \check{\bb}=\bSigma_{\by \by}^{-1/2}\check{\mathbf r}_1=\check{\bb}_1.
$$
The full set of population CC weight vectors is given by
$$
\check{\ba}_j =\bSigma_{\bx \bx}^{-1/2}\check{\bq}_j \qquad \text{and} \qquad
\check{\bb}_j=\bSigma_{\by \by}^{-1/2}\check{\mathbf r}_1, \qquad , j=1, \ldots , t,
$$
and the $j$th population CC coefficient is given by $\check{\xi}_j$.

## Invariance/equivariance properties of CCA

Suppose we apply  orthogonal transformations and translations to the $\bx_i$ and the $\by_i$ of the form
\begin{equation}
{\mathbf h}_i={\mathbf T}\bx_i + {\pmb \mu} \qquad \text{and} \qquad {\mathbf k}_i={\mathbf V}\by_i +{\pmb \eta},
\qquad i=1,\ldots , n,
(\#eq:transformations)
\end{equation}
where $\mathbf T$ ($p \times p$) and $\mathbf V$ ($q \times q$) are orthogonal matrices, and $\pmb \mu$ ($p \times 1$) and
$\pmb \eta$ ($q \times 1$) are fixed vectors.

How do these transformations affect the CC analysis?

First of all, since the CCA depends only on sample covariance matrices, it follows that the translation vectors $\pmb \mu$ and $\pmb \eta$ have no effect on the analysis, so we can ignore $\pmb \mu$ and $\pmb \eta$, and without loss of generality we shall set each to be the zero vector.

As seen in the previous section, the CCA in the original coordinates depends on
\begin{equation}
\bA \equiv \bA_{\bx \by}=\bS_{\bx \bx}^{-1/2}\bS_{\bx \by}\bS_{\by \by }^{-1/2}.
(\#eq:Axy)
\end{equation}
In the new coordinates we have
$$
\tilde{\bS}_{\mathbf h h}={\mathbf T} \bS_{\bx \bx}{\mathbf T}^\top, \qquad \tilde{\bS}_{\mathbf kk}={\mathbf V}\bS_{\by \by}{\mathbf V}^\top,
$$
$$
\tilde{\bS}_{\mathbf hk}={\mathbf T}\bS_{\bx \by}{\mathbf V}^\top \qquad \text{and} \qquad
\tilde{\bS}_{\mathbf kh}={\mathbf V}\bS_{\by \bx}{\mathbf T}^\top=\bS_{\mathbf h k}^\top,
$$
where here and below, a tilde above a symbol is used to indicate that the corresponding term is defined in terms of the new $\bh$, $\bk$ coordinates, rather
than the old $\bx$, $\by$ coordinates.
Moreover, due to the fact that $\mathbf T$ and $\mathbf V$ are orthogonal,
$$
\tilde{\bS}_{\mathbf hh}^{ 1/2}={\mathbf T}\bS_{\bx \bx}^{ 1/2}{\mathbf T}^\top, \qquad
\tilde{\bS}_{\mathbf hh}^{ -1/2}={\mathbf T}\bS_{\bx \bx}^{ -1/2}{\mathbf T}^\top
$$
$$
\tilde{\bS}_{\mathbf kk}^{ 1/2}={\mathbf V}\bS_{\by \by}^{ 1/2}{\mathbf V}^\top \qquad \text{and} \qquad
\tilde{\bS}_{\mathbf kk}^{ -1/2}={\mathbf V}\bS_{\by \by}^{- 1/2}{\mathbf V}^\top.
$$
The analogue of \@ref(eq:Axy) in the new coordinates is given by
\begin{align*}
\tilde{\bA}_{\mathbf h k}&=\tilde{\bS}_{\mathbf hh}^{-1/2}\tilde{\bS}_{\mathbf h k}\tilde{\bS}_{\mathbf kk}^{-1/2}\\
&={\mathbf T} \bS_{\bx \bx}^{-1/2}{\mathbf T}^\top {\mathbf T}\bS_{\bx \by}{\mathbf V}^\top {\mathbf V}\bS_{\by \by}^{-1/2}{\mathbf V}^\top\\
&={\mathbf T}\bS_{\bx \bx}^{-1/2}\bS_{\bx \by}\bS_{\by \by }^{-1/2}{\mathbf V}^\top\\
&={\mathbf T} \bA_{\bx \by}{\mathbf V}^\top.
\end{align*}
So, again using the fact that  $\mathbf T$ and $\mathbf V$ are orthogonal matrices, if $\bA_{\bx \by}$ has SVD $\sum_{j=1}^t \xi_j {\mathbf q}_j {\mathbf r}_j^\top$, then $\tilde{\bA}_{\mathbf hk}$ has SVD
\begin{align*}
\tilde{\bA}_{\mathbf hk}&={\mathbf T }\bA_{\bx \by}{\mathbf V}^\top
={\mathbf T} \left ( \sum_{j=1}^t \xi_j {\mathbf q}_j {\mathbf r}_j^\top \right){\mathbf V}^\top\\
&=\sum_{j=1}^t \xi_j {\mathbf T}{\mathbf q}_j {\mathbf r}_j^\top {\mathbf V}^\top
=\sum_{j=1}^t \xi_j \left ( {\mathbf T} {\mathbf q}_j \right )\left ({\mathbf V}{\mathbf r}_j  \right )^\top
=\sum_{j=1}^t \xi_j \tilde{\bq}_j \tilde{\mathbf r}_j^\top,
\end{align*}
where, for $j=1, \ldots,t$, the  $\tilde{\bq}_j={\mathbf T}\bq_j$ are mutually orthogonal unit vectors,
and the $\tilde{\mathbf r}_j={\mathbf V}{\mathbf r}_j$ are also mutually orthogonal unit vectors.

Consequently, $\tilde{\bA}_{\mathbf h k}$ has the same singular values as $\bA_{\bx \by}$, namely $\xi_1, \ldots , \xi_t$ in both cases, and so the  canonical correlation coefficients are invariant with respect to the transformations \@ref(eq:transformations).  Moreover, since the optimal linear combinations  for the $j$th CC in the original coordinates are given by $\ba_j =\bS_{\bx \bx}^{-1/2}{\mathbf q}_j$ and $\bb_j=\bS_{\by \by}^{-1/2}{\mathbf r}_j$, the optimal linear combinations in the new coordinates are given by
\begin{align*}
\tilde{\ba}_{j}&=\bS_{\mathbf hh}^{-1/2}{\mathbf T}{\mathbf q}_j\\
&={\mathbf T}\bS_{\bx \bx}^{-1/2}{\mathbf T}^\top {\mathbf T}{\mathbf q}_j\\
&={\mathbf T}\bS_{\bx \bx}^{-1/2}{\mathbf q}_j \\
&={\mathbf T}\ba_{j},
\end{align*}
and a similar argument shows that $\tilde{\bb}_{j}={\mathbf V}\bb_{j}$.  So under transformations \@ref(eq:transformations),
the optimal vectors $\ba_{j}$ and $\bb_{j}$ transform in an equivariant manner to $\tilde{\ba}_{j}$ and $\tilde{\bb}_{j}$, respectively, $j=1, \ldots , t$.

If either of $\mathbf T$ or $\mathbf V$ in \@ref(eq:transformations) is not an orthogonal matrix then the singular values are not invariant and the cc vectors do not transform in an equivariant manner.

## Testing for zero canonical correlation coefficients

So far in Part II of this module we have not considered formal statistical inference (e.g. hypothesis testing, construction of confidence regions).  Inference in various multivariate settings is considered in Part III.  However, before moving on, we briefly explain how to perform tests for zero correlations in the CCA setting, under the assumption that the $\bz_i = (\bx_i^\top , \by_i^\top)^\top$ are IID multivariate normal.

As previously, suppose that the $\bx_i$ are $p \times 1$ vectors and the $\by_i$ are $q \times 1$ vectors and the sample size, i.e.  the number of $\bz_i$ vectors, is $n$.  Let $\bSigma_{\bx \by} =\text{Cov}(\bx,\by)$ denote the population cross-covariance matrix as before and consider the null hypothesis
$$
H_0: \, \bSigma_{\bx \by}={\mathbf 0}_{p,q},
$$
i.e. $\bSigma_{\bx \by}$ is the $p \times q$ matrix of zeros.  Let $H_A$ denote the general alternative
$$
H_A:\, \bSigma_{\bx \by} \quad *unrestricted*.
$$

Then the large-sample log-likelihood ratio test statistic for testing $H_0$ versus $H_A$ is as follows:
$$
W_0=-\left \{n-\frac{1}{2}(p+q+3)  \right \}\sum_{j=1}^{\min(p,q)} \log (1-\xi_j^2),
$$
where $\xi_1\geq \xi_2 \cdots \geq \xi_{\min(p,q)} \geq 0$ are the sample canonical correlations.
Moreover, when $n$ is large, $W_0$ is approximately $\chi_{pq}^2$ under $H_0$,  and $H_0$ should be rejected
when $W_0$ is sufficiently large.

We now consider a test concerning the rank of $\bSigma_{\bx \by}$.  For $0 \leq t <\min(p,q)$, consider  the hypothesis:
$$
H_t: \,   \text{at most $t$ of the CC coefficients are non-zero}.
$$
It turns out there is a similar statistic to $W_0$ above, for testing $H_t$ against $H_A$, defined by
$$
W_t=-\left \{n-\frac{1}{2}(p+q+3)  \right \}\sum_{j=t+1}^{\min(p,q)} \log (1-\xi_j^2),
$$
where, under $H_t$ with $n$ large, $W_t$ is approximately $\chi_{(p-t)(q-t)}^2$.  Also, we reject
$H_t$ when $W_t$ is sufficiently large.


**Example \@ref(exm:prem) (continued)**.
  Here $p=q=2$, $n=20$ and $\xi_1=0.974$ and $\xi_2=0.508$.
So we should refer $W_0$ to $\chi_4^2$ and refer $W_1$ to $\chi_1^2$.  Here, $W_0=53.92$ and $W_1=4.92$.  So hypothesis $H_0$ is strongly rejected, with $p$-value $<\,<0.001$.  In contrast, $H_1$ is rejected at the $0.05$ level but is not rejected at the $0.01$ level.  So there is only moderate evidence that the 2nd CC coefficient is non-zero.

