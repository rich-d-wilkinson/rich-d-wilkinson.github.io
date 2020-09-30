# PART II: Dimension reduction methods {-}


In applications of statistics in many different fields it is common to measure several (or even a large number) of variables on each experimental unit under study.  For example,  experimental units could be individual people and variables could be measurements obtained in a general health check-up (e.g. age, blood pressure, cholesterol level, lung function measurements, BMI and other variables).

When analysing data of moderate or high dimension, it is often desirable to seeks ways to restructure the data and reduce its dimension *in such a way that we retain the most important information within the data*.  In reduced dimensions it is often much easier to understand and appreciate the most important features of a dataset.

In Part II of this module we investigate three different methods for dimension reduction: Principal Components Analysis (PCA) in Chapter \@ref(pca); Canonical Correlation Analysis (CCA) in Chapter \@ref(cca); and Multidimensional Scaling (MDS) in Chapter \@ref(mds).

Matrix algebra (see Chapter \@ref(linalg-prelim)) plays a key role in all three of these techniques.


# Principal component analysis {#pca}

Although it is very common to collect multivariate data, we often want to reduce the dimension of such data *in a sensible way*.

For example, exam marks across different modules are
averaged to produce a single overall mark for each
student.  Similarly, in a football league table we convert the
numbers of wins, draws and losses to a single measure of
points.

Mathematically, we can express these examples of
dimension reduction as a linear combination of the
original variables, $y = \bu^\top \bx$.  For the exam mark
example, suppose each student sits $p=4$ modules
with marks, $x_1,x_2,x_3,x_4$.  Then, writing $\bx=(x_1, x_2 , x_3, x_4)^\top$ and choosing $\bu = \lb \frac{1}{4}, \frac{1}{4}, \frac{1}{4}, \frac{1}{4} \rb ^\top$
 gives an overall average,
$$ y =\bu^\top \bx= \begin{pmatrix} \frac{1}{4} & \frac{1}{4} & \frac{1}{4} & \frac{1}{4} \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \frac{x_1}{4} + \frac{x_2}{4} + \frac{x_3}{4} + \frac{x_4}{4}.$$
For the football league table, if $w$ is the number of wins, $d$ is the number of draws and $l$ is the number of losses then, writing
 ${\mathbf r}=(w,d,l)^\top$, we choose $\bu = \lb 3,1,0 \rb^\top$ to get the points score
$$ y = \bu^\top {\mathbf r}=\begin{pmatrix} 3 & 1 & 0 \end{pmatrix} \begin{pmatrix} w \\ d \\ l \end{pmatrix} = 3w + 1d + 0l=3w+d.$$

In these examples we use $\bu$ to convert our original variables, the components of $\bx$, to a new variable, $y$.  These choices of $\bu$ are fairly standard for these types of data.  However,we should ask whether we can do better.  In a more general setting, how should we choose $\bu$?

 A key objective of principal component analysis (PCA): to find the linear combination of the original variables that **maximises the variability** in the new variable.  Intuitively, this seems sensible for the exam mark data because a large variance in $y$ would separate out the better students from the weaker students, making it easier to rank them.

## Principal component vectors and scores

Let $\bx_1,\ldots,\bx_n$ be $p \times 1$ vectors of measurements on $n$ experimental units with sample mean $\bar{\bx} = \frac{1}{n} \sum_{i=1}^n \bx_i$ and sample covariance matrix\\
 $\bS = \frac{1}{n} \sum_{i=1}^n (\bx_i - \bar{\bx}) (\bx_i - \bar{\bx})^\top$.

We wish to project the data onto a lower-dimensional subspace in which the data
displays *maximal variation*, using appropriate scalar products of the observation
vectors.

Let $\bu$ be a unit vector (i.e. $\| \bu \| = 1$ or $\bu^\top \bu=1$) and define
$$y_i= \bu^\top (\bx_i - \bar{\bx})$$
for $i=1,\ldots,n$.

Now
$$ \sum_{i=1}^n y_i = \sum_{i=1}^n \bu^\top (\bx_i - \bar{\bx})
= \bu^\top \sum_{i=1}^n (\bx_i - \bar{\bx})
= \bu^\top (n \bar{\bx} - n \bar{\bx}) = 0,$$
by the definition of $\bar{\bx}$, so $\bar{y} = \frac{1}{n} \sum_{i=1}^n y_i = 0$.

The sample variance of the $y_i$'s is
\begin{eqnarray*}
s^2[\bu] &=& \frac{1}{n} \sum_{i=1}^n (y_i - \bar{y})^2 = \frac{1}{n} \sum_{i=1}^n y_i^2 \\
&=& \frac{1}{n} \sum_{i=1}^n \lsb \bu^\top (\bx_i - \bar{\bx}) \rsb \lsb (\bx_i - \bar{\bx})^\top \bu \rsb \\
&=& \bu^\top \lsb \frac{1}{n} \sum_{i=1}^n (\bx_i - \bar{\bx})(\bx_i - \bar{\bx})^\top \rsb \bu \\
&=& \bu^\top \bS \bu.
\end{eqnarray*}

We would like to find the $\bu$ which maximises the sample variance, $s^2[\bu] = \bu^\top \bS \bu$ over unit vectors $\bu$.

Since $\bS$ is symmetric, then by the spectral decomposition theorem we can write
$$\bS = \bQ \bLambda \bQ^\top = \sum_{j=1}^p \lambda_j \bq_j \bq_j^\top $$
with $\bQ = [ \bq_1,  \ldots , \bq_p ]$ an orthogonal matrix (so $\bQ \bQ^\top = \bQ^\top \bQ = \bI_p$) and $\bLambda = \tdiag  \{ \lambda_1, \ldots, \lambda_p \}$ where we may assume $\lambda_1  \geq \cdots \geq \lambda_p$ and, since $\bS$ is a covariance matrix and therefore non-negative definite, $\lambda_p \geq 0$. Note that $\lambda_j$ and $\bq_j$, $j=1,\ldots,p$, are eigenvalues and eigenvectors, respectively, of $\bS$.

Then,
\begin{eqnarray*}
s^2[\bu] &=& \bu^\top \bS \bu = \bu^\top \bQ \bLambda \bQ^\top \bu
= \bu^\top \lb \sum_{j=1}^p \lambda_j \bq_j \bq_j^\top \rb \bu \\
&=& \sum_{j=1}^p \lambda_j (\bu^\top \bq_j) (\bq_j^\top \bu)
= \sum_{j=1}^p \lambda_j (\bu^\top \bq_j)^2 \\
&\leq& \sum_{j=1}^p \lambda_1 (\bu^\top \bq_j)^2
\end{eqnarray*}
since $\lambda_1 \geq \lambda_j, j=1,\ldots,p$.  Therefore, using Proposition \@ref(prp:two1),
$$ s^2[\bu] \leq \lambda_1 \sum_{j=1}^p (\bu^\top \bq_j)^2
= \lambda_1 \bu^\top \left ( \sum_{j=1}^p \bq_j \bq_j^\top  \right) \bu
= \lambda_1 \bu^\top \bu =\lambda_1,$$
since, by assumption, $\| \bu \| = 1$.

Therefore, the maximum $s^2[\bu]$ is at most $\lambda_1$, where $\lambda_1$ is the largest eigenvalue of $\bS$.

Recall that
$$ \bq_i^\top \bq_j = \left\{ \begin{array}{ll} 0 & \text{if } j \neq i,\\ 1 & \text{if } j=i. \end{array} \right.$$
because eigenvectors are orthogonal to each other, so if we take $\bu = \bq_1$ then
\begin{eqnarray*}
\bq_1^\top \bS \bq_1 &=& \bq_1^\top \lb \sum_{j=1}^p \lambda_j \bq_j \bq_j^\top \rb \bq_1
= \sum_{j=1}^p \lambda_j (\bq_1^\top \bq_j) (\bq_j^\top \bq_1) \\
&=& \sum_{j=1}^p \lambda_j (\bq_1^\top \bq_j)^2
= \lambda_1 (\bq_1^\top \bq_1)^2
= \lambda_1
\end{eqnarray*}


So $s^2[\bu] = \bu^\top \bS \bu$ is maximised over unit vectors $\bu$ when $\bu = \bq_1$ where $\bq_1$ is the unit eigenvector corresponding to the largest eigenvalue, $\lambda_1$.  By maximising $\bu^\top \bS \bu$ over unit vectors $\bu$, we are in effect choosing a projection onto a 1-dimensional subspace which captures as much of the sample variation as possible.

We can repeat this procedure and look for the largest sample variance of the $y_i$'s, when
$\bu$ is chosen to be orthogonal to $\bq_1$ (i.e. restrict attention to those $\bu$ such that $\bu^\top \bq_1 = 0$).  Similar reasoning shows that this constrained maximum occurs when $\bu = \bq_2$, where $\bq_2$ is
the eigenvector corresponding to the second largest eigenvalue, $\lambda_2$; and the corresponding maximum of $\bu^\top \bS \bu$ is $\lambda_2$.

We can repeat the process for $j=1,\ldots,p$ to define $p$ new variables.  In general, to find PC $j$, we solve the following optimisation problem:
\begin{equation}
\max_{\bu: \, \vert \vert \bu \vert \vert =1}\bu^\top \bS \bu
(\#eq:pcmaxgen)
\end{equation}
subject to
\begin{equation}
\bq_k^\top \bu =0, \qquad k=1, \ldots , j-1.
(\#eq:pccongen)
\end{equation}
It turns out that the maximum of \@ref(eq:pcmaxgen)
subject to \@ref(eq:pccongen) is equal to $\lambda_j$ and is obtained when $\bu=\bq_j$.

The 1st PC scores are $y_{i1} = \bq_1^\top (\bx_i - \bar{\bx}), \quad i=1,\ldots,n$. \\
The 2nd PC scores are $y_{i2} = \bq_2^\top (\bx_i - \bar{\bx}), \quad i=1,\ldots,n$.
$$ \vdots $$
The $p$th PC scores are $y_{ip} = \bq_p^\top (\bx_i - \bar{\bx}), \quad i=1,\ldots,n$.

We summarise these findings in the following result.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:pca1"><strong>(\#prp:pca1) </strong></span>Let $\bx_1, \ldots , \bx_n$ denote a sample of vectors in $\mathbb{R}^p$ with sample mean vector $\bar{\bx}$ and sample covariance matrix $\bS$.  Suppose $\bS$ has spectral decomposition (see Proposition \@ref(prp:spectraldecomp))
$$
\bS=\bQ \bLambda \bQ^\top = \sum_{j=1}^p  \lambda_j \bq_j \bq_j^\top,
$$
where $\bQ$ is orthogonal, $\bLambda=\text{diag}\{\lambda_1, \ldots, \lambda_p\}$ and $\lambda_1 \geq \lambda_2 \geq \lambda_p \geq 0$.  Then the following holds:



1. The maximum of \@ref(eq:pcmaxgen)
subject to \@ref(eq:pccongen) is equal to $\lambda_j$ and is obtained when $\bu=\bq_j$.

2.  For $j=1, \ldots , p$, the scores of the $j$th principal component (PC)  are given  by
$$
y_{ij}=\bq_j^\top(\bx_i - \bar{\bx}), \qquad i=1, \ldots , n,
$$
where $\bq_j$ is the vector of *loadings* for the $j$th PC.  Moreover,
$$
\by_i=( y_{i1}, y_{i2}, \ldots , y_{ip})^\top = \bQ^\top (\bx_i -\bar{\bx}), \qquad i=1, \ldots ,n.
$$

3.  In matrix form, the full set of PC scores is given in the matrix
$$
\bY = [\by_1 , \ldots , \by_n]^\top =\bH \bX \bQ,
$$
where $\stackrel{n \times n}{\bH}$ is the centering matrix and $\bX=[\bx_1, \ldots , \bx_n]^\top$ is the original  data matrix.

4.  The sample mean vector of $\by_1, \ldots , \by_n$ is the zero vector ${\mathbf 0}_p$ and the sample covariance matrix is $\bLambda$.

\EndKnitrBlock{proposition}


\BeginKnitrBlock{example}
<span class="example" id="exm:unnamed-chunk-1"><strong>(\#exm:unnamed-chunk-1) </strong></span>We consider the marks of $n=10$ students who studied G11PRB and G11STA.

\EndKnitrBlock{example}

```
## Warning: package 'kableExtra' was built under R version 3.6.2
```

\begin{table}[H]
\centering
\begin{tabular}{rrr}
\toprule
student & PRB & SMM\\
\midrule
1 & 81 & 75\\
2 & 79 & 73\\
3 & 66 & 79\\
4 & 53 & 55\\
5 & 43 & 53\\
\addlinespace
6 & 59 & 49\\
7 & 62 & 72\\
8 & 79 & 92\\
9 & 49 & 58\\
10 & 55 & 56\\
\bottomrule
\end{tabular}
\end{table}


The sample mean vector  and sample covariance matrix are
$$
\bar{\bx} = \begin{pmatrix} 62.6 \\ 66.2 \end{pmatrix}\qquad \text{and} \qquad \bS = \begin{pmatrix} 162.04 & 135.38 \\ 135.38 & 175.36 \end{pmatrix}.
$$

```r
library(dplyr)
```

```
## Warning: package 'dplyr' was built under R version 3.6.2
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following object is masked from 'package:kableExtra':
## 
##     group_rows
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
secondyr %>% select(2:3) %>% colMeans -> xbar
secondyr %>% select(2:3) %>% cov(use="everything")*9/10 -> S
eigs = eigen(S)
```

DELETE THIS - ASSUME THEY CAN DO IT, OR DO ON A COMPUTER.

To find the eigenvalues we need to solve $|\bS - \lambda \bI| = 0$, where
\begin{eqnarray*}
|\bS - \lambda \bI_2| &=& (162.04-\lambda)(175.36-\lambda) - 135.38^2 \\
&=& \lambda^2 - 337.4 \lambda + 10887.59.
\end{eqnarray*}



Using the quadratic equation formula we find,
$$ \lambda = \frac{337.4 \pm \sqrt{337.4^2 - 4(10087.59)}}{2} = \frac{337.4 \pm \sqrt{73488.4}}{2}. $$
So $\lambda_1 = 304.24$ and $\lambda_2 = 33.16$.

To find the first eigenvector we solve $(\bS - \lambda_1 \bI_2) \bq_1 = \bzero$.  To simplify, we use row operations:
\begin{eqnarray*}
\bS - \lambda_1 \bI_2 = \begin{pmatrix} -142.20 & 135.38 \\ 135.38 & -128.88 \end{pmatrix} &\rightarrow& \begin{pmatrix} 1 & -0.952 \\ 135.38 & -128.88 \end{pmatrix} \\
&\rightarrow& \begin{pmatrix} 1 & -0.952 \\ 0 & 0 \end{pmatrix}.
\end{eqnarray*}
If we let $\bq_1 = (q_{11}, q_{21})^\top$ then solving $(\bS - \lambda_1 \bI_2) \bq_1 = \bzero$ is equivalent to solving
$$ \begin{pmatrix} 1 & -0.952 \\ 0 & 0 \end{pmatrix} \begin{pmatrix} q_{11} \\ q_{21} \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}. $$
So $q_{11} = 0.952 q_{21}$ and the eigenvectors are of the form $t\begin{pmatrix} 0.952 \\ 1 \end{pmatrix}$ where $t \ne 0$ is a constant.  We choose $t$ such that $\| \bq \| = 1$, so
$$
{\ds t = \pm \frac{1}{\sqrt{0.952^2 + 1^2}} = \pm 0.724}.
$$
Therefore,
$$ \bq_1 = 0.724 \begin{pmatrix} 0.952 \\ 1 \end{pmatrix} = \begin{pmatrix} 0.690 \\ 0.724 \end{pmatrix}.$$

To find the second eigenvector we use the same method to solve $(\bS - \lambda_2 \bI_2) \bq_2 = \bzero$ and find that $\bq_2 = \begin{pmatrix} -0.724 \\ 0.690 \end{pmatrix}$.

The plot below shows the original data.  The two lines, centred on $\bar{\bx}$, have the direction of the eigenvectors, and their lengths are $2 \sqrt{\lambda_j}$, $j=1,2$.

![](03-pca_files/figure-latex/unnamed-chunk-4-1.pdf)<!-- --> 


We can now compute the PC scores using
\begin{eqnarray*}
y_{i1} &=& \bq_1^\top (\bx_i - \bar{\bx}) = 0.690 (x_{1i} - \bar{x}_1) + 0.724 (x_{2i} - \bar{x}_2) \\
y_{i2} &=& \bq_2^\top (\bx_i - \bar{\bx}) = -0.724 (x_{1i} - \bar{x}_1) + 0.690 (x_{2i} - \bar{x}_2),
\end{eqnarray*}
which gives


FIX FIX
\begin{tabular}{l|rrrrrrrrrr}
Student & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & 10 \\ \hline
$y_{[1]}^\top$ & 19.1 & 16.2 & 11.6 & -14.7 & -23.1 & -14.9 & 3.8 & 30.0 & -15.3 & -12.6 \\
$y_{[2]}^\top$ & -7.3 & -7.2 & 6.4 & -0.8 & 5.1 & -9.3 & 4.4 & 5.9 & 4.2 & -1.5
\end{tabular}

Note that these new variables have sample mean $\bar{\by}=\bzero$ and sample covariance matrix (see part 4. of Proposition \@ref(prp:pca1))
$$
\bLambda = \tdiag(\lambda_1,\lambda_2) =  \begin{pmatrix} 304.24 & 0 \\ 0 & 33.16 \end{pmatrix}.
$$
The plot below shows the PC scores $(y_{i1},y_{i2})^\top$.  The two lines shown have lengths $2\sqrt{\lambda_j}$, $j=1,2$.  Note that $\sqrt{\lambda_j}$ is the standard deviation of the $j$th PC.
\begin{center}
\includegraphics[width=12cm,angle=0]{figs/prbsta_pca2.pdf}
\end{center}



Sometimes the new variables have an obvious interpretation.  Note that the first PC gives positive, roughly equal, weight to PRB and STA and thus represents some form of ``average'' mark.  For example, a student that has a high mark on PRB and STA will have a high value for $y_1$.  The second PC, meanwhile, represents a contrast between PRB and STA.  For example, a large positive value for $y_2$ implies the student did much better on STA than PRB, and a large negative value implies the opposite.

Note that we could have chosen $t=-0.724$ instead of $t=+0.724$.  The only difference would be that the first eigenvector was $\bq_1^\ast = -\bq_1$.  In this case, a student who scored a high mark on PRB and STA would have a low value for $y_1$.  This is perfectly legitimate but makes the interpretation less intuitive.  One can always change the sign of the eigenvectors if it makes interpretation easier.
```



## Properties of principal components

Let  $\bx_1, \ldots, \bx_n$ have sample mean $\bar{\bx}$ and sample covariance matrix $\bS$, with spectral decomposition $\bS=\bQ\bLambda \bQ^\top$ where
 $\bQ=[\bq_1, \ldots , \bq_p]$ is orthogonal and $\bLambda=\text{diag}\{\lambda_1, \ldots , \lambda_p\}$. The transformed variables have some important properties.
\vskip 0.2truein

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:pca2"><strong>(\#prp:pca2) </strong></span>For $j,k=1, \ldots , p$, the following results hold.

1. $\bar{y}_{+j} = n^{-1} \sum_{i=1}^n y_{ij}=n^{-1}\sum_{i=1}^n \bq_j^\top (\bx_i-\bar{\bx})=0$;
2. $\bq_j^\top \bS\bq_j = \lambda_j$;
3. $\bq_j^\top \bS \bq_k = 0$ for $j \neq k$;
4. $\bq_1^\top \bS \bq_1 \geq \bq_2^\top \bS \bq_2 \geq \ldots \geq \bq_p^\top \bS \bq_p\geq 0$;
5. $\sum_{j=1}^p \bq_j^\top \bS \bq_j = \sum_{j=1}^p \lambda_j = \ttr(\bS)$;
6. $\prod_{j=1}^p \bq_j^\top \bS \bq_j = \prod_{j=1}^p \lambda_j = |\bS|$.

\EndKnitrBlock{proposition}

In words:

- part 1. tells us that the sample mean of $y_{1j}, \ldots , y_{nj}$ for each fixed $j$ is $0$;

- part 2. tells us that, for each fixed $j$, the sample variance of
the $y_{ij},\, i=1, \ldots , n$ is $\lambda_j$;

- part 3. states that the sample covariance of the pairs $(y_{ij}, y_{ik})$, $i=1, \ldots , n$, is $0$ if $j \neq k$;

- part 4. states that the sample variance of $y_{ij}, \, i=1, \ldots , n$, is not less than the sample variance of $y_{ik}, \, i=1, \ldots , n$, if $j\leq k$;

-  part 5. states that the sum of the sample variances is equal to the trace of $\bS$;

- and part 6. states that the product of the sample variances is equal to the determinant of $\bS$.


From these properties we say that a proportion
$$\frac{\lambda_j}{\lambda_1 + \ldots + \lambda_p}$$
of the variability in the sample is `explained' by the $j$th PC.

For the G11PRB and G11STA data above,
$$\frac{\lambda_1}{\lambda_1 + \lambda_2} = \frac{304.24}{304.24+33.16} = 0.90,$$
so 90\% of the variability in the sample is explained by the 1st PC.

```{Example} We can apply PCA to a football league table where $W$, $D$, $L$ are the number of matches won, drawn and lost and $F$ and $A$ are the goals scored for and against.  An extract of the table for a recent Premiership season is:
FIX FIX


\begin{tabular}{lrrrrr}
\toprule
Team & W & D & L & F & A\\
\midrule
Chelsea & 27 & 5 & 6 & 103 & 32\\
Manchester United & 27 & 4 & 7 & 86 & 28\\
Arsenal & 23 & 6 & 9 & 83 & 41\\
Tottenham Hotspur & 21 & 7 & 10 & 67 & 41\\
Manchester City & 18 & 13 & 7 & 73 & 45\\
\bottomrule
\end{tabular}




The sample mean vector is

$$\bar{\bx} =\begin{pmatrix}14.2 \\9.6 \\14.2 \\52.6 \\52.6 \\\end{pmatrix}$$

<!--$$\bar{\bx} = \begin{pmatrix} 14.20 \\ 9.60 \\ 14.20 \\ 52.65 \\ 52.65 \end{pmatrix},$$-->

and the sample covariance matrix is

\begin{equation}
\bS= \begin{pmatrix}39.4&-8.27&-31.1&116&-81.9 \\-8.27&8.14&0.13&-29.4&6.01 \\-31.1&0.13&31&-86.3&75.9 \\116&-29.4&-86.3&392&-209 \\-81.9&6.01&75.9&-209&231 \\\end{pmatrix}
(\#eq:PLES)
\end{equation}

<!--\begin{equation}
\bS = \begin{pmatrix} 39.4 & -8.3 & -31.1 & 115.7 & -81.9 \\
-8.3 & 8.1 & 0.1 & -29.4 & 6.0 \\
-31.1 & 0.1 & 31.0 & -86.3 & 75.9 \\
115.7 & -29.4 & -86.3 & 392.2 & -208.7 \\
-81.9 & 6.0 & 75.9 & -208.7 & 230.9 \end{pmatrix}.
(\#eq:PLES)
\end{equation}-->

The eigenvalues of $\bS$ are
$$\bLambda = \tdiag \begin{pmatrix}631&96.7&8.83&2.44&-4.97e-14 \\\end{pmatrix}$$

<!--$$\bLambda = \tdiag(599.06, 91.85, 8.39, 2.32, 0.00).$$-->

Note that we have a zero eigenvalue because one of our variables is a linear combination of the other variables, $L = 38 - W - D$.  The corresponding eigenvectors are
$$\bQ = [\bq_1 \ldots \bq_5] =\begin{pmatrix}0.251&-0.0133&-0.116&0.768&0.577 \\-0.0477&-0.146&0.74&-0.309&0.577 \\-0.204&0.16&-0.624&-0.459&0.577 \\0.776&0.582&0.0674&-0.234&-2e-15 \\-0.539&0.784&0.213&0.222&1.83e-15 \\\end{pmatrix}$$


<!--$$\bQ = [\bq_1 \ldots \bq_5] =
\begin{pmatrix} 0.25 & 0.01 & -0.12 & 0.77 & 0.58 \\
-0.05 & 0.15 & 0.74 & -0.31 & 0.58 \\
-0.20 & -0.16 & -0.62 & -0.46 & 0.58 \\
0.78 & -0.58 & 0.07 & -0.23 & 0.00 \\
-0.54 & -0.78 & 0.21 & 0.22 & 0.00 \end{pmatrix}.$$-->

The proportion of variability explained by each of the PCs is:
$$
\begin{pmatrix}0.854&0.131&0.012&0.0033&-6.73e-17 \\\end{pmatrix}
$$

<!--$$
%(0.854, 0.131, 0.012, 0.003, 0.000).
$$-->

There is no point computing the scores for PC 5 because PC5 does not explain any of the variability in the data.  Similarly, there is little value in computing the scores for PCs 3 \& 4 because they only account for 1.5\% of the variability in the data.

We can, therefore, choose to compute only the first two PC scores.  We are reducing the dimension of our data set from $p=5$ to $p=2$ while still retaining 98.5\% of the variability.  The first PC is given by:
\begin{align*}
y_{i1} &= 0.25(W_i-\bar{W}) +-0.05(D_i-\bar{D}) +-0.2(L_i-\bar{L})\\
& \qquad +0.78(F_i-\bar{F}) +-0.54(A_i-\bar{A}),
\end{align*}
and similarly for PC 2.

<!--\begin{align*}
y_{i1} &= 0.25(W_i-\bar{W}) -0.05(D_i-\bar{D}) -0.20(L_i-\bar{L})\\
& \qquad +0.78(F_i-\bar{F}) -0.54(A_i-\bar{A}),
\end{align*}-->

The first five rows of our revised ``league table'' are now
\begin{table}[H]
\centering
\begin{tabular}{lrr}
\toprule
Team & PC1 & PC2\\
\midrule
Chelsea & 55.3 & 12.3\\
Manchester United & 44.1 & -0.4\\
Arsenal & 33.3 & 8.1\\
Tottenham Hotspur & 20.1 & -1.2\\
Manchester City & 22.2 & 4.1\\
\bottomrule
\end{tabular}
\end{table}



<!--FIX FIX
\begin{tabular}{l|rr}
Team &		PC 1 & PC 2 \\ \hline
Chelsea &	55.3 & -12.3 \\
Man Utd &	44.1 & 0.4 \\
$\vdots$ && $\vdots$ \\
Portsmouth &	-25.4 & -1.7
\end{tabular}-->



Now that we have reduced the dimension to $p=2$, we can visualise the differences between the teams.

![](03-pca_files/figure-latex/unnamed-chunk-8-1.pdf)<!-- --> 


<!--\begin{center}
\includegraphics[width=12cm,angle=0 ]{figs/premleag_pca1.pdf}
\end{center}-->

We might interpret the PCs as follows.  The first PC seems to measure overall performance.  It rewards teams with 0.78 for every goal they score and 0.25 for every match they win, while penalising them by 0.54 for every goal they concede, 0.2 for every match they lose and 0.05 for every match they draw.

We could, therefore, rank teams by PC 1 and compare this with the rankings using 3 points for a win and 1 point for a draw.  The rankings are the same for the top three teams but differ below that.  Under our system Wigan would be relegated in place of Portsmouth.

The second PC has a strong negative loading for both goals for and against.  A team with a large negative PC 2 score was, therefore, involved in matches with lots of goals.  We could, therefore, interpret PC 2 as an ``entertainment'' measure, ranking teams according to their involvement in high-scoring games.

The above example raises the question of how many PCs should we use in practice.  If we reduce the dimension to $p=1$ then we can rank observations and analyse our new variable with univariate statistics.  If we reduce the dimension to $p=2$ then it is still easy to visualise the data.  However, reducing the dimension to $p=1$ or $p=2$ may involve losing lots of information and a sensible answer should depend on the objectives of the analysis and the data itself.

One tool for looking at the contributions of each PC is to look at the **scree graph** which plots the percentage of variance explained by PC $j$ against $j$.  The scree graph for the football example is:

![](03-pca_files/figure-latex/unnamed-chunk-9-1.pdf)<!-- --> 


<!--\begin{center}
\includegraphics[width=12cm,angle=0]{figs/premleag_pca2.pdf}
\end{center}-->

Possible methods for choosing the number of PCs include:


- retain enough PCs to explain, say, 90\% of the total variation;
- retain PCs where the eigenvalue is above the average.

For the football example, the first method would retain 2 PCs whereas the second method would only retain 1 PC.

```


## Population PCA

So far we have considered sample PCA based on the sample covariance matrix
$$
\bS=\frac{1}{n}\sum_{i=1}^n (\bx_i-\bar{\bx})(\bx_i-\bar{\bx})^\top.
$$

We note now that there is a *population* analogue of PCA based on the population
covariance matrix $\bSigma$.  Although the population version of PCA is not of as much direct practical
relevance as sample PCA, it is nevertheless of conceptual importance.

Let $\bx$ denote a $p \times 1$ random vector with $E(\bx)={\pmb \mu}$ and $\text{Var}(\bx)={\pmb \Sigma}$.  As defined,
$\pmb \mu$ is the population mean vector and $\pmb \Sigma$ is the population covariance matrix.

Since $\pmb \Sigma$ is symmetric, the spectral decomposition theorem tells us that
$$
{\pmb \Sigma}=\sum_{j=1}^p \check{\lambda}_j \check{\bq}_j \check{\bq}_j^\top=\check{\bQ} \check{\bLambda}\check{\bQ}^\top
$$
where the `check' symbol  $\quad \check{} \quad$   is used to distinguish population quantities from their sample analogues.

Then:

- the first population PC is defined by $Y_1=\check{\bq}_1^\top (\bx-{\pmb \mu})$;
-the second population PC is defined by $Y_2=\check{\bq}_2^\top (\bx-{\pmb \mu})$;
-  \$ldots$
- the $p$th population PC is defined by $Y_p=\check{\bq}_p^\top (\bx-{\pmb \mu})$.


The $Y_1, \ldots , Y_p$ are random variables, unlike the sample PCA case, where the $y_{ij}$ are observed quantities.
In the sample PCA case, the $y_{ij}$ can often be regarded as the observed values of random variables.

In matrix form, the above definitions can be summarised by writing
$$
\by =\begin{pmatrix} Y_1 \\ Y_2 \\ ... \\...\\Y_p   \end{pmatrix} = \check{\bQ}^\top (\bx-{\pmb \mu}).
$$

The population PCA analogues of the 6 sample PCA properties listed in Proposition \@ref(prp:pca2)  are now given.  Note that the
$Y_j$'s are random variables as opposed to observed values of random variables.


\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:pca3"><strong>(\#prp:pca3) </strong></span>The following results hold for the random variables $Y_1, \ldots , Y_p$ defined above.

1.  $E(Y_j)=0$ for $j=1, \ldots , p$;\\

2.   $\text{Var}(Y_j)=\check{\lambda}_j$ for $j=1,\ldots, p$;\\

3.  $\text{Cov}(Y_j,Y_k)=0$ if $j \neq k$;\\

4.  $\text{Var}(Y_1) \geq \text{Var}(Y_2) \geq \cdots \geq \text{Var}(Y_p) \geq 0$;\\

5.  $\sum_{j=1}^p \text{Var}(Y_j)=\sum_{j=1}^p \check{\lambda}_j=\text{tr}(\bSigma)$;\\

6.  $\prod_{j=1}^p \text{Var}(Y_j)=\prod_{j=1}^p \check{\lambda}_j=\vert \bSigma \vert$.
 
\EndKnitrBlock{proposition}

Note that, defining $\by=(Y_1, \ldots , Y_p)^\top$ as before,  part 1. implies that $E(\by)={\mathbf 0}_p$ and parts 2. and 3. together imply that
$$
\text{Var}(\by)=\bLambda \equiv \text{diag}(\check{\lambda}_1, \ldots , \check{\lambda}_p).
$$

\BeginKnitrBlock{example}
<span class="example" id="exm:unnamed-chunk-10"><strong>(\#exm:unnamed-chunk-10) </strong></span>Suppose
$$
\bSigma=\bI_p+ \delta {\mathbf 1}_p {\mathbf 1}_p^\top,
$$
where $\delta>0$.  What is the proportion of variability explained by the first PC?  Since $\delta>0$, the largest eigenvalue is $\lambda_1=1+p\delta$ which is achieved when $\check{\bq}_1$ is the unit vector $p^{-1/2}{\mathbf 1 }_p$.  This and related examples are dealt with in more detail in the example sheets.
\EndKnitrBlock{example}

Consider now a repeated sampling framework in which we assume that $\bx_1, \ldots , \bx_n$ are IID random vectors from a population
with mean vector $\pmb \mu$ and covariance matrix $\bSigma$.

What is the relationship between the sample PCA based on the sample of observed vectors $\bx_1, \ldots , \bx_n$, and the population PCA based on the unobserved random vector $\bx$,
from the same population?

Assuming $n$ is large, and the elements of $\bSigma$ are all finite, the elements of the sample covariance matrix $\bS$ will be close with high probability to the corresponding elements
of the population covariance matrix $\bSigma$.  Justification of this statement comes from the weak law of large numbers applied to the components of $\Sigma$ (details omitted).

Consequently, when $n$ is large, sample PCA and the corresponding population PCA may be expected to give similar results.

## An Alternative Derivation of PCA

Consider a sample $\bx_1, \ldots ,  \bx_n \in \mathbb{R}^p$.

Recall from \S2.8 that any line in $\mathbb{R}^p$ may be written in the form
$\{\ba +u \bb: \, u \in \mathbb{R}\}$ where $\ba , \bb \in \mathbb{R}^p$ are
fixed.

Here we consider the following problem: *find the best-fitting line to the sample*
$\bx_1, \ldots , \bx_n$.

We first formulate this problem more precisely.  Define the function
\begin{align*}
F(\ba , \bb; u_1, \ldots , u_n)&=\sum_{i=1}^n \vert \vert \bx_i - \ba - u_i \bb \vert \vert^2\\
& =\sum_{i=1}^n (\bx_i - \ba -u_i \bb)^\top (\bx_i - \ba -u_i \bb).
\end{align*}

We wish to solve the following problem:
\begin{align}
&\textit{minimise $F(\ba, \bb; u_1, \ldots , u_n)$ subject to the}\nonumber \\
&\textit{constraints that $\ba$ and $\bb$ are orthogonal, i.e. $\ba^\top \bb=0$,}(\#eq:linemanager)\\
&\textit{and $\bb$ is a unit vector, i.e. $\vert \vert \bb \vert \vert =1$.} \nonumber
\end{align}

\BeginKnitrBlock{theorem}
<span class="theorem" id="thm:unnamed-chunk-11"><strong>(\#thm:unnamed-chunk-11) </strong></span>The solution to optimisation problem \@ref(eq:linemanager)
is given by
\begin{equation}
\hat{\ba}=\left ({\mathbf I}_p-\bq_1 \bq_1^\top\right)\bar{\bx}, \quad \hat{\bb}= \bq_1 \quad
\hbox{and} \quad \hat{u}_i=\bx_i^\top \bq_1, \quad i=1,\ldots , n,
(\#eq:hatitems)
\end{equation}
where $\bar{\bx}=n^{-1}\sum_{i=1}^n \bx_i$ is the sample mean and the unit vector $\bq_1$ is the direction of the first sample PC.
\EndKnitrBlock{theorem}

Note that  the quantities $u_i -\bar{u}=\bq_1^\top (\bx_i-\bar{\bx})$ are the PC scores associated with the first PC.

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}The proof is broken into two steps.

**Step 1**.   In Step 1, we  want to minimise $F(\ba, \bb; u_1,\ldots, u_n)$ subject to the constraint $\ba^\top \bb=0$, with $\bb$ an arbitrary *fixed* unit vector in $\mathbb{R}^p$.  So we introduce a Lagrangian term for the constraint $\ba^\top \bb=0$ and minimise
$$
\bar{F}(\ba; u_1,\ldots ,  u_n; \gamma)
 \equiv \left \{\sum_{i=1}^n (\bx_i -\ba-u_i\bb)^\top
(\bx_i-\ba-u_i\bb)\right \}+ \gamma \ba^\top \bb
$$
over $\ba$, $u_1, \ldots , u_n$ and $\gamma$.  Then, for $i=1, \ldots , n$,
\begin{equation}
\frac{\partial \bar{F}}{\partial u_i}=-2\bb^\top (\bx_i-\ba -u_i\bb);
(\#eq:barF1)
\end{equation}
\begin{align}
\frac{\partial \bar{F}}{\partial \ba}&=-2\left \{\sum_{i=1}^n (\bx_i -\ba -u_i\bb)\right \}+\gamma \bb\nonumber\\
&=-2n\{\bar{\bx}-\ba-(\bar{u}+\gamma/(2n))\bb\},
(\#eq:barF2)
\end{align}
where $\bar{u}=n^{-1}\sum_{i=1}^n u_i$; and
\begin{equation}
\frac{\partial \bar{F}}{\partial \gamma}=\ba^\top \bb.
(\#eq:barF3)
\end{equation}
Setting the partial derivatives \@ref(eq:barF1), \@ref(eq:barF2) and \@ref(eq:barF3}) to zero,
$$
\frac{\partial \bar{F}}{\partial \gamma} =0 \implies \hat{\ba}^\top \bb=0;
$$
\begin{equation}
\frac{\partial \bar{F}}{\partial u_i}=0 \implies \hat{u}_i=\bb^\top \bx_i,
(\#eq:uhat)
\end{equation}
and therefore
\begin{equation}
\bar{\hat{u}}\equiv n^{-1} \sum_{i=1}^n \hat{u}_i=\bb^\top \bar{\bx};
(\#eq:barlambda)
\end{equation}
and
$$
\frac{\partial \bar{F}}{\partial \ba}={\mathbf 0}_p \implies \hat{\ba}=\bar{\bx}-\{\bar{\hat{u}}+\hat{\gamma}/(2n)\}\bb.
$$
Using \@ref(eq:barlambda) and the fact that $\bb^\top \hat{\ba}=0$, it follows that
$$
0=\bb^\top \hat{\ba} =\bb^\top [\bar{\bx}-\{\hat{\bar{u}}+\hat{\gamma}/(2n)\}\bb]=\bb^\top \bar{\bx}-\bb^\top \bar{\bx}
+\hat{\gamma}/(2n),
$$
which implies that $\hat{\gamma}=0$.  Consequently,
\begin{equation}
\hat{\ba}=\bar{\bx}-\bar{\hat{u}}\bb=\bar{\bx}-\bb \bb^\top \bar{\bx}=\left ( \bI_p-\bb \bb^\top \right ) \bar{\bx};
(\#eq:bahat)
\end{equation}
and so
\begin{align}
&\bar{F}(\hat{\ba}; \hat{u}_1, \ldots, \hat{u}_n; \hat{\gamma})\\
&=\sum_{i=1}^n (\bx_i-\hat{\ba}-\hat{u_i} \bb)^\top (\bx_i -\hat{\ba}-\hat{u}_i\bb)\nonumber \\
&=\sum_{i=1}^n \left \{\bx_i-({\mathbf I}_p-\bb \bb^\top)\bar{\bx}- \bb \bb^\top \bx_i \right \}^\top \left \{\bx_i -({\mathbf I}_p-\bb \bb^\top)\bar{\bx} -\bb \bb^\top \bx_i\right \}\nonumber \\
&\sum_{i=1}^n (\bx_i -\bar{\bx})^\top ({\mathbf I}_p -\bb \bb^\top)^2 (\bx_i -\bar{\bx})\nonumber \\
&= n \text{tr}\left \{ ({\mathbf I}_p-\bb \bb^\top )\bS \right \} \nonumber \\
&=n\left \{\text{tr}(\bS)-\bb^\top \bS \bb \right\},
(\#eq:object37)
\end{align}
where $\bS$ is the sample covariance of the $\bx_i$.

**Step 2.**  We now minimise \@ref(eq:object37) over unit vectors $\bb \in \mathbb{R}^p$.  But minimising \@ref(eq:object37) is equivalent to maximising $\bb^\top \bS \bb$, so
from Proposition \@ref(prp:two8), $\hat{\bb}=  \bq_1$, and so $\hat{\ba}=(\bI_p-\bq_1 \bq_1^\top)\bar{\bx}$ from \@ref(eq:bahat), and
from \@ref(eq:uhat), $\hat{u}_i=\bq_1^\top \bx_i$, $i=1,\ldots , n$, all of which agrees with the expressions in \@ref(eq:hatitems).  
\EndKnitrBlock{proof}



## PCA under transformations of variables

Let us return to the example of $n=10$ students who studied G11PRB and G11STA.  Earlier, we calculated the sample mean, sample variance matrix and the eigenvalues/vectors of $\bS$,
\begin{eqnarray*}
\bar{\bx} = \begin{pmatrix} 62.6 \\ 66.2 \end{pmatrix}, &\quad&
\bS = \begin{pmatrix} 162.04 & 135.38 \\ 135.38 & 175.36 \end{pmatrix} \\
\bLambda = \begin{pmatrix} 304.24 & 0 \\ 0 & 33.16 \end{pmatrix}, &\quad&
\bQ = \begin{pmatrix} 0.690 & -0.724 \\ 0.724 & 0.690 \end{pmatrix}
\end{eqnarray*}
with PC 1 scores
$$y_i = \bq_1^\top (\bx_i - \bar{\bx}) = 0.690 (x_{1i} - \bar{x}_1) + 0.724 (x_{2i} - \bar{x}_2).$$

We now consider what happens to the above quantities under various transformations of the $\bx_i$, the $2 \times 1$
response vectors.

**Addition transformation**

Firstly, we consider the transformation of addition where, for example, the G11PRB lecturer decides to add 5 marks for all the students.  We can write this transformation as $\bz_i = \bx_i + \bc$, where $\bc$ is a fixed vector.  Under this transformation the sample mean changes, $\bar{\bz} = \bar{\bx} + \bc$, but the sample variance remains $\bS$.  Consequently, the eigenvalues and eigenvectors of $\bS$ remain the same and, therefore, so does the PC 1 score,
$$y_i = \bq_1^\top (\bz_i - \bar{\bz}) = \bq_1^\top(\bx_i + \bc - (\bar{\bx} + \bc)) = \bq_1^\top (\bx_i - \bar{\bx}).$$
We say that the principal components are **invariant** under the addition transformation.  An important special case is to choose $\bc = -\bar{\bx}$ so that the PC 1 score is simply $y_i = \bq_1^\top \bz_i$.

**Scale transformation**

Secondly, we consider the scale transformation where, for example, the G11PRB lecturer decides to double the marks for all students.  A scale transformation occurs more naturally when we convert units of measurement from, say, metres to kilometres.  We can write this transformation as $\bz_i = \bD \bx_i$, where $\bD$ is a diagonal matrix with positive elements.  Under this transformation the sample mean changes from $\bar{\bx}$ to $\bar{\bz} = \bD \bar{\bx}$, and the sample covariance matrix changes from $\bS$ to $\bD \bS \bD$.  Consequently, the principal components also change.

This lack of scale-invariance is undesirable.  One solution is to choose
$$
\bD = \tdiag(s_{11}^{-1/2}, \ldots , s_{pp}^{-1/2}),
 $$
 where $s_{ii}$ is the $i$th diagonal element of $\bS$.  In effect, we have standardised all the new variables to have variance 1.  In this case the sample covariance matrix of the $\bz_i$'s is simply the sample correlation matrix of the original variables, $\bx_i$.  Therefore, we can carry out PCA on the sample correlation matrix, $\bR$, which is invariant to changes of scale.

In summary: $\bR$ is scale-invariant while $\bS$ is not.

\BeginKnitrBlock{example}
<span class="example" id="exm:unnamed-chunk-13"><strong>(\#exm:unnamed-chunk-13) </strong></span>For the G11PRB/G11STA data, we choose
$$
\bD = \tdiag(162.04,175.36)^{-1/2} = \tdiag(0.079,0.076)
$$
so that $\bz_i = \bD \bx_i$.
 The sample correlation matrix is then
\begin{align*}
\bR &= \bD \bS \bD\\
 &= \begin{pmatrix} 0.079 & 0 \\ 0 & 0.076 \end{pmatrix}
\begin{pmatrix} 162.04 & 135.38 \\ 135.38 & 175.36 \end{pmatrix}
\begin{pmatrix} 0.079 & 0 \\ 0 & 0.076 \end{pmatrix} \\
&= \begin{pmatrix} 1.000 & 0.803 \\ 0.803 & 1.000 \end{pmatrix}.
\end{align*}
The eigenvalues and eigenvectors of $\bR$ are then
$$\bLambda = \begin{pmatrix} 1.803 & 0 \\ 0 & 0.197 \end{pmatrix}, \qquad
\bQ = \begin{pmatrix} 0.707 & 0.707 \\ 0.707 & -0.707 \end{pmatrix},$$
and the PC 1 score is
\begin{eqnarray*}
y_i &=& \bq_1^\top (\bz_i - \bar{\bz}) = \bq_1^\top \bD (\bx_i - \bar{\bx}) \\
&=& 0.707 \times 0.079 (x_{1i} - \bar{x}_1) + 0.707 \times 0.076 (x_{2i} - \bar{x}_2).
\end{eqnarray*}

In the example above, there is little difference between using $\bS$ and $\bR$ for the PCA because the variances for G11PRB and G11STA are similar.  In other cases, particularly when the variables are measured on wildly different scales, the difference will be notable.  For example, in the football data the sample variances of $F$ and $A$ are much larger than the sample variances of $W$, $D$ and $L$.
\EndKnitrBlock{example}

**Orthogonal transformation**

Thirdly, we consider a transformation by an orthogonal matrix, $\stackrel{p \times p}{\bA}$, such that $\bA \bA^\top = \bA^\top \bA = \bI_p$, and write $\bz_i = \bA \bx_i$.  This is equivalent to rotating and/or reflecting the original data.

Let $\bS$ be the sample covariance matrix of the $\bx_i$ and let $\bT$ be the sample covariance matrix of the $\bz_i$.  Under this transformation the sample mean changes from $\bar{\bx}$ to $\bar{\bz} = \bA \bar{\bx}$, and the sample covariance matrix $\bS$ changes from $\bS$ to $\bT = \bA \bS \bA^\top$.

However, if we write $\bS$ in terms of its spectral decomposition $\bS = \bQ \bLambda \bQ^\top$, then $\bT = \bA \bQ \bLambda \bQ^\top \bA^\top = \bB \bLambda \bB^\top$ where $\bB = \bA \bQ$ is also orthogonal.  It is therefore apparent that the eigenvalues of $\bT$ are the same as those of  $\bS$; and the eigenvectors of $\bT$ are given by $\bb_j$ where $\bb_j = \bA \bq_j$, $j=1,\ldots,p$.  The PC 1 scores of the transformed variables are
$$ y_i = \bb_1^\top (\bz_i - \bar{\bz}) = \bq_1^\top \bA^\top \bA (\bx_i - \bar{\bx}) = \bq_1^\top (\bx_i - \bar{\bx}),$$
and so they are identical to the PC 1 scores of the original variables.

Therefore, under an orthogonal transformation the eigenvalues and PC scores are unchanged and the PCs are orthogonal transformations of the original PCs.  We say that the principal components are **equivariant** with respect to orthogonal transformations.


\BeginKnitrBlock{example}
<span class="example" id="exm:unnamed-chunk-14"><strong>(\#exm:unnamed-chunk-14) </strong></span>Suppose we rotate the G11PRB/G11STA data by the matrix $\bA = \begin{pmatrix} 0.866 & -0.500 \\ 0.500 & 0.866 \end{pmatrix}$.  The sample covariance matrix of the rotated data is
\begin{align*}
\bT &= \bA \bS \bA^\top\\
&= \begin{pmatrix} 0.866 & -0.500 \\ 0.500 & 0.866 \end{pmatrix}
\begin{pmatrix} 162.04 & 135.38 \\ 135.38 & 175.36 \end{pmatrix}
\begin{pmatrix} 0.866 & 0.500 \\ -0.500 & 0.866 \end{pmatrix} \\
&= \begin{pmatrix} 48.13 & 61.92 \\ 61.92 & 289.27 \end{pmatrix}.
\end{align*}
The eigenvalues of $\bT$ are $304.24$ and $33.16$ (same as for $\bS$).  The eigenvectors of $\bT$ are then
\begin{eqnarray*}
\bB &=& \bA \bQ = \begin{pmatrix} 0.866 & -0.500 \\ 0.500 & 0.866 \end{pmatrix} \begin{pmatrix} 0.690 & -0.724 \\ 0.724 & 0.690 \end{pmatrix} \\
&=& \begin{pmatrix} 0.235 & -0.972 \\ 0.972 & 0.235 \end{pmatrix}
\end{eqnarray*}
and the PC 1 scores are unchanged.
\EndKnitrBlock{example}


## PCA based on $\bS$ versus PCA based on $\bR$ 

Recall the distinction between the sample covariance matrix $\bS$ and the sample correlation matrix $\bR$.

Note that all correlation matrices are also covariance matrices, but not all covariance matrices are correlation matrices.

So in practice we have a choice of using $\bS$ or $\bR$ for PCA.  As we have seen, PCA based on $\bR$ is scale invariant, but PCA based on $\bS$ is not; while PCA based on $\bS$ is invariant (eigenvalues and PC scores) and equivariant (eigenvectors) under orthogonal transformation, whereas $\bR$ is not.

This raises the important practical question: for a given dataset, should we use PCA based on $\bS$ or $\bR$?

If the $p$ variables represent very different types of quantity or show marked differences in variances, then it will usually be better to use $\bR$ rather than $\bS$.  However, in some circumstances, we may wish to use $\bS$, such as when the $p$ variables are measuring similar entities and the sample variances are not too different.

Bearing in mind that the required numerical calculations are so easy to perform in R, we might wish to do it both ways and see if it makes much difference.


