
# PART I:  Prerequisites {-}


In Chapter \@ref(stat-prelim) we explain what we mean by multivariate analysis and give some examples of multivariate data.  We also introduce basic definitions and concepts such as the sample covariance matrix, the sample correlation matrix and graphical techniques. We also briefly discuss random vectors and random matrices and derive some of their elementary properties.



In Chapter \@ref(linalg-prelim) we summarise the definitions, ideas and results from matrix algebra that will be needed later in the module.

\renewcommand{\bY}{\boldsymbol Y}
\renewcommand{\bx}{\boldsymbol x}
\renewcommand{\bX}{\boldsymbol X}
\renewcommand{\bH}{\boldsymbol H}
\renewcommand{\by}{\boldsymbol y}
\renewcommand{\bz}{\boldsymbol z}
\renewcommand{\bS}{\boldsymbol S}
\renewcommand{\bR}{\boldsymbol R}
\renewcommand{\bI}{\boldsymbol I}
\renewcommand{\bmu}{\boldsymbol \mu}
\renewcommand{\bSigma}{\boldsymbol \Sigma}
\renewcommand{\bLambda}{\boldsymbol \Lambda}
\renewcommand{\bgamma}{\boldsymbol \gamma}
\renewcommand{\blambda}{\boldsymbol \lambda}
\renewcommand{\bA}{\boldsymbol A}
\renewcommand{\bB}{\boldsymbol B}
\renewcommand{\bD}{\boldsymbol D}
\renewcommand{\bM}{\boldsymbol M}
\renewcommand{\bP}{\boldsymbol P}
\renewcommand{\bQ}{\boldsymbol Q}
\renewcommand{\bT}{\boldsymbol T}
\renewcommand{\bW}{\boldsymbol W}
\renewcommand{\ba}{\boldsymbol a}
\renewcommand{\bb}{\boldsymbol b}
\renewcommand{\bc}{\boldsymbol c}
\renewcommand{\bd}{\boldsymbol d}
\renewcommand{\bh}{\boldsymbol h}
\renewcommand{\bp}{\boldsymbol p}
\renewcommand{\bq}{\boldsymbol q}
\renewcommand{\bu}{\boldsymbol u}
\renewcommand{\bzero}{\boldsymbol 0}
\renewcommand{\mR}{\mathbb R}
\renewcommand{\cR}{\mathcal R}

\renewcommand{\bs}{\boldsymbol}
\renewcommand{\ds}{\displaystyle}
\renewcommand{\tdiag}{\text{diag}}
\renewcommand{\ttr}{\text{tr}}
\renewcommand{\tmin}{\text{min}}
\renewcommand{\tmax}{\text{max}}
\renewcommand{\tdet}{\text{det}}

\renewcommand{\tcov}{\text{cov}}
\renewcommand{\texp}{\text{exp}}
\renewcommand{\lb}{\left(}
\renewcommand{\rb}{\right)}
\renewcommand{\lsb}{\left[}
\renewcommand{\rsb}{\right]}

<!-- NEED TO OPEN PROJECT TO MAKE IT WORK PROPERLY! -->




#  Statistical Preliminaries {#stat-prelim}

In this chapter blah blah


## Multivariate data

What is multivariate analysis (MVA)?  Analysis of data where two or more response variables are measured on each object under study.  If we measure $p$ variables on $n$ objects then the data can be presented in a $n \times p$ **data matrix**.

We shall often write the data matrix as $\mathbf X$ ($n \times p$) where
$$
{\mathbf X}=\left[ \begin{array}{c}
\bx_1^\top\\
\bx_2^\top\\
..\\
..\\
..\\
\bx_n^\top
\end{array}\right ]= [ \bx_1, \ldots , \bx_n]^\top.
$$

In words: the rows of $\mathbf X$ are $\bx_1^\top, \ldots , \bx_n^\top$ and the columns of $\bX^\top$ are
$\bx_1, \ldots , \bx_n$.

In this setup, we think of the $\bx_1, \ldots , \bx_n$ are being the observation vectors, and the $p$ columns of $\mathbf X$
correspond to the $p$ variables being measured.


**Important remark on notation:**  Throughout the module we shall use non-bold letters, whether upper or lower case, to indicate scalar (i.e. real-valued) quantities; lower-case letters in bold to signify column vectors; and upper case letters in bold to signify matrices.  This convention for bold letters will also apply to random quantities.  So, in particular, for a random vector we always use (bold) lower case, and for a random matrix we always use bold upper-case, regardless of whether we are referring to (i) the unobserved random quantity or (ii) its observed value.  It should always be clear from the context which of these two interpretations (i) or (ii)  is appropriate.


\BeginKnitrBlock{example}
<span class="example" id="exm:unnamed-chunk-1"><strong>(\#exm:unnamed-chunk-1) </strong></span>Football league table where $W =$ number of wins, $D =$ number of draws, $F =$ number of goals scored and $A =$ number of goals conceded for four teams. 
In this example we have $p=4$ variables $(W, D, F, A)^\top$ measured on $n=4$ objects (teams).
\EndKnitrBlock{example}

\begin{table}[H]
\centering
\begin{tabular}{lrrrr}
\toprule
Team & W & D & F & A\\
\midrule
USA & 1 & 2 & 4 & 3\\
England & 1 & 2 & 2 & 1\\
Slovenia & 1 & 1 & 3 & 3\\
Algeria & 0 & 1 & 0 & 2\\
\bottomrule
\end{tabular}
\end{table}


\BeginKnitrBlock{example}
<span class="example" id="exm:unnamed-chunk-3"><strong>(\#exm:unnamed-chunk-3) </strong></span>Exam marks for a set of $n$ students where $P =$ mark in probability and $S =$ mark in statistics. 
Note that $x_{ij}$ denotes the $j$th variable measured on the $i$th subject.
\EndKnitrBlock{example}

\begin{table}[H]
\centering
\begin{tabular}{lll}
\toprule
Student & P & S\\
\midrule
1 & $x_{11}$ & $x_{12}$\\
2 & $x_{21}$ & $x_{22}$\\
$\vdots$ & $\vdots$ & $\vdots$\\
n & $x_{n1}$ & $x_{n2}$\\
\bottomrule
\end{tabular}
\end{table}





In MVA we attempt to answer questions such as:

- What is the joint distribution of marks?
- How can we visualise the data?
- Can we simplify the data?  For example, we rank football teams using $3W+D$ and we rank students by their average module mark but is this fair?  Can we reduce the dimension in a better way?
- Can we use the data to discriminate, for example, between male and female students?

We could just apply standard univariate techniques to each variable in turn but this ignores possible dependencies between the variables which we must represent to draw valid conclusions.

Finally, before moving on, we ask the question: what is the difference between MVA and standard linear regression? Answer: in standard linear regression we have a scalar response variable, $y$ say, and a vector of covariates, $\bx$, say.  The focus of interest is on how knowledge of $\bx$ influences the distribution of $y$ (in particular, the mean of $y$).  In contrast, with MVA the focus of interest is a response vector $\by$, in which all the components of $\by$ are viewed as responses rather than covariates.  However,  there are also situations where the response is a vector $\by$ but we also have covariate information   $\bx$.  This leads to study of the multivariate linear model, which we will investigate later on in Chapter ???.




## Summary statistics

In univariate statistics we define the sample mean and sample variance as
$$ \bar{x} = \frac{1}{n} \sum_{i=1}^n x_i \quad \text{and} \quad s^2 = \frac{1}{n} \sum_{i=1}^n (x_i - \bar{x})^2. $$

The analogous multivariate **sample mean** and **sample covariance matrix** are direct extensions:
$$ \bar{\bx} = \frac{1}{n} \sum_{i=1}^n \bx_i \quad \text{and} \quad \bS = \frac{1}{n} \sum_{i=1}^n (\bx_i - \bar{\bx}) (\bx_i - \bar{\bx})^\top, $$
where $\bx_i$ is the $p$ dimensional vector denoting the $p$ observations on the $i$th object.

Note that the $j$th entry in $\bar{\bx}$ is simply the (univariate) sample mean of the $j$th variable.  Similarly, if $s_{ab}$ is the $(a,b)$th entry of $\bS$, then $s_{jj}$ is the (univariate) sample variance of the $j$th variable and $s_{ab}$ is the sample covariance of variables $a$ and $b$.  Note that $\bS$ is symmetric since $s_{ab}=s_{ba}$.

An equivalent alternative formula for $\bS$ is
$$\bS = \frac{1}{n} \lb \sum_{i=1}^n \bx_i \bx_i^\top \rb - \bar{\bx} \bar{\bx}^\top.$$

Similarly, let $\bR$ be the **sample correlation matrix** where the $(a,b)$th entry of $\bR$ is the sample correlation between variables $a$ and $b$, that is
$$ r_{ab} = s_{ab}/\sqrt{s_{aa}s_{bb}}. $$
Note that
$$ \bR = \bD^{-1} \bS \bD^{-1} $$
where $\bD = \tdiag(\sqrt{s_{11}}, \dots, \sqrt{s_{pp}})$.  Note that $\bR$ is symmetric, the diagonal entries are always exactly 1 (each variable is perfectly correlated with itself) and that $|r_{ab}| \leq 1$.

Note that if we change the unit of measurement for the $\bx_i$'s then $\bS$ will change but $\bR$ will not.

The **total variation** in the data set is usually measured by $\ttr(\bS)$ where $\ttr()$ is the trace function that sums the diagonal elements of the matrix.  That is,
$$\ttr(\bS) = s_{11} + s_{22} + \ldots + s_{pp}.$$

In MVA it is much easier to work with vector and matrix notation.

\BeginKnitrBlock{example}
<span class="example" id="exm:unnamed-chunk-5"><strong>(\#exm:unnamed-chunk-5) </strong></span>The table below  shows the module marks for 5 students on the modules G11PRB ($P$) and G11STA ($S$).
\EndKnitrBlock{example}




\begin{table}[H]
\centering
\begin{tabular}{lrr}
\toprule
Student & P & S\\
\midrule
A & 41 & 63\\
B & 72 & 82\\
C & 46 & 38\\
D & 77 & 57\\
E & 59 & 85\\
\bottomrule
\end{tabular}
\end{table}


Calculate the sample mean, sample covariance, sample correlation and total variation.

The sample mean is $\bar{\bx} = \begin{pmatrix} 59 \\ 65 \end{pmatrix}$.

The sample covariance matrix is $\bS = \begin{pmatrix} 197.2 & 92.8 \\ 92.8 & 297.2 \end{pmatrix}$.

The sample correlation matrix is
\begin{align*}
\bR &= \bD^{-1} \bS \bD^{-1}\\
&= \begin{pmatrix} 14.0 & 0 \\ 0 & 17.2 \end{pmatrix}^{-1} \begin{pmatrix} 197.2 & 92.8 \\ 92.8 & 297.2 \end{pmatrix} \begin{pmatrix} 14.0 & 0 \\ 0 & 17.2 \end{pmatrix}^{-1} \\
&= \begin{pmatrix} 0.071 & 0 \\ 0 & 0.058 \end{pmatrix} \begin{pmatrix} 197.2 & 92.8 \\ 92.8 & 297.2 \end{pmatrix} \begin{pmatrix} 0.071 & 0 \\ 0 & 0.058 \end{pmatrix} \\
&= \begin{pmatrix} 1.000 & 0.383 \\ 0.383 & 1.000 \end{pmatrix}.
\end{align*}

The total variation is $\ttr(\bS) = 197.2 + 297.2 = 494.4$.

To calculate these in R use, 'colMeans', 'cov', and 'cor'. These assume each column is a different variable, and each row a different observation.

```r
library(dplyr)
Ex1 <- data.frame(
  Student=LETTERS[1:5],
  P = c(41,72,46,77,59),
  S = c(63,82,38,57,85)
  )

Ex1 %>% select_if(is.numeric) %>% colMeans
```

```
##  P  S 
## 59 65
```

```r
Ex1 %>% select_if(is.numeric) %>% cov
```

```
##       P     S
## P 246.5 116.0
## S 116.0 371.5
```

```r
Ex1 %>% select_if(is.numeric) %>% cor
```

```
##           P         S
## P 1.0000000 0.3833276
## S 0.3833276 1.0000000
```

NOTE R USES 1/n-1 whereas hand calculation si 1/n.....


## Graphical techniques

We can draw histograms and scatter plots to view the distribution when $p=1$ and $p=2$ respectively.  For $p \geq 3$ the task is much harder.  One solution is a matrix of pair-wise scatter plots using the `pairs` command in R.  The graph below shows the relationship between goals scored (F), goals against (A) and points (PT) for 20 teams during a recent Premiership season.

![Scatter plots](figs/pairs.png)

WOULD BE BETTER TO CREATE THIS FIGURE HERE

You can also use the `plot3d` command in the `rgl` library to create an interactive 3D plot of the data.  The difficulty of displaying multivariate data is further motivation for developing a method for reducing the number of dimensions in the data.


## Random Vectors and Matrices


The *population mean vector* of the random vector $\bx$ is
$$\bmu = E(\bx).$$

The *population covariance matrix* of $\bx$ is
$$ \bSigma = \text{Var}(\bx) = E \lb (\bx-\bmu)(\bx-\bmu)^\top \rb.$$

The covariance between $\bx$ ($p \times 1$) and $\by$ ($q \times 1$) is
$$ \text{Cov}(\bx,\by) = E \lb (\bx - E(\bx))(\by - E(\by))^\top \rb. $$

Let $\bA$ ($q \times p$)  denote a constant matrix and let $\bb$ ($q \times 1$) denote a constant vector.  Then the following properties hold:

- $E(\bA \bx + \bb) = \bA E(\bx) + \bb=\bA \bmu +\bb$.
- $\bSigma = E(\bx \bx^\top) - \bmu \bmu^\top$.
- $\text{Var}(\bA \bx + \bb) = \bA \bSigma \bA^\top$, where $\bA$ is a $q \times p$ constant matrix.
- A covariance matrix $\bSigma$ is always non-negative
definite.  Moreover,  $\bSigma$ is positive definite if and only
if all its eigenvalues are positive, in which case its
determinant is positive and $\bSigma$ is non-singular.
- $\text{Cov}(\bx,\by) = E(\bx \by^\top) - E(\bx) E(\by)^\top$.
- $\text{Cov}(\bx,\bx) = \bSigma$.
 $\text{Cov}(\bx,\by) = \text{Cov}(\by,\bx)^\top$.
- If $\bx$ and $\by$ are independent then $\text{Cov}(\bx,\by) = {\mathbf 0}_{p,q}$.
- If $p=q$ then
$$
\text{Var}(\bx + \by) = \text{Var}(\bx) + \text{Var}(\by) + \text{Cov}(\bx,\by) + \text{Cov}(\by,\bx).
$$




## Unbiased Estimators
Recall from univariate statistics that an estimator $\hat{\theta}$ of a parameter $\theta$ is unbiased if $E(\hat{\theta}) = \theta$ for all $\theta$.  This concept readily transfers to the multivariate context.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-9"><strong>(\#prp:unnamed-chunk-9) </strong></span>Let $\bx_1, \ldots \bx_n$ be independent and identically distributed (i.i.d.), sampled from a population with mean $\bmu$ and covariance matrix $\bSigma$.  If $\bar{\bx}$ and $\bS$ are the sample mean and covariance matrix, respectively, then

1.  $E(\bar{\bx}) = \bmu$.

2.  $\text{Var}(\bar{\bx}) = {\ds \frac{1}{n}} \bSigma$.

3. $E(\bS) = {\ds \frac{n-1}{n}} \bSigma$.

\EndKnitrBlock{proposition}

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}**Part 1**  Since expectations behave linearly,
$$E(\bar{\bx}) = E \lb \frac{1}{n} \sum_{i=1}^n \bx_i \rb = \frac{1}{n} \sum_{i=1}^n E(\bx_i) = \frac{1}{n} n \bmu = \bmu.$$
**Part 2**
We have 
\begin{eqnarray*}
\text{Var}(\bar{\bx}) &=& \text{Var} \lb \frac{1}{n} \sum_{i=1}^n \bx_i \rb \\
&=& \sum_{i,j=1}^n \text{Cov}\left (\frac{1}{n}\bx_i, \frac{1}{n}\bx_j \right)\\
&=& \sum_{i=1}^n \text{Var} \lb \frac{1}{n} \bx_i \rb + \sum_{i \neq j} \text{Cov} \lb \frac{1}{n} \bx_i, \frac{1}{n} \bx_j \rb \\
&=& \frac{1}{n^2} \lb \sum_{i=1}^n \text{Var}(\bx_i) + \sum_{i \neq j} \text{Cov}(\bx_i,\bx_j) \rb \\
&=& \frac{1}{n^2} \lb \sum_{i=1}^n \text{Var}(\bx_i) \rb \\
&=& \frac{1}{n^2} n \bSigma \\
&=& \frac{1}{n} \bSigma.
\end{eqnarray*}

**Part 3**  From the definition of the sample covariance,

\begin{eqnarray*}
\bS &=& \frac{1}{n} \sum_{i=1}^n  (\bx_i -\bar{\bx})( \bx_i - \bar{\bx})^\top \\
&=& \frac{1}{n} \sum_{i=1}^n (\bx_i - \bmu+ \bmu -\bar{\bx})(\bx_i - \bmu +\bmu-\bar{\bx})^\top  \\
&=& \frac{1}{n} \sum_{i=1}^n \bigg \{(\bx_i - \bmu)(\bx_i - \bmu)^\top +(\bar{\bx}-\bmu)(\bar{\bx}-\bmu)^\top \\
&&  \qquad -(\bx_i-\bmu)(\bar{\bx}-\bmu)^\top -(\bar{\bx}-\bmu)(\bx_i-\bmu)^\top \bigg \}\\
&=& \frac{1}{n} \left \{\sum_{i=1}^n (\bx_i - \bmu)(\bx_i - \bmu)^\top \right \}- (\bar{\bx}-\bmu)(\bar{\bx}-\bmu)^\top
\end{eqnarray*}

Since $E(\bar{\bx})=\bmu$ and $\text{Var}(\bar{\bx})=n^{-1}\bSigma$, it follows that

\begin{eqnarray*}
E(\bS) &=& E\left \{\frac{1}{n}\sum_{i=1}^n  (\bx_i - \bmu)(\bx_i - \bmu)^\top \right \}\\
&& \qquad -E\left \{ (\bar{\bx}-\bmu)(\bar{\bx}-\bmu)^\top \right \}\\
&=&\left \{\frac{1}{n}\sum_{i=1}^n E \left [(\bx_i - \bmu)(\bx_i - \bmu)^\top \right ]\right \} -\text{Var}(\bar{\bx})\\
&=&E \left \{(\bx_1 - \bmu)(\bx_1  - \bmu)^\top \right \} -\text{Var}(\bar{\bx})\\
&=& \text{Var}(\bx_1)-\text{Var}(\bar{\bx}) \\
&=& \bSigma - \frac{1}{n}\bSigma \\
&=& \frac{n-1}{n}  \bSigma,
\end{eqnarray*}

which completes the proof.  
\EndKnitrBlock{proof}


  
 An implication of this theorem is that $\bar{\bx}$ is an unbiased estimator for $\bmu$ but that $\bS$ is a biased estimator of $\bSigma$.  Note, however, that ${\ds \frac{n}{n-1}} \bS$ is an unbiased estimator of $\bSigma$, i.e.
 $$
 \frac{n}{n-1}E[\bS]=\bSigma.
 $$


