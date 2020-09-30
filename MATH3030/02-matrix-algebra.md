# Review of Matrix algebra  {#linalg-prelim}

In this chapter blah blah


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


## Basic definitions

The matrix ${\mathbf A}$ will be referred to in the following equivalent ways:
\begin{eqnarray*}
{\mathbf A}=\stackrel{n\times p}{\mathbf A} &=& \left[\begin{array}{cccc}
a_{11}&a_{12}&\dots&a_{1p}\\
a_{21}&a_{22}&\dots&a_{2p}\\
\vdots&\vdots&&\vdots\\
a_{n1}&a_{n2}&\dots&a_{np}
\end{array} \right] \\
&=& [a_{ij}: i=1,\dots ,n;\ j=1,\dots ,p\,]=(a_{ij}),
\end{eqnarray*}
where the $a_{ij}$ are real numbers.

A matrix of order $1\times 1$ is called a *scalar*.

A matrix of order $n\times 1$ is called a *(column) vector*.

A matrix of order $1\times p$ is called a *(row) vector*.

e.g. $\stackrel{n\times 1}{\mathbf a}=\left(
\begin{array}{c}
a_1\\\vdots\\a_n
\end{array}
\right)$\quad is a column vector.



The $n\times n$ *identity matrix* ${\mathbf I}_n$ has diagonal elements equal to 1
and off-diagonal elements equal to zero.

A *diagonal* matrix is an $n \times n$ matrix whose
off-diagonal elements are zero.  Sometimes we denote a diagonal
matrix by $\tdiag\{a_1,\ldots, a_n\}$.

## Elementary matrix operations 


1. *Addition/Subtraction*.  If $\stackrel{n\times p}{\mathbf A}=[a_{ij}]$ and $\stackrel{n\times p}{\mathbf B}=[b_{ij}]$ are
given matrices then
$$ {\mathbf A}+{\mathbf B}=[a_{ij}+b_{ij}] \qquad \text{and} \qquad {\mathbf A}-{\mathbf B}=[a_{ij}-b_{ij}].$$

2. *Scalar Multiplication*.  If $\lambda$ is a scalar and ${\mathbf A}=[a_{ij}]$ then
$$\lambda {\mathbf A}=[\lambda a_{ij}].$$

3. *Matrix Multiplication*.  If $\stackrel{n\times p}{\mathbf A}$ and $\stackrel{p\times q}{\mathbf B}$ are
given matrices then ${\mathbf AB}=\stackrel{n\times q}{\mathbf C}=[c_{ij}]$ where
$$c_{ij}=\sum _{k=1}^p a_{ik}b_{kj}, \qquad i=1,\dots,n, \qquad j=1,\dots ,q.$$

4. *Matrix Transpose*.  If $\stackrel{m \times n}{\bA}=[a_{ij}: i=1, \ldots , m; j=1, \ldots , n]$, then the transpose of $\bA$, written
$\bA^\top$, is given by the $n \times m$ matrix
$$
\bA^\top =[a_{ji}: j=1, \ldots , n; i=1, \ldots, m].
$$
Note from the definitions that $({\mathbf AB})^\top={\mathbf B}^\top {\mathbf A}^\top$.
In the special case in which ${\mathbf A}=\stackrel{n\times 1}{\mathbf a}$ and
${\mathbf B}=\stackrel{n\times 1}{\mathbf b}$, the quantity ${\mathbf A}^\top {\mathbf B}={\mathbf a}^\top {\mathbf b}$ is real-valued and is
known as the *scalar product of $\mathbf a$ and $\mathbf b.$\newline <br>
Note that ${\mathbf a}^\top {\mathbf a}=\sum _{i=1}^n a_i^2\geq 0$ with equality iff
${\mathbf a}={\mathbf 0}_n$ where $\stackrel{n\times 1}{\mathbf 0}_n=(0,0,\dots ,0)^\top$.  We may
interpret ${\mathbf a}^\top {\mathbf a}$ as the (length)$^2$ of the vector ${\mathbf a}$.
The *norm* of a vector ${\mathbf a}$ is defined by
$$||{\mathbf a}||=({\mathbf a}^\top {\mathbf a})^{1/2}.$$
The scalar product has an alternative (but equivalent) representation:
$$
\ba^\top \bb = \vert \vert \ba \vert \vert \cdot \vert \vert \bb \vert \vert \cos(\theta),
$$
where $\theta$ is the angle (in radians) between the vectors $\ba$ and $\bb$.\newline <br>
A *unit vector* $\mathbf a$ is a vector satisfying $||{\mathbf a}||^2=
{\mathbf a}^\top {\mathbf a}=1$.\\
A matrix $\stackrel{n\times n}{\mathbf Q}$ satisfying ${\mathbf QQ}^\top = {\mathbf Q}^\top {\mathbf Q}={\mathbf I}_n$ is called an
*orthogonal matrix*.  Equivalently, a matrix $\mathbf Q$ is orthogonal iff ${\mathbf Q}^{-1}={\mathbf Q}^\top$.\newline <br>
If ${\mathbf Q}=[\bq_1,\ldots, \bq_n]$ is an orthogonal matrix, then the columns $\bq_1, \ldots, \bq_n$ are mutually orthogonal unit vectors, i.e.
$$
\bq_j^\top \bq_k=\begin{cases} 1 &\hbox{ if }  j=k\\
0 &\hbox{ if }   j \neq k \\
\end{cases}
$$

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:two1"><strong>(\#prp:two1) </strong></span>If $\bq_1, \ldots , \bq_n$ are mutually orthogonal $n \times 1$ unit vectors then
$$
\sum_{i=1}^n \bq_i \bq_i^\top = {\mathbf I}_n,
$$
the $n \times n$ identity matrix.
\EndKnitrBlock{proposition}

5. *Matrix Inverse*.  The inverse of a matrix $\stackrel{n\times n}{\mathbf A}$ (if it exists) is a
matrix $\stackrel{n\times n}{\mathbf B}$ such that ${\mathbf AB}={\mathbf BA}={\mathbf I}_n.$  We denote
the inverse by ${\mathbf A}^{-1}$.  Note that if ${\mathbf A}_1$ and ${\mathbf A}_2$ are both invertible,
then $({\mathbf A}_1 {\mathbf A}_2)^{-1}={\mathbf A}_2^{-1}{\mathbf A}_1^{-1}$.

6. *Trace*.  The trace of a matrix $\stackrel{n\times n}{\mathbf A}$ is given by
$$ \text{tr}({\mathbf A})=\sum _{i=1}^n a_{ii}.$$

\BeginKnitrBlock{lemma}
<span class="lemma" id="lem:unnamed-chunk-1"><strong>(\#lem:unnamed-chunk-1) </strong></span>For any matrices $\bA$ ($ n \times m$) and $\bB$ ($m \times n$),
$$
\text{tr}(\bA \bB) = \text{tr}(\bB \bA).
$$
\EndKnitrBlock{lemma}



## Linear independence and determinants

\noindent Vectors $\stackrel{n\times 1}{\mathbf x}_1 ,\dots , \stackrel{n\times 1}{\mathbf x}_p$
are said to be *linearly dependent* if there exist scalars
$\lambda _1, \dots ,\lambda _p$ not all zero such that
$$ \lambda _1 {\mathbf x}_1+\lambda _2 {\mathbf x}_2+ \dots + \lambda _p {\mathbf x}_p={\mathbf 0}.$$
Otherwise, these vectors are said to be *linearly independent*.

The *rank* of a matrix is equal to the maximum number of linearly independent rows (equivalently, columns).

The *determinant* of a square matrix $\stackrel{n\times n}{\mathbf A}$ is
defined as
$$ \tdet ({\mathbf A})=\sum (-1)^{|\tau |} a_{1\tau(1)}\dots a_{n\tau (n)} $$
where the summation is taken over all permutations $\tau$ of $\{1,2,\dots ,n\}$,
and we define $|\tau |=0$ or $1$ depending on whether $\tau$ can be written as an even or
odd number of transpositions.

E.g. If ${\mathbf A}=\left[
\begin{array}{cc}
a_{11}&a_{12}\\
a_{21}&a_{22}
\end{array}
\right]$,
then $\tdet ({\mathbf A})=a_{11}a_{22}-a_{12}a_{21}$.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-2"><strong>(\#prp:unnamed-chunk-2) </strong></span>For any matrices $\stackrel{n\times n}{\mathbf A}$,
$\stackrel{n\times n}{\mathbf B}$, $\stackrel{n\times n}{\mathbf C}$ such that ${\mathbf C}={\mathbf{AB}}$,
$$ \tdet ({\mathbf C})=\tdet ({\mathbf A}) \cdot \tdet ({\mathbf B}).$$
\EndKnitrBlock{proposition}

## Eigenvalues and eigenvectors 

If $\stackrel{n\times n}{\mathbf A}$ is a given matrix then
$R(\lambda )=\tdet ({\mathbf A}-\lambda {\mathbf I}_n)$ is an $n^{\text{th}}$ order polynomial in $\lambda$.
The $n$ roots $\lambda _1, \dots , \lambda _n$ of $R(\lambda )$ (possibly complex
numbers) are called *eigenvalues* of $\mathbf A$.

For any eigenvalue $\lambda$ of $\stackrel{n\times n}{\mathbf A}$, there exists
a non-zero vector $\stackrel{n\times 1}{\mathbf x}$, called an *eigenvector*, such
that ${\mathbf A} {\mathbf x} = \lambda {\mathbf x}$.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-3"><strong>(\#prp:unnamed-chunk-3) </strong></span>If $\mathbf A$ is symmetric (i.e.\ ${\mathbf A}^\top ={\mathbf A}$) then the
eigenvalues of $\mathbf A$ are all *real* and all the eigenvectors of $\mathbf A$ have
*real* components.
\EndKnitrBlock{proposition}

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-4"><strong>(\#prp:unnamed-chunk-4) </strong></span>The rank of a symmetric matrix is equal to the number of non-zero eigenvalues.
\EndKnitrBlock{proposition}

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:spectraldecomp"><strong>(\#prp:spectraldecomp) </strong></span>**(Spectral decomposition theorem)**.  Any symmetric matrix $\stackrel{n\times n}{\mathbf A}$ can
be written
$$ {\mathbf A}={\mathbf Q \Lambda Q}^\top = \sum _{i=1}^{n} \lambda _i {\mathbf q}_i {\mathbf q}_i^\top ,$$
where $\stackrel{n\times n}{\mathbf \Lambda}=\tdiag \{ \lambda _1, \dots , \lambda _n \}$ is
a diagonal matrix consisting of the eigenvalues of $\mathbf A$ and $\stackrel{n\times n}{\mathbf Q}$ is
an orthogonal matrix whose columns are unit eigenvectors
${\mathbf q}_1, \dots , {\mathbf q}_n$ of $\mathbf A$.
\EndKnitrBlock{proposition}


\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-5"><strong>(\#prp:unnamed-chunk-5) </strong></span>If $\stackrel{n\times n}{\mathbf A}$ is a symmetric matrix
then its determinant is the product of its eigenvalues, i.e.\ $\tdet ({\mathbf A})=\lambda _1 \dots
\lambda _n$.
\EndKnitrBlock{proposition}

Let $\stackrel{n\times n}{\mathbf A}$ be a symmetric matrix
with (necessarily real) eigenvalues $\lambda _1 \geq \lambda _2 \geq
\dots \geq \lambda _n$.  Then $\mathbf A$ is said to be *positive definite*
if and only if  $\lambda _n >0$, and $\bA$ is said to be *non-negative definite* if and only if $\lambda _n\geq 0$.


For a given symmetric $\stackrel{n\times n}{\mathbf A}$, define the
*quadratic form*
$\mathcal{Q}({\mathbf x})={\mathbf x}^\top {\mathbf A} {\mathbf x}$.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:two8"><strong>(\#prp:two8) </strong></span>In the above notation,
$$\displaystyle{\tmax_{{\mathbf x}:{\mathbf x}^\top {\mathbf x}=1}} \mathcal{Q}({\mathbf x})=\lambda_1, $$
where the maximum occurs at $\bx=\pm \bq_1$.
\EndKnitrBlock{proposition}

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-6"><strong>(\#prp:unnamed-chunk-6) </strong></span>In the above notation,
$$\displaystyle{\tmin_{{\mathbf x}:{\mathbf x}^\top {\mathbf x}=1}} \mathcal{Q}({\mathbf x})=\lambda_n,$$
where the minimum occurs at $\bx = \pm \bq_n$.
\EndKnitrBlock{proposition}

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-7"><strong>(\#prp:unnamed-chunk-7) </strong></span>We have (i)  $\mathcal{Q}(\bx)>0$ for all
$\bx\neq {\mathbf 0}_n$ if and only if $\mathbf A$ is positive  definite; and
(ii) $\mathcal{Q}(\bx)\geq 0$ for all $\bx$ if and only if $\bA$ is non-negative definite.
\EndKnitrBlock{proposition}

A matrix $\stackrel{n \times n}{\bP}$ is a *projection*
matrix if it is symmetric, i.e. $\bP^\top =\bP$, and
$$
\bP^2 =\bP.
$$


\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-8"><strong>(\#prp:unnamed-chunk-8) </strong></span>The eigenvalues of a projection matrix $\bP$ are all $0$ or $1$.
\EndKnitrBlock{proposition}

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-9"><strong>(\#prp:unnamed-chunk-9) </strong></span>If $\stackrel{n \times n}{\bP}$ is a projection matrix then ${\mathbf I}_n - \bP$ is also
a projection matrix.
\EndKnitrBlock{proposition}


## Matrix square roots

From time to time we shall need to consider square roots of  symmetric non-negative definite matrices.  From Proposition \@ref(prp:spectraldecomp), a symmetric
matrix $\mathbf A$ may be written as
 $\mathbf A=Q \Lambda Q^\top$ where $\mathbf \Lambda$ is a diagonal matrix and $\mathbf Q$ is an orthogonal matrix.  Moreover, $\bA$ is non-negative definite if and only if the diagonal elements of $\mathbf \Lambda$ (the
eigenvalues of $\mathbf A$) are all non-negative.

For such a matrix we define ${\bA}^{1/2}$, a matrix square root of $\bA$, by ${\bA}^{1/2}=\bQ \bLambda^{1/2} \bQ^\top$ where $\bLambda^{1/2}=\text{diag}\{\lambda_1^{1/2}, \ldots , \lambda_n^{1/2}\}$.  This definition makes sense because
\begin{align*}
\bA^{1/2}\bA^{1/2}&=\bQ \bLambda^{1/2}\bQ^\top \bQ \bLambda^{1/2} \bQ^\top\\
&=\bQ \bLambda^{1/2}\bLambda^{1/2}\bQ^\top\\
&=\bQ \bLambda \bQ^\top\\
&=\bA,
\end{align*}
where $\bQ^\top \bQ=\bI_n$ and $\bLambda^{1/2}\bLambda^{1/2}=\bLambda$.  The matrix $\bA^{1/2}$ is not the only matrix square root of $\bA$, but it *is* the only symmetric, non-negative definite square root of $\bA$.

If $\bA$ is positive definite (as opposed to just non-negative definite), then all the $\lambda_i$ are positive and so we can also define $\bA^{-1/2}=\bQ \bLambda^{-1/2}\bQ^\top$ where $\bLambda^{-1/2}=\text{diag}\{\lambda_1^{-1/2},\ldots , \lambda_n^{-1/2}\}$.  Note that
$$
\bA^{-1/2}\bA^{-1/2}=\bQ\bLambda^{-1/2}\bQ^\top \bQ \bLambda^{-1/2}\bQ^\top =\bQ \bLambda^{-1}\bQ^\top=\bA^{-1},
$$
so that, as defined above,  $\bA^{-1/2}$ the matrix square root of $\bA^{-1}$.  Furthermore, similar calculations show that
$$
\bA^{1/2}\bA^{-1/2}=\bA^{-1/2}\bA^{1/2}=\bI_n,
$$
so that $\bA^{-1/2}$, as defined above, is the matrix inverse of $\bA^{1/2}$.

## Singular Value Decomposition 

The spectral decomposition theorem (Proposition \@ref(prp:spectraldecomp)) gives a decomposition of any symmetric matrix. We now give a generalisation of this result which applies to *all* matrices.  We will need this extra generality in Chapter 4 on Canonical Correlation Analysis.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-10"><strong>(\#prp:unnamed-chunk-10) </strong></span>**(Singular value decomposition)**.
Let $\bA$ be a $p \times q$ matrix of rank $t$, where $1 \leq t \leq \min(p,q)$.  Then there exists a $p \times t$ matrix $\bQ=[\bq_1,\ldots , \bq_t]$, a $q \times t$ matrix $\bR=[{\mathbf r}_1,\ldots ,{ \mathbf r}_t]$,  and a $t \times t$  diagonal matrix ${\mathbf \Xi}=\text{diag}\{\xi_1,\ldots , \xi_t\}$ such that
$$
\bA=\bQ {\mathbf \Xi} \bR^\top =\sum_{i=1}^t \xi_i \bq_i {\mathbf r}_i^\top,
$$
where $\bQ^\top \bQ = \bI_t = \bR^\top \bR$ and the $\xi_1 \geq \ldots \geq \xi_t >0$.
\EndKnitrBlock{proposition}


Note that the $\bq_i$ and the ${\mathbf r}_i$ are necessarily unit vectors.

The scalars $\xi_1, \ldots , \xi_t$ are called *singular values*.

When $\bA$ is symmetric, we take ${\mathbf R}=\bQ$, and  the spectral decomposition theorem is recovered, and in this case (but not in general) the singular values of $\bA$ are in fact eigenvalues of $\bA$.

The following result relates  $\mathbf Q$, $\mathbf \Xi$ and $\mathbf R$ to certain eigenvalues and eigenvectors.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-11"><strong>(\#prp:unnamed-chunk-11) </strong></span>Let $\bA$ be any matrix of rank $t$.  Then the non-zero eigenvalues of both $\bA \bA^\top$ and $\bA^\top \bA$ are given by $\xi_1^2, \ldots , \xi_t^2$; the corresponding unit eigenvectors of $\bA \bA^\top$ are given by the columns of $\mathbf Q$; and the corresponding unit eigenvectors of $\bA^\top \bA$ are given by the columns of $\mathbf R$.
\EndKnitrBlock{proposition}

The following result is important in Canonical Correlation Analysis.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-12"><strong>(\#prp:unnamed-chunk-12) </strong></span>For any matrix $\bA$ of rank $t$ with singular values $\xi_1 \geq \xi_2 \geq \ldots \geq \xi_t >0$, we have
$$
\max_{\bx, \by:\, \vert \vert \bx \vert \vert=\vert \vert \by \vert \vert =1} \bx^\top \bA \by =\xi_1.
$$
\EndKnitrBlock{proposition}

## The Centering Matrix

From time to time in this module an important role will be played by the *centering matrix*
\begin{equation}
\bH=\bI_n - \frac{1}{n} {\mathbf 1}_n {\mathbf 1}_n^\top.
(\#eq:Hcentre)
\end{equation}
Note that, in the above, $\bI_n$ is the $n \times n$ identity matrix, while ${\mathbf 1}_n$ is an $n \times 1$ column vector of ones.

The reason for the terminology *centering* will become clear below.

Important properties of the matrix $\bH$ in Equation \@ref(eq:Hcentre) are now listed.  These properties are proved in the example sheets.


1. The matrix $\bH$ is a projection matrix, i.e. $\bH^\top =\bH$ and $\bH^2=\bH$.
2. Writing ${\mathbf 0}_n$ for the $n \times 1$ vector of zeros, we have
$\bH {\mathbf 1}_n={\mathbf 0}_n$ and  ${\mathbf 1}_n^\top \bH={\mathbf 0}_n^\top.$
3. If $\bx=(x_1, \ldots , x_n)^\top$, then $\bH\bx = \bx - \bar{x}{\mathbf 1}_n$ where $\bar{x}=n^{-1}\sum_{i=1}^n x_i$.
4.  With $\bx$ as in (iii), we have
$$
\bx^\top \bH \bx = \sum_{i=1}^n (x_i-\bar{x})^2,
$$
and so
$$
\frac{1}{n}\bx^\top \bH \bx=\frac{1}{n}\sum_{i=1}^n (x_i-\bar{x})^2 = \hat{\sigma}^2,
$$
where $\hat{\sigma}^2$ is the sample variance.
5. If $\bX=[\bx_1, \ldots , \bx_n]^\top$ is an $n \times p$ data matrix then
$$
\bH\bX=\left[ \begin{array}{c}
(\bx_1-\bar{\bx})^\top\\
(\bx_2 -\bar{\bx})^\top\\
..\\
..\\
..\\
(\bx_n - \bar{\bx})^\top
\end{array}\right ]= [ \bx_1 -\bar{\bx}, \ldots , \bx_n-\bar{\bx}]^\top
$$
6. With $\bX$ as in (v),
$$
\frac{1}{n}\bX^\top \bH \bX =\frac{1}{n} \sum_{i=1}^n (\bx_i -\bar{\bx})(\bx_i -\bar{\bx})^\top =\bS,
$$
where $\bS$ is the sample covariance matrix.
7.  If $\bA=(a_{ij})_{i,j=1}^n$ is a symmetric  $n \times n$ matrix, then
$$
\bB=\bH \bA\bH= \bA - {\mathbf 1}_n \bar{\ba}_+^\top -\bar{\ba}_+{\mathbf 1}_n^\top +\bar{a}_{++}{\mathbf 1}_n {\mathbf 1}_n^\top,
$$
or, equivalently,
$$
b_{ij}=a_{ij}-\bar{a}_{i+}-\bar{a}_{+j}+\bar{a}_{++}, \qquad i,j=1, \ldots , n,
$$
where
$$
\bar{\ba}_{+}\equiv (\bar{a}_{1+}, \ldots , \bar{a}_{n+})^\top=\frac{1}{n}\bA{\mathbf 1}_n,
$$
$\bar{a}_{+j}=\bar{a}_{j+}$ \, for \, $j=1, \ldots , n$,\, and \, $\bar{a}_{++}=n^{-2}\sum_{i,j=1}^n a_{ij}$.

Note that Property 3. is a special case of Property 5., and Property 4. is a special case of Property 6.
However, it is useful to see these results in the simpler scalar case before moving onto the the general matrix case.

## Quadratic forms and ellipses

A standard ellipse in $\mathbb{R}^2$ is given by the equation
$$
\frac{x^2}{a^2}+\frac{y^2}{b^2}=1 \quad (a>b>0).
$$
The interior (the shaded region) is given by
\begin{equation}
\frac{x^2}{a^2}+\frac{y^2}{b^2}\leq 1.   (\#eq:ellipse)
\end{equation}
Note that a standard ellipse has axes of symmetry given by the $x$-axis and $y$-axis
(if $a>b$, the former is the major axis and the latter the minor axis).

If we define
${\mathbf A}=\left[
\begin{array}{cc}
a^2&0\\
0&b^2
\end{array}
\right]$
then Equation \@ref(eq:ellipse) can be written in the form
$$ \binom{x}{y}^\top {\mathbf A}^{-1} \binom{x}{y}\leq 1. $$
If we write ${\mathbf x}=\binom{x_1}{x_2}$ and generalise to an arbitrary symmetric
positive definite matrix $\stackrel{2 \times 2}{\mathbf A}$, what is the set
$$ \{ {\mathbf x} \in \mathbb{R}^2 : {\mathbf x}^\top {\mathbf A}^{-1} {\mathbf x} \leq 1\} ? $$
We get a rotated ellipse with axes of symmetry given by the eigenvectors of $\mathbf A$,
with the major axis determined by the eigenvector corresponding to the larger
eigenvalue of $\mathbf A$, and the minor axis determined by the eigenvector corresponding
to the smaller eigenvalue of $\mathbf A$.

Note that, for $c>0$,
$${\mathbf x}^\top {\mathbf A}^{-1} {\mathbf x}\leq c \qquad \Leftrightarrow \qquad {\mathbf x}^\top (c{\mathbf A})^{-1}{\mathbf x}\leq 1 ,$$
where $c{\mathbf A}$ is a scalar multiple of $\mathbf A$.

If ${\mathbf m}$ is a fixed 2-vector, then what is the set
$$
\{ {\mathbf x} \in \mathbb{R}^2 : ({\mathbf x}-{\mathbf m})^\top {\mathbf A}^{-1}({\mathbf x}-{\mathbf m})\leq 1\} ?
$$
Since
$$
\{ {\mathbf x} : ({\mathbf x}-{\mathbf m})^\top {\mathbf A}^{-1}({\mathbf x}-{\mathbf m})\leq 1 \}=\{ {\mathbf z}+
{\mathbf m} : {\mathbf z}^\top {\mathbf A}^{-1} {\mathbf z}\leq 1 \} ,
$$
it follows that
$$
\{ {\mathbf x} : ({\mathbf x}-{\mathbf m})^\top {\mathbf A}^{-1}({\mathbf x}-{\mathbf m})\leq 1\}
$$
is just the ellipse $\{ {\mathbf z}:{\mathbf z}^\top {\mathbf A}^{-1}{\mathbf z}\leq 1\}$ translated by
${\mathbf m}$.

Analogous results for ellipsoids and quadratic forms hold in three and higher dimensions.

## Lines and Hyperplanes in $\mathbb{R}^p$

For any $\ba, \bb \in \mathbb{R}^p$, the set
\begin{equation}
\mathcal{L}=\mathcal{L}(\ba, \bb)=\{\ba+\gamma \bb: \gamma \in \mathbb{R}\}
(\#eq:straightl) 
\end{equation}
is a *straight line* in $\mathbb{R}^p$.

If $\ba^\top \bb=0$, i.e. $\ba$ and $\bb$ are orthogonal, then $\ba$ is the perpendicular from the origin ${\mathbf 0}_p$
to the line $\mathcal{L}(\ba,\bb)$.

For fixed $\ba \in \mathbb{R}^p$ and $\gamma \in \mathbb{R}$,
$$
\mathcal{H}=\mathcal{H}(\ba, \gamma) =\{\bx \in \mathbb{R}^p:\, \ba^\top \bx =\gamma\}
$$
is a hyperplane of dimension $p-1$ in $\mathbb{R}^p$.  The vector $\ba$ is the perpendicular from the origin ${\mathbf 0}_p$ to the hyperplane $\mathcal{H}(\ba, \gamma)$.

There is an alternative way to define hyperplanes in $\mathbb{R}^P$.  Suppose that, for $1 \leq r <p$, $\stackrel{p \times 1}{\ba}_1,
\ldots , \stackrel{p \times 1}{\ba}_r, \stackrel{p \times 1}{\ba}_{r+1}$ are linearly independent.  Then
$$
\mathcal{H}=\left \{ \sum_{j=1}^{r+1} \gamma_j \ba_j: \, \sum_{j=1}^{r+1}\gamma_j =1  \right \}
$$
is an $r$-dimensional hyperplane in $\mathbb{R}^p$.

When $r=1$, using the fact that $\gamma_1+\gamma_2=1$, we may write
$$
\gamma_1 \ba_1 + \gamma_2 \ba_2=(1-\gamma_2)\ba_1 + \gamma_2 \ba_2 = \ba_1 +\gamma_2(\ba_2-\ba_1),
$$
which agrees with $\ba+\gamma \bb$ in \@ref(eq:straightl) when $\ba = \ba_1$, $\bb = \ba_2 -\ba_1$ and $\gamma=\gamma_2$.  So we have shown that the two definitions agree in the case of a straight line.

## Vector Differentiation

 Consider a real-valued function $f: \mathbb{R}^p \rightarrow \mathbb{R}$ of a vector variable $\bx=(x_1, \ldots , x_p)^\top$.  Sometimes we will want to differentiate $f$.    We define the partial derivative of $f(\bx)$ with respect to $\bx$ to be
the vector of partial derivatives, i.e.
\begin{equation}
\frac{\partial f}{\partial \bx}(\bx)=\left [ \begin{array}{c} \frac{\partial f}{\partial x_1}(\bx)\\
 ..\\
 ..\\
 ..\\
 \frac{\partial f}{\partial x_p}(\bx)
\end{array} \right ]
(\#eq:derivx)
\end{equation}
The following examples can be worked out directly from the definition \@ref(eq:derivx), using the chain rule in some cases.

\BeginKnitrBlock{example}
<span class="example" id="exm:unnamed-chunk-13"><strong>(\#exm:unnamed-chunk-13) </strong></span>If $f(\bx)=\ba^\top \bx$ where $\ba \in \mathbb{R}^p$ is a constant vector, then
$$
\frac{\partial f}{\partial \bx}(\bx)=\ba.
$$
\EndKnitrBlock{example}

\BeginKnitrBlock{example}
<span class="example" id="exm:unnamed-chunk-14"><strong>(\#exm:unnamed-chunk-14) </strong></span>If $f(\bx)=(\bx-\ba)^\top \bA (\bx-\ba)$ for a fixed vector $\ba \in \mathbb{R}^p$
and $\bA$ is a constant symmetry $p \times p$ matrix, then
$$
\frac{\partial f}{\partial \bx}(\bx)=2\bA (\bx-\ba).
$$
\EndKnitrBlock{example}

\BeginKnitrBlock{example}
<span class="example" id="exm:unnamed-chunk-15"><strong>(\#exm:unnamed-chunk-15) </strong></span>Suppose that $g: \, \mathbb{R} \rightarrow \mathbb{R}$ is a differentiable function with derivative $g^\prime$.  Then, using the chain rule for partial derivatives,
$$
\frac{\partial g(\ba^\top \bx)}{\partial \bx}=g^{\prime}(\ba^\top\bx)\frac{\partial}{\partial \bx}\left \{\ba^\top \bx \right \}=g^{\prime}(\ba^\top\bx) \ba.
$$
\EndKnitrBlock{example}

\BeginKnitrBlock{example}
<span class="example" id="exm:unnamed-chunk-16"><strong>(\#exm:unnamed-chunk-16) </strong></span>If $f$ is defined as in Example 2 and $g$ is as in Example 3 then, using the chain rule again,
$$
\frac{\partial }{\partial \bx} g\{f(\bx)\}=g^{\prime} \{f(\bx)\}\frac{\partial f}{\partial \bx}(\bx)
=2 g^{\prime}\{(\bx - \ba)^\top \bA (\bx - \ba)\}\bA(\bx -\ba).
$$
\EndKnitrBlock{example}


If we wish to find a maximum or minimum of $f(\bx)$ we should search for stationary points of $f$, i.e.
solutions to the system of equations
$$
\frac{\partial f}{\partial \bx}(\bx)\equiv \left [ \begin{array}{c} \frac{\partial f}{\partial x_1}(\bx)\\
 ..\\
 ..\\
 ..\\
 \frac{\partial f}{\partial x_p}(\bx)
\end{array} \right ]={\mathbf 0}_p.
$$
The nature of a stationary point is determined by the *Hessian*, i.e. the matrix of second derivatives.
The Hessian is the $p \times p$ matrix
$$
\frac{\partial^2f}{\partial \bx \partial \bx^\top}(\bx) =\left \{ \partial^2 f(\bx)/\partial x_j \partial x_k\right \}_{j,k=1}^p.
$$

If the Hessian is positive (negative) definite at a stationary point $\bx$, then the stationary point is a minimum (maximum).

If the Hessian has both positive an negative eigenvalues at $\bx$ then the stationary point will be a *saddle point*.

