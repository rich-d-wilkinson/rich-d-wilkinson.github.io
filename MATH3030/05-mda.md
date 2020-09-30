
# Multidimensional Scaling {#mds}


In this chapter our starting point is somewhat different.  Suppose we have a sample of $n$ experimental units
and we have a way to measure `distance' or `dissimilarity' between any pair of experimental units $i$ and $j$, leading to a measure of distance or dissimilarity $d_{ij}$, $i,j=1, \ldots , n$.  The starting point for Multidimensional Scaling (MDS) is a distance matrix $\bD=(d_{ij}: \, i,j=1, \ldots , n)$.    A key goal in MDS is to determine coordinates of  a set of points in a low-dimensional Euclidean space, e.g. $\mathbb{R}$ or $\mathbb{R}^2$, whose inter-point distances (or dissimilarities) are approximately equal to the $d_{ij}$.  Using this approximate approach we are able to perform a statistical study of the original experimental units in a lower-dimensional space than the original one.  We shall also see that there is a close connection between MDS and PCA.

## Multidimensional Scaling

We call an $n \times n$ matrix $\bD=(d_{ij})_{i,j=1}^n$ a **distance matrix** or, equivalently, a **dissimilarity matrix**, if the following properties are satisfied:

1. For $i=1, \ldots , n$, $d_{ii}=0$.
2.  Symmetry: $d_{ij}=d_{ji} \geq 0$ for all $i,j=1,\ldots, n$.
3.  Definiteness:  $d_{ij}=0$ implies $i=j$.


A comment on our terminology.  We do not require distances necessarily to satisfy the triangle inequality
\begin{equation}
d_{ik} \leq d_{ij}+d_{jk}.
(\#eq:triangle)
\end{equation}
A distance function which always satisfies the triangle inequality is called a **metric distance**
or just a **metric**, and a distance function which does not always satisfy the triangle inequality is called
**non-metric** distance.

Suppose $\bx_1,\ldots , \bx_n$ are points in $\mathbb{R}^p$.  If the $d_{ij}$ are of the form
$$
d_{ij}=\vert \vert \bx_i -\bx_j \vert \vert =\sqrt{(\bx_i-\bx_j)^\top (\bx_i-\bx_j)}.
$$
Then each $d_{ij}$ is called a **Euclidean distance** and, in this case, $\bD$ is called a **Euclidean distance matrix**.  Since Euclidean distances satisfy the triangle
inequality \@ref(eq:triangle), it follows that Euclidean distance is a metric distance.

Given a distance matrix ${\mathbf D}=\{d_{ij}\}_{i,j=1}^n$, define the matrix
\begin{equation}
\bA=\{a_{ij}\}_{i,j=1}^n,  \quad \text{where} \qquad a_{ij}=-\frac{1}{2}d_{ij}^2.
(\#eq:defA)
\end{equation}
Note that, for $i=1,\ldots , n$, $a_{ii}=-d_{ii}^2/2=0$.

Now define the matrix
\begin{equation}
{\mathbf B}={\mathbf H} \bA {\mathbf H},
(\#eq:defB)
\end{equation}
where
\begin{equation}
{\mathbf H}={\mathbf I}_n -n^{-1}{\mathbf 1}_n {\mathbf 1}_n^\top
(\#eq:defH)
\end{equation}
is the $n \times n$  **centering matrix**; see \S 2.7. For reasons that will soon become clear, $\bB$ defined by \@ref(eq:defB) is known as a centred inner-product matrix.

Let $\bx_1, \ldots , \bx_n$ denote $n$ points in $\mathbb{R}^p$.  Then the $n \times p$ matrix
$\mathbf X=[\bx_1, \ldots , \bx_n]^\top$ is the data matrix, as before.

We now present the key result for classical MDS.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:five1"><strong>(\#prp:five1) </strong></span>Let $\bD$ denote an $n \times n$ distance matrix and suppose $\bA$, $\bB$ and $\mathbf H$ be as defined in \@ref(eq:defA), \@ref(eq:defB) and \@ref(eq:defH), respectively.

1. The matrix $\bD$ is a Euclidean distance matrix if and only if $\bB$ is a non-negative definite matrix.

2. If $\bD$ is a Euclidean distance matrix for the sample of $n$ vectors $\bx_1,\ldots , \bx_n$,  then
\begin{equation}
b_{ij}=(\bx_i-\bar{\bx})^\top (\bx_j - \bar{\bx}), \qquad i,j=1,\ldots , n,
(\#eq:bijB)
\end{equation}
where $\bar{\mathbf x}=n^{-1}\sum_{i=1}^n \bx_i$ is the sample mean vector.  Equivalently, we may write
$$
\bB = ({\mathbf H} {\mathbf X})({\mathbf H} {\mathbf X})^\top,
$$
where ${\mathbf X}=[\bx_1,\ldots , \bx_n]^\top$ is the data matrix, and $\bH$ is the $n \times n$ centering matrix.  Consequently,
$\bB$ is non-negative definite.

3. Suppose $\bB$ is non-negative definite with positive eigenvalues $\lambda_1 \geq \lambda_2 \cdots \geq \lambda_k$ and spectral decomposition $\bB={\mathbf Q} {\pmb \Lambda}{\mathbf Q}^\top$, where ${\pmb \Lambda}=\text{diag}\{\lambda_1 \ldots \lambda_k\}$ and $\mathbf Q$ is $n \times k$ and satisfies ${\mathbf Q}^\top {\mathbf Q}={\mathbf I}_k$.  Then ${\mathbf X}=[\bx_1, \ldots , \bx_n]^\top={\mathbf Q}{\pmb \Lambda}^{1/2}$ is an $n \times k$ data matrix for points $\bx_1, \ldots , \bx_n$ in $\mathbb{R}^k$, which have inter-point distances given by $\bD=(d_{ij})$. Moreover, for this data matrix $\bar{\bx}={\mathbf 0}_k$ and $\bB$ represents the inner product matrix with elements given by \@ref(eq:bijB).

\EndKnitrBlock{proposition}

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}Part 1. is a direct consequence of parts 2. and 3.  Parts 2. and 3. are proved in the example sheets.
\EndKnitrBlock{proof}


**Important Point**: Proposition \@ref(prp:five1) may be useful even if ${\mathbf D}$ is not a Euclidean distance matrix, in which case $\bB$ has some negative eigenvalues.  What we can do is to replace $\bB$ by its positive part.  If $\bB$ has spectral decomposition $\sum_{j=1}^p \lambda_j \bq_j \bq_j^\top$, then its positive definite part is defined by
$$
\bB_{\text{pos}}=\sum_{j: \, \lambda_j>0} \lambda_j \bq_j \bq_j^\top.
$$
In other words, we sum over those $j$ such that $\lambda_j$ is positive.
Then $\bB_{\text{pos}}$ is non-negative definite and so we can use  Theorem 5.1(iii) to determine a Euclidean configuration which has centred inner-product matrix $\bB_{\text{pos}}$.  Then, provided the negative eigenvalues are small in absolute value relative to the positive eigenvalues, the inter-point distances of the new points in Euclidean space should provide a good approximation to the original inter-point distances $(d_{ij})$.


\BeginKnitrBlock{example}
<span class="example" id="exm:mdsm1"><strong>(\#exm:mdsm1) </strong></span>Consider the five point in $\mathbb{R}^2$:
$$
\bx_1=(0,0)^\top,  \bx_2 =(1,0)^\top, \quad \bx_3 =(0,1)^\top
$$
$$
\bx_4 =(-1,0)^\top \quad \text{and} \quad \bx_5=(0,-1)^\top.
$$
The resulting distance matrix is
$$
\bD=\left [ \begin{array}{ccccc}
0&1&1&1&1\\
1&0&\sqrt{2}&2&\sqrt{2}\\
1&\sqrt{2}&0&\sqrt{2}&2\\
1&2&\sqrt{2}&0&\sqrt{2}\\
1&\sqrt{2}&2&\sqrt{2}&0
\end{array} \right ].
$$
Using \@ref(eq:defA) first to calculate $\bA$, and then using \@ref(eq:defB) to calculate $\bB$, we find that
$$
\bA=-\left [ \begin{array}{ccccc}
0&0.5&0.5&0.5&0.5\\
0.5&0&1&2&1\\
0.5&1&0&1&2\\
0.5&2&1&0&1\\
0.5&1&2&1&0
\end{array} \right ]
$$
and
$$
\bB=\left [ \begin{array}{ccccc}
 0& 0&0&0&0\\
0&1&0&-1&0\\
0&0&1&0&-1\\
0&-1&0&1&0\\
0&0&-1&0&1
\end{array} \right ].
$$
Further numerical calculations using R show that the eigenvalues of $\bB$ are
$$
\lambda_1=\lambda_2=2 \qquad \text{and} \qquad \lambda_3=\lambda_4=\lambda_5=0.
$$
Note that, as expected from Proposition \@ref(prp:five1),  $\bB$ is non-negative definite because it is a Euclidean
distance matrix.

The following mutually orthogonal unit eigenvectors corresponding to the repeated eigenvalue $2$
are produced by R:
$$
\bq_1= \begin{pmatrix}0 \\ -0.439 \\ -0.554 \\ 0.439 \\ 0.554 \end{pmatrix} \qquad
\text{and} \qquad \bq_2 =\begin{pmatrix}0 \\ 0.554 \\ -0.439 \\ -0.554\\ 0.439 \end{pmatrix}.
$$
So the coordinates of five points in $\mathbb{R}^2$ which have the same inter-point distance matrix, $\bD$, as the original five points in $\mathbb{R}^2$, are given by the rows of the matrix
$$
\bQ \bLambda^{1/2}=\sqrt{2}[\bq_1 , \bq_2]=
\begin{pmatrix}
0&0\\
-0.621 & 0.784\\
-0.784 & -0.621\\
0.621 & -0.784\\
0.784 & 0.621
\end{pmatrix}.
$$
In the example sheets you asked to verify that there is an orthogonal transformation which maps the original five points onto the new five points.
\EndKnitrBlock{example}

## Principal Coordinates

Starting with a distance matrix $\bD$, and using the matrix $\bB$, we now show how to calculate exact or approximate Euclidean coordinates for the $n$ objects under study.  We already know from Proposition \@ref(prp:five1)  how to do this when the distance matrix $\bD$ is Euclidean, but we will see now that this construction works more generally.  Moreover, there is a very close connection with principal components analysis.


- **Step 1**: Given a distance matrix $\stackrel{n \times n}{\bD}$, calculate $\bA$ according to \@ref(eq:defA).

- **Step 2**:  Calculate $\bB=(b_{ij})_{i,j=1}^n$ in \@ref(eq:defB) using
$$
b_{ij}=a_{ij}-\bar{a}_{i+}-\bar{a}_{+j}+\bar{a}_{++}, \qquad i,j=1, \ldots ,n,
$$
where
$$
\bar{a}_{i+}=n^{-1}\sum_{j=1}^n a_{ij}, \quad \bar{a}_{+j}=n^{-1}\sum_{i=1}^n a_{ij}\quad
\text{and} \quad \bar{a}_{++}=n^{-2}\sum_{i,j=1}^n a_{ij}.
$$

- **Step 3**:  Assume that the $k$ largest eigenvalues of $\bB=(b_{ij})_{i,j=1}^n$, $\lambda_1 > \lambda_2 > \cdots > \lambda_k$ are all positive and have associated unit eigenvectors $\bv_1, \ldots , \bv_k$.

- **Step 4**: Define $\bV=[\bv_1 , \ldots , \bv_k]$ and
$$
\bX \equiv [\bx_1, \ldots , \bx_n]^\top = \bV\bLambda^{1/2}=[\sqrt{\lambda_1}\bv_1, \ldots, \sqrt{\lambda_k}\bv_k].
$$

Then $\bx_i \in \mathbb{R}^k$, $i=1, \ldots, n$, are the principal coordinates of the $n$ points in $k$ dimensions.

It turns out that there is a very close connection between principal coordinate and principal components.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:mds0"><strong>(\#prp:mds0) </strong></span> Let $\bX$ be an $n \times p$ data matrix with associated Euclidean distance matrix
$$
d_{ij}^2 = (\bx_i -\bx_j)^\top(\bx_i -\bx_j),
$$
where $\bx_1^\top, \ldots , \bx_n^\top$ are the rows of $\bX$.  Then the centred PC scores based on the first $k$ principal components are  principal coordinates of the $n$ points in $k$ dimensions based on the distance matrix $\bD$.
\EndKnitrBlock{proposition}


## Similarity measures

Recap: so far in this chapter we have considered distances matrices $\bD=(d_{ij})_{i,j=1}^n$ with distances $d_{ij}$.  In this setting, the larger $d_{ij}$ is, the more distant, or dissimilar, object $i$ is from object $j$.

Recall that we have distinguished between metric distances (``metrics''), which satisfy the triangle inequality \@ref(eq:triangle), and non-metric distances, or dissimilarities, which need not satisfy \@ref(eq:triangle).

In this section, we now consider the analysis of measures of *similarity* as opposed to measures of dissimilarity.

A *similarity* matrix is defined to be an $n \times n$ matrix $\mathbf=(f_{ij})_{i,j=1}^n$ with the following properties:

1.  Symmetry, i.e. $f_{ij} =f_{ji}$, $i,j=1, \ldots , n$.
2. For all $i,j=1, \ldots , n$, $f_{ij} \leq f_{ii}$.

Note that when working with similarities $f_{ij}$, the larger $f_{ij}$ is, the more similar objects $i$ and $j$ are.

Condition 1. implies that object $i$ is as similar to object $j$ as object $j$ is to object $i$ (symmetry).

Condition 2. implies that an object is at least as similar to itself as it is to any other object.

One important class of problems is when the similarity between any two objects is measured by the number of common attributes.  We illustrate this through two examples.

\BeginKnitrBlock{example}
<span class="example" id="exm:unnamed-chunk-2"><strong>(\#exm:unnamed-chunk-2) </strong></span>Suppose there are 4 attributes we wish to consider.

1.  Attribute 1: Carnivore? If yes, put $a_1=1$; if no, put $a_1=0$.
2. Attribute 2: Mammal?  If yes, put $a_2=1$; if no, put $a_2=0$.
3. Attribute 3: Natural habitat in Africa?  If yes, put $a_3=1$; if no, put $a_3=0$.
4. Attribute 4: Can climb trees?  If yes, put $a_4=1$; if no, put $a_4=0$.

Consider a lion.  Each of the attributes is present so $a_1=a_2=a_3=a_4=1$.

A tiger?  In this case, 3 of the attributes are present (1, 2 and 4) but 3 is absent.
So for a tiger, $a_1=a_2=a_4=1$ and $a_3=0$.

How might we measure the similarity of lions and tigers based on the presence or absence of these four attributes?

First form a $2 \times 2$ table as follows.
$$
\begin{array}{cccc}
 &1 &0\\
1& a & b\\
0& c & d
\end{array}
$$
Here $a$ counts the number of attributes common to both  lion and tiger; $b$ counts the number of attributes the lion has but the tiger does not have; $c$ counts the number of attributes the tigher has that the lion does not have; and $d$ counts the number of attributes which neither the lion nor the tiger has.

In the above, $a=3$, $b=1$ and $c=d=0$.

How might we make use of the information in the $2 \times 2$ table to construct a measure of similarity?

The simplest measure of similarity is the proportion of the attributes which are shared.
$$
\frac{a}{a+b+c+d},
$$
which gives $0.75$ in this example.
A second similarity measure, which gives the same value in this example but not in general, is known as the *similarity matching coefficient* and is given by
\begin{equation}
\frac{a+d}{a+b+c+d}.
(\#eq:smc)
\end{equation}
There are many other possibilities, e.g. we could consider weighted versions of the above if we wish to weight different attributes differently.
\EndKnitrBlock{example}

\BeginKnitrBlock{example}
<span class="example" id="exm:mds1"><strong>(\#exm:mds1) </strong></span>Let us now consider a similar but more complex example with 6 unspecified attributes (not the same attributes as in Example 1) and 5 types of living creature, with the following data matrix, consisting of zeros and ones.
$$
\begin{array}{lcccccc}
&1&2&3&4&5&6\\
Lion&1&1&0&0&1&1\\
Giraffe&1&1&1&0&0&1\\
Cow&1&0&0&1&0&1\\
Sheep&1&0&0&1&0&1\\
Human&0&0&0&0&1&0
\end{array}
$$
Suppose we decide to use the similarity matching coefficient \@ref(eq:smc) to measure similarity.  Then the following similarity matrix is obtained.
$$
\bF=\begin{array}{lccccc}
&\text{Lion}&\text{Giraffe}&\text{Cow}&\text{Sheep}&\text{Human}\\
Lion&1&2/3&1/2&1/2&1/2\\
Giraffe&2/3&1&1/2&1/2&1/6\\
Cow&1/2&1/2&1&1&1/3\\
Sheep&1/2&1/2&1&1&1/3\\
Human&1/2&1/6&1/3&1/3&1
\end{array}
$$
It is easily checked from the definition that $\mathbf=(f_{ij})_{i,j=1}^5$ is a similarity matrix.
\vskip 0.2truein
We now return to the general case.  What should we do once we have calculated a similarity matrix?  It turns out there is a nice transformation from a similarity
matrix to a distance matrix $\bD=(d_{ij})_{i,j=1}^n$ defined by
\begin{equation}
d_{ij}=\left ( f_{ii}+f_{jj} -2f_{ij} \right )^{1/2}, \qquad i,j=1, \ldots , n.
(\#eq:defD)
\end{equation}
Note that, provided $\bF$ is a similarity matrix, the $d_{ij}$ are well-defined (i.e. real, not imaginary) because
$f_{ii}+f_{jj}-2f_{ij}\geq 0$ by condition 2., so the bracket is non-negative.
\EndKnitrBlock{example}

We have the following result.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:mds2"><strong>(\#prp:mds2) </strong></span>Suppose that $\bF$ is a similarity matrix.  If, in addition, $\bF$ is non-negative definite, then $\bD$ defined in \@ref(eq:defD) is Euclidean with centred inner product matrix
\begin{equation}
\bB=\bH\bF \bH,
(\#eq:BHFH)
\end{equation}
where $\bH=\bI_n - n^{-1}{\mathbf 1}_n {\mathbf 1}_n^\top$ is the centering matrix.
\EndKnitrBlock{proposition}

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}Since $\bF$ is non-negative definite by assumption, and $\bH^\top =\bH$ by definition of $\bH$, it follows that $\bH \bF \bH$ must also be non-negative definite.  So by Result 5.1, we just need to show that \@ref(eq:BHFH) holds, where $\bB$ is given by  $\bB=\bH\bA\bH$ and $\bA$ is defined as in \@ref(eq:defA), and the $d_{ij}$ are defined by \@ref(eq:defD).  Then
$$
a_{ij}=-\frac{1}{2}d_{ij}^2 =f_{ij}-\frac{1}{2}(f_{ii}+f_{jj}).
$$
Define
$$
t=n^{-1}\sum_{i=1}^n f_{ii}.
$$
Then, summing over $j=1, \ldots , n$ for fixed $i$,
$$
\bar{a}_{i+}=n^{-1}\sum_{j=1}^n a_{ij} = \bar{f}_{i+}-\frac{1}{2}(f_{ii}+t);
$$
similarly,
$$
\bar{a}_{+j}=n^{-1}\sum_{i=1}^n a_{ij}=\bar{f}_{+j}-\frac{1}{2}(f_{jj}+t),
$$
and also
$$
\bar{a}_{++}=n^{-2}\sum_{i,j=1}^n a_{ij}=\bar{f}_{++}-\frac{1}{2}(t+t).
$$
So, using part (vii) of section 7 of Chapter 2 (FIX FIX),
\begin{align*} 
b_{ij}&=a_{ij}-\bar{a}_{i+}-\bar{a}_{+j}+\bar{a}_{++}\\
&=f_{ij}-\frac{1}{2}(f_{ii}+f_{jj})-\bar{f}_{i+}+\frac{1}{2}(f_{ii}+t)\\
& \qquad -\bar{f}_{+j}+\frac{1}{2}(f_{jj}+t) +\bar{f}_{++}-t\\
& =\qquad f_{ij}-\bar{f}_{i+}-\bar{f}_{+j}+\bar{f}_{++}.
\end{align*}
Consequently, $\bB=\bH\bF\bH$, using part (vii) of \S 2.7 again, and the result is proved.
\hfill$\square$
\EndKnitrBlock{proof}
