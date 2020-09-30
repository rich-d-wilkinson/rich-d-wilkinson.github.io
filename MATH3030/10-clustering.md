
# Cluster Analysis {#cluster}


A key difference between Discriminant Analysis, covered in \S 9, and Cluster Analysis is that in the former, a training sample is available, whereas in the latter, we do not have access to a training sample.  In Machine Learning and Data Mining, Cluster Analysis is typically referred to as Unsupervised Learning (''Unsupervised'' here refers to the fact that no training sample is available).

In this chapter, we shall limit ourselves to two topics within Cluster Analysis:

-  likelihood-based clustering which, as we will see, includes the widely-used method $K$-means clustering as a special case; and
- hierarchical clustering methods.


## Likelihood-based clustering 

Suppose we have a sample of $n$ random vectors $\bx_1, \ldots , \bx_n$, assumed independent.  Suppose that each $\bx_i$ comes from one of $g$ sub-populations, where the $j$th sub-population has
probability density function $f_j(\bx; \btheta_j)$, $j=1,\ldots , g$. For simplicity it is assumed in this section  that $g$ is given.  However, in many applications, $g$ is unknown, and we do not assume that $g$ is known later in this chapter.  The key point is that we do not know which sub-population each $\bx_i$ comes from, i.e. we do not know how to allocate the observations to sub-populations.

**The goal of cluster analysis**:  to allocate each of $n$ observational vectors to one of $g$ clusters, in such a way that observation vectors within a cluster tend to be more similar, in a suitable sense, than observations in different clusters.

In principle we can estimate the optimal allocation, along with the (usually) unknown parameter vectors $\btheta_1, \ldots , \btheta_g$, by maximum likelihood, as we shall see.

It will be convenient to introduce two equivalent ways to describe an arbitrary allocation of each $\bx_i$ to one of the $g$ clusters.  Write $\bdelta=(\delta_1,\ldots ,  \delta_n)^\top$.  Then consider the equivalence
\begin{equation}
\delta_i = j \iff \bx_i \in \mathcal{C}_j, \qquad j=1,\ldots , g,
(\#eq:two-way)
\end{equation}
where
$$
\bigcup_{j=1}^g \mathcal{C}_j = \{\bx_1, \ldots , \bx_n\}, \qquad \text{and} \qquad \mathcal{C}_j \cap \mathcal{C}_{k} =\emptyset,\qquad j \neq k.
$$
Then the likelihood for $\btheta_1, \ldots , \btheta_g$ and the allocation $\bdelta$  may be written
$$
L(\btheta_1,\ldots , \btheta_g; \bdelta)=\left \{\prod_{\bx \in \mathcal{C}_1} f(\bx; \btheta_1)\right \}\cdots \left \{\prod_{\bx \in \mathcal{C}_g} f(\bx; \btheta_g)\right \}.
$$
The sets $\mathcal{C}_1, \ldots , \mathcal{C}_g$ are called **clusters**.

Let $\hat{\btheta}_1, \ldots , \hat{\btheta}_g$ and $\hat{\bdelta}$ denote the maximum likelihood estimators of $\btheta_1, \ldots , \btheta_g$ and the unknown allocation $\bdelta$.
Also, let $\hat{\mathcal{C}}_1, \ldots, \hat{\mathcal{C}}_g$ denote the maximum likelihood clusters, which are determined by $\hat{\bdelta}$ via \@ref(eq:two-way).

Then we have the following result.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:ten1"><strong>(\#prp:ten1) </strong></span>If $\bx \in \hat{\mathcal{C}}_j$ then
$$
\sup_{1 \leq k \leq g}f(\bx; \hat{\btheta}_k) = f(\bx; \hat{\btheta}_j).
$$
\EndKnitrBlock{proposition}

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Proof. </em></span>  \fi{}If we move an $\bx$ from $\hat{\mathcal{C}}_j$ to $\hat{\mathcal{C}}_k$, then the likelihood changes to
$$
L(\hat{\btheta}_1, \ldots , \hat{\btheta}_g; \bdelta)f(\bx; \hat{\btheta}_k)/f(\bx; \hat{\btheta}_j).
$$
But by definition of $L(\hat{\btheta}_1, \ldots , \hat{\btheta}_g; \hat{\bdelta})$,
$$
L(\hat{\btheta}_1, \ldots , \hat{\btheta}_g; \hat{\bdelta})f(\bx; \hat{\btheta}_k)/f(\bx; \hat{\btheta}_j) \leq L(\hat{\btheta}_1, \ldots , \hat{\btheta}_g; \hat{\bdelta}),
$$
from which we conclude that Proposition \@ref(prp:ten1) holds.
\EndKnitrBlock{proof}

Note that this result is closely related to the sample ML discriminant rule considered in \S 9.3.

We now consider the case where the sub-populations are multivariate Gaussian, i.e. $f(\bx; \btheta_j)$ is the density of $N_p(\bmu_j, \bSigma_j)$ for $j=1, \ldots , g$.


In the general case, when the mean vector and covariance matrix are different for each sub-population, we know how to maximise the likelihood when the allocation $\bdelta$ is given.
Here,  $\btheta_j$ consists of $\bmu_j$ and $\bSigma_j$ for each $j=1, \ldots , g$, and from \S 8.4, we know that the maximised log-likelihood for cluster $j$ is
$$
\ell(\hat{\bmu}_j[\bdelta], \hat{\bSigma}_j[ \bdelta])=-\frac{n_j[\bdelta]}{2}\log (\vert \hat{\bSigma}_j[\bdelta]\vert)-\frac{n_j[\bdelta]}{2}p(1+\log 2\pi),
$$
where $n_j[\bdelta]$ is the number of elements in cluster $j$, $\mathcal{C}_j[\bdelta]$, for the given allocation $\bdelta$.  Similarly, $\hat{\bmu}_j[\bdelta]$ and $\hat{\bSigma}_j[\bdelta]$ are the MLEs of $\bmu_j$ and $\bSigma_j$ for the given allocation $\bdelta$, i.e. the sample mean and covariance matrix
$$
\hat{\bmu}_j[\bdelta]=\frac{1}{n_j[\bdelta]}\sum_{\bx \in \mathcal{C}_j[\bdelta]} \bx =\bar{\bx}_j[\bdelta]
$$
and
$$
\hat{\bSigma}_j[\bdelta]=\frac{1}{n_j[\bdelta]}\sum_{\bx \in \mathcal{C}_j[\bdelta]}(\bx - \bar{\bx}_j[\bdelta])(\bx - \bar{\bx}_j[\bdelta])^\top.
$$

It follows that the MLE of $\bdelta$ is the choice of $\bdelta$ which maximises the log-likelihood
$$
-\frac{1}{2}\sum_{j=1}^g n_j[\bdelta] \log (\vert\hat{\bSigma}_j[\bdelta]\vert)\, + \, \text{constant}
$$
over $\bdelta$.

Under the assumption that the population covariance matrices are the same, i.e. $\bSigma_1=\cdots = \bSigma_g$, the maximised log-likelihood for a given allocation $\bdelta$ is
given by
\begin{equation}
-\frac{n}{2}\log (\vert \bW[\bdelta]\vert)\, + \,\text{constant},
(\#eq:Wdelta)
\end{equation}
where $n=\sum_{j=1}^g n_j$ and
$$
\bW[\bdelta]=\sum_{j=1}^g \sum_{x \in \mathcal{C}_j[{\pmb \delta}]} (\bx - \hat{\bmu}_j[\bdelta])(\bx - \hat{\bmu}_j[\bdelta])^\top
$$
is the ''within'' sum of squares and products matrix for the given allocation.  So the maximum likelihood allocation $\hat{\bdelta}$ is the $\bdelta$ which maximises \@ref(eq:Wdelta).

The final case we consider is where $\bSigma_1=\cdots = \bSigma_g=\sigma^2 \bI_p$, i.e. the common covariance matrix is a scalar multiple of the $p \times p$ identity matrix.
This is a version of the so-called \text{$k$-means clustering} approach, and here the maximum likelihood allocation $\hat{\bdelta}$ is obtained by the $\bdelta$ which minimises
$$
\sum_{j=1}^g \sum_{\bx \in \mathcal{C}_j} \vert \vert \bx - \bar{\bx}_j\vert \vert^2.
$$

Although clustering based on the likelihood function has a certain intuitive appeal, in most practical situation it is not feasible to find the global maximum due to the
computational explosion in the number of possible allocations when $n$ is even of moderate size, e.g. $n=100$.  A further problem is that there may be a large number of local maxima of the likelihood function.  However, despite these challenges, likelihood-based clustering, such as  $k$-means clustering, is widely used and can lead to useful, even if sub-optimal, solutions to the clustering problem.

## Hierarchical clustering methods

 The adjective *hierarchical* applies to clustering methods which have the following property: the arrangement of the experimental units into $g$ clusters and into $g+1$ clusters have the properties that $g-1$ of the clusters are identical, and the remaining single cluster in the $g$ clusters is split into $2$ clusters in the $g+1$ clusters.

 A hierarchical clustering method is usually of one of two types:
 
1. an *agglomerative* clustering method progressively combines clusters, usually starting with $n$ singleton clusters; and

2. a *divisive* clustering method progressively splits, or divides, clusters, usually starting with a single cluster containing $n$ elements.

 Many clustering methods a based on an $n \times n$ matrix of inter-point distances $\bD=(d_{ij})_{i,j=1}^n$ of the type we considered in \S 5.  However, the goal here is somewhat different to that in Multidimensional Scaling.

 Given $\bD$, two common clustering procedures are:

1. the *single linkage* method, sometimes called the  *nearest neighbour* method; and

2. the *complete linkage* method, sometimes called the *furthest neigbour* method.

 We now explain in more detail how these methods are applied.  First, we consider single linkage.  It is assumed that the set of distances is ordered, so that
 \begin{equation}
 d_{a_1,b_1}\leq d_{a_2,b_2} \leq \cdots \leq d_{a_{N-1}, b_{N-1}} \leq d_{a_N,b_N},
(\#eq:orderdist)
 \end{equation}
 where $N=n(n-1)/2$, and  we adopt the following  conventions: (i) $a_t < b_t$; and (ii) to break a tie such as  $d_{a_{t},b_{t}} = d_{a_{t+1},b_{t+1}}$,  we write
 $$
 \cdots \leq d_{a_{t},b_{t}} \leq d_{a_{t+1},b_{t+1}} \leq \cdots
 $$
  if $a_{t} < a_{t+1}$, or $a_{t}=a_{t+1}$ and $b_t<b_{t+1}$; and we write
 $$
 \cdots \leq d_{a_{t+1},b_{t+1}} \leq d_{a_{t},b_{t}} \leq \cdots
 $$
 if $a_{t+1}<a_t$, or $a_{t+1}=a_t$ and $b_{t+1}<b_t$.


1. Start with singleton clusters $\mathcal{C}_1, \ldots , \mathcal{C}_n$, i.e. $\mathcal{C}_i =\{i\}$.

2. Combine experimental units $a_1$ and $b_1$ into a single new cluster, so that we now have one doubleton cluster and $n-2$ singleton clusters.

3. Next, consider $a_2$ and $b_2$, where $d_{a_2,b_2}$ is the second smaller distance.  If both $a_2 \notin \{a_1,b_1\}$ and $b_2 \notin \{a_1,b_1\}$, then we combine $a_2$ and $b_2$ into a second doubleton cluster,
     $\{a_2, b_2\}$, leading to $n-2$ clusters altogether ($2$ doubleton clusters and $n-4$ singleton clusters).  If, on the other hand, $a_2 \in \{a_1,b_1\}$, then necessarily $b_2 \notin \{a_1,b_1\}$, and so we form the trippleton cluster $\{a_1,b_1,b_2\}$; while if $b_2 \in \{a_1, b_1\}$, then necessarily $a_2 \notin \{a_1, b_1\}$, and we form the trippleton cluster $\{a_1, a_2, b_1\}$.  Either way, in the latter two cases, we are left with $n-2$ clusters altogether (one trippleton and $n-3$ singleton clusters).

4. We continue this process as we pass through the $N$ inter-point distances.  However, sometimes we may wish to terminate this process at some threshold, $T$; i.e. we stop the process at the smallest $t$ such that
     $d_{a_t, b_t}>T$.


 In the single linkage approach, in effect we are defining the distance between two clusters $\mathcal{C}_u$ and $\mathcal{C}_v$ by
 $$
 d_S^\ast(\mathcal{C}_u, \mathcal{C}_v)=\min_{i \in \mathcal{C}_u, \, j \in \mathcal{C}_v} d_{ij}.
 $$

In contrast, with the complete linkage method, in effect we are defining  the distance between two clusters $\mathcal{C}_u$ and $\mathcal{C}_v$ by
$$
d_C^\ast(\mathcal{C}_u, \mathcal{C}_v)=\max_{i \in \mathcal{C}_u, \, j \in \mathcal{C}_v} d_{ij}.
$$

A convenient graphical way to present the output from either a single linkage procedure or a complete linkage procedure is to plot a **dendrogram**.

\BeginKnitrBlock{example}
<span class="example" id="exm:exten1"><strong>(\#exm:exten1) </strong></span>The following distance matrix was based on the relative gene frequencies for the four blood-group systems $A_1$, $A_2$, $B$ and $O$ for large samples from four human populations: (1) Inuit, (2) African, (3) English and (4) Korean.  The inter-point distances were determined using the Mahalanobis distance between $\bmu_i$ and $\bmu_j$ defined by
$$
d_{ij}=\sqrt{(\bmu_i -\bmu_j)^\top \bSigma^{-1}(\bmu_i - \bmu_j)}.
$$
The matrix of inter-point distances is given by:
$$
\begin{pmatrix}
0&23.26&16.34 & 16.87\\
&0&9.85 & 20.43 \\
&&0&19.60\\
&&&0
\end{pmatrix}.
$$

**Hypothesis of interest**:  *there there is a natural clustering*
$$
\{\text{African}=2, \text{English}=3\}\qquad \text{and}\qquad
\{\text{Inuit}=1, \text{Korean}=4\}.
$$

Let us first of all look at single linkage.  The ordering of the distances is given by
$$
d_{23} < d_{13} <d_{14}<d_{34}<d_{24}<d_{12}.
$$

*Single Linkage*
  
- At Stage $0$ the clusters are $\{1\}$, $\{2\}$, $\{3\}$ and $\{4\}$.
- At Stage $1$ the clusters are $\{1\}$, $\{2,3\}$ and  $\{4\}$.
- At Stage $2$ the clusters are $\{1,2,3\}$ and $\{4\}$.
- At Stage $3$ we have a single cluster $\{1,2,3,4\}$.

Now let us look at complete linkage.  Here, Stage $0^\prime$ and Stage $1^\prime$ are the same at Stage $0$ and Stage $1$ of single linkage, but Stage $2^\prime$ is different.

*Complete Linkage*
- At Stage $0^\prime$ the clusters are $\{1\}$, $\{2\}$, $\{3\}$ and $\{4\}$.
- At Stage $1^\prime$ the clusters are $\{1\}$, $\{2,3\}$ and  $\{4\}$.
- At Stage $2^\prime$ the clusters are $\{1,4\}$ and $\{2,3\}$.
- At Stage $3^\prime$ we have a single cluster $\{1,2,3,4\}$.

The reason that we combine $\{1,4\}$ at Stage $2^\prime$ is because
$$
d_{1,4}=16.87<d^\ast_C(\{1\}, \{2,3\})=\max \{d_{12}, d_{13}\}=23.26
$$
and
$$
d_{1,4}=16.87<d^\ast_C(\{4\}, \{2,3\})=\max \{d_{24}, d_{34}\}=20.43.
$$
So with complete linkage we should combine $\{1\}$ and $\{4\}$ before combining either $\{1\}$ or $\{4\}$ with $\{2,3\}$.

In this example, single linkage and complete linkage lead to different conclusions: complete linkage supports the hypothesis of interest, while single linkage does not.  Some people have argued that single linkage is more appropriate here.


FIX THESE PLOTS
\vfill\eject
\centerline{Sketch of Dendrogram: Single Linkage}

\vskip 3.5truein

\centerline{Sketch of Dendrogram: Complete Linkage}

\vfill\eject
\EndKnitrBlock{example}

## Further Points



In this chapter we have focused on just two approaches to Cluster Analysis: likelihood-based methods, especially those based on the multivariate Gaussian model; and hierarchical clustering methods, especially single linkage and complete linkage approaches.

There are many algorithms for performing Cluster Analysis, but there is no generally accepted ''best'' method.  Moreover, different algorithms (or even the same algorithm with a different initialisation) do not necessarily produce the same results on a given dataset, and there is often a fairly large subjective element in the assessment of any particular method.

One way to test a clustering algorithm is to apply it on data with a known group structure.  Experience suggests that this will only produce good results when the groups are very distinct. When, on the other hand,  there is a lot of overlap between groups, clustering algorithms are not likely to perform particularly well.

However, despite these cautionary remarks, clustering algorithms are often useful in practice, but it is an area where usually the most one can hope for is to find a good, but sub-optimal, solution.

