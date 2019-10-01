# [Revisiting Frank-Wolfe: Projection-Free Sparse Convex Optimization](http://proceedings.mlr.press/v28/jaggi13-supp.pdf)

This is a modern treatment of the conditional gradient or the Frank-Wolfe method, which was originally described in 1956. General constrained convex optimisation problems of the form
$$
\text{min}_{x \in \mathcal{D}} f(x),
$$
are considered. Here, $f$ is convex and continuously differentiable, and the domain $\mathcal{D}$ is a compact convex subset of any vector space.

The conditional gradient method amounts to linearising the objective function at a position $x$ and moving towards a minimiser of the linear function (over the same domain $\mathcal{D}$). At each step of the algorithm, at most one point is added to the solution, thus making it possible to obtain sparse solutions












































<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/MathJax.js?config=TeX-MML-AM_CHTML">
</script>
