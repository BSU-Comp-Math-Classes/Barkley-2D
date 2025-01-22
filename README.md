# Barkley-2D
Barkley model [1] is an example of diffusion-reaction equations, which have wide applications in chemistry, ecology and life sciences. The principal idea is that we model time change of the spatial distribution of two species $u(x,y)$ and $v(x,y)$, which can react with each other and are also diffused in the process. The model is given by

 $$ \frac{\partial u}{\partial t} = f(u,v) + \nabla^2 u   $$
 $$ \frac{\partial v}{\partial t} = g(u,v) $$

where $f(u,v)$ and $g(u,v)$ are the production (reaction) terms of the two chemical species. In the Barkley model, we choose

$$f(u,v) = \frac1{\epsilon}u(1-u)\left(u-\frac{v+b}{a}\right)$$
$$g(u,v) = u - v$$,

where $\epsilon$, $a$, $b$ are constant parameters. For more discussion on the physical meaning of those parameters and the model see [1]. 

[1]  Dwight Barkley (2008), Scholarpedia, 3(11):1877, http://www.scholarpedia.org/article/Barkley_model

## Model set-up
We will solve the Barkley model on a square domain $(x,y) \in [-L,\ L] \times [-L,\ L]$ and use  the no-flux boundary conditions on all boundaries (i.e. the species cannot escape the square box). The no-flux condition means that on the left and right boundaries we have:

$$ u(x,y) = \begin{cases}  1& \text{if } y 
\geq 0 \\ 0 & \text{if } y<0 \end{cases} $$

