# Barkley-2D
Barkley model [1] is an example of diffusion-reaction equations, which have wide applications in chemistry, ecology and life sciences. The principal idea is that we model time change of the spatial distribution of two species $u(x,y)$ and $v(x,y)$, which can react with each other and are also diffused in the process. The model is given by

 $$ \frac{\partial u}{\partial t} = f(u,v) + \nabla^2 u   $$
 $$ \frac{\partial v}{\partial t} = g(u,v) $$

where $f(u,v)$ and $g(u,v)$ are the production (reaction) terms of the two chemical species. In the Barkley model, we choose

$$f(u,v) = \frac1{\epsilon}u(1-u)\left(u-\frac{v+b}{a}\right)$$

$$g(u,v) = u - v$$

where $\epsilon$, $a$, $b$ are constant parameters. For more discussion on the physical meaning of those parameters and the model see [1]. 

[1]  Dwight Barkley (2008), Scholarpedia, 3(11):1877, http://www.scholarpedia.org/article/Barkley_model

## Model set-up
We will solve the Barkley model on a square domain $(x,y) \in [-L,\ L] \times [-L,\ L]$ and use  the no-flux boundary conditions on all boundaries (i.e. the species cannot escape the square box). The no-flux condition means that on the left and right boundaries we have:

$$\frac{\partial u}{\partial x} = 0$$

and on the top and bottom boundaries we prescribe

$$\frac{\partial u}{\partial y} = 0$$

We will begin the simulation with the two species initially distributed only in the top ($u$) and only left ($v$) half of the domain:

$$ u(x,y) = \begin{cases}  1& \text{if } y 
\geq 0 \newline 0 & \text{if } y<0 \end{cases} $$

$$v(x,y) = \begin{cases}  1& \text{if } x 
\geq 0 \newline 0 & \text{if } x<0 \end{cases}$$

By default, we will use the parameter values $\epsilon = 0.02$, $a = 0.75$, $b=0.01$ and run the simulation until the final time of 40 s. THe result you should see is quite unexpected, as for this particular choice of parameters, the $u$ species concentration forms a spiral

<img width="366" alt="image" src="https://github.com/user-attachments/assets/631f17a2-9e59-4348-9647-f416fc48b9ad" />

## Numerical discretization

We solve the system of equations by first dividing the domain into a set of $N \times N$ elements

<img width="832" alt="image" src="https://github.com/user-attachments/assets/b330a24a-294b-45f7-a33f-6380e2429045" />

and locate the solution variables $u_{i,j}$, $v_{ij}$ inside each of the elements. We then discretize the diffusion operator

$$\nabla^2 u_{i,j} \approx \frac1{\Delta x^2}\left( u_{i-1,j} + u_{i+1,j} + u_{i,j-1} + u_{i,j+1} - 4u_{i,j} \right)$$

where $\Delta x$ is the size of each element.

We use forward Euler explicit time integration (not a good idea in general, but will do for our project), which leads to the following:

$$u_{i,j}^{n+1} = u_{i,j}^{n} + \Delta t \left(f(u_{i,j}^n,\ v_{i,j}^n) + \frac{u_{i-1,j}^n + u_{i+1,j}^n + u_{i,j-1}^n + u_{i,j+1}^n - 4u_{i,j}^n}{\Delta x^2}   \right)$$

$$v_{i,j}^{n+1} = v_{i,j}^{n} + \Delta t\ g(u_{i,j}^n,\ v_{i,j}^n) $$

where superscript $^n$ represents the time at which the solution is evaluated. 





