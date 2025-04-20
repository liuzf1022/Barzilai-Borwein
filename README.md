This project compares the Barzilai-Borwein (BB) method and the Steepest Descent method for solving unconstrained optimization problems, particularly in the context of least squares.

Project Includes:

1. Barzilai-Borwein method with adaptive step size

2. Classical Steepest Descent method

3. Comparison of convergence behavior (iteration paths, gradient norms, and objective values)

4. 2D contour plots and semilog plots for visualization

Background：
The Barzilai-Borwein (BB) method is a gradient-based optimization algorithm proposed in 1988. It accelerates convergence by smartly choosing step sizes based on previous iterations, inspired by quasi-Newton techniques. It is especially effective for large-scale or ill-conditioned problems.
In contrast, the Steepest Descent method uses exact line search but often suffers from slow convergence, especially in narrow valleys.
