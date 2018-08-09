function f = dottheta(omega, lambda, r, psi, theta)
f = omega + lambda * r * sin(psi - theta);
end