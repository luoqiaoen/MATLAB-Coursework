function scores = rastrigins(pop)
scores = 10.0 * size(pop,2) + sum((pop./10) .^2 - 10.0*cos(2*pi.*(pop./10)),2);
  


