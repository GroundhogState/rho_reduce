L = 5;
vec = rand(2^L,1);
rho = vec'.*vec;
rho = rho/trace(rho);
A=rho_to_graph(rho)
    