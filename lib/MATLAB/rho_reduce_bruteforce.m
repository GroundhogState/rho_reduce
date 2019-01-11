function A = rho_reduce_bruteforce(rho)
% INPUT: One density matrix. 
%OUTPUT: One Aleph. Diagonal is 2x VN entropy, off-diag is MI.
% SUPERSEDES vec_to_graph
        
        rho = rho/trace(rho);
        L = log2(length(rho));
        A = cell(L,L);
        ent_list = zeros(L,1);
        dims = 2*ones(1,L);
        systems = 1:L;

        for ii=1:L 
            for jj = ii+1:L
               trace_pair = systems;
               trace_pair(ii) =[];
               trace_pair(jj-1) = [];
               rho_red = TrX(rho,trace_pair,dims);
               A{ii,jj} =  rho_red;
            end
            if ii==L
                A{ii,ii} = TrX(rho_red,1,[2,2]);
            else
                A{ii,ii} = TrX(rho_red,2,[2,2]);
            end
        end

end