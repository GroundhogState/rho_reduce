% fprintf('======Testing with product state=======')
% 
Nmax = 14;
rec_times = zeros(1,Nmax);
full_times = zeros(1,Nmax);
tic
profile on
for num_qubits = 2:Nmax
    fprintf('num_qubits = %u\n',num_qubits)
    Q = cell(1,num_qubits);
    for ii=1:num_qubits
        Q{ii} = rand_qubit();
    end
    rho = Tensor(Q);
    rec_test = @() rho_reduce(rho,2*ones(1,length(Q)),length(Q));
    full_test = @() rho_reduce_bruteforce(rho);

    rec_times(num_qubits) = timeit(rec_test);
    full_times(num_qubits) = timeit(full_test);
end
toc
profile off
profile viewer

%%
figure(1)
clf
subplot(1,2,1)
plot(1:Nmax, rec_times,'bx')
hold on
plot(1:Nmax,full_times,'ko')
set(gca,'Yscale','log')
legend('Recursive','Full','location','southeast')
xlabel('num qubits')
ylabel('Walltime [sec]')
title('Recursive vs naive 2-body trace')
hold off
subplot(1,2,2)
plot(full_times(1:13)./rec_times(1:13),'kx')
title('Relative speedup')
xlabel('num qubits')



%% Below: Correctness checks (be more rigorous...)
% num_qubits = 2
% Q = cell(1,num_qubits);
% for ii=1:num_qubits
%     Q{ii} = rand_qubit();
% end
% rho = Tensor(Q);
% delta = rho_reduce(rho,2*ones(1,length(Q)),length(Q));
% for ii=1:length(Q)
%     fprintf('OBDM %u error %f\n',ii,norm(delta{ii,ii} - Q{ii}))
% end
% for ii = 1:length(Q)
%     for jj = ii+1:length(Q)
%         fprintf('2BDM (%u,%u) error %f\n',ii,jj,norm(delta{ii,jj} - Tensor(Q{ii},Q{jj})))
%     end
% end
% % 
% % % Todo: test with entangled states
% % 
% 
% 
