% Plot results
   % 2d plot/colour 
% fit asymptotic complexity
% Profile findings




%The max number of qubits in the system
figure(1)
clf;
Nmax = 11;
profile on
for N=2:Nmax
    dim = 2*ones(1,N); % all qubits for now
    for num_sys = 1:N-1
        sys = randperm(N);
        sys = sys(1:num_sys);

        % Generate a big state with real coefficients (extend)
        psi = randn(2^N,1);
        rho = toDM(psi); %toDM normalizes input

        TrX_call = @() TrX(rho,sys,dim);
        TrX_time = timeit(TrX_call);
        ptrace_call = @() ptrace(rho,sys,dim);
        ptrace_time = timeit(ptrace_call);
        
        
        subplot(1,2,1)
        plot(N,TrX_time,'kx')
        hold on
        plot(N,ptrace_time,'rx')
        subplot(1,2,2)
        plot(N,TrX_time/ptrace_time,'bx')
        hold on
    end
end

profile off
profile viewer
hold off

%%
subplot(1,2,1)
set(gca,'Yscale','log')
legend('TrX','ptrace')
xlabel('N, with subsystems up to N-1')
ylabel('Walltime (s)')
title('TrX vs ptrace: MATLAB')

subplot(1,2,2)
set(gca,'Yscale','log')
xlabel('N, various subsys sizes')
ylabel('sec/sec')
title('t(TrX)/t(ptrace)')

saveas(gcf,'trace_time_matlab.png')
saveas(gcf,'trace_time_matlab.fig')
