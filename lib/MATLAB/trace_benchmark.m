%% Benchmarking trace algorithms

% the QETLab PartialTrace is pretty trash, and very inconsistent.
%  TrX really is quick, even for big things. Maybe I can win by
   % vectorizing?


% Case 1: Product of qudits
% Case 2: Uniform products of non-uniform systems
% Case 3: Entangled states of qudits
% Case 4: Entangled states of non-uniform systems

clear all
close all

%% Case 1
% Fix test parameters
Nmax = 11;
TrX_times = zeros(Nmax);
QET_times = zeros(Nmax);
ptr_times = zeros(Nmax);
profile on
for N = 2:Nmax
    N
    for M = 1:N-1

    %Generate and Tensor together N random qubits
    sys_cell = cell(N,1);
    for n = 1:N
        sys_cell{n} = rand_qubit;
    end
    rho = Tensor(sys_cell);
    dim = 2*ones(N,1);

    % Randomly identify M subsystems to trace out
    sys = sort(rand_sel(M,N,[]));

    % Consistency test
%     logic_label = {'FALSE' 'TRUE'};
%     fprintf('|TrX-QET|    %u\n', sum(sum(abs(PartialTrace(rho,sys,dim)-TrX(rho,sys,dim)))));
%     fprintf('|TrX-ptr|    %u\n', sum(sum(abs(ptrace(rho,sys,dim)-TrX(rho,sys,dim)))));


    TrX_test = @() TrX(rho,sys,dim);
%     QET_test = @() PartialTrace(rho, sys,dim);
    ptr_test = @() ptrace(rho,sys,dim);
    % Test the functions
    TrX_times(N,M) = timeit(TrX_test);
%     QET_times(N,M) = timeit(QET_test);
    ptr_times(N,M) = timeit(ptr_test);
    % Print outputs for now
%     fprintf('==Results==\n')
%     fprintf('TrX        %f\n', TrX_times(N,M))
%     fprintf('QET        %f\n', QET_times(N,M))
%     fprintf('ptr        %f\n', ptr_times(N,M))
    end
end
profile off
profile viewer
figure()
plot(TrX_times,'b.')
hold on
% plot(QET_times,'ro')
plot(ptr_times,'kx')
xlabel('N qubits')
ylabel('Wall time varying subsys')
title('Toby (blue) v Jacob(black)')
set(gca,'Yscale','log')

%% TrX: 
%    RHO = TrX(PSI,SYS,DIM) traces out the subsystems specified in
%    vector SYS of state PSI (a state vector or densitry matrix) whose
%    subsystem dimensions are specified by the vector DIM.

%


% PartialTrace
%   XPT = PartialTrace(X) is the partial trace of the matrix X,
%   where it is assumed that length(X) is a perfect squares and both
%   subsystems have equal dimension. The trace is taken over the second
%   subsystem.
%
%   This function has three optional arguments:
%     SYS (default 2)
%     DIM (default has all subsystems of equal dimension)
%     MODE (default -1)
%

% ptrace
% Accepts one density matrix or state vector,  an array of dimensions of the subsystems, 
%       and an array of subsystems to trace out.

function sys = rand_sel(m,N,sys)
    if length(sys) >= m
        sys = sys(1:m);
    elseif length(sys) >= N
        sys = sys(1:m);
    else
        i = randi([1,N]);
        if any(sys ==i)
            sys = rand_sel(m,N,sys);
        else
            sys = [sys, i];
            sys = rand_sel(m,N,sys);
        end
    end
end
        