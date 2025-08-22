% ItsSchurCohn: checks Schur–Cohn stability and (optionally) plots eigenvalues
% Inputs:
%   A, B    - system matrices
%   K1      - state feedback gain
%   do_plot - if true, plots eigenvalues and unit circle
%   r       - (optional) radius for a secondary circle overlay
% Outputs:
%   none (error if open-loop with feedback is unstable)

function ItsSchurCohn(A, B, K1, do_plot, r)

    if nargin < 4
        do_plot = 0;
    end

    [n, ~] = size(B);
    angles = linspace(0, 2*pi, 100);
    unit_x = cos(angles);
    unit_y = sin(angles);

    eig_A   = eig(A);
    eig_Acl = eig(A + B*K1);

    stable = all(abs(eig_Acl) <= 1);
    if ~stable
        do_plot = 0;
        error("Closed-loop system without saturation is unstable");
    end

    if do_plot
        figure;
        for i = 1:n
            h1 = plot(real(eig_A(i)),   imag(eig_A(i)),   'r*', 'MarkerSize', 10); hold on;
            h2 = plot(real(eig_Acl(i)), imag(eig_Acl(i)), 'b.', 'MarkerSize', 10); hold on;
        end
        h3 = plot(unit_x, unit_y, 'Color', [0.8 0.8 0.8]);

        if nargin >= 5
            circ_x = r*cos(angles);
            circ_y = r*sin(angles);
            h4 = plot(circ_x, circ_y, 'm--'); hold on;
            text(-0.12, 0.1 + r, sprintf('r = %.2f', r), 'Color', 'm'); hold on;
            legend([h1 h2 h3 h4], {'$\lambda(A)$', '$\lambda(A+BK_1)$', 'Unit circle', ...
                   sprintf('$\\text{Circle }(r=%.2f)$', r)}, 'Interpreter','latex');
        else
            legend([h1 h2 h3], {'$\lambda(A)$', '$\lambda(A+BK_1)$', 'Unit circle'}, 'Interpreter','latex');
        end

        axis equal;
        xlabel('Real'); ylabel('Imaginary');
        title('Eigenvalues - Schur–Cohn');
    end
end

