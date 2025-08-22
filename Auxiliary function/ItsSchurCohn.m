function ItsSchurCohn(A, B, K1, do_plot,r)

    if nargin < 4
        do_plot = 0;
    end

    [n,~] = size(B);

    angles = linspace(0, 2*pi, 100);
    x = cos(angles);
    y = sin(angles);

    eigvalues1 = eig(A);
    eigvalues2 = eig(A + B*K1);

    stable = all(abs(eigvalues2) <= 1);

    if ~stable
        do_plot=0;
        error("Sistema sem saturação não é estavel");
    end

    if do_plot
        figure;

        for i = 1:n
            h1=plot(real(eigvalues1(i)), imag(eigvalues1(i)), 'r*', 'MarkerSize', 10); hold on;
            h2=plot(real(eigvalues2(i)), imag(eigvalues2(i)), 'b.', 'MarkerSize', 10); hold on;
        end

        h3=plot(x, y, 'Color', [0.8 0.8 0.8]);

        if ~(nargin<5)
            x1 = r*cos(angles);
            y1 = r*sin(angles);
            h4=plot(x1, y1, 'm--'); hold on;
            text(-0.12, 0.1+r, sprintf('r = %.2f', r), 'Color', 'm'); hold on;
            label=sprintf('$Circle (r = %.2f)$', r);
            legend([h1 h2 h3 h4], {'$\lambda (A)$', '$\lambda (A+BK_1)$', 'Unit circle',label}, 'Interpreter','latex');
        else
            legend([h1 h2 h3], {'$\lambda (A)$', '$\lambda (A+BK_1)$', 'Unit circle'}, 'Interpreter','latex');
        end

        axis equal;
        xlabel('Real'); ylabel('Imaginary');
        title("Eigvalues - Schur Cohn");
    end
end
