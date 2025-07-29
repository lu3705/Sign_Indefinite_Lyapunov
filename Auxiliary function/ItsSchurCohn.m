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
            plot(real(eigvalues1(i)), imag(eigvalues1(i)), 'r*', 'MarkerSize', 10); hold on;
            plot(real(eigvalues2(i)), imag(eigvalues2(i)), 'b.', 'MarkerSize', 10); hold on;
        end

        plot(x, y, 'Color', [0.8 0.8 0.8]);

        if ~(nargin<5)
            x1 = r*cos(angles);
            y1 = r*sin(angles);
            plot(x1, y1, 'm--'); hold on;
            text(-0.12, 0.1+r, sprintf('r = %.2f', r), 'Color', 'm'); hold on;
            legend('eig(A)', 'eig(A+BK1)', 'unit circle');
        else
            legend('eig(A)', 'eig(A+BK1)');
        end

        axis equal;
        xlabel('Real'); ylabel('Imaginary');
        title("Eigvalues - Schur Cohn");
    end
end
