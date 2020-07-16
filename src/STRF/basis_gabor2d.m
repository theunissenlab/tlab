function basis = basis_gabor2d()

    basis = struct;
    basis.type = 'gabor2d';
    basis.ndim = 2;
    basis.func = @(x, y, x0, y0, params) gabor2d(x, y, x0, y0, [0 1 1 1 0]);
    basis.numParams = 0;
    