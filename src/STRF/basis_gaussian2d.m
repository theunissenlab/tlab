function basis = basis_gaussian2d()
    
    basis = struct;
    basis.type = 'gaussian2d';
    basis.ndim = 2;
    basis.func = @(x, y, x0, y0, params) gaussian2d(x, y, x0, y0, [0.8 0.8]);
    basis.numParams = 0;
    