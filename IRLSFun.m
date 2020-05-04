function [ varargout ] = IRLSFun( params, Z, Model, maxIter )
N = size(Z, 2);
%%%%%%%%% optimize Aft + dX - Ref alternatively.
% W = ones(1, N);
% params.W = W / sum(W); 
Aft = params.Aft;
XArray = [];
rstTf = eye(4); 
dR = eye(3); 
dT = zeros(3, 1); 
W = ones(1, N);
for i = 1 : 1 : maxIter
    %% update W.
    params.Aft = Loc2Glo(params.Aft, dR', dT);
    params.V   = params.Aft - params.Ref;
    r = compute_residualFun(params);
    % [Model, ~, Z, ~]  = VBMoEFun(params.Err, params.VBParams, 0);
    r( r < 1e-10) = 1e-10;  % regularition for numerical stable.
    W = zeros(1, N);
    for k = 1 : 1 : size(Z, 1)
        p = Model(k).prior;
        W = W + Model(k).para * Z(k, :) .* r .^(p-2);
    end
    %%%%%%% remember to push W into the params.
    params.W = W/sum(W);
    %% update transformation.
    params.W = W/sum(W);
    if strcmpi(params.mode, 'point2point')
        [dR, dT] = reg_point2pointFun(params);
    end
    if strcmpi(params.mode, 'point2plane')
        [dR, dT] = reg_point2planeFun(params);
    end
    if strcmpi(params.mode, 'plane2plane')
        [dR, dT] = reg_plane2planeFun(params);
    end
    tmp = [];
    tmp.Tf = [dR dT];
    % XArray = [XArray tmp];
    rstTf = [dR dT; 0 0 0 1] * rstTf;
    if norm( dT ) <= 1e-6
        break;
    end
end
params.Err = r; 
if nargout == 1
    % varargout{1} =  XArray;
end
if nargout == 2
    varargout{1} = rstTf(1:3, 1:3);
    varargout{2} = rstTf(1:3, end);
end
if nargout == 3
    varargout{1} = rstTf(1:3, 1:3);
    varargout{2} = rstTf(1:3, end);
    if isrow(r)
        r = r'; 
    end
    varargout{3} = r; 
end
end