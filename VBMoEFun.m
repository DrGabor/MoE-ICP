%%%%%%%%% fit a mixture model via Variational inference. 
function [out_model, prob_model, Z, L] = VBMoEFun(X, VBParams, IS_SHOW)
if isrow(X)
    X = X'; 
end
X0 = X;
%%%%%% sub-samplling for efficiency. 
if VBParams.reduceNum < length(X0)
    Idx = randperm(length(X0), VBParams.reduceNum);
    X = X0(Idx);
end
global prior; 
prior = VBParams; 
model = init(X, length(VBParams.SPV)); 
ModelArray = model; 
L = -inf; 
tol = VBParams.tol; 
nCnt = 0; 
for nIter = 1 : 1 : VBParams.maxIter
    model = ModelArray(end); 
    Z = expectation(X, model); 
    model = maximation(X, model, Z); 
    L(end+1) = ELBO(model, Z);
    if abs(L(end)-L(end-1)) <= tol*abs(L(end))
        nCnt = nCnt + 1;
    end
    if nCnt >= 3
        break;
    end
    ModelArray = [ModelArray model]; 
end
%%%% get the full Z and ELBO.
Z = expectation(X0, model);
L(end+1) = ELBO(model, Z);
model.Z  = Z; 
[prob_model, Z] = extractMode(model);
out_model = []; 
for id = 1 : 1 : length(prob_model.SPV)
    tmp = []; 
    tmp.para  = prob_model.mode(id);
    tmp.prior = prob_model.SPV(id);
    tmp.coeff = prob_model.weight(id); 
    out_model = [out_model tmp]; 
end
if IS_SHOW
    figure;
    hold on;
    grid on;
    plot(L, 'b.-');
    title(sprintf( 'Convergence Curve, Iter = %02d', length(L)) );
    if exist('outParams')
        Center = outParams.Center;
        Edges  = outParams.Edges;
        gtProb = outParams.gtProb;
    else
        [~, Edges] = histcounts(X, 100); % linspace(min(X), max(X), 50);
        Center = (Edges(1:1:end-1) + Edges(2:1:end))/2;
    end
    etProb = zeros(1, length(Center)); 
    SPV = prob_model.SPV;
    absEdges = abs(Edges); 
    for id = 1 : 1 : length(prob_model.weight)
        tao = prob_model.mode(id);
        s    = SPV(id);
        H = zeros(1, length(Edges)); 
        if VBParams.EP_Plus
            H = gammainc(tao * absEdges.^s, 1/s);
        else
            EffIdx = find(Edges >= 0.0);
            H(EffIdx) = 0.5 + 0.5*gammainc(tao * absEdges(EffIdx).^s, 1/s);
            NffIdx = find(Edges < 0);
            H(NffIdx) = 0.5 - 0.5*gammainc(tao * absEdges(NffIdx).^s, 1/s);
        end
        dVal = H(2:end) - H(1:end-1);
        etProb = etProb + prob_model.weight(id) * dVal;
    end
    figure;
    hold on;
    grid on;
    [val, ~] = histcounts(X, Edges);
    val = val / length(X);
    bar(Center,  val, 'r');
    hold on; 
    plot(Center, etProb, 'b.-', 'linewidth', 2);
    if exist('outParams')
        plot(Center, gtProb, 'g.-', 'linewidth', 2);
    end
    legend({'Empirical', 'Estimated', 'GT'}, 'box', 'off'); 
    if isfield(VBParams, 'K0')
        str = sprintf('VB-MoE, K_o = %02d, reduces to K = %02d', VBParams.K0, length(prob_model.weight));
    else
        str = sprintf('VB-MoE, K = %02d', length(prob_model.weight));
    end
    title(str);
    bTest = 1; 
end
end

function [prob_model, Z] = extractMode(model0)
EffIdx = find(model0.alpha >= 10);
%%%%%% if size(SPV) = 1, just keep it.  
if length(model0.alpha) == 1
    EffIdx = 1; 
end
Z = model0.Z(EffIdx, :); 
prob_model = [];
prob_model.alpha = model0.alpha(EffIdx);
prob_model.mu    = zeros(1, length(EffIdx));
prob_model.a     = model0.a(EffIdx);
prob_model.b     = model0.b(EffIdx);
prob_model.SPV   = model0.SPV(EffIdx); 
prob_model.NkArray    = model0.NkArray(EffIdx); 
prob_model.LamdaArray = model0.LamdaArray(EffIdx); 

prob_model.mode = [];
for id = 1 : 1 : length(prob_model.alpha)
    a = prob_model.a(id);
    b = prob_model.b(id);
    prob_model.mode(end+1) = a / b; % a / b
end
weight = prob_model.alpha / sum(prob_model.alpha);
if ~isrow(weight)
    weight = weight'; 
end
prob_model.weight = weight; 
prob_model = rmfield( prob_model, {'alpha', 'a', 'b'});
end

function model = init(X, K)
global prior; 
model = []; 
model.alpha  = repelem(prior.alpha, K)'; 
model.a      = repelem(prior.a,  K)'; 
model.b      = repelem(prior.b,  K)';
model.SPV    = prior.SPV; 
model.NkArray       = []; 
model.LamdaArray    = []; 
label = ceil(K*rand(1,length(X)));
% [label, ~] = kmeans(X, K);
N = length(X); 
Z = full(sparse(label, 1:N,1,K,N,N));
model  = maximation(X, model, Z);
end
%%%%%%%%
function Z = expectation(X, model)
if length(model.SPV) == 1
    Z = ones(1, length(X)); 
    return; 
end
K = length(model.alpha); 
Z = zeros(K, length(X));
expTao   = expLogGamFun(model.a, model.b); 
expAlpha = expLogDirFun(model.alpha);
SPV = model.SPV; 
logCArray = log(SPV)  - gammaln(1.0 ./ SPV) - log(2); 
absX = abs(X); 
for nK = 1 : 1 : K
    a     = model.a(nK);
    b     = model.b(nK); 
    logC  = logCArray(nK);
    spv   = SPV(nK); 
    val = logC  + expTao(nK)/spv - absX .^ spv * a/b + expAlpha(nK); 
    if ~isrow(val)
        val = val'; 
    end
    Z(nK, :) = val;  
end 
Z = Z - repmat(max(Z), K, 1); 
expZ = exp(Z); 
Z = expZ ./ sum(expZ); 
% for id = 1 : 1 : length(X)
%     val = Z(:, id); 
%     val = val - max(val); 
%     val = exp(val); 
%     val = val / sum(val); 
%     Z(:, id) = val; 
% end
end
%%%%%%%%
function [model]  = maximation(X, model0, Z)
global prior; 
model = compute_statistics(model0, X, Z); 
K = length(model.alpha); 
alpha0 = prior.alpha; 
a0     = prior.a; 
b0     = prior.b; 
NkArray    = model.NkArray; 
LamdaArray = model.LamdaArray; 
SPV        = model.SPV; 
for nK = 1 : 1 : K
    Nk    = NkArray(nK); 
    lamda  = LamdaArray(nK); 
    spv    = SPV(nK); 
    model.alpha(nK) = alpha0 + Nk; 
    model.a(nK)     = a0 + Nk/spv; 
    model.b(nK)     = b0 + lamda; 
end
model.Z = Z; 
end
%%%%%%%% Expectation on Gamma distribution. 
function expTao = expLogGamFun(A, B)
expTao = []; 
for id = 1 : 1 : length(A)
    a = A(id); 
    b = B(id); 
    expTao(end+1) = psi(a) - log(b); 
end
end
%%%%% Expectation on Direclete distribution. 
function expAlpha = expLogDirFun(Alpha)
expAlpha = []; 
SumAlpha = sum(Alpha); 
for id = 1 : 1 : length(Alpha)
    expAlpha(end+1) = psi(Alpha(id)) - psi(SumAlpha); 
end
end

function [model] = compute_statistics(model0, X, Z)
   model = model0; 
   K = size(Z, 1); 
   if ~isrow(X)
       X = X'; 
   end
   XX = repmat(X, K, 1);
   NkArray = sum(Z, 2);  
   LamdaArray = []; 
   SPV    = model.SPV; 
   absX = abs(X); 
   for nK = 1 : 1 : K 
       spv = SPV(nK);  
       LamdaArray(end+1)    = sum(Z(nK, :) .* absX .^ spv ); 
   end
   model.NkArray    = NkArray; 
   model.LamdaArray = LamdaArray; 
end
function L = ELBO(model, Z)
L = 0.0; 
global prior; 
K = length(model.a); 
expTaoArray   = expLogGamFun(model.a, model.b); 
expAlphaArray = expLogDirFun(model.alpha);
alpha0 = prior.alpha; 
a0     = prior.a; 
b0     = prior.b; 
AlphaArray = model.alpha;
AArray     = model.a;
BArray     = model.b;
SPV        = model.SPV; 
NkArray    = model.NkArray; 
LamdaArray = model.LamdaArray; 
logCArray = log(SPV)  - gammaln(1.0 ./ SPV) - log(2); 
logCAlpha0 = gammaln(K*alpha0) - K * gammaln(alpha0); 
logCAlpha = gammaln(sum(AlphaArray)) - sum(gammaln(AlphaArray)); 
L = logCAlpha0 - logCAlpha; 
logB0 = a0 * log(b0) - gammaln(a0); 
logBArray = AArray .* log(BArray) - gammaln(AArray); 
for nK = 1 : 1 : K 
    Nk       = NkArray(nK); 
    lamda    = LamdaArray(nK);
    expTao   = expTaoArray(nK); 
    expAlpha = expAlphaArray(nK); 
    alpha = AlphaArray(nK);  
    a     = AArray(nK); 
    b     = BArray(nK); 
    spv   = SPV(nK); 
    %%%%% expectation on likelihood term. 
    expLikehood = Nk * logCArray(nK) + Nk * expTao / spv - lamda * a / b; 
    %%%% expectation on latent variable.
    r = Z(nK, :);
    r = r( r >= eps);
    expZ = Nk * expAlpha - sum( r.* log(r) );
    %%%%% expectation on weight. 
    expWeight = (alpha0 - alpha) * expAlpha; 
    %%%%% expectation on Gaussian Gamma
    expPrecision = -(b0 - b)*a/b + ( a0 - a) * expTao + logB0 - logBArray(nK); 
    L = L + expLikehood + expZ + expWeight + expPrecision;  
end
end

