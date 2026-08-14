function [x, fval, exitflag, output] = simplexSolve(c, A, b, varargin)
%SIMPLEXSOLVE Solve a linear program with the two-phase revised simplex method.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = SIMPLEXSOLVE(C,A,B) solves
%
%       maximize    C' * X
%       subject to  A * X <= B
%                   X >= 0.
%
%   Name-value options:
%       'ObjectiveSense'  - 'max' (default) or 'min'
%       'ConstraintSense' - M-by-1 array containing '<=', '>=', or '='
%       'FreeVariables'   - logical N-vector or indices of unrestricted X
%       'Tolerance'       - numerical tolerance (default 1e-9)
%       'MaxIterations'   - limit for each phase (default 1000)
%
%   EXITFLAG is 1 for an optimum, 0 for the iteration limit, -2 for an
%   infeasible problem, -3 for an unbounded problem, and -4 for a
%   numerical failure.

validateattributes(c, {'numeric'}, {'vector', 'real', 'finite'}, mfilename, 'c', 1);
validateattributes(A, {'numeric'}, {'2d', 'real', 'finite'}, mfilename, 'A', 2);
validateattributes(b, {'numeric'}, {'vector', 'real', 'finite'}, mfilename, 'b', 3);

c = c(:);
b = b(:);
[constraintCount, variableCount] = size(A);

if numel(c) ~= variableCount
    error('simplexSolve:DimensionMismatch', ...
        'The number of elements in c must equal the number of columns in A.');
end
if numel(b) ~= constraintCount
    error('simplexSolve:DimensionMismatch', ...
        'The number of elements in b must equal the number of rows in A.');
end

parser = inputParser;
parser.FunctionName = mfilename;
addParameter(parser, 'ObjectiveSense', 'max', ...
    @(value) ischar(value) || (isstring(value) && isscalar(value)));
addParameter(parser, 'ConstraintSense', repmat("<=", constraintCount, 1), ...
    @(value) ischar(value) || isstring(value) || iscellstr(value));
addParameter(parser, 'FreeVariables', false(variableCount, 1), ...
    @(value) islogical(value) || isnumeric(value));
addParameter(parser, 'Tolerance', 1e-9, ...
    @(value) isnumeric(value) && isscalar(value) && isfinite(value) && value > 0);
addParameter(parser, 'MaxIterations', 1000, ...
    @(value) isnumeric(value) && isscalar(value) && isfinite(value) ...
    && value >= 1 && value == floor(value));
parse(parser, varargin{:});

objectiveSense = lower(string(parser.Results.ObjectiveSense));
if ~ismember(objectiveSense, ["max", "min"])
    error('simplexSolve:InvalidObjectiveSense', ...
        'ObjectiveSense must be ''max'' or ''min''.');
end

constraintSense = normalizeConstraintSense( ...
    parser.Results.ConstraintSense, constraintCount);
freeVariables = normalizeFreeVariables( ...
    parser.Results.FreeVariables, variableCount);
tolerance = parser.Results.Tolerance;
maxIterations = parser.Results.MaxIterations;

% Replace each unrestricted variable x_j by x_j^+ - x_j^-.
[variableTransform, transformedA] = transformFreeVariables(A, freeVariables);
transformedC = variableTransform' * c;
transformedVariableCount = size(transformedA, 2);

if objectiveSense == "min"
    maximizationC = -transformedC;
else
    maximizationC = transformedC;
end

% A nonnegative right-hand side provides an obvious initial slack or
% artificial-variable basis. Flip the corresponding constraint when needed.
negativeRows = b < -tolerance;
transformedA(negativeRows, :) = -transformedA(negativeRows, :);
b(negativeRows) = -b(negativeRows);
constraintSense(negativeRows) = flipConstraintSense(constraintSense(negativeRows));
b(abs(b) <= tolerance) = 0;

if constraintCount == 0
    [standardSolution, basis, status, iterations, message, reducedCosts] = ...
        revisedSimplex(transformedA, b, maximizationC, zeros(0, 1), ...
        tolerance, maxIterations);
    phase1Iterations = 0;
    phase2Iterations = iterations;
    removedRows = zeros(0, 1);
else
    [standardA, phase2C, basis, artificialColumns] = ...
        buildStandardForm(transformedA, maximizationC, constraintSense);

    phase1Iterations = 0;
    removedRows = zeros(0, 1);

    if ~isempty(artificialColumns)
        phase1C = zeros(size(standardA, 2), 1);
        phase1C(artificialColumns) = -1;

        [phase1Solution, basis, phase1Status, phase1Iterations, ...
            phase1Message] = revisedSimplex(standardA, b, phase1C, basis, ...
            tolerance, maxIterations);

        if phase1Status ~= 1
            x = nan(variableCount, 1);
            fval = NaN;
            exitflag = phase1Status;
            output = makeOutput(phase1Message, phase1Iterations, 0, basis, ...
                zeros(0, 1), zeros(0, 1), removedRows);
            return
        end

        phase1Objective = phase1C' * phase1Solution;
        if phase1Objective < -10 * tolerance
            x = nan(variableCount, 1);
            fval = NaN;
            exitflag = -2;
            output = makeOutput('The constraints are infeasible.', ...
                phase1Iterations, 0, basis, zeros(0, 1), zeros(0, 1), removedRows);
            return
        end

        [standardA, b, basis, phase2C, removedRows, cleanupSucceeded] = ...
            removeArtificialVariables(standardA, b, basis, phase2C, ...
            artificialColumns, tolerance);

        if ~cleanupSucceeded
            x = nan(variableCount, 1);
            fval = NaN;
            exitflag = -4;
            output = makeOutput('The Phase I basis could not be cleaned safely.', ...
                phase1Iterations, 0, basis, zeros(0, 1), zeros(0, 1), removedRows);
            return
        end
    end

    [standardSolution, basis, status, phase2Iterations, message, reducedCosts] = ...
        revisedSimplex(standardA, b, phase2C, basis, tolerance, maxIterations);
end

exitflag = status;
if exitflag == 1
    transformedSolution = standardSolution(1:transformedVariableCount);
    x = variableTransform * transformedSolution;
    x(abs(x) <= tolerance) = 0;
    fval = c' * x;
elseif exitflag == -3
    x = nan(variableCount, 1);
    if objectiveSense == "max"
        fval = Inf;
    else
        fval = -Inf;
    end
else
    x = nan(variableCount, 1);
    fval = NaN;
end

output = makeOutput(message, phase1Iterations, phase2Iterations, basis, ...
    standardSolution, reducedCosts, removedRows);
end


function sense = normalizeConstraintSense(value, constraintCount)
sense = string(value);
sense = strip(sense(:));

if numel(sense) ~= constraintCount
    error('simplexSolve:DimensionMismatch', ...
        'ConstraintSense must contain one entry for each row of A.');
end

sense(sense == "<" | sense == "=<" | sense == "≤") = "<=";
sense(sense == ">" | sense == "=>" | sense == "≥") = ">=";
sense(sense == "==") = "=";

if any(~ismember(sense, ["<=", ">=", "="]))
    error('simplexSolve:InvalidConstraintSense', ...
        'Each ConstraintSense entry must be ''<='', ''>='', or ''=''.');
end
end


function freeVariables = normalizeFreeVariables(value, variableCount)
if islogical(value)
    freeVariables = value(:);
    if numel(freeVariables) ~= variableCount
        error('simplexSolve:DimensionMismatch', ...
            'A logical FreeVariables value must have one entry per variable.');
    end
    return
end

validateattributes(value, {'numeric'}, {'vector', 'integer', 'positive'});
indices = value(:);
if any(indices > variableCount)
    error('simplexSolve:InvalidFreeVariableIndex', ...
        'FreeVariables contains an index larger than the number of variables.');
end
freeVariables = false(variableCount, 1);
freeVariables(indices) = true;
end


function [transform, transformedA] = transformFreeVariables(A, freeVariables)
variableCount = size(A, 2);
transformedVariableCount = variableCount + nnz(freeVariables);
transform = zeros(variableCount, transformedVariableCount);

nextColumn = 0;
for variable = 1:variableCount
    nextColumn = nextColumn + 1;
    transform(variable, nextColumn) = 1;
    if freeVariables(variable)
        nextColumn = nextColumn + 1;
        transform(variable, nextColumn) = -1;
    end
end

transformedA = A * transform;
end


function flipped = flipConstraintSense(sense)
flipped = sense;
flipped(sense == "<=") = ">=";
flipped(sense == ">=") = "<=";
end


function [standardA, phase2C, basis, artificialColumns] = ...
        buildStandardForm(A, c, constraintSense)
[constraintCount, variableCount] = size(A);
extraColumnCount = nnz(constraintSense == "<=") ...
    + 2 * nnz(constraintSense == ">=") ...
    + nnz(constraintSense == "=");

extraColumns = zeros(constraintCount, extraColumnCount);
isArtificialExtra = false(extraColumnCount, 1);
basis = zeros(constraintCount, 1);
nextExtraColumn = 0;

for row = 1:constraintCount
    switch constraintSense(row)
        case "<="
            nextExtraColumn = nextExtraColumn + 1;
            extraColumns(row, nextExtraColumn) = 1;
            basis(row) = variableCount + nextExtraColumn;
        case ">="
            nextExtraColumn = nextExtraColumn + 1;
            extraColumns(row, nextExtraColumn) = -1;
            nextExtraColumn = nextExtraColumn + 1;
            extraColumns(row, nextExtraColumn) = 1;
            isArtificialExtra(nextExtraColumn) = true;
            basis(row) = variableCount + nextExtraColumn;
        case "="
            nextExtraColumn = nextExtraColumn + 1;
            extraColumns(row, nextExtraColumn) = 1;
            isArtificialExtra(nextExtraColumn) = true;
            basis(row) = variableCount + nextExtraColumn;
    end
end

standardA = [A, extraColumns];
phase2C = [c; zeros(extraColumnCount, 1)];
artificialColumns = variableCount + find(isArtificialExtra);
end


function [x, basis, status, iterations, message, reducedCosts] = ...
        revisedSimplex(A, b, c, basis, tolerance, maxIterations)
[constraintCount, variableCount] = size(A);
iterations = 0;
reducedCosts = c;

if constraintCount == 0
    x = zeros(variableCount, 1);
    if any(c > tolerance)
        status = -3;
        message = 'The objective is unbounded.';
    else
        status = 1;
        message = 'An optimal solution was found.';
    end
    return
end

while iterations < maxIterations
    basisMatrix = A(:, basis);
    if rcond(basisMatrix) < eps(class(basisMatrix))
        x = nan(variableCount, 1);
        status = -4;
        message = 'The basis matrix is numerically singular.';
        return
    end

    basicValues = basisMatrix \ b;
    if any(basicValues < -100 * tolerance) || any(~isfinite(basicValues))
        x = nan(variableCount, 1);
        status = -4;
        message = 'The current basis lost primal feasibility.';
        return
    end
    basicValues(abs(basicValues) <= tolerance) = 0;

    dualPrices = basisMatrix' \ c(basis);
    reducedCosts = c - A' * dualPrices;
    reducedCosts(basis) = 0;

    isNonbasic = true(variableCount, 1);
    isNonbasic(basis) = false;
    enteringCandidates = find(isNonbasic & reducedCosts > tolerance);

    if isempty(enteringCandidates)
        x = zeros(variableCount, 1);
        x(basis) = basicValues;
        x(abs(x) <= tolerance) = 0;
        status = 1;
        message = 'An optimal solution was found.';
        return
    end

    % Bland's rule selects the lowest-index eligible entering variable.
    enteringVariable = enteringCandidates(1);
    direction = basisMatrix \ A(:, enteringVariable);
    eligibleRows = direction > tolerance;

    if ~any(eligibleRows)
        x = nan(variableCount, 1);
        status = -3;
        message = 'The objective is unbounded.';
        return
    end

    ratios = inf(constraintCount, 1);
    ratios(eligibleRows) = basicValues(eligibleRows) ./ direction(eligibleRows);
    minimumRatio = min(ratios);
    tiedRows = find(abs(ratios - minimumRatio) ...
        <= tolerance * max(1, abs(minimumRatio)));

    % Bland's tie break selects the lowest-index leaving basic variable.
    [~, tiedPosition] = min(basis(tiedRows));
    leavingRow = tiedRows(tiedPosition);
    basis(leavingRow) = enteringVariable;
    iterations = iterations + 1;
end

basisMatrix = A(:, basis);
basicValues = basisMatrix \ b;
x = zeros(variableCount, 1);
x(basis) = basicValues;
status = 0;
message = 'The iteration limit was reached.';
end


function [A, b, basis, c, removedRows, succeeded] = ...
        removeArtificialVariables(A, b, basis, c, artificialColumns, tolerance)
artificialMask = false(size(A, 2), 1);
artificialMask(artificialColumns) = true;
rowIdentifiers = (1:size(A, 1))';
removedRows = zeros(0, 1);
succeeded = true;

while true
    artificialBasicRows = find(artificialMask(basis));
    if isempty(artificialBasicRows)
        break
    end

    row = artificialBasicRows(1);
    basisMatrix = A(:, basis);
    if rcond(basisMatrix) < eps(class(basisMatrix))
        succeeded = false;
        return
    end

    transformedRows = basisMatrix \ A;
    isNonbasic = true(size(A, 2), 1);
    isNonbasic(basis) = false;
    candidates = find(~artificialMask & isNonbasic ...
        & abs(transformedRows(row, :)') > tolerance);

    if ~isempty(candidates)
        basis(row) = candidates(1);
    else
        basicValues = basisMatrix \ b;
        if abs(basicValues(row)) > 100 * tolerance
            succeeded = false;
            return
        end
        removedRows(end + 1, 1) = rowIdentifiers(row); %#ok<AGROW>
        A(row, :) = [];
        b(row) = [];
        basis(row) = [];
        rowIdentifiers(row) = [];
    end
end

keepColumns = ~artificialMask;
columnMap = zeros(size(A, 2), 1);
columnMap(keepColumns) = 1:nnz(keepColumns);
basis = columnMap(basis);
A = A(:, keepColumns);
c = c(keepColumns);
end


function output = makeOutput(message, phase1Iterations, phase2Iterations, ...
        basis, standardSolution, reducedCosts, removedRows)
output = struct( ...
    'algorithm', 'two-phase revised simplex', ...
    'message', message, ...
    'iterations', phase1Iterations + phase2Iterations, ...
    'phase1Iterations', phase1Iterations, ...
    'phase2Iterations', phase2Iterations, ...
    'basis', basis, ...
    'standardSolution', standardSolution, ...
    'reducedCosts', reducedCosts, ...
    'removedRedundantRows', removedRows);
end
