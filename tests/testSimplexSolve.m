function tests = testSimplexSolve
%TESTSIMPLEXSOLVE Regression tests for the revised simplex implementation.
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
testDirectory = fileparts(mfilename('fullpath'));
sourceDirectory = fullfile(fileparts(testDirectory), 'src');
addpath(sourceDirectory);
testCase.TestData.sourceDirectory = sourceDirectory;
end


function teardownOnce(testCase)
rmpath(testCase.TestData.sourceDirectory);
end


function testProductionMaximum(testCase)
c = [50; 40];
A = [2, 4; 3, 1];
b = [80; 60];

[x, fval, exitflag] = simplexSolve(c, A, b);

verifyEqual(testCase, exitflag, 1);
verifyEqual(testCase, x, [16; 12], 'AbsTol', 1e-8);
verifyEqual(testCase, fval, 1280, 'AbsTol', 1e-8);
end


function testGreaterThanMinimizationUsesPhaseOne(testCase)
c = [1; 1];
A = [1, 2];
b = 4;

[x, fval, exitflag, output] = simplexSolve(c, A, b, ...
    'ObjectiveSense', 'min', ...
    'ConstraintSense', ">=");

verifyEqual(testCase, exitflag, 1);
verifyEqual(testCase, x, [0; 2], 'AbsTol', 1e-8);
verifyEqual(testCase, fval, 2, 'AbsTol', 1e-8);
verifyGreaterThan(testCase, output.phase1Iterations, 0);
end


function testFreeVariableAndNegativeRightHandSide(testCase)
c = 1;
A = [1; 1];
b = [2; -3];
sense = ["<="; ">="];

[x, fval, exitflag] = simplexSolve(c, A, b, ...
    'ConstraintSense', sense, ...
    'FreeVariables', true);

verifyEqual(testCase, exitflag, 1);
verifyEqual(testCase, x, 2, 'AbsTol', 1e-8);
verifyEqual(testCase, fval, 2, 'AbsTol', 1e-8);
end


function testInfeasibleProblem(testCase)
c = 1;
A = [1; 1];
b = [1; 2];
sense = ["<="; ">="];

[x, fval, exitflag] = simplexSolve(c, A, b, ...
    'ConstraintSense', sense);

verifyEqual(testCase, exitflag, -2);
verifyTrue(testCase, all(isnan(x)));
verifyTrue(testCase, isnan(fval));
end


function testUnboundedProblem(testCase)
c = [3; -2; 0];
A = [1, 3, 2; 2, 1, 1];
b = [7; 4];

[x, fval, exitflag] = simplexSolve(c, A, b, ...
    'ObjectiveSense', 'min', ...
    'ConstraintSense', [">="; ">="]);

verifyEqual(testCase, exitflag, -3);
verifyTrue(testCase, all(isnan(x)));
verifyEqual(testCase, fval, -Inf);
end


function testRedundantEqualityIsRemoved(testCase)
c = [1; 2; 1];
A = [3, 1, -1; 8, 4, -1; 2, 2, 1];
b = [15; 50; 20];

[x, fval, exitflag] = simplexSolve(c, A, b, ...
    'ConstraintSense', ["="; "="; "="]);

verifyEqual(testCase, exitflag, 1);
verifyEqual(testCase, x, [2.5; 7.5; 0], 'AbsTol', 1e-8);
verifyEqual(testCase, fval, 17.5, 'AbsTol', 1e-8);
end


function testCorrectedBasisRecoveryExample(testCase)
c = [1; 3; 2; 0];
A = [1, -1, -1, 2; 0, 1, 1, -1; 0, 0, 1, -1];
b = [2; 3; 1];

[x, fval, exitflag] = simplexSolve(c, A, b, ...
    'ConstraintSense', ["="; "="; "="]);

verifyEqual(testCase, exitflag, 1);
verifyEqual(testCase, x, [0; 2; 6; 5], 'AbsTol', 1e-8);
verifyEqual(testCase, fval, 18, 'AbsTol', 1e-8);
end
