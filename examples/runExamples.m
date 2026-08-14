function results = runExamples()
%RUNEXAMPLES Solve every problem preserved from the original project.

exampleDirectory = fileparts(mfilename('fullpath'));
projectDirectory = fileparts(exampleDirectory);
sourceDirectory = fullfile(projectDirectory, 'src');
addpath(sourceDirectory);
pathCleanup = onCleanup(@() rmpath(sourceDirectory));

problems = problemCatalog();
results = repmat(struct( ...
    'source', '', ...
    'x', [], ...
    'fval', NaN, ...
    'exitflag', NaN, ...
    'message', ''), numel(problems), 1);

fprintf('%-12s %-10s %-14s %s\n', 'Example', 'Exit flag', 'Objective', 'Solution');
fprintf('%s\n', repmat('-', 1, 72));

for index = 1:numel(problems)
    problem = problems(index);
    [x, fval, exitflag, output] = simplexSolve( ...
        problem.c, problem.A, problem.b, ...
        'ObjectiveSense', problem.objectiveSense, ...
        'ConstraintSense', problem.constraintSense, ...
        'FreeVariables', problem.freeVariables);

    results(index).source = problem.source;
    results(index).x = x;
    results(index).fval = fval;
    results(index).exitflag = exitflag;
    results(index).message = output.message;

    if exitflag == 1
        solutionText = mat2str(x', 6);
    else
        solutionText = output.message;
    end
    fprintf('%-12s %-10d %-14g %s\n', ...
        problem.source, exitflag, fval, solutionText);
end
end
