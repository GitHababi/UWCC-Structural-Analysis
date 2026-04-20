% Do not call this script. This is intended to hide implementation details
% from Main.mlx

% anon function to generate box of where water is for a given waterline.
waterbox =  @(waterline, padding) polyshape( ...
    [-ymax_all - padding,  ymax_all + padding, ymax_all + padding, -ymax_all - padding], ...
    [zmin_all - padding, zmin_all - padding, waterline, waterline]);
% turns buoyantload func. into just func. of the guess trim and draft.
buoyantdist = @(guess) buoyantload(guess(1), guess(2), ...
    internalShapes,xUnique,waterbox,1);

for i = 1:numel(caseNames)
    guessInitial = [draftGuess + zmin_all, 0];
    if isfield(LoadCases.(caseNames{i}),"NotFloating") % for special cases, not in water
        LoadCases.RaceStands.BuoyantDistribution = zeros(1,INTERPOLATION_RATE); 
        continue
    end
    paddlers = LoadCases.(caseNames{i}).PointLoads;
    
    % net force and moment as function of guess
    netforcefromguess = @(guess) trapz(xInterp,selfWeightInterp) + trapz(xUnique,buoyantdist(guess)) + sum(paddlers(2,:));
    netmomentfromguess = @(guess) trapz(xInterp,xInterp.*selfWeightInterp) + trapz(xUnique,xUnique.*buoyantdist(guess)) + sum(paddlers(1,:) .* paddlers(2,:));
    objective = @(guess) (netforcefromguess(guess)^2 + netmomentfromguess(guess)^2)^2 ; % squaring turns 0 into a local min
    % do the optimization!
    options = optimset('PlotFcns','optimplotfval','Display','off');
    [bestGuess, ~] = fminsearch(objective, guessInitial, options);
    disp("for loadcase " + caseNames{i});
    if zmax_all - bestGuess(1) < minFreeboard
        warning('Freeboard is less than required amount.');
    end
    disp(zmax_all - bestGuess(1) + " (in.). freeboard at angle " + bestGuess(2) +  "(rad.)");
    % DEBUG:
    % disp("Sanity Check net force: " + netforcefromguess(bestGuess) + "(lbf) moment: " + ...
    % netmomentfromguess(bestGuess) + "(lbf*in.));
    LoadCases.(caseNames{i}).Draft = bestGuess(1);
    LoadCases.(caseNames{i}).Trim  = bestGuess(2);
    buoyancy = buoyantdist(bestGuess);
    buoyantInterp = interp1(xUnique,buoyancy,xInterp,"linear");
    LoadCases.(caseNames{i}).BuoyantDistribution = buoyantInterp;
end
