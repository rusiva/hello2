function positive = biGaussFilter(data,suspects,options)

%% Input validation
arguments
    data {mustBeNumeric}
    suspects {mustBeVector,mustBeNumeric}
    options.to_exclude {mustBeNumericOrLogical} = [];
    options.min_separation {mustBeNumeric,mustBeScalarOrEmpty} = 0;
    options.display {mustBeNumericOrLogical} = false;
end

if ndims(data) < 2 || ndims(data) > 3
    error('Size of "data" must be 2 or 3')
else
    frame_size = [];
    if ndims(data) == 3
        frame_size = size(data,[1 2]);
        data = reshape(data,[],size(data,3));
    end
end

if ismatrix(options.to_exclude)
    if isempty(frame_size)
        frame_size = size(options.to_exclude);
    end
    options.to_exclude = find(options.to_exclude);
end
%% Auxiliary functions
% addpath('I:\Programs\Ruslans\Auxiliary')
dualGauss = @(x,a1,b1,a2,b2,c) a1*exp(-(x-b1).^2/c^2) + a2*exp(-(x-b2).^2/c^2);
objFn = @(p, x, y) sum( (y - dualGauss(x, p(1), p(2), p(3), p(4), p(5))) .^ 2 );
fit_options = optimset('Display','off');
%%
% Remove NOPs
suspects( ismember(suspects,options.to_exclude) ) = [];

% Initialize 2-gaussian fitting
coeffs = zeros(length(suspects),5);
fval = zeros(length(suspects),1);
fitexitflag = zeros(length(suspects),1);

fprintf('Fitting started: Total of %u cells\n',length(suspects))
tic
for i = 1:length(suspects)
    if mod(i,1e3) == 0, fprintf('%u\n',i);  end
    [y,x] = histcounts(data(suspects(i),:),round(length(data(suspects(i),:))/2.5));
    x = x(1:end-1);

    fit_model = @(params) objFn(params,x,y);
    Start_params = [max(y)/2, x(1)+0.15*(x(end)-x(1)), max(y)/2, x(end)-0.25*(x(end)-x(1)), 100];

    [coeffs(i,:),fval(i),fitexitflag(i),~] = fminsearchbnd( fit_model, Start_params,...
        [0 x(1) 0 x(1) 10], [10*max(y) x(end) 10*max(y) x(end) 5e3], fit_options );
end
toc
deltaX = abs(coeffs(:,4)-coeffs(:,2));
widthX = coeffs(:,5);

% Rayleigh threshold for discriminating two closely-spaced peaks
condition = widthX < deltaX/(2*sqrt(2*log(1/0.4)));

deltaX_sorted = sort(deltaX(condition));
cumulative_distribution = arrayfun(@(thresh) nnz(deltaX_sorted>thresh),unique(deltaX_sorted));

condition = condition & deltaX >= options.min_separation;

%% Output
% positive = struct([]);
positive.linds = suspects(condition);
if ~isempty(frame_size)
    positive.map = false(frame_size);
    positive.map( suspects(condition) ) = true;
end
positive.deltaX = deltaX(condition);
positive.widthX = widthX(condition);

positive.deltaX_sorted = unique(deltaX_sorted);
positive.cumcount = cumulative_distribution;

positive.data = data(positive.linds,:)';
positive.unconfirmed_data = data(suspects(~condition),:)';

% rmpath('I:\Programs\Ruslans\Auxiliary')

%% Display
if options.display
    % Plot the fitting results for all pixels as a scatter plot
    tmp = condition;
    coloredScatterPlot(deltaX(tmp),widthX(tmp),...
        [],[],...
        'Noise vs inter-peak distance',...
        'on',0,[]);
    xlabel('Interpeak width (dig. u.)'); ylabel('Standard deviation of individual levels (dig. u.)')
    hold on; plot( linspace(0,max(deltaX(tmp)),100) , linspace(0,max(deltaX(tmp)),100)/(2*sqrt(2*log(1/0.4))) , 'k--' )
    text( max(deltaX(tmp))*0.7 , max(widthX(tmp))*0.3 ,...
    sprintf('%d pixels',nnz(condition & deltaX>options.min_separation)) )

    figure;
    plot(unique(deltaX_sorted)/options.min_separation,cumulative_distribution)
    xlabel(sprintf('Flickering amplitude (\\timesNoise)')); ylabel('Number of RTS candidates'); grid on

end

end