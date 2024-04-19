%% Run samples of the Inventory simulation
%
% Collect statistics and plot histograms along the way.

%% Set up

% Set-up and administrative cost for each batch requested.
K = 25.00;

% Per-unit production cost.
c = 3.00;

% Lead time for production requests.
L = 2;

% Holding cost per unit per day.
h = 0.05/7;

% Reorder point.
ROP = 141;

% Batch size.
Q = 757;

% How many samples of the simulation to run.
NumSamples = 100;


% Run each sample for this many days.
MaxTime = 1000;

%% Run simulation samples


% Make this reproducible
rng("default");

% Samples are stored in this cell array of Inventory objects
InventorySamples = cell([NumSamples, 1]);

% Run samples of the simulation.
% Log entries are recorded at the end of every day
for SampleNum = 1:NumSamples
    fprintf("Working on %d\n", SampleNum);
    inventory = Inventory( ...
        RequestCostPerBatch=K, ...
        RequestCostPerUnit=c, ...
        RequestLeadTime=L, ...
        HoldingCostPerUnitPerDay=h, ...
        ReorderPoint=ROP, ...
        OnHand=Q, ...
        RequestBatchSize=Q);
    run_until(inventory, MaxTime);
    InventorySamples{SampleNum} = inventory;
end

%% Collect statistics

% Pull the RunningCost from each complete sample.
TotalCosts = cellfun(@(i) i.RunningCost, InventorySamples);


% Express it as cost per day and compute the mean, so that we get a number
% that doesn't depend directly on how many time steps the samples run for.
meanDailyCost = mean(TotalCosts/MaxTime);
fprintf("Mean daily cost: %f\n", meanDailyCost);



%% Make pictures

% record days with backlog
BacklogDaysPerSample = zeros(1,NumSamples);

for k = 1:NumSamples
    inventory = InventorySamples{k};
    backlogDays = 0;

    Temp = inventory.Log.Backlog;

    for i = 1:MaxTime
        u = Temp(i);
        if  u > 0
            backlogDays = backlogDays + 1;
        end
    end
    BacklogDaysPerSample(1,k) = backlogDays;
end

% Print the mean fraction of days with a non zero backlog
meanBacklogDays = mean(BacklogDaysPerSample/MaxTime);
fprintf("Mean fraction of days with a non-zero backlog: %f\n", meanBacklogDays);

fig1 = figure();
t1 = tiledlayout(fig1,1,1);
ax1 = nexttile(t1);

histogram(ax1, BacklogDaysPerSample/MaxTime, Normalization="probability")
title("Fraction of Days With Non-Zero Backlog")
xlabel(ax1,"Fraction of Days That Have Backlog")
ylabel(ax1,"Probability")


%backlog order fraction
backlogfracvec= [];
for i =1:length(InventorySamples)
    backlogcount = 0;
    totalorders = length(InventorySamples{i,1}.Fulfilled) + length(InventorySamples{i,1}.Backlog);
    for j = 1:(length(InventorySamples{i,1}.Fulfilled))
        if InventorySamples{i,1}.Fulfilled{1,j}.OriginalTime ~= InventorySamples{i,1}.Fulfilled{1,j}.Time
            backlogcount = backlogcount + 1;
        end
    end
    backlogcount = backlogcount + length(InventorySamples{i,1}.Backlog);
    backlogfrac = backlogcount/totalorders;
    backlogfracvec(end+1) = backlogfrac;
end

% print the mean fraction of orders that get backloged
meanBacklogFrac = mean(backlogfracvec);
fprintf("Mean fraction of orders that get backlogged: %f\n", meanBacklogFrac);

fig2 = figure();
t2 = tiledlayout(fig2,1,1);
ax2 = nexttile(t2);

histogram(ax2, backlogfracvec, Normalization="probability")
title("Fraction of Orders the Get Backlogged Histogram")
xlabel(ax2,"Fraction of Orders That Get Backloged")
ylabel(ax2,"Probability")


%For days that experience a backlog, the total backlog amount
totalbacklogvec = [];
for i =1:length(InventorySamples)
    for j = 1:MaxTime+1
        comp = InventorySamples{i,1}.Log.Backlog(j);
        if comp > 0
            totalbacklogvec(end+1) = InventorySamples{i,1}.Log.Backlog(j);
        end
    end
end

% print mean total backlog for days with a non zero backlog
meanBacklog = mean(totalbacklogvec);
fprintf("Mean total backlog on days with a non-zero backlog: %f\n", meanBacklog);

fig3 = figure();
t3 = tiledlayout(fig3,1,1);
ax3 = nexttile(t3);


histogram(ax3, totalbacklogvec, Normalization="probability")
title("Total Backlog Amount For Days With Non-Zero Backlog")
xlabel(ax3,"Backlog Amount")
ylabel(ax3,"Probability")



%Delay Time
delayvec = [];
for i =1:length(InventorySamples)
    for j = 1:(length(InventorySamples{i,1}.Fulfilled))
        if InventorySamples{i,1}.Fulfilled{1,j}.OriginalTime ~= InventorySamples{i,1}.Fulfilled{1,j}.Time
            delaytime = InventorySamples{i,1}.Fulfilled{1,j}.Time - InventorySamples{i,1}.Fulfilled{1,j}.OriginalTime;
            delayvec(end+1) = delaytime;
        end
    end
end

% print mean delay time
meanDelay = mean(delayvec);
fprintf("Mean delay time for orders that get backlogged: %f\n", meanDelay);

fig4 = figure();
t4 = tiledlayout(fig4,1,1);
ax4 = nexttile(t4);

histogram(ax4, delayvec, Normalization="probability")
title("Delay Time of Orders That Get Backlogged")
xlabel(ax4,"Delay")
ylabel(ax4,"Probability")


% create hustogram for days with non-zero backlog

%f = figure('Name','Backlog days per sample');
%ft = tiledlayout(f,1,1);
%axf = nexttile(ft);


%b = histogram(axf, BacklogDaysPerSample/MaxTime, Normalization="probability", BinWidth=.05);

%xlabel(axf, "Fraction of days with non-zero backlog");
%ylabel(axf, "Probability");


% Make a figure with one set of axes.
fig = figure();
t = tiledlayout(fig,1,1);
ax = nexttile(t);


  
% Histogram of the cost per day.
h = histogram(ax, TotalCosts/MaxTime, Normalization="probability", ...
    BinWidth=5);

% Add title and axis labels
title(ax, "Daily total cost");
xlabel(ax, "Dollars");
ylabel(ax, "Probability");

% Fix the axis ranges
ylim(ax, [0, 0.5]);
xlim(ax, [240, 290]);

% Wait for MATLAB to catch up.
pause(2);

% Save figure as a PDF file
exportgraphics(fig, "Daily cost histogram.pdf");
