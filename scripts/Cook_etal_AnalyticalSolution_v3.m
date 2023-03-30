%% Basic solution and learning how to manipulate it

% This code shows how you can use a hemi-analytical approach in which you
% ask MATLAB to generate the numerical values of the equilibria for a set
% of parameter values.


% First, define the symbols
syms C I M S N

% Second, define the parameters

phiC = 0.001; %coral open recruitment
phiI = 0.0001; %macroalgal open recruitment
gC = 0.1; % growth of coral over turf
gI = 0.6; % growth of immature mac over turf
gM = 0.6; % growth of mature mac over turf ## NOTE: SHOULD THIS GO INTO I? ##
dC = 0.05; % coral mortality
dI = 0.5; % immature mac mortality
dM = 0.3; % mature mac mortality
d = 0.001; % fish non-density-dependent mortality
r = 0.5; % recruitment of immature mac from reproducing local adults
gamma = 0.5; % scaling of mac overgrowth onto coral
omega = 2; % maturation rate of immature mac into mature mac
alpha = 10; % grazing pressure
K = 0.2; % grazer carrying capacity per reef
eff = 0.02; % grazer conversion efficiency; note can't use 'e' because of exponential
zeta = 0.2; % scaling constant reducing herbivory on mature mac
rho = 4; % parrotfish growth relative to unicornfish
f = 0.02662; % fishing rate
sigma = 0; % harvester preference



% Third, input the system of five differential equations. Note that here
% we've made the substitution that T = 1 - C - I - M.
eq1 = phiC*(1-C-I-M) + gC*C*(1-C-I-M) -dC*C-gamma*gM*M*C == 0;
eq2 = phiI*(1-C-I-M) + r*M*(1-C-I-M) + gI*I*(1-C-I-M)-dI*I-omega*I-alpha*N*I-alpha*S*I==0;
eq3 = omega*I+gM*M*(1-C-I-M) + gamma*gM*M*C - dM*M - zeta*alpha*N*M == 0;
eq4 = rho*eff*alpha*((1-C-I-M)+I)*S*(1-S/K)-d*S-f*sigma*S == 0;
eq5 = eff*alpha*((1-C-I-M)+I+zeta*M)*N*(1-N/K) - d*N - f*(1-sigma)*N == 0;

% Use a numerical solver to find the exact values of the equilibria
sol = vpasolve([eq1, eq2, eq3, eq4, eq5],[C,I,M,S,N]); % VPASOLVE finds equilibria numerically
sol.C % Shows, as a vertical vector, the equilibrium values of live coral that are possible

% Extract the relevant results; here are just some examples of how you can
% examine them. Note that we are only interested in real, positive
% equilibria, so you'll see that we subset for these below.
res = [sol.C,sol.I,sol.M,sol.S,sol.N]'; % stores the equilibria in a matrix. We use the transpose here so that the rows are state variables, and the columns are equilibria

res = res(:,abs(sum(imag(res)))==0); % removes imaginary equilibria
posres = res(:,min(res)>=0); % extracts non-negative equilibria
relevantres = posres(:,posres(4,:)>0); % extracts equilibria w/ positive Scaridae
if size(relevantres,2) > 1 % If there are multiple equilibria output...
    relevantres = relevantres(:,relevantres(5,:)>0); % extracts equilibria w/ positive Naso
end
t_res  = relevantres'; % create a transpose of the equilibria so that we can sort by row
sorted_equil = sort(t_res,1) % sorts equilibria from lowest coral to highest coral
    

%% Determining stability

% One thing that complicates this model is that there are two fish species.
% Thus, there could be many more possible equilibria:
%   only Scaridae
%   only Naso
%   neither
%   both

% Examination in R suggests that the two fish species can invade the system
% and can invade one another.... so we should focus on scenarios with
% positive fish abundance first?

% Let's think through the logic of going from low fishing pressure to high
% fishing pressure.

% When fishing pressure = 0, for this model parameterization, both fish
% will be present at the stable equilibrium.

% Eventually, you could fish hard enough to eliminate all of the fish,
% but it looks like with the max effort set at 0.15, this won't happen.
% However, if f = 0.15 and sigma = 0, then unicornfish can be eliminated.

% Thus we can just count equilibria based on the ones that have
% positive densities of Scaridae.


posres = res(:,min(res)>=0); % extracts non-negative equilibria
relevantres = posres(:,posres(4,:)>0); % extracts equilibria w/ positive Scaridae
if size(relevantres,2) > 1 % If there are multiple equilibria output...
    relevantres = relevantres(:,relevantres(5,:)>0); % extracts equilibria w/ positive Naso
end

t_res  = relevantres'; % create a transpose of the equilibria so that we can sort by row
sorted_equil = sort(t_res,1) % sort from lowest coral to highest coral


%% Large simulation to generate surface
% Vectors of parameter values to iterate over
fset = linspace(0,0.15,100); % vector of fishing efforts
sigmaset = linspace(0,1,101); % vector of harvest selectivities

% Set up efficient data collection by creating placeholder matrices
highstab = NaN(length(fset),length(sigmaset),5); % Placeholder for high stable equilibrium
unstab = highstab; % Placeholder for unstable equilibrium
lowstab = highstab; % Placeholder for low unstable equilibrium

% These placeholder arrays have # of rows = # of fishing efforts, # of
% columns = # of selectivities, and a "depth"/"layer" for each of the five
% state variables.



%%

% Collect data on equilibria
tic % sets a timer; useful to understand the performance of this loop.
for i = 1:length(fset) % Loop over every fishing effort
    f = fset(i) % update herbivory value
    eqctr = 0; % Re-set the equilibrium counter
    for j = 1:length(sigmaset) % Loop over every refuge strength
        sigma = sigmaset(j); % update refuge value
        
        % Use the routine explored above to find the real, nonnegative
        % equilibria
        eq1 = phiC*(1-C-I-M) + gC*C*(1-C-I-M) -dC*C-gamma*gM*M*C == 0;
        eq2 = phiI*(1-C-I-M) + r*M*(1-C-I-M) + gI*I*(1-C-I-M)-dI*I-omega*I-alpha*N*I-alpha*S*I==0;
        eq3 = omega*I+gM*M*(1-C-I-M) + gamma*gM*M*C - dM*M - zeta*alpha*N*M == 0;
        eq4 = rho*eff*alpha*((1-C-I-M)+I)*S*(1-S/K)-d*S-f*sigma*S == 0;
        eq5 = eff*alpha*((1-C-I-M)+I+zeta*M)*N*(1-N/K) - d*N - f*(1-sigma)*N == 0;

        sol = vpasolve([eq1, eq2, eq3, eq4, eq5],[C,I,M,S,N]); % VPASOLVE finds equilibria numerically

        res = [sol.C,sol.I,sol.M,sol.S,sol.N]'; % stores the equilibria in a matrix. We use the transpose here so that the rows are state variables, and the columns are equilibria
        res = res(:,abs(sum(imag(res)))==0); % removes imaginary equilibria
        posres = res(:,min(res)>=0); % extracts non-negative equilibria
        relevantres = posres(:,posres(4,:)>0); % extracts equilibria w/ positive Scaridae
        if size(relevantres,2) > 1 % If there are multiple equilibria output...
            relevantres = relevantres(:,relevantres(5,:)>0); % extracts equilibria w/ positive Naso
        end
        t_res  = relevantres'; % create a transpose of the equilibria so that we can sort by row
        sorted_equil = sort(t_res,1); % sorts equilibria from lowest coral to highest coral

        % Save the equilibria in their appropriate holding matrices
        if size(sorted_equil,1)>1 % If there are multiple equilibria, we fill them from low to high. 
            % Note that this represents an ASSUMPTION that the MIDDLE
            % equilibrium is UNSTABLE.
            if size(sorted_equil,1) == 2 % if there are exactly 2 equilibria one's going to be semi-stable (on the boundary w/ the unstable equilibrium)
                lowstab(i,j,:) = sorted_equil(1,:);
                highstab(i,j,:) = sorted_equil(2,:);
                %if eqctr == 0 % if we're just starting to enter the region of bistability
                %    unstab(j,i) = CLstar(2); % then the unstable equilibrium is the same as the high-coral stable one
                %else % Otherwise, it's the same as the low-coral stable one
                %    unstab(j,i) = CLstar(1);
                %end
            else
                eqctr = 1; % Once we get to the 3-equilibrium state, set equilibrium counter to 1
                highstab(i,j,:) = sorted_equil(3,:);
                unstab(i,j,:) = sorted_equil(2,:);
                lowstab(i,j,:) = sorted_equil(1,:);
            end
        else % If there's only one equilibrium, we have to decide if it's the 'high' one or the 'low' one
            if eqctr == 0
                % Remember how we set eqctr == 0 outside of this for loop?
                % Within THIS for loop (the one indexed by j), we're
                % iterating over a range of fishing values. When the
                % fishing = 0, we expect to be at the 'high' coral
                % equilibrium (or, at least, the highest value we can be
                % at). So we'll save the sole equilibrium here.
                if sorted_equil(1,1) > 0.2
                    highstab(i,j,:) = sorted_equil(1,:);
                else
                    lowstab(i,j,:) = sorted_equil(1,:);
                end
            else
                % Otherwise, when we're 'past' the region of bistability
                % (that is, eqctr has now been set to 1 because we've gone
                % through a region with more than 1 equilibrium), we're
                % going to save the equilibrium as the "low coral" state.
                if sorted_equil(1,1) > 0.2
                    highstab(i,j,:) = sorted_equil(1,:);
                else
                    lowstab(i,j,:) = sorted_equil(1,:);
                end
            end
        
        end
        
    end
    toc % Reports out the time that iterating through this value went
    save('Run1_100x101') 
    % ^Highly recommend saving at each iteration of the loop, as this is a 
    % computationally intensive analysis, and it'd be a shame to lose the
    % data. But note that you want to be very careful to not over-write
    % your matrices with zeros, which is why this piece of code is in a
    % separate chunk.
end


%% Plot the surfaces

figure(2)
close(2)
figure(2)
coords = get(gcf,'Position');
set(gcf,'Position',[coords(1),coords(2),coords(3)*2,coords(4)*1.1])

subplot(1,2,1)
surf(sigmaset,fset,highstab(:,:,1),'EdgeColor','n','FaceLighting','flat')
hold on
surf(sigmaset,fset,lowstab(:,:,1))
surf(sigmaset,fset,unstab(:,:,1),'FaceColor','red','EdgeColor','n')
xlabel('Sigma')
ylabel('Fishing effort')
zlabel('Coral Cover')
meshgrid off
shading interp
view(-53,28)
%view(33,30)
%view(70,35)
lighting none


subplot(1,2,2)
%surf(sigmaset,fset,highstab(:,:,2)+highstab(:,:,3),'EdgeColor','n','FaceLighting','flat')
hold on
surf(sigmaset,fset,lowstab(:,:,2)+lowstab(:,:,3))
surf(sigmaset,fset,unstab(:,:,2)+unstab(:,:,3),'FaceColor','red','EdgeColor','n')
xlabel('Sigma')
ylabel('Fishing effort')
zlabel('Macroalgal Cover')
meshgrid off
shading interp
view(-53,28)
%view(33,30)
%view(70,35)
lighting none





%% The surface is pretty ugly so let's try another approach in which we make many lines

% First, make a single bifurcation diagram

% plot hysteresis line
hystchoice = round(length(sigmaset)/10*3,0);

% select lower edge
lowhyst = lowstab(:,hystchoice,1);
lowhyst = lowhyst(isfinite(lowhyst));
lowhystx = fset(isfinite(lowstab(:,hystchoice,1)));

% select unstable edge
unhyst = unstab(:,hystchoice,1);
unhyst = unhyst(isfinite(unhyst));
unhystx = fset(isfinite(unstab(:,hystchoice,1)));

% select upper edge
highhyst = highstab(:,hystchoice,1);
highhyst = highhyst(isfinite(highhyst));
highhystx = fset(isfinite(highstab(:,hystchoice,1)));

% concatenate and plot
xset = [highhystx,wrev(unhystx),lowhystx];
yset = sigmaset(hystchoice)*ones(1,length(xset));
zset = [highhyst',wrev(unhyst'),lowhyst'];
%plot3(xset,yset,zset,'LineWidth',3,'Color',[0.6350 0.0780 0.1840])
plot3(xset,yset,zset,'LineWidth',3,'Color','k')
plot3(xset,yset,zset,'LineWidth',3)%,'Color','k')

view(0,0)
cla
patch([xset nan],[yset nan],[zset nan],[zset nan],'EdgeColor','interp','FaceColor','none','LineWidth',5)


%% Use a for-loop to plot all the lines
hold on
close(1)
figure(1)
clf(1)
for i = 1:length(sigmaset)
% plot hysteresis line
hystchoice = i;

% select lower edge
lowhyst = lowstab(:,hystchoice,1);
lowhyst = lowhyst(isfinite(lowhyst));
lowhystx = fset(isfinite(lowstab(:,hystchoice,1)));

% select unstable edge
unhyst = unstab(:,hystchoice,1);
unhyst = unhyst(isfinite(unhyst));
unhystx = fset(isfinite(unstab(:,hystchoice,1)));

% select upper edge
highhyst = highstab(:,hystchoice,1);
highhyst = highhyst(isfinite(highhyst));
highhystx = fset(isfinite(highstab(:,hystchoice,1)));

% concatenate and plot
xset = [highhystx,wrev(unhystx),lowhystx];
yset = sigmaset(hystchoice)*ones(1,length(xset));
zset = [highhyst',wrev(unhyst'),lowhyst'];
%plot3(xset,yset,zset,'LineWidth',3,'Color',[0.6350 0.0780 0.1840])
%plot3(xset,yset,zset,'LineWidth',3,'Color','k')
%plot3(xset,yset,zset,'LineWidth',3)%,'Color','k')

view(-34,47)
view(117,28)
%cla
patch([xset nan],[yset nan],[zset nan],[zset nan],'EdgeColor','interp','FaceColor','none','LineWidth',10)
xlabel('Fishing effort, f')
ylabel('Fisher preference, \sigma')
zlabel('Coral cover, C*')
xlim([0,.145])

end




%% 
% I am not so happy with the 3D plot but I think an "operating diagram"
% that shows the regions in parameter space with coral-only, mac-only, and
% bistability would be easier to look at and interpret.

% The easiest way to make this is to choose a range of sigma values, and
% find crit_M and crit_C for each one.

% For any sigma value, when f = 0, we will have a coral-dominated state.
% As f increases, we'll move into bistability. The transition point is
% crit_M
% Eventually, we move from bistability into a mac-only state. This
% transition point is crit_C

% The larger sigma is, the greater the critical values of f.


% Let's implement a binary search algorithm to find crit_M and crit_C.

sigma_set = linspace(0,1,100);
critC = NaN*sigma_set;
critM = NaN*sigma_set;


iters = 10; % number of iterations of the binary search algorithm at each step


% We'll start at the lower bound of sigma

sigma = 0;
f = 0.01;
f_low = f;
f_up = -1;

f_holder = f;

while f_up < 0
    
    eq1 = phiC*(1-C-I-M) + gC*C*(1-C-I-M) -dC*C-gamma*gM*M*C == 0;
    eq2 = phiI*(1-C-I-M) + r*M*(1-C-I-M) + gI*I*(1-C-I-M)-dI*I-omega*I-alpha*N*I-alpha*S*I==0;
    eq3 = omega*I+gM*M*(1-C-I-M) + gamma*gM*M*C - dM*M - zeta*alpha*N*M == 0;
    eq4 = rho*eff*alpha*((1-C-I-M)+I)*S*(1-S/K)-d*S-f*sigma*S == 0;
    eq5 = eff*alpha*((1-C-I-M)+I+zeta*M)*N*(1-N/K) - d*N - f*(1-sigma)*N == 0;
    
    sol = vpasolve([eq1, eq2, eq3, eq4, eq5],[C,I,M,S,N]); % VPASOLVE finds equilibria numerically
    
    res = [sol.C,sol.I,sol.M,sol.S,sol.N]'; % stores the equilibria in a matrix. We use the transpose here so that the rows are state variables, and the columns are equilibria
    res = res(:,abs(sum(imag(res)))==0); % removes imaginary equilibria
    posres = res(:,min(res)>=0); % extracts non-negative equilibria
    relevantres = posres(:,posres(4,:)>0); % extracts equilibria w/ positive Scaridae
    if size(relevantres,2) > 1 % If there are multiple equilibria output...
        relevantres = relevantres(:,relevantres(5,:)>0); % extracts equilibria w/ positive Naso
    end
    t_res  = relevantres'; % create a transpose of the equilibria so that we can sort by row
    sorted_equil = sort(t_res,1); % sorts equilibria from lowest coral to highest coral
    
    if size(sorted_equil,1) > 1
        % If there's more than one equilibrium now
        f_up = f;
        f = (f_low + f_up)/2;
    else
        % If we still have only one equilibrium
        if sorted_equil(1,1) < 0.2 % if we've jumped the region of bistability into low coral
            f_up = f; % use the current f value as the upper bound
            f = (f_low + f_up)/2; % test a new f in the middle of the region
        else
            f_low = f; % use the current f value as the lower bound
            f = 1.2*f; % increase f by a small amount and see if that puts us into the new region
        end
    end
    
    f_holder = [f_holder,f];
    
end


f_holder2 = f_holder;

% Now perform binary search algorithm to refine guess.

for i = 1:iters
    eq1 = phiC*(1-C-I-M) + gC*C*(1-C-I-M) -dC*C-gamma*gM*M*C == 0;
    eq2 = phiI*(1-C-I-M) + r*M*(1-C-I-M) + gI*I*(1-C-I-M)-dI*I-omega*I-alpha*N*I-alpha*S*I==0;
    eq3 = omega*I+gM*M*(1-C-I-M) + gamma*gM*M*C - dM*M - zeta*alpha*N*M == 0;
    eq4 = rho*eff*alpha*((1-C-I-M)+I)*S*(1-S/K)-d*S-f*sigma*S == 0;
    eq5 = eff*alpha*((1-C-I-M)+I+zeta*M)*N*(1-N/K) - d*N - f*(1-sigma)*N == 0;
    
    sol = vpasolve([eq1, eq2, eq3, eq4, eq5],[C,I,M,S,N]); % VPASOLVE finds equilibria numerically
    
    res = [sol.C,sol.I,sol.M,sol.S,sol.N]'; % stores the equilibria in a matrix. We use the transpose here so that the rows are state variables, and the columns are equilibria
    res = res(:,abs(sum(imag(res)))==0); % removes imaginary equilibria
    posres = res(:,min(res)>=0); % extracts non-negative equilibria
    relevantres = posres(:,posres(4,:)>0); % extracts equilibria w/ positive Scaridae
    if size(relevantres,2) > 1 % If there are multiple equilibria output...
        relevantres = relevantres(:,relevantres(5,:)>0); % extracts equilibria w/ positive Naso
    end
    t_res  = relevantres'; % create a transpose of the equilibria so that we can sort by row
    sorted_equil = sort(t_res,1); % sorts equilibria from lowest coral to highest coral
    
    if size(sorted_equil,1) > 1
        f_up = f;
    else
        f_low = f;
    end
    
    f = (f_up+f_low)/2;
    
    
    
    f_holder2 = [f_holder2,f];
end

f_holder

f_holder2

critM(1) = f_low;

f_max = 0.2;

for j = 2:length(sigma_set)
    sigma = sigma_set(j);
    % Use the last effort level as your lower bound
    f_low = critM(j-1);
    f = 1.2*f_low;
    f_up = -1;
    
    % First, find a maximum value for f
    while f_up < 0
        if f_low < f_max % if you exceed the max bound, stop computing
            
            eq1 = phiC*(1-C-I-M) + gC*C*(1-C-I-M) -dC*C-gamma*gM*M*C == 0;
            eq2 = phiI*(1-C-I-M) + r*M*(1-C-I-M) + gI*I*(1-C-I-M)-dI*I-omega*I-alpha*N*I-alpha*S*I==0;
            eq3 = omega*I+gM*M*(1-C-I-M) + gamma*gM*M*C - dM*M - zeta*alpha*N*M == 0;
            eq4 = rho*eff*alpha*((1-C-I-M)+I)*S*(1-S/K)-d*S-f*sigma*S == 0;
            eq5 = eff*alpha*((1-C-I-M)+I+zeta*M)*N*(1-N/K) - d*N - f*(1-sigma)*N == 0;
            
            sol = vpasolve([eq1, eq2, eq3, eq4, eq5],[C,I,M,S,N]); % VPASOLVE finds equilibria numerically
            
            res = [sol.C,sol.I,sol.M,sol.S,sol.N]'; % stores the equilibria in a matrix. We use the transpose here so that the rows are state variables, and the columns are equilibria
            res = res(:,abs(sum(imag(res)))==0); % removes imaginary equilibria
            posres = res(:,min(res)>=0); % extracts non-negative equilibria
            relevantres = posres(:,posres(4,:)>0); % extracts equilibria w/ positive Scaridae
            if size(relevantres,2) > 1 % If there are multiple equilibria output...
                relevantres = relevantres(:,relevantres(5,:)>0); % extracts equilibria w/ positive Naso
            end
            t_res  = relevantres'; % create a transpose of the equilibria so that we can sort by row
            sorted_equil = sort(t_res,1); % sorts equilibria from lowest coral to highest coral
            
            
            
            if size(sorted_equil,1) > 1
                % If there's more than one equilibrium now
                f_up = f;
                f = (f_low + f_up)/2;
            else
                % If we still have only one equilibrium
                if sorted_equil(1,1) < 0.2 % if we've jumped the region of bistability into low coral
                    f_up = f; % use the current f value as the upper bound
                    f = (f_low + f_up)/2; % test a new f in the middle of the region
                else
                    f_low = f; % use the current f value as the lower bound
                    f = 1.2*f; % increase f by a small amount and see if that puts us into the new region
                end
            end
        else
            f_up = f_max;
        end
    end
    
    % Second, run binary search algorithm for refinment
    if f_low < f_max % if you're still below the critical threshold
        for i = 1:iters
            eq1 = phiC*(1-C-I-M) + gC*C*(1-C-I-M) -dC*C-gamma*gM*M*C == 0;
            eq2 = phiI*(1-C-I-M) + r*M*(1-C-I-M) + gI*I*(1-C-I-M)-dI*I-omega*I-alpha*N*I-alpha*S*I==0;
            eq3 = omega*I+gM*M*(1-C-I-M) + gamma*gM*M*C - dM*M - zeta*alpha*N*M == 0;
            eq4 = rho*eff*alpha*((1-C-I-M)+I)*S*(1-S/K)-d*S-f*sigma*S == 0;
            eq5 = eff*alpha*((1-C-I-M)+I+zeta*M)*N*(1-N/K) - d*N - f*(1-sigma)*N == 0;
            
            sol = vpasolve([eq1, eq2, eq3, eq4, eq5],[C,I,M,S,N]); % VPASOLVE finds equilibria numerically
            
            res = [sol.C,sol.I,sol.M,sol.S,sol.N]'; % stores the equilibria in a matrix. We use the transpose here so that the rows are state variables, and the columns are equilibria
            res = res(:,abs(sum(imag(res)))==0); % removes imaginary equilibria
            posres = res(:,min(res)>=0); % extracts non-negative equilibria
            relevantres = posres(:,posres(4,:)>0); % extracts equilibria w/ positive Scaridae
            if size(relevantres,2) > 1 % If there are multiple equilibria output...
                relevantres = relevantres(:,relevantres(5,:)>0); % extracts equilibria w/ positive Naso
            end
            t_res  = relevantres'; % create a transpose of the equilibria so that we can sort by row
            sorted_equil = sort(t_res,1); % sorts equilibria from lowest coral to highest coral
            
            if size(sorted_equil,1) > 1
                f_up = f;
            else
                f_low = f;
            end
            
            f = (f_up+f_low)/2;
        end
        
        
        % Third, save the result.
        critM(j) = f_low;
        
    else
        break % If you're "off the charts", stop computing
    end
    save('critical_thresh_v1')
    
end

%%
figure(2)
plot(sigma_set,critM)

%% 
% OK but seriously there's a much better way to do this which is just to
% extract the biologically relevant equilibria and ask whether the max
% coral value is above or below the 0.2 threshold.... duh. Let's try that
% for critC.

for j = 1:length(sigma_set)
    sigma = sigma_set(j);
    
    f_up = f_max;
    f_low = 0;
    f = (f_up+f_low)/2;
    
    if f_low < f_max % if you're still below the critical threshold
        for i = 1:iters
            eq1 = phiC*(1-C-I-M) + gC*C*(1-C-I-M) -dC*C-gamma*gM*M*C == 0;
            eq2 = phiI*(1-C-I-M) + r*M*(1-C-I-M) + gI*I*(1-C-I-M)-dI*I-omega*I-alpha*N*I-alpha*S*I==0;
            eq3 = omega*I+gM*M*(1-C-I-M) + gamma*gM*M*C - dM*M - zeta*alpha*N*M == 0;
            eq4 = rho*eff*alpha*((1-C-I-M)+I)*S*(1-S/K)-d*S-f*sigma*S == 0;
            eq5 = eff*alpha*((1-C-I-M)+I+zeta*M)*N*(1-N/K) - d*N - f*(1-sigma)*N == 0;
            
            sol = vpasolve([eq1, eq2, eq3, eq4, eq5],[C,I,M,S,N]); % VPASOLVE finds equilibria numerically
            
            res = [sol.C,sol.I,sol.M,sol.S,sol.N]'; % stores the equilibria in a matrix. We use the transpose here so that the rows are state variables, and the columns are equilibria
            res = res(:,abs(sum(imag(res)))==0); % removes imaginary equilibria
            posres = res(:,min(res)>=0); % extracts non-negative equilibria
            relevantres = posres(:,posres(4,:)>0); % extracts equilibria w/ positive Scaridae
            if size(relevantres,2) > 1 % If there are multiple equilibria output...
                relevantres = relevantres(:,relevantres(5,:)>0); % extracts equilibria w/ positive Naso
            end
            t_res  = relevantres'; % create a transpose of the equilibria so that we can sort by row
            sorted_equil = sort(t_res,1); % sorts equilibria from lowest coral to highest coral
            
            if max(sorted_equil(:,1)) > 0.2 % if the high-coral equilibrium is still present
                f_low = f;
            else
                f_up = f;
            end
            
            f = (f_up+f_low)/2;
        end
        
        
        % Third, save the result.
        critC(j) = f_low;
        
    else
        break % If you're "off the charts", stop computing
    end
    
    
end
save('critical_thresh_v1')
%%
load('critical_thresh_v1')
figure(2)
clf(2)
hold on
xlabel('Fisher selectivity, \sigma')
ylabel('Fishing effort, f')
ylim([0,.15])
xlim([0,1])
critM_noNaN = critM;
critM_noNaN(isnan(critM_noNaN)) = f_max;
area([0 1],[1 1],'FaceColor',[0.9290 0.6940 0.1250],'Facealpha',.5,'edgealpha',0)
area(sigma_set,critC,'FaceColor',[0 0.4470 0.7410],'Facealpha',.5,'Edgealpha',0)
area([sigma_set 1],[critM_noNaN 1],'FaceColor','white','Edgealpha',0)
area([sigma_set 1],[critM_noNaN 1],'FaceColor',[0 0.4470 0.7410],'Facealpha',.5,'Edgealpha',0)
plot(sigma_set,critM,'LineWidth',2,'Color','black','LineStyle','--')
plot(sigma_set,critC,'LineWidth',2,'Color','black','LineStyle','-')
box on
text(.6,.02,'Coral-dominated','FontSize',18)
text(.1,.12,'Macroalgae-dominated','FontSize',18)
text(.15,.041,'Bistable','FontSize',18)
set(gca,'FontSize',18)



%% Figure 4 - bifurcation diagrams in f
% In the main text figure, we include panels that show the bifurcation
% diagrams along a specific transect. Here, we generate those model outputs
% for plotting.

    
phiC = 0.001; %coral open recruitment
phiI = 0.0001; %macroalgal open recruitment
gC = 0.1; % growth of coral over turf
gI = 0.6; % growth of immature mac over turf
gM = 0.6; % growth of mature mac over turf ## NOTE: SHOULD THIS GO INTO I? ##
dC = 0.05; % coral mortality
dI = 0.5; % immature mac mortality
dM = 0.3; % mature mac mortality
d = 0.001; % fish mortality
r = 0.5; % recruitment of immature mac from reproducing local adults
gamma = 0.5; % scaling of mac overgrowth onto coral
omega = 2; % maturation rate of immature mac into mature mac
alpha = 10; % grazing pressure
K = 0.2; % grazer carrying capacity per reef
eff = 0.02; % grazer conversion efficiency; note can't use 'e' because of exponential
zeta = 0.2; % scaling constant reducing herbivory on mature mac
rho = 4; % parrotfish growth relative to unicornfish



sigma_set = [0.25, 0.5, .75]; % harvester preference
f_set = linspace(0,0.15,100); % fishing rates

Cstar_high = NaN(length(f_set),length(sigma_set));
Cstar_low = Cstar_high;
Cstar_un = Cstar_high;

Mstar_high = Cstar_high;
Mstar_low = Cstar_high;
Mstar_un = Cstar_high;

for k = 1:length(sigma_set)
    
    sigma = sigma_set(k); % update sigma value
    eqctr = 0;
    
for j = 1:length(f_set)
    f = f_set(j); % update fishing value
    
    % Use the routine explored above to find the real, nonnegative
    % equilibria
    eq1 = phiC*(1-C-I-M) + gC*C*(1-C-I-M) -dC*C-gamma*gM*M*C == 0;
    eq2 = phiI*(1-C-I-M) + r*M*(1-C-I-M) + gI*I*(1-C-I-M)-dI*I-omega*I-alpha*N*I-alpha*S*I==0;
    eq3 = omega*I+gM*M*(1-C-I-M) + gamma*gM*M*C - dM*M - zeta*alpha*N*M == 0;
    eq4 = rho*eff*alpha*((1-C-I-M)+I)*S*(1-S/K)-d*S-f*sigma*S == 0;
    eq5 = eff*alpha*((1-C-I-M)+I+zeta*M)*N*(1-N/K) - d*N - f*(1-sigma)*N == 0;
    
    sol = vpasolve([eq1, eq2, eq3, eq4, eq5],[C,I,M,S,N]); % VPASOLVE finds equilibria numerically
    
    res = [sol.C,sol.I,sol.M,sol.S,sol.N]'; % stores the equilibria in a matrix. We use the transpose here so that the rows are state variables, and the columns are equilibria
    res = res(:,abs(sum(imag(res)))==0); % removes imaginary equilibria
    posres = res(:,min(res)>=0); % extracts non-negative equilibria
    relevantres = posres(:,posres(4,:)>0); % extracts equilibria w/ positive Scaridae
    if size(relevantres,2) > 1 % If there are multiple equilibria output...
        relevantres = relevantres(:,relevantres(5,:)>0); % extracts equilibria w/ positive Naso
    end
    t_res  = relevantres'; % create a transpose of the equilibria so that we can sort by row
    sorted_equil = sort(t_res,1); % sorts equilibria from lowest coral to highest coral
    
    % Save the equilibria in their appropriate holding matrices
    if size(sorted_equil,1)>1 % If there are multiple equilibria, we fill them from low to high.
        % Note that this represents an ASSUMPTION that the MIDDLE
        % equilibrium is UNSTABLE.
        if size(sorted_equil,1) == 2 % if there are exactly 2 equilibria one's going to be semi-stable (on the boundary w/ the unstable equilibrium)
            Cstar_high(j,k) = sorted_equil(2,1);
            Cstar_low(j,k) = sorted_equil(1,:);
            Mstar_high(j,k) = sorted_equil(1,2)+sorted_equil(1,3);
            Mstar_low(j,k) = sorted_equil(2,2)+sorted_equil(2,3);
        else
            eqctr = 1; % Once we get to the 3-equilibrium state, set equilibrium counter to 1
            Cstar_high(j,k) = sorted_equil(3,1);
            Cstar_un(j,k) = sorted_equil(2,1);
            Cstar_low(j,k) = sorted_equil(1,1);
            Mstar_high(j,k) = sorted_equil(3,2)+sorted_equil(3,3);
            Mstar_un(j,k) = sorted_equil(2,2)+sorted_equil(2,3);
            Mstar_low(j,k) = sorted_equil(1,2)+sorted_equil(1,3);
        end
    else % If there's only one equilibrium, we have to decide if it's the 'high' one or the 'low' one
        if sorted_equil(1,1) > 0.2
            Cstar_high(j,k) = sorted_equil(1,1);
            Mstar_low(j,k) = sorted_equil(1,2)+sorted_equil(1,3);
        else
            Cstar_low(j,k) = sorted_equil(1,1);
            Mstar_high(j,k) = sorted_equil(1,2)+sorted_equil(1,3);
        end
    end
    
    
    
    
    
end

    

 

end


%% Make the figure

% Create the figure and size it
figure(3)
close(3)
figure(3)
coords = get(gcf,'Position');
set(gcf,'Position',[coords(1),coords(2),coords(3)*1.3,coords(4)*2])


% Set up subpanels
rowsperpanel = 3;
colsperpanel = 3;
subplotpanelsC = [1 2 3 7 8 9 13 14 15];
subplotpanelsM = [4 5 6 10 11 12 16 17 18];

% Set up the color scheme
coralcol = [70/255,129/255,184/255];
maccol = [117/255,90/255,69/255];

% Make the plots
for k = 1:length(sigma_set)
    
subplot(length(sigma_set)*rowsperpanel,2*colsperpanel,subplotpanelsC+(k-1)*2*rowsperpanel*colsperpanel)
plot(f_set,Cstar_low(:,k),'LineWidth',3,'Color',coralcol)
hold on
plot(f_set,Cstar_un(:,k),'LineStyle','--','LineWidth',3,'Color',coralcol)
plot(f_set,Cstar_high(:,k),'LineWidth',3,'Color',coralcol)
ylim([0,.7]) % set common y-axis limit from 0 to 0.7
set(gca,'FontSize',15)

% Except for bottom panel, suppress x-axis ticks
if k ~= length(sigma_set)
    set(gca,'Xtick',[])
end
% Add bounding box around bistability region
test = isnan(Cstar_un(:,k));
lowbound = f_set(min(find(test==0)));
highbound = f_set(max(find(test==0)));
if length(lowbound) > 0
area([lowbound highbound highbound lowbound],[-1 -1 1 1],'FaceColor','k','Edgealpha',0,'FaceAlpha',.20)
end
% Add title if it's the top plot
if k == 1
    title('Coral Equilibrium, C*')
end
% Add y-axis label
ylabel('Proportion Benthic Cover')

% Add x-axis label if it's the bottom plot
if k == length(sigma_set)
    xlabel('Fishing Pressure, f')
end

subplot(length(sigma_set)*rowsperpanel,2*colsperpanel,subplotpanelsM+(k-1)*2*rowsperpanel*colsperpanel)
plot(f_set,Mstar_low(:,k),'LineWidth',3,'Color',maccol)
hold on
plot(f_set,Mstar_un(:,k),'LineStyle','--','LineWidth',3,'Color',maccol)
plot(f_set,Mstar_high(:,k),'LineWidth',3,'Color',maccol)
ylim([0,.7]) % set common y-axis limit from 0 to 0.7
set(gca,'FontSize',15)

% Except for bottom panel, suppress x-axis ticks
if k ~= length(sigma_set)
    set(gca,'Xtick',[])
end
% Suppress y-axis ticks
set(gca,'YTick',[])

% Add bounding box around bistability region
if length(lowbound) > 0
area([lowbound highbound highbound lowbound],[-1 -1 1 1],'FaceColor','k','Edgealpha',0,'FaceAlpha',.20)
end

% Add title if it's the top plot
if k == 1
    title('Macroalgal Equilibrium, M*')
end
% Add x-axis label if it's the bottom plot
if k == length(sigma_set)
    xlabel('Fishing Pressure, f')
end


end
   