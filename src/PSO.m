classdef PSO < handle
    % ParticleSwarmOptimizer Class definition for particle swarm optimization.
    %
    % This class implements a particle swarm optimization algorithm.
    % It allows for customization of various parameters and tracks the optimization process.

	properties
    	costFunction       % Handle to the function that will be minimized.
    	numParticles       % Number of particles in the swarm.
    	numIterations      % Maximum number of iterations to run the optimization.
    	lowerBound         % Lower bounds for the position of particles.
    	upperBound         % Upper bounds for the position of particles.
    	kappa = 1          % Constriction coefficient (influence of the current velocity).
    	phi1 = 2.05        % Cognitive coefficient, personal best influence.
    	phi2 = 2.05        % Social coefficient, global best influence.
    	wDamp = 0.985      % Damping ratio of the inertia coefficient.
    	minW = 0.1         % Minimum value of the inertia coefficient.
    	recordPart = 0     % Number of particles to keep detailed records of.
    	minError = -Inf    % Minimum acceptable error to terminate optimization early.
    	particles          % Array of structures representing the swarm particles.
    	globalBest         % Structure holding the global best position and its cost.
    	bestCosts          % Array to track the best cost found in each iteration.
    	bestPositions      % Matrix to track the best positions found in each iteration.
    	particleRecord     % Array to record detailed information about particles.
    	w                  % Inertia weight of particles.
    	c1                 % Acceleration coefficient to personal best.
    	c2                 % Acceleration coefficient to global best.
    	maxVel             % Maximum velocity a particle can have.
		minVel	           % Minimum velocity a particle can have.
		nVar               % Number of varibles
		varSize            % 1 by 2 array of elements [1, nVar]
	end
    
    methods
		function obj = PSO(costFun, nParticles, nIter, lowerBound, upperBound, varargin)
    		% Constructor for the PSO class initializes properties and processes
    		% optional parameters.
    		obj.costFunction = costFun;
    		obj.numParticles = nParticles;
    		obj.numIterations = nIter;
    		obj.lowerBound = lowerBound(:)';  % Ensure it's a row vector
    		obj.upperBound = upperBound(:)';  % Ensure it's a row vector
		
    		% Create an input parser instance to handle parameter-value pairs.
    		p = inputParser;
		
    		addParameter(p, 'kappa', obj.kappa, @(x) isnumeric(x) && isscalar(x) && x > 0);
			addParameter(p, 'phi1', obj.phi1, @(x) isnumeric(x) && isscalar(x) && x > 0);
			addParameter(p, 'phi2', obj.phi2, @(x) isnumeric(x) && isscalar(x) && x > 0);
			addParameter(p, 'wDamp', obj.wDamp, @(x) isnumeric(x) && isscalar(x) && x > 0);
			addParameter(p, 'minW', obj.minW, @(x) isnumeric(x) && isscalar(x) && x > 0);
			addParameter(p, 'recordPart', obj.recordPart, @(x) isnumeric(x) && isscalar(x) && x > 0);
			addParameter(p, 'minError', obj.minError, @(x) isnumeric(x) && isscalar(x) && x > 0);
		
    		% Parse input parameters and overwrite default if provided by user.
    		parse(p, varargin{:});
    		obj.kappa = p.Results.kappa;
    		obj.phi1 = p.Results.phi1;
    		obj.phi2 = p.Results.phi2;
    		obj.wDamp = p.Results.wDamp;
    		obj.minW = p.Results.minW;
    		obj.recordPart = p.Results.recordPart;
    		obj.minError = p.Results.minError;

			obj.nVar = length(upperBound);
			obj.varSize = [1, obj.nVar];

			phi = obj.phi1 + obj.phi2;
			chi = 2*obj.kappa/abs(2-phi-sqrt(phi^2-4*phi));
			obj.w = chi;
			obj.c1 = chi*obj.phi1;
			obj.c2 = chi*obj.phi2;

			obj.maxVel = 0.2*(upperBound - lowerBound);
			obj.minVel = -obj.maxVel;

    		% Initialize particles after all parameters are set.
    		obj.initializeParticles();
		end
        
        function initializeParticles(obj)
            obj.nVar = length(obj.upperBound);
			emptyParticle = struct('position', [], 'velocity', [], 'cost', Inf, 'bestPosition', [], 'bestCost', Inf);  % Create Empty Particle
            obj.particles = repmat(emptyParticle, obj.numParticles, 1);  % Create Population Array
			obj.particleRecord = repmat(emptyParticle, obj.recordPart, obj.numIterations);  % Create Population Array to be Recorded

            obj.globalBest.cost = Inf;  % Initialize Glabal Best

            for ii = 1:obj.numParticles
                obj.particles(ii).position = obj.unifrandom(obj.lowerBound, obj.upperBound, obj.nVar);  % Generate Random Solution
                obj.particles(ii).velocity = zeros(1, obj.nVar);  % Initialize Velocity
                obj.particles(ii).cost = obj.costFunction(obj.particles(ii).position);  % Evaluation
                obj.particles(ii).bestPosition = obj.particles(ii).position;  % Update Personal Best Posistion
                obj.particles(ii).bestCost = obj.particles(ii).cost;  % Update Personal Best Velocity

                if obj.particles(ii).bestCost < obj.globalBest.cost  % Update Global Best
                    obj.globalBest.position = obj.particles(ii).position;
                    obj.globalBest.cost = obj.particles(ii).cost;
                end
			end
			obj.bestCosts = zeros(obj.numIterations, 1);  % Array to Hold Best Cost
			obj.bestPositions = zeros(obj.numIterations, obj.nVar);  % Matrix to Hold the Best Position at every Iteration
        end
        
        function updateParticles(obj, iter)
			for ii = 1:obj.numParticles
				obj.particles(ii).velocity = obj.w*obj.particles(ii).velocity...  % Update Velocity
					+ obj.c1*rand(obj.varSize).*(obj.particles(ii).bestPosition - obj.particles(ii).position)...
					+ obj.c2*rand(obj.varSize).*(obj.globalBest.position - obj.particles(ii).position);
				
				obj.particles(ii).position = obj.particles(ii).position + obj.particles(ii).velocity;  % Update Posistion
				obj.particles(ii).velocity = max(obj.particles(ii).velocity, obj.minVel);  % Apply Velocity Lower Bound
				obj.particles(ii).velocity = min(obj.particles(ii).velocity, obj.maxVel);  % Apply Velocity Upper Bound
				obj.particles(ii).position = max(obj.particles(ii).position, obj.lowerBound);  % Apply Position Lower Bound
				obj.particles(ii).position = min(obj.particles(ii).position, obj.upperBound);  % Apply Position Upper Bound
				
				obj.particles(ii).cost = obj.costFunction(obj.particles(ii).position);  % Evaluation
				if obj.particles(ii).cost < obj.particles(ii).bestCost
					obj.particles(ii).bestPosition = obj.particles(ii).position;  % Update Personal Best Posistion
					obj.particles(ii).bestCost = obj.particles(ii).cost;  % Update Personal Best Velocity
					if obj.particles(ii).bestCost < obj.globalBest.cost  % Update Global Best
						obj.globalBest.position = obj.particles(ii).bestPosition;
						obj.globalBest.cost = obj.particles(ii).bestCost;
					end
				end
				if ii <= obj.recordPart
					obj.particleRecord(ii,iter) = obj.particles(ii);
				end
			end
		end
        
        function run(obj)
            obj.bestCosts = zeros(obj.numIterations, 1);
            obj.bestPositions = zeros(obj.numIterations, length(obj.upperBound));
            iter = 1;
			while iter <= obj.numIterations && obj.globalBest.cost > obj.minError
                obj.updateParticles(iter);
                obj.bestCosts(iter) = obj.globalBest.cost;
                obj.bestPositions(iter, :) = obj.globalBest.position;
				obj.w = max(obj.w * obj.wDamp, obj.minW); % Sets w
				iter = iter + 1;
            end
        end
        
        %%% Plotting 

        function drawLandscape2d(obj, ax)
            assert(obj.nVar == 2, 'Number of variables must be 2');
            
            size = 150;
            landscape = zeros(size);
            xSpace = linspace(obj.lowerBound(1), obj.upperBound(1), size);
            ySpace = linspace(obj.lowerBound(2), obj.upperBound(2), size);
            
            [X, Y] = meshgrid(xSpace, ySpace);
            
            for ii = 1:size
                for jj = 1:size
                    xVal = X(ii, jj);
                    yVal = Y(ii, jj);
                    landscape(ii, jj) = obj.costFunction([xVal, yVal]);
                end
            end
            
            % Plot the filled contour
            contourf(ax, X, Y, landscape, 50, 'LineStyle', 'none');
            hold(ax, "on")
            plot(ax, obj.globalBest.position(1), obj.globalBest.position(2), 'r.', 'MarkerSize', 20)
        end

        function ax = animateParticles(obj, numAnimateParticles)
            if obj.recordPart == 0
                disp('No particles are recorded.');
                return;
            end
            
            % Draw the landscape
            ax = obj.drawLandscape2d();
            hold on;

            % Get the recorded particles
            recordedParticles = obj.particleRecord;
        
            % Set up the plot for the particles
            h = plot(zeros(numAnimateParticles, 1), zeros(numAnimateParticles, 1), 'r.', 'MarkerFaceColor', 'r');
            
            % Animate the particles
            for iter = 1:obj.numIterations
                if iter > size(recordedParticles, 2) || isempty(recordedParticles(1, iter).position)
                    break;
                end
                
                for ii = 1:numAnimateParticles
                    h.XData(ii) = recordedParticles(ii, iter).position(1);
                    h.YData(ii) = recordedParticles(ii, iter).position(2);
                end
                
                drawnow;
                pause(0.1);  % Adjust the pause duration for animation speed
            end
            
            hold off;
        end

        function drawCostCurve(obj)
            figure();
            plot(obj.bestCosts, 'k-', 'LineWidth', 2)
            xlabel("Iteration Number")
            ylabel("Cost")
            grid on
        end
    end

    methods (Static)
        function out = unifrandom(a, b, size)
	        out = a + (b-a).*rand(1, size);
        end
    end
end
