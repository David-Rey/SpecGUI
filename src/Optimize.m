
classdef Optimize < handle
    properties
        minScale = 0;
        maxScale;
        minPhaseShift = -0.8;
        maxPhaseShift = 0.8;
        minY = 1E4;
        maxY = 3E4;
        centerWavelength;
        readRaman;
        raman;
        realWavelength;
        realIntensity;
        genWavelength;
        genIntensity;
        priorityFunction;
        psoOptinos;
		numParticlesRecord = 0; % number of particles to record
        lowerBound
        upperBound
        pso  % PSO object
        costFun  % function to pass into PSO
        bestShift
        bestScale
        yOffset
    end

    methods
		function obj = Optimize(dataPath, instrumentPath, B_ev, environment, centerWavelength, psoOptinos)
            % Constructor: Initializes the optimization setup with provided parameters and paths.
            obj.centerWavelength = centerWavelength;
            obj.readRaman = ReadSPE(dataPath, true);
            obj.raman = RamanSpec(instrumentPath, B_ev, environment, centerWavelength);
            obj.raman.wavelengthDelta = obj.readRaman.wavelengthDelta * 1E-9;
            obj.realWavelength = obj.readRaman.wavelength;
            obj.realIntensity = obj.readRaman.intensity;
            obj.genWavelength = linspace(obj.readRaman.minWavelength + obj.minPhaseShift, ...
                                          obj.readRaman.maxWavelength + obj.maxPhaseShift, 1000);
            obj.psoOptinos = psoOptinos;
            obj.raman.generateRamanSpectrum(1000);
            obj.genIntensity = obj.raman.intensityArray;
            obj.priorityFunction = obj.generatePriorityFunction();

            % Find max scale
            maxGenIntensity = max(obj.genIntensity);
            maxRealIntensity = max(obj.realIntensity);
            obj.maxScale = maxRealIntensity / maxGenIntensity;
            
            % Set bounds on PSO
            obj.lowerBound = [obj.minPhaseShift, obj.minScale];
            obj.upperBound = [obj.maxPhaseShift, obj.maxScale];
        end

        function priorityFunction = generatePriorityFunction(obj)
            % Generates priority function that multiplies with the error.
            % This prioritizes the peaks of the Raman function.
            a1 = 150;
            a2 = 0.15;
            priorityFunction = zeros(size(obj.realWavelength));
            for i = 1:length(obj.realWavelength)
                for l = obj.raman.peakWavelengths
                    wavelength = obj.realWavelength(i);
                    k1 = wavelength - (l * 1E9);
                    priorityFunction(i) = priorityFunction(i) + exp(-a1 * k1 ^ 2);
                end
            end

            for i = 1:length(obj.realWavelength)
                wavelength = obj.realWavelength(i);
                k2 = wavelength - obj.centerWavelength;
                priorityFunction(i) = priorityFunction(i) * exp(-a2 * k2 ^ 2);
            end
        end

        function mse = errorFunc(obj, shift, scale, y)
            % Calculates the mean squared error between the generated and real intensity data after applying scale and shift.
            testWavelength = obj.genWavelength + shift;
            testIntensity = (obj.genIntensity * scale) + double(y);
            compIntensity = interp1(testWavelength, testIntensity, obj.realWavelength, 'linear', 'extrap');
            diff = compIntensity - double(obj.realIntensity.');
            diffAdj = diff .* obj.priorityFunction;
            mse = mean(diffAdj .^ 2);
        end

        function optimize(obj)
            % Conducts 2D optimization using PSO.
            sortedIntensity = sort(obj.readRaman.intensity);
            numPoints = length(sortedIntensity);
            index99percent = ceil(numPoints * 0.01);
            obj.yOffset = sortedIntensity(index99percent);
    		obj.costFun = @(x) obj.errorFunc(x(1), x(2), obj.yOffset);
            numParticles = obj.psoOptinos.numParticles;
            numIterations = obj.psoOptinos.numIterations;
            obj.pso = PSO(obj.costFun, numParticles, numIterations, obj.lowerBound, obj.upperBound);
			obj.pso.run()
            obj.bestShift = obj.pso.globalBest.position(1);
            obj.bestScale = obj.pso.globalBest.position(2);
        end

        function drawOverlay(obj, ax)
            testWavelength = obj.genWavelength + obj.bestShift;
            testIntensity = (obj.genIntensity * obj.bestScale) + double(obj.yOffset);
            compIntensity = interp1(testWavelength, testIntensity, obj.realWavelength, 'linear', 'extrap');
            plot(ax, obj.realWavelength, obj.realIntensity, 'DisplayName', 'Real', 'Color', [0 0.3470 0.7410]);
            hold(ax, "on");
            plot(ax, obj.realWavelength, compIntensity, 'DisplayName', 'Generated', 'Color', [0.8500 0.3250 0.0980]);
            plot(ax, [obj.realWavelength(1), obj.realWavelength(end)], [obj.yOffset, obj.yOffset], 'k--');
			xlim(ax, [obj.realWavelength(1), obj.realWavelength(end)])
        end
        
        function drawLandscape(obj, ax)
            obj.pso.drawLandscape2d(ax);
        end

        function animateParticles(obj, varargin)
            numAnimateParticles = obj.pso.recordPart;
            if ~isempty(varargin)
                numAnimateParticles = min(varargin{1}, obj.pso.recordPart);
            end
            ax = obj.pso.animateParticles(numAnimateParticles);
            xlabel(ax, "Phase")
            ylabel(ax, "Scale")
        end
    end
end

