classdef Thomson  < handle
	properties
		filenameSignal  % filename of signal spe (this includes background)
		filenameBackground  % filename of background spe
		specSignalReader  % ReadSPE object that gets the data from the spe
		specBackgroundReader  % ReadSPE object that gets the data from the spe

        psoOptinos  % PSO Options
        pso  % PSO Object

		wavelength  % raw wavelength
		signalIntensity  % signal intensity
		signalIntensityGaussian  % signal intensity smoothed
		gaussianSmoothFactor = 100

		centerWavelength  % center wavlength according to local min of gaussian
		centerIntensity  % intensity at center wavlength

		wavelegthFilteredLeft  % raw wavlengths with center part removed
		intensityFilteredLeft  % raw intensity with center part removed

		wavelegthFilteredRight  % raw wavlengths with center part removed
		intensityFilteredRight  % raw intensity with center part removed

		wavelegthFiltered  % raw wavlengths with center part removed
		intensityFiltered  % raw intensity with center part removed

		gaussianWavelengths  % wavelengths of fitted gaussian
		intensityWavelengths  % intensity of fitted gaussian

		area  % area under the curve in nm*counts
		area_SI  % area under the curve m*counts

        FWHM  % full width at half maximum in nm
        FWHM_SI  % full width at half maximum in m

        temperature  % temperature of thomson in eV

        drift  % election drift in m/s
	end

	methods
		function obj = Thomson(filenameSignal, filenameBackground, psoOptinos)
			obj.filenameSignal = filenameSignal;
			obj.filenameBackground = filenameBackground;
			obj.specSignalReader = ReadSPE(obj.filenameSignal, false);  % create specSignalReader object
			obj.specBackgroundReader = ReadSPE(obj.filenameBackground, false);  % create specBackgroundReader object

            % PSO Options
            obj.psoOptinos = psoOptinos;

			obj.wavelength = obj.specSignalReader.wavelength;  % set the wavelenth array (can be either wavlengths)

			 % subtract the differance specSignalReader and specBackgroundReader
			 % to get raw signal
			obj.signalIntensity = (obj.specSignalReader.intensity - obj.specBackgroundReader.intensity).';

			% apply a smoothing factor to signal intensity
			obj.signalIntensityGaussian = smoothdata(obj.signalIntensity, "gaussian", obj.gaussianSmoothFactor);

			obj.removeFilter()
			obj.curveFit()
		end

		function removeFilter(obj)
			apprxCenter = (obj.wavelength(end) + obj.wavelength(1)) / 2;
			[localMinX, localMinY, localMaxX, ~] = obj.findLocalExtrema(obj.wavelength, obj.signalIntensityGaussian);
			
			% add check to see if therer are local mins

			diff = abs(apprxCenter - localMinX);
			[~, minIndex] = min(diff);
			obj.centerWavelength = localMinX(minIndex);
			obj.centerIntensity = localMinY(minIndex);
			
			localMaxXRightHalf = localMaxX(localMaxX > obj.centerWavelength);
			% add check if localMaxXRightHalf is empty
			
			rightCutoffWavelength = localMaxXRightHalf(1);
			rightCutoffIndex= find(obj.wavelength == rightCutoffWavelength);

			localMaxXLefttHalf = localMaxX(localMaxX < obj.centerWavelength);
			% add check if localMaxXLefttHalf is empty
			
			leftCutoffWavelength = localMaxXLefttHalf(end);
			leftCutoffIndex = find(obj.wavelength == leftCutoffWavelength);

			obj.wavelegthFilteredLeft = obj.wavelength(1:leftCutoffIndex);
			obj.intensityFilteredLeft = obj.signalIntensity(1:leftCutoffIndex);

			obj.wavelegthFilteredRight = obj.wavelength(rightCutoffIndex:end);
			obj.intensityFilteredRight = obj.signalIntensity(rightCutoffIndex:end);

			obj.wavelegthFiltered = [obj.wavelegthFilteredLeft, obj.wavelegthFilteredRight];
			obj.intensityFiltered = [obj.intensityFilteredLeft, obj.intensityFilteredRight];
		end

		function curveFit(obj)
			gaussian = @(x, a, mean, sigma, offset) a .* exp(-((x - mean) .^ 2) / (2 * sigma ^ 2)) + offset;
			x_data = obj.wavelegthFiltered;
			y_data = obj.intensityFiltered;
			objective_function = @(params) sum((gaussian(x_data, params(1), params(2), params(3), params(4)) - y_data) .^ 2);

            delta = max(obj.signalIntensity) - min(obj.signalIntensity);
            min_std = 0.1;
            max_std = 5;

            numParticles = obj.psoOptinos.numParticles;
            numIterations = obj.psoOptinos.numIterations;

            lowerBound = [0, obj.wavelength(1), min_std, min(obj.signalIntensity)];
            upperBound = [delta, obj.wavelength(end), max_std, max(obj.signalIntensity)];

            obj.pso = PSO(objective_function, numParticles, numIterations, lowerBound, upperBound);
            obj.pso.run()

            a = obj.pso.globalBest.position(1);
            mean_nm = obj.pso.globalBest.position(2);
            std = obj.pso.globalBest.position(3);
            offset = obj.pso.globalBest.position(4);

            mean = mean_nm / 1E9;
			
			obj.gaussianWavelengths = linspace(obj.wavelength(1), obj.wavelength(end));
			obj.intensityWavelengths = gaussian(obj.gaussianWavelengths, a, mean_nm, std, offset);
			
            % calculate area / electron density
			obj.area = a * std * sqrt(2*pi);
			obj.area_SI = obj.area / 1E9;
            
            shift_nm = abs(532 - mean_nm);
            shift = shift_nm / 1E9;
            obj.drift = obj.cal_election_drift(shift, mean);

            %obj.FWHM = 2*sqrt(2*log(2)) * std;
            %obj.FWHM_SI = obj.FWHM / 1E9;
            obj.temperature = obj.cal_temperature(std / 1E9, mean);
            obj.saveThomson()
		end
		
		function drawThomson(obj, ax)
			plot(ax, obj.wavelength, obj.signalIntensity, 'Color', [0.3 0.5 0.9 0.3])
			hold(ax, "on")
			plot(ax, obj.wavelegthFilteredLeft, obj.intensityFilteredLeft, 'Color', [0 0 0.9])
			plot(ax, obj.wavelegthFilteredRight, obj.intensityFilteredRight, 'Color', [0 0 0.9])
			%plot(ax, obj.wavelength, obj.signalIntensityGaussian, 'Color', [1 0 0])
			plot(ax, obj.gaussianWavelengths, obj.intensityWavelengths, 'k-')
            xlim(ax, [obj.wavelength(1), obj.wavelength(end)])
            ylim(ax, [min(obj.signalIntensity), max(obj.signalIntensity)])
        end

        function saveThomson(obj)
            % Define the filename
            filename = 'ThomsonData.csv';
            
            % Combine the wavelength and signal intensity into a matrix
            data = [obj.wavelength(:), obj.signalIntensity(:)];
            
            % Write the matrix to a CSV file
            writematrix(data, filename);
        end
	end

	methods (Static)
		function [localMinX, localMinY, localMaxX, localMaxY] = findLocalExtrema(x, y)
    		% Initialize arrays to hold the x and y values of local minimums and maximums
    		localMinX = [];
    		localMinY = [];
    		localMaxX = [];
    		localMaxY = [];
    		
    		% Loop through the array y to find local minimums and maximums
    		for i = 2:length(y)-1
        		if y(i) < y(i-1) && y(i) < y(i+1)
            		localMinX = [localMinX, x(i)];
            		localMinY = [localMinY, y(i)];
        		elseif y(i) > y(i-1) && y(i) > y(i+1)
            		localMaxX = [localMaxX, x(i)];
            		localMaxY = [localMaxY, y(i)];
        		end
    		end
        end

        function T_e_K = cal_temperature(std, mean)
            m_e = 9.10938356e-31;  % Electron mass in kg
            k_B = 1.380649e-23;  % Boltzmann's constant in J/K
            c = 3e8;  % speed of light
            theta = 90;  % incident laser wavelength and scattering angle
            
            lambda_1e = sqrt(2) * std * 1E9;

            % Calculate electron temperature in K
            %T_e = (m_e * c^2) / (8 * k_B * sind(theta / 2)^2) * (spectrum_width / mean)^2;
            T_e = 5238 * lambda_1e^2;

            % Convert temperature from J to eV (1 eV = 1.60218e-19 J)
            T_e_K = T_e / 11604;
        end

        function v_d = cal_election_drift(shift, mean)
            c = 3e8;  % speed of light
            theta = 90;  % incident laser wavelength and scattering angle
            v_d = (shift * c) / (2 * mean * sind(theta / 2));
        end
	end
end
