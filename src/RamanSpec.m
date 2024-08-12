classdef RamanSpec < handle
    %
    % A class to analyze Raman spectra data, calculate Raman intensities, and visualize the spectral data.
    %
    properties
        dataPath  % Instrument function data path
        instrumentTable  % Table of instrument function from path
        wavelengths  % Wavelengths of instrument function
        intensities  % Intensity of instrument function
        centerWavelengthNm  % Center wavelength in nanometers
        centerWavelengthM  % Center wavelength in meters

        maxJ  % Maximum rotational energy state
        B  % Rotational constant
        environment  % environment struct that includes temperature, pressure and volume

        h = 6.62607015e-34  % Planck's constant in J*s
        c = 299792458  % Speed of light in m/s
        k = 1.380649e-23  % Boltzmann's constant in J/K
        gammaSquared = 0.505E-48  % Gamma squared value
        wavelengthDelta  % Difference between wavelengths of read data
        partitionFunctionValue  % Partition function value
        centerInstrumentFunction  % Center wavelength of the instrument function
        maxInstrumentFunction  % Maximum value of the instrument function
        n  % Number density of gas (atoms/m^3)

        peakWavelengths  % Wavelength array corresponding to peaks of Raman spectra
        peakIntensities  % Intensities array corresponding to peaks of Raman spectra
        wavelengthArray  % Wavelengths of Raman spectra
        intensityArray  % Intensities of Raman spectra
    end

    methods
        function obj = RamanSpec(dataPath, rotationalConstantEv, environment, centerWavelength)
            % Constructor method to initialize the RamanSpec object with data path, rotational constant, temperature, and center wavelength.
            obj.dataPath = dataPath;  % Path to instrument function
            obj.instrumentTable = readtable(dataPath);  % Table of instrument function from path
            obj.wavelengths = obj.instrumentTable.Wavelength;  % Wavelengths of instrument function
            obj.intensities = obj.instrumentTable.Intensity;  % Intensity of instrument function
            obj.centerWavelengthNm = centerWavelength;  % Center wavelength in nanometers
            obj.centerWavelengthM = centerWavelength * 1E-9;  % Center wavelength in meters
            obj.maxJ = 30;  % Maximum rotational energy state
            obj.B = rotationalConstantEv * 1.60218e-19;  % Conversion from eV to Joules
            obj.environment = environment;  % Temperature
            obj.partitionFunctionValue = obj.calculatePartitionFunction();  % Partition function value
            [obj.centerInstrumentFunction, obj.maxInstrumentFunction] = obj.calculateCenterInstrumentFunction();
			obj.n = obj.getGasDensity();
        end

        function [center, maxValue] = calculateCenterInstrumentFunction(obj)
            % Function to calculate the center and max value of the instrument function.
            [maxValue, maxIndex] = max(obj.intensities);
            center = obj.wavelengths(maxIndex) * 1E-9;
		end

        function intensity = interpolateIntensity(obj, targetWavelength)
            % Interpolates intensity values for a specific target wavelength given known wavelengths and their corresponding intensities.
            intensity = interp1(obj.wavelengths, obj.intensities, targetWavelength * 1e9, 'linear', 'extrap');
        end

        function partitionFunctionValue = calculatePartitionFunction(obj)
            % Computes the partition function to normalize the population distribution.
            partitionFunctionValue = 0;
            for J = 0:obj.maxJ-1
                E_J = obj.B * J * (J + 1);
                g_J = 6 * (mod(J, 2) == 0) + 3 * (mod(J, 2) == 1);
                partitionFunctionValue = partitionFunctionValue + g_J * (2 * J + 1) * exp(-E_J / (obj.k * obj.environment.temperature));
            end
        end

        function n_J = calculatePopulation(obj, J)
            % Calculate the population of rotational level J.
            E_J = obj.B * J * (J + 1);
            g_J = 6 * (mod(J, 2) == 0) + 3 * (mod(J, 2) == 1);
            n_J = obj.n / obj.partitionFunctionValue * g_J * (2 * J + 1) * exp(-E_J / (obj.k * obj.environment.temperature));
        end

        function generatePeakWavelengths(obj)
            % Generates peaks of Raman spectrum.
            obj.peakWavelengths = zeros(1, 2 * obj.maxJ - 2);  % Create array for peak wavelengths
            obj.peakIntensities = zeros(1, 2 * obj.maxJ - 2);  % Create array for peak intensities
            for J = 1:obj.maxJ  % Sum over J from 0 to maxJ
                lambda_J_P2 = obj.centerWavelengthM + (obj.centerWavelengthM ^ 2 / (obj.h * obj.c)) * obj.B * (4 * J + 6);
                lambda_J_M2 = obj.centerWavelengthM - (obj.centerWavelengthM ^ 2 / (obj.h * obj.c)) * obj.B * (4 * J - 2);
                obj.peakWavelengths(2*J-1:2*J) = [lambda_J_P2, lambda_J_M2];
                obj.peakIntensities(2*J-1:2*J) = [obj.calculateRamanIntensity(lambda_J_P2), obj.calculateRamanIntensity(lambda_J_M2)];
            end
        end

        function intensity = calculateRamanIntensity(obj, wavelength)
            % Calculates the Raman intensity for a given wavelength, accounting for both Stokes and anti-Stokes contributions.
            intensity = 0;
            for J = 1:obj.maxJ  % Sum over J from 0 to maxJ
                lambda_J_P2 = obj.centerWavelengthM + (obj.centerWavelengthM ^ 2 / (obj.h * obj.c)) * obj.B * (4 * J + 6);
                lambda_J_M2 = obj.centerWavelengthM - (obj.centerWavelengthM ^ 2 / (obj.h * obj.c)) * obj.B * (4 * J - 2);
                b_J_P2 = (3 * (J + 1) * (J + 2)) / (2 * (2 * J + 1) * (2 * J + 3));
                b_J_M2 = (3 * J * (J - 1)) / (2 * (2 * J + 1) * (2 * J - 1));
                n_J = obj.calculatePopulation(J);

                lambdaArray = [lambda_J_P2, lambda_J_M2];
                b_jArray = [b_J_P2, b_J_M2];
                for i = 1:2
                    l = lambdaArray(i);  % Wavelength
                    b_J = b_jArray(i);
                    omega = 1 / (100 * l);
                    d_sigma_d_omega = (64 * pi^4 * b_J * omega^4 * obj.gammaSquared) / 1E4;
                    wavelengthJ = wavelength - l;
                    interpolatedIntensity = obj.interpolateIntensity(wavelengthJ + obj.centerInstrumentFunction);
                    intensity = intensity + n_J * obj.wavelengthDelta * obj.environment.power * d_sigma_d_omega * interpolatedIntensity;
                end
            end
        end

        function generateRamanSpectrum(obj, numPoints)
            % Generates the Raman spectrum with a specified number of points.
            specWidth = 4;  % Width of the spectrum
            obj.generatePeakWavelengths();  % Generate peak wavelengths
            obj.wavelengthArray = linspace(obj.centerWavelengthNm - specWidth, obj.centerWavelengthNm + specWidth, numPoints) * 1e-9;
            obj.intensityArray = arrayfun(@(w) obj.calculateRamanIntensity(w), obj.wavelengthArray);
        end

        function drawRamanSpectra(obj)
            % Draws the Raman spectra.
            figure;
            set(gcf, 'Units', 'inches', 'Position', [1, 1, 8, 6]);
            plot(obj.wavelengthArray * 1e9, obj.intensityArray, 'DisplayName', 'Raman Intensity');
            xlabel('Wavelength (nm)');
            ylabel('Intensity');
            title('Raman Intensity vs. Wavelength');
            legend show;
            grid on;
        end

        function atoms = getGasDensity(obj)
		    P = obj.environment.pressure;  % Pa
		    T = obj.environment.temperature;  % kelvin
            V = obj.environment.volume;  % m^3
    
		    R = 8.314;  % J / mol
		    Av = 6.022E23;  % atoms/mol
		    
		    moles = (P*V) / (R*T);
		    atoms = moles * Av;
        end
    end
end
