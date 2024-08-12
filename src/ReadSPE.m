

classdef ReadSPE < handle
	properties
    	filepath
    	wavelength
    	intensity
	
    	maxIntensity
    	minIntensity
    	maxWavelength
    	minWavelength
    	wavelengthDelta
	end
	
	methods
		function obj = ReadSPE(filepath, shift)
        	obj.filepath = filepath;
        	[obj.wavelength, obj.intensity] = obj.spectra_from_spe();
            if shift
                obj.wavelength = obj.wavelength + 2;
            end
        	obj.maxIntensity = max(obj.intensity);
        	obj.minIntensity = min(obj.intensity);
        	obj.maxWavelength = max(obj.wavelength);
        	obj.minWavelength = min(obj.wavelength);
        	obj.wavelengthDelta = obj.wavelength(2) - obj.wavelength(1);
		end

        function [wavelengthArray, intensity] = spectra_from_spe(obj)
            [width, height, frame, wavelengthArray] = obj.openXMLline();
	        data = obj.readBinaryData(width, height, frame);
        
	        intensity = double(data{1});
			
        end

        function data = readBinaryData(obj, Width, Height, Frame)
            % Open the file in 'read' mode
            fileID = fopen(obj.filepath, 'r');
        
            % Read the entire file into a uint8 array
            bytes = fread(fileID, '*uint8');
        
            % Close the file
            fclose(fileID);

            dataTypes = {'uint32', 'int32', 'int16', 'uint16', [], 'double', 'uint8', [], 'uint32'};
            np_type = dataTypes(bytes(109)+1);
            np_type = np_type{1};
        
            % Get the item size in bytes for the MATLAB data type
            itemsize = numel(typecast(cast(0, np_type), 'uint8'));
        
            % Calculate the total number of data points per frame
            Count = Width * Height;
        
            % Data extraction from the buffer
            data = cell(1, Frame);
            for i = 0:Frame-1
                offset = 4100 + i * Count * itemsize;
                data{i+1} = typecast(bytes(offset + (1:Count*itemsize)), np_type);
            end
        end

        function [width, height, frame, wavelengthArray] = openXMLline(obj)
            % Read the file
            fid = fopen(obj.filepath, 'rb');
            if fid == -1
                error('Failed to open file: %s', obj.filepath);
            end
            lines = textscan(fid, '%s', 'Delimiter', '\n');
            lines = lines{1};
            fclose(fid);
        
            XMLline = lines{end};
        
	        framePos = strfind(XMLline, "Frame");
	        frame = XMLline(framePos(1)+14:end); 
	        framePosEnd = strfind(frame, '"');
	        frame = str2double(frame(1:framePosEnd(1)-1));
	        
        
	        widthPos = strfind(XMLline, "width");
	        width = XMLline(widthPos(1)+7:end);
	        widthPosEnd = strfind(width, '"');
	        width = str2double(width(1:widthPosEnd(1)-1));
        
	        heightPos = strfind(XMLline, "height");
	        height = XMLline(heightPos(1)+8:end);
	        heightPosEnd = strfind(height, '"');
	        height = str2double(height(1:heightPosEnd(1)-1));
        
	        %dataStartPos = strfind(XMLline, "<Wavelength ");
	        %dataEndPos = strfind(XMLline, "</Wavelength>");
	        %waveData = XMLline(dataStartPos(1): dataEndPos(1));
        
	        % Extract the data between the tags
	        dataStart = strfind(XMLline, '<Wavelength ') + 1;
	        dataEnd = strfind(XMLline, '</Wavelength') - 1;
	        dataString = XMLline(dataStart(1):dataEnd(1));
	        waveDataStart = strfind(dataString, '>');
	        waveData = dataString(waveDataStart(1)+1:end);
	        
	        % Split the string by commas and convert to an array of doubles
	        wavelengthArray = str2double(strsplit(waveData, ','));
	        wavelengthArray = wavelengthArray(1:width);
        end

        function drawRaman(obj)
			figure;
            set(gcf, 'Units', 'inches', 'Position', [1, 1, 8, 5]);
            plot(obj.wavelength, obj.intensity, 'DisplayName', 'Raman Intensity')
            xlabel('Wavelength (nm)')
            ylabel('Intensity')
            xlim([obj.wavelength(1), obj.wavelength(end)])

			titleStr = {'Raman Intensity vs. Wavelength', obj.filepath};
            title(titleStr, 'Interpreter', 'none')
            grid on
            legend show
        end
    end
end

