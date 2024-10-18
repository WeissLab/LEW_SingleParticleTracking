%% Simulate 3D movies of diffusion with EMCCDs
% version 1.0
% (c) Lucien Weiss, 2024
% Use with citaiton to the Github page


% Output settings
display_images = true;
write_images = false;

%% Emitter movement
photon_rate = 150000;
imagenoise_per_pixel = 3;
D = 4; % um2/s
wavelength = 0.6; % microns

%% Camera settings
exposuretime = 0.0005;
gap_between_exposures = 0.0008;
Camera_gain = 100;
Gain_registers = 1000;

%% Microscope
NA = 1.45; % Numerical Aperture
refractiveIndex = 1.33; % Immersion oil
pixelSize = 0.1; % microns
focal_position = 0; % microns
xySize = 5; % microns: total length of the FOV

%% Simulation_parameters
fine_time_steps = 10; % This number simulates shorter frames, then adds them up to longer frames.
number_of_frames = 10;
number_of_molecules = 10;
initialization_rangeX = 1; % All emitters appear in a -N to N position
initialization_rangeY = 1; % All emitters appear in a -N to N position
initialization_rangeZ = 0.25; % All emitters appear in a Zero to N position

%% Start
empty_image = zeros(xySize/pixelSize+1);
for molecule_number = 1:number_of_molecules
    % Initialize position
    emitterPos = [2*initialization_rangeX*rand()-initialization_rangeX, 2*initialization_rangeY*rand()-initialization_rangeY, initialization_rangeZ*rand()]; % Position of the emitter in microns
    
    % Choose a number of photons
    n_photons = poissrnd(photon_rate*exposuretime/fine_time_steps);
    
    % Molecule loop
    for i = 1:(number_of_frames)
        psf = empty_image;

        % Start with a gap, then do each exposure.
        emitterPos = emitterPos + sqrt(2*1*D*gap_between_exposures)*randn(1,3); % This is the gap
        for finestep = 1:fine_time_steps
            emitterPos = emitterPos + sqrt(2*1*D*exposuretime/fine_time_steps)*randn(1,3);
            psf_finestep = calculate3DPSF(wavelength, NA, refractiveIndex, pixelSize, focal_position, xySize, emitterPos,n_photons);
            psf = psf + psf_finestep;
        end

        % Add background noise
        bg_noise = poissrnd(imagenoise_per_pixel,size(psf));
        % Simulate read noise
        rd_noise = 10*randn(size(psf));
        % Simulate Gain
        Gained_image = gamrnd((psf+bg_noise)*Camera_gain,Gain_registers/Camera_gain,size(psf)) + rd_noise;

        % Display the PSF
        if display_images
            figure(1);
            imagesc(Gained_image); 
            shading interp;
            colormap gray;
            colorbar;
            title({['Molecule ' num2str(molecule_number)];num2str(emitterPos);['frame ' num2str(i)]});
        end

        % Write images
        if write_images
            imwrite(uint16(Gained_image),'examplemovie_3.tif','WriteMode','append');
        end
    end
end

%% calculate3DPSF 4
% Lucien Weiss 2024
function binned_psf = calculate3DPSF(wavelength, NA, refractiveIndex, pixelSize, zPosition, xySize, emitterPos, n_photons)
    % This function calculates a 3D PSF given microscope parameters and emitter position.
    % Inputs:
    % - Emitter wavelength in microns
    % - Numerical Aperture of the objective
    % - Refractive index of the medium (water = 1.33)
    % - Pixel size in the xy plane in microns (we will simulate at 1/2 of this and then rebin pixels, it works better)
    % - zPosition of the focus - This is a simple model, so moving the focus is equivalent to moving the emitter
    % - xySize: Size of the camera in microns
    % - emitterPos: Position of the emitter [x_0, y_0, z_0] in microns
    % - Number of photons to distribute

    % This gets used for binning the small pixels later on.
    fun = @(block_struct) sum(block_struct.data(:));

    % Setup
    k = 2 * pi / wavelength; % Wavenumber
    z = zPosition;
    x = -xySize/2 : pixelSize/2 : xySize/2;
    y = x;
    [X, Y, Z] = meshgrid(x, y, z);
    
    % Unpack the emitter position
    x0 = emitterPos(1);
    y0 = emitterPos(2);
    z0 = emitterPos(3);
    
    % Radial distance in the XY plane relative to the emitter
    r = sqrt((X - x0).^2 + (Y - y0).^2);
    
    % Calculate the lateral (XY) and axial (Z) components of the PSF
    % using the diffraction pattern equations
    alpha = asin(NA / refractiveIndex); % Aperture angle
    rho = (2 * pi * NA / wavelength) * r; % Scaled radial distance
    
    % Introduce a Z-dependent term for the lateral component that blurs with defocus
    defocusTerm = (1 + ((Z - z0) / (wavelength / NA)).^2); % Defocus causes lateral spreading
    rhoZ = rho ./ sqrt(defocusTerm); % Lateral spread increases with Z
    
    % Calculate the Airy disk (sinc pattern) in 3D
    J = besselj(1, rhoZ); % Bessel function of the first kind (order 1)
    
    % Handle the division by zero at the center of the Airy disk (rho=0)
    psfXY = (2 * J ./ rhoZ).^2; % Lateral PSF in XY
    psfXY(r == 0) = 1; % Set the central maximum (for r = 0 where J(0)/0 is undefined)
    
    % Calculate axial intensity profile using a Gaussian approximation
    sigmaZ = (wavelength * refractiveIndex) / (NA^2); % Approximate axial PSF width
    psfZ = exp(-((Z - z0).^2) / (2 * sigmaZ^2)); % Gaussian axial intensity profile centered at z0
    
    % Calculate the 3D PSF by multiplying the XY and Z components
    psf = psfXY .* psfZ;
    
    % Normalize the PSF
    psf = psf / sum(psf(:)) * n_photons;

    binned_psf = poissrnd(blockproc(psf, [2 2], fun));
end
