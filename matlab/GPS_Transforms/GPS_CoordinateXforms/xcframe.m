function [transPos, transVel, transAcc] = xcframe(time, orgPos, orgVel, ...
  orgAcc, sourceFrame, destFrame, input, output) 
% [destpos, destvel, destacc] = XCFRAME(time, pos, vel, acc, source, 
%    dest, input, output)
%
% Transforms position from one coordinate frame to another. For a GEI/J2000/GEO
% <-> GEI/J2000/GEO transformation velocity and acceleration is also supported.
% The position might be given, and returned, in Cartesian (x,y,z) or spherical
% (azimuth, elevation, distance) coordinates. For the GEO frame, geodetic
% coordinates (latitude,longitude,altitude) is also supported.
%
% Currently supported coordinate frames are,
% 
% Geo (Earth) centered systems:
% GEO        Geographic, also known as Greenwich Rotating Coordinates (GRC),
%            Earth-fixed Greenwich (EFG) or Earth Centered Earth Fixed (ECEF).
% GEI        Geocentric Equatorial Inertial, also known as True Equator and 
%            True Equinox of Date, True of Date (TOD), ECI, GCI.
% J2000      Geocentric Equatorial Inertial for epoch J2000.0 (GEI2000), also
%            known as Mean Equator and Mean Equinox of J2000.0
% MAG        Geomagnetic
% GSE        Geocentric Solar Ecliptic
% GSM        Geocentric Solar Magnetospheric
% SM         Solar Magnetic
% RTN        Radial Tangential Normal (Earth-centered)*
% GSEQ       Geocentric Solar Equatorial
% * = this system has not been implemented/tested and should be treated as such
%
% Helio (Sun) centered systems:
% HEE        Heliocentric Earth Ecliptic
% HAE        Heliocentric Aries Ecliptic
% HEEQ       Heliocentric Earth Equatorial
%
%
% INPUT
% time     - Required if "source" is different from "dest": time(s) in the 
%            source coordinate frame. Should be a date vector/matrix where each
%            row is [year month day hours minutes seconds] in UTC Gregorian
%            time. Use the MATLAB functions "datetime" and "datavec" to convert
%            to this format from other time formats. For example, use,
%            time  = datevec(datetime(jtime, 'ConvertFrom', 'juliandate'));
%            to convert from time given in Julian days in "jtime".
%
% pos      - Required: position(s) in the source coordinate frame. Should be a
%            three-element vector or a matrix with three columns where each row
%            is a position vector. See also the "input" parameter.
%
% vel      - Optional (set to [] if not used): velocity vector(s) in the
%            source coordinate frame in Cartesian (vx,vy,vz) coordinates.
%            Should be a three-element vector or a matrix with three columns
%            where each row is a velocity vector. The unit must be position
%            length unit/s. This input is only supported for GEI/J2000/GEO <->
%            GEI/J2000/GEO transformations.
%
% acc      - Optional (set to [] if not used): acceleration vector(s) in the
%            source coordinate frame in Cartesian (ax,ay,az) coordinates.
%            Should be a three-element vector or a matrix with three columns
%            where each row is an acceleration vector. The unit must be position
%            length unit/s^2. This input is only supported for GEI/J2000/GEO <->
%            GEI/J2000/GEO transformations and you must then also supply the
%            velocity in "vel" above.
%
% source   - Required: source coordinate frame as one of the strings above.
%
% dest     - Required: the coordinate frame to convert as one of the strings 
%            above.
%
% input    - Optional (omit or set to '' to use default): string specifying the
%            input position coordinate representation. Available options are:
%            'cart'   - Cartesian coordinates (default)
%            'sph'    - Spherical coordinates with azimuth/elevation in degrees.
%            'geod72' - Geodetic coordinates with latitude/longitude in degrees
%                       and elevation in meters. Uses the WGS72 ellipsoid. This
%                       option is only available in the source frame is 'GEO'.
%            'geod84' - As above but uses the WGS84 ellipsoid.
%
% output   - Optional (omit or set to '' to use default): string specifying the
%            output position coordinate representation, same options as above.
%            Note that the 'geod72' and 'geod84' options are only available in
%            the destination frame is 'GEO' and the input position length unit
%            must then be meters.
%
%
% OUTPUT
% transPos - The position in the new coordinate frame.
%
% transVel - Velocity in the new coordinate frame, in position length unit/s, if
%            velocity in the source frame was supplied. Only supported for
%            GEI/J2000/GEO <-> GEI/J2000/GEO transformations.
%
% transAcc - Acceleration in the new coordinate frame, in position length
%            unit/s^2, if acceleration in the source frame was supplied. Only
%            supported for GEI/J2000/GEO <-> GEI/J2000/GEO transformations.
%
%
% NOTE
% You might specify the same source and destination frame, to do coordinate 
% representation transformations. for example: 
% geodPos = xcframe([], cartPos, [], [], 'GEO', 'GEO', 'cart', 'geod84')
% transforms cartesian (x,y,z) coordinates in the GEO frame to geodetic
% (latitude,longitude,altitude) coordinates in the same frame.
%
%
% This implementation is based on the CXFORM C code written by Ed Santiago
% (LANL) and Ryan Boller (NASA/GSFC) that can be found here,
% <a href="https://spdf.sci.gsfc.nasa.gov/pub/software/old/selected_software_from_nssdc/coordinate_transform/#Misc">CXFORM: Coordinate transformation package</a> 
% The transformations are based on Mike Hapgood's work in, 
% Planet. Space Sci. Vol. 40, No. 5. pp. 71l-717, 1992.
%
% VERSION 1.0, 2017 
% By Patrik Forssén (SatStar Ltd & Karlstad University) 
%
% See also: datetime datevec

% Default input
if (nargin < 7 || isempty(input)) , input  = 'cart'; end
if (nargin < 8 || isempty(output)), output = 'cart'; end

% Default output
transPos = [];
transVel = [];
transAcc = [];


%---------------------------------SANITY CHECKS---------------------------------
% Check number of input parameters
if (nargin < 6)
  error('xcframe:tooFewParameters', ['xcframe requires at least 6 input ', ...
    'parameters'])
end

% Check frames
frames = {'J2000', 'GEI', 'GEO', 'MAG', 'GSE', 'GSM', 'SM', 'RTN', 'GSEQ', ...
  'HEE', 'HAE', 'HEEQ'};
if (~ischar(sourceFrame) && ~ismember(upper(sourceFrame), frames))
  error('xcframe:wrongSourceFrame', ['xcframe parameter "source" ', ...
    'should be a string with a supported coordinate frame'])
end
if (~ischar(destFrame) && ~ismember(upper(destFrame), frames))
   error('xcframe:wrongDestFrame', ['xcframe parameter "dest" ', ...
    'should be a string with a supported coordinate frame'])
end

% Check orgPos
errID  = 'xcframe:wrongInputPosition';
errStr = ['xcframe parameter "pos" should be a position vector with ', ...
  'three elements or a position vector with three columns'];
if (isempty(orgPos) || ~isnumeric(orgPos) || ~isreal(orgPos) || ...
    ~ismatrix(orgPos))
  error(errID, errStr)
end
sz = size(orgPos);
if (sz(2) ~= 3 && sz(1) == 3)
  % Assume that the input was given transposed...
  orgPos = transpose(orgPos);
  sz     = size(orgPos);
end
if (sz(2) ~= 3)
  error(errID, errStr)
end
nPos     = size(orgPos, 1);
transPos = orgPos;  % Default output

% Check time input
errID  = 'xcframe:wrongInputTime';
if (nPos == 1)
  errStr = ['xcframe parameter "time" should be a time vector with ', ...
    'six real scalars'];
else
  errStr = ['xcframe parameter "time" should be a real time matrix ', ...
    'with ', num2str(nPos), ' rows and six columns'];
end
if (~strcmpi(sourceFrame, destFrame))
  if (isempty(time) || ~isnumeric(time) || ~isreal(time) || ~ismatrix(time))
    error(errID, errStr)
  end
  sz = size(time);
  if (sz(2) ~= 6 && sz(1) == 6)
    % Assume that the input was given transposed...
    time = transpose(time);
    sz   = size(time);
  end
  if (sz(2) ~= 6), error(errID, errStr), end
  if (max(time(:, 2)) > 12 || min(time(:, 2)) < 1 || ...
      max(time(:, 3)) > 31 || min(time(:, 3)) < 1 || ...
      max(time(:, 4)) > 23 || min(time(:, 4)) < 0 || ...
      max(time(:, 5)) > 59 || min(time(:, 5)) < 0 || ...
      max(time(:, 6)) > 60 || min(time(:, 6)) < 0)
    error(errID, errStr)
  end
else
  % Not used
  time = [];
end

velAccFrames = {'J2000', 'GEI', 'GEO'};
% Check orgVel
errID  = 'xcframe:wrongInputVelocity';
if ((ismember(sourceFrame, velAccFrames) && ...
    ismember(destFrame, velAccFrames)) || strcmpi(sourceFrame, destFrame))
  if (~isempty(orgVel))
    if (nPos == 1)
      errStr = ['xcframe parameter "vel" should be a velocity ', ...
        'vector with three real scalars'];
    else
      errStr = ['xcframe parameter "vel" should be a real velocity ', ...
        'vector matrix with ', num2str(nPos), ' rows and three columns'];
    end
    if (~isnumeric(orgVel) || ~isreal(orgVel) || ~ismatrix(orgVel))
      error(errID, errStr)
    end
    sz = size(orgVel);
    if (sz(2) ~= 3 && sz(1) == 3)
      % Assume that the input was given transposed...
      orgVel = transpose(orgVel);
      sz     = size(orgVel);
    end
    if (sz(2) ~= 3 || sz(1) ~= nPos)
      error(errID, errStr)
    end
  end
else
  % Not used
  if (~isempty(orgVel))
    % Issue a warning
    warning('xcframe:velocityNotSupported', ['Velocity transformation ', ...
      'from frame ''', sourceFrame, ''' to ''', destFrame, ''' not supported'])
  end
  orgVel = [];  
end
if (nargout >= 2)
  transVel = orgVel;  % Default output
else
  % Not requsted, no need to calculate
  orgVel = [];
end

% Check orgAcc
errID  = 'xcframe:wrongInputAcceleration';
if ((ismember(sourceFrame, velAccFrames) && ...
    ismember(destFrame, velAccFrames)) || strcmpi(sourceFrame, destFrame))
  if (~isempty(orgAcc))
    if (nPos == 1)
      errStr = ['xcframe parameter "acc" should be an acceleration ', ...
        'vector with three real scalars'];
    else
      errStr = ['xcframe parameter "vel" should be a real ', ...
        'acceleration vector matrix with ', num2str(nPos), ...
        ' rows and three columns'];
    end
    if (~isnumeric(orgAcc) || ~isreal(orgAcc) || ~ismatrix(orgAcc))
      error(errID, errStr)
    end
    sz = size(orgAcc);
    if (sz(2) ~= 3 && sz(1) == 3)
      % Assume that the input was given transposed...
      orgAcc = transpose(orgAcc);
      sz     = size(orgAcc);
    end
    if (sz(2) ~= 3 || sz(1) ~= nPos)
      error(errID, errStr)
    end
  end
else
  % Not used
  if (~isempty(orgAcc))
    % Issue a warning
    warning('xcframe:accelerationNotSupported', ['Acceleration ', ...
      'transformation from frame ''', sourceFrame, ''' to ''', destFrame, ...
      ''' not supported'])
  end
  orgAcc = [];  
end
if (~isempty(orgAcc) && ~strcmpi(sourceFrame, destFrame) && isempty(orgVel))
  % Velocity must also be supplied
  error('xcframe:velocityNotSupplied', ['You must supply velocity in ', ...
    'order to transform acceleration'])
end
if (nargout >= 3)
  transAcc = orgAcc;  % Default output
else
  % Not requsted, no need to calculate
  orgAcc = [];
end

% Check "input" parameter
errStr = ['xcframe parameter "input" should be a string ', ...
  '''cart'', ''sph'', ''geod72'' or ''geod84'''];
errID  = 'xcframe:wrongInputSpecification';
if (~ischar(input)), error(errID, errStr), end
% Dont be to hard, one letter is sufficient...
if (input(1) == 'c' || input(1) == 'C')
  input = 'cart'; 
elseif (input(1) == 's' || input(1) == 'S')
  input = 'sph';
elseif (input(1) == 'g' || input(1) == 'G')
  wgsType = str2double(regexprep(input, '[^\d]', ''));
  if (isnan(wgsType)) 
    input = 'geod84'; 
  elseif (wgsType == 72)
    input = 'geod72';
  elseif (wgsType == 84)
    input = 'geod84';
  else
    error(errID, errStr)
  end
else
  error(errID, errStr)
end
if (~strcmpi(sourceFrame, 'GEO') && (strcmpi(input, 'geod72') || ...
    strcmpi(input, 'geod84')))
  error(errID, ['xcframe input coordinate system ''', input, ''' not ', ...
    'supported for source frame ''', sourceFrame, ''''])
end

% Check "output" parameter
errStr = ['xcframe parameter "output" should be a string ', ...
  '''cart'', ''sph'', ''geod72'' or ''geod84'''];
errID  = 'xcframe:wrongOutputSpecification';
if (~ischar(output)), error(errID, errStr), end
% Dont be to hard, one letter is sufficient...
if (output(1) == 'c' || output(1) == 'C')
  output = 'cart'; 
elseif (output(1) == 's' || output(1) == 'S')
  output = 'sph';
elseif (output(1) == 'g' || output(1) == 'G')
  wgsType = str2double(regexprep(output, '[^\d]', ''));
  if (isnan(wgsType)) 
    output = 'geod84'; 
  elseif (wgsType == 72)
    output = 'geod72';
  elseif (wgsType == 84)
    output = 'geod84';
  else
    error(errID, errStr)
  end
else
  error(errID, errStr)
end
if (~strcmpi(destFrame, 'GEO') && (strcmpi(output, 'geod72') || ...
    strcmpi(output, 'geod84')))
  error(errID, ['xcframe output coordinate system ''', output, ''' not ', ...
    'supported for destination frame ''', destFrame, ''''])
end


%-------------------------------CALCULATIONS------------------------------------
if (~isempty(time))
  % Convert times to "ephemeris seconds past J2000 (1 Jan 2000 12:00)"
  et = (juliandate(datetime(time)) - 2451545)*86400;
  
  % Check IGRF constants and possibly issue warnings
  fracYearIndex = (min(et) + 3155803200)/157788000;
  fracYear      = rem(fracYearIndex, 1);
  [~, warnFlag] = calcG01(fracYearIndex, fracYear);
  if (warnFlag == -1)
    warning('xcframe:igrfExtrapolationBefore', ['Extrapolation of IGRF ', ...
      'constants will be done before first available value'])
  end
  fracYearIndex = (max(et) + 3155803200)/157788000;
  fracYear      = rem(fracYearIndex, 1);
  [~, warnFlag] = calcG01(fracYearIndex, fracYear);
  if (warnFlag == 1)
    warning('xcframe:igrfExtrapolationAfter', ['Extrapolation of IGRF ', ...
      'constants will be done after last available value'])
  end
  etNow = (juliandate(datetime(datevec(now))) - 2451545)*86400;
  fracYearIndex = (etNow + 3155803200)/157788000;
  fracYear      = rem(fracYearIndex, 1);
  [~, warnFlag] = calcG01(fracYearIndex, fracYear);
  if (warnFlag == 1)
    warning('xcframe:igrfOutdated', ['IGRF constants are outdated, ', ...
      'please update the code'])
  end
end

% Possibly have to convert input to Cartesian coordinates
if (strcmpi(input, 'sph'))
  [x, y, z] = sph2cart(orgPos(:, 1)*pi/180, orgPos(:, 2)*pi/180, orgPos(:, 3));
  orgPos    = [x(:), y(:), z(:)];
elseif (strcmpi(input(1:3), 'geo'))
  orgPos    = geod2cart(orgPos, input);
end

% Get conversion cell
convCell = getConvCell(sourceFrame, destFrame);

% Do conversion
if (~isempty(convCell))
  % Source and destination frame are different!
  transPos = zeros(nPos, 3);
  nSteps   = length(convCell);
  for posNo = 1 : nPos
    oldPos = orgPos(posNo, :);
    for stepNo = 1 : nSteps
      stepCell  = convCell{stepNo};
      funHandle = stepCell{1};
      dirStr    = stepCell{2};
      newPos    = funHandle(et(posNo), oldPos, dirStr);
      oldPos    = newPos;
    end
    transPos(posNo, :) = newPos(:)';
  end
  
  % Possibly calculate velocity and acceleration
  if (~isempty(orgVel))
    [transVel, transAcc] = calcVelAcc(time, orgPos, orgVel, orgAcc, ...
      transPos, sourceFrame, destFrame);
  end
else
  % No need to transform
  transPos = orgPos;
end

% Possibly convert result
if (strcmpi(output, 'sph'))
  [azim, elev, r] = cart2sph(transPos(:, 1), transPos(:, 2), ...
    transPos(:, 3));
  transPos(:, 1)  = azim*180/pi;
  transPos(:, 2)  = elev*180/pi;
  transPos(:, 3)  = r;
elseif (strcmpi(output(1:3), 'geo'))
  transPos = cart2geod(transPos, output);
end

end



%-----------------------------CALCULATION ROUTINES------------------------------

function cartPos = geod2cart(geodPos, ellipsoid)

% Ellipsoid data
switch ellipsoid
  case 'geod72'
    a = 6378135;          % Semimajor axis [m]
    f = 1/298.26;         % Flattening
  case 'geod84'
    a = 6378137;          % Semimajor axis [m]
    f = 1/298.257223563;  % Flattening
end

% Input data
lat = geodPos(:, 1)*pi/180; lng = geodPos(:, 2)*pi/180; alt = geodPos(:, 3);

% Calculations
sinlat = sin(lat);
e2     = f*(2 - f);
n      = a./sqrt(1 - e2*sinlat.^2);
rho    = (n + alt).*cos(lat);
z      = (n*(1 - e2) + alt).*sinlat;

[x, y]  = pol2cart(lng, rho);

% Output data
cartPos = [x(:), y(:), z(:)];

end



function geodPos = cart2geod(cartPos, ellipsoid)

% Ellipsoid data
switch ellipsoid
  case 'geod72'
    a = 6378135;          % Semimajor axis [m]
    f = 1/298.26;         % Flattening
  case 'geod84'
    a = 6378137;          % Semimajor axis [m]
    f = 1/298.257223563;  % Flattening
end

% Input data
x = cartPos(:, 1); y = cartPos(:, 2); z = cartPos(:, 3);

% Longitude
[lng, rho] = cart2pol(x, y);

% Ellipsoid data
b   = (1 - f)*a;     % Semiminor axis
e2  = f*(2 - f);     % Squared first eccentricity
ep2 = e2/(1 - e2);   % Squared second eccentricity

% Initial values for beta and latitude
beta = atan2(z, (1 - f)*rho);
lat  = atan2(z + b*ep2*sin(beta).^3, rho - a*e2*cos(beta).^3);

% Iterate using Bowring's formula to get latitude
betaNew = atan2((1 - f)*sin(lat), cos(lat));
iter    = 0;
while (any(beta(:) ~= betaNew(:)) && iter < 5)
  beta    = betaNew;
  lat     = atan2(z + b*ep2*sin(beta).^3, rho - a*e2*cos(beta).^3);
  betaNew = atan2((1 - f)*sin(lat), cos(lat));
  iter    = iter + 1;
end

% Altitude
sinphi  = sin(lat);
n       = a./sqrt(1 - e2*sinphi.^2);
alt     = rho.*cos(lat) + (z + e2*n.*sinphi).*sinphi - n;

% Output
geodPos = [lat(:)*180/pi, lng(:)*180/pi, alt(:)];

end



function theta = gmst(et)
% Caculates the Greenwich mean sidereal time in degrees given ephemeris seconds
% past J2000

% T0 is Julian centuries from J2000
temp  = -6.2e-6*T0(et).^3 + 0.093104*T0(et).^2 + ...
  (876600*3600 + 8640184.812866)*T0(et) + 67310.54841;
 % Sidereal time in degrees
theta = temp/240;

% Sidereal time should be 0 - 360 deg, but as it's only used in trig functions
% here is no need to ensure this.
% theta = rem(theta, 360);

end



function [transVelOut, transAccOut] = calcVelAcc(time, orgPos, orgVel, ...
  orgAcc, transPos, sourceFrame, destFrame)

transVelOut = [];
transAccOut = [];
nPos        = size(time, 1);

% Systems with analytical conversion and flag for earth rotation correction
sysConv = {...
  'GEI'  , 'J2000', 0; ...
  'GEI'  , 'GEO'  , 1; ...
  'J2000', 'GEO'  , 1};
convDir = 0;
% Check if analytic conversion is available, forward
tmpInd = find(ismember(sysConv(:, 1), sourceFrame));
if (~isempty(tmpInd))
  tmpConv = sysConv(tmpInd, :);
  tmpInd  = find(ismember(tmpConv(:, 2), destFrame));
  if (~isempty(tmpInd))
    if (tmpConv{tmpInd, 3}), convDir = 1; end
  end
end
% Check if analytic conversion is available, backward
tmpInd = find(ismember(sysConv(:, 2), sourceFrame));
if (~isempty(tmpInd))
  tmpConv = sysConv(tmpInd, :);
  tmpInd  = find(ismember(tmpConv(:, 1), destFrame));
  if (~isempty(tmpInd))
    if (tmpConv{tmpInd, 3}), convDir = -1; end
  end
end
        
% Earth rotation
omegaearth = [0; 0; 7.29211514670698e-05];
if (convDir == 1)
  % Earth rotation must be adjusted for
  
  % Rotation + possible translatation of the velocities and acceleration
  transVel   = zeros(nPos, 3);
  transAcc   = zeros(nPos, 3);
  tmpVel     = xcframe(time, orgVel, [], [], sourceFrame, destFrame);
  if (~isempty(orgAcc))
    tmpAcc   = xcframe(time, orgAcc, [], [], sourceFrame, destFrame);
  end
 
  for i = 1 : nPos
    temp             = cross(omegaearth, transPos(i, :));
    transVel(i, :)   = tmpVel(i, :) - temp;
    if (~isempty(orgAcc))
      transAcc(i, :) = tmpAcc(i, :) - cross(omegaearth, temp) - ...
        2*cross(omegaearth, transVel(i, :));
    end
  end
  
elseif (convDir == -1)
  % Earth rotation must be adjusted for
  
  tmpVel = zeros(nPos, 3);
  tmpAcc = zeros(nPos, 3);
  for i = 1 : nPos
    tmpVel(i, :) = orgVel(i, :) + cross(omegaearth, orgPos(i, :));
    temp         = cross(omegaearth, orgPos(i, :));
    if (~isempty(orgAcc))
      tmpAcc(i, :) = orgAcc(i, :) + cross(omegaearth, temp) + ...
        2*cross(omegaearth, orgVel(i, :));
    end
  end
  transVel   = xcframe(time, tmpVel, [], [], sourceFrame, destFrame);
  if (~isempty(orgAcc))
    transAcc = xcframe(time, tmpAcc, [], [], sourceFrame, destFrame);
  end
 
else  
  
  transVel   = xcframe(time, orgVel, [], [], sourceFrame, destFrame);
  if (~isempty(orgAcc))
    transAcc = xcframe(time, orgAcc, [], [], sourceFrame, destFrame);
  end
end

% Output
transVelOut = transVel;
if (~isempty(orgAcc))
  transAccOut = transAcc;
end
  
end



function convCell = getConvCell(sourceFrame, destFrame)
% A more compact notation could probably be used, but this works. A smart
% automated system, as used for the old C code, to generate this would be nice
% if other transformations is added.

% Default output
convCell = [];

if (strcmpi(sourceFrame, destFrame)), return, end

invFlag = 0;
if (strcmpi(sourceFrame, 'GEI') || strcmpi(destFrame, 'GEI'))
  sourceConv = {...
    'J2000', {{@j2000_twixt_gei, 'b'}}; ...
    'GEO'  , {{@gei_twixt_geo  , 'f'}}; ...
    'MAG'  , {{@gei_twixt_geo  , 'f'}, {@geo_twixt_mag , 'f'}}; ...
    'GSE'  , {{@gei_twixt_gse  , 'f'}}
    'GSM'  , {{@gei_twixt_gse  , 'f'}, {@gse_twixt_gsm , 'f'}}; ...
    'SM'   , {{@gei_twixt_gse  , 'f'}, {@gse_twixt_gsm , 'f'} , ...
      {@gsm_twixt_sm , 'f'}}; ...
    'RTN'  , {{@gei_twixt_gse  , 'f'}, {@gse_twixt_rtn , 'f'}}; ...
    'GSEQ' , {{@gei_twixt_gse  , 'f'}, {@gse_twixt_gseq, 'f'}}; ...
    'HEE'  , {{@gei_twixt_gse  , 'f'}, {@gse_twixt_hee , 'f'}}; ...
    'HAE'  , {{@gei_twixt_gse  , 'f'}, {@gse_twixt_hee , 'f'},  ...
      {@hae_twixt_hee, 'b'}}; ...
    'HEEQ' , {{@gei_twixt_gse  , 'f'}, {@gse_twixt_hee , 'f'},  ...
      {@hae_twixt_hee, 'b'}, {@hae_twixt_heeq, 'f'}}};
  % If you need to check the whole transformation table, uncomment below!  
  % checkSourceConvCell('GEI', '', sourceConv, 'table')
  if (strcmpi(sourceFrame, 'GEI'))
    convCell = sourceConv{ismember(sourceConv(:, 1), destFrame), 2};
  else
    % Inverse transformation
    convCell = sourceConv{ismember(sourceConv(:, 1), sourceFrame), 2};
    invFlag  = 1;
  end

elseif (strcmpi(sourceFrame, 'J2000') || strcmpi(destFrame, 'J2000'))
  sourceConv = {...
    'GEO' , {{@j2000_twixt_gei, 'f'}, {@gei_twixt_geo, 'f'}}; ...
    'MAG' , {{@j2000_twixt_gei, 'f'}, {@gei_twixt_geo, 'f'}, ...
      {@geo_twixt_mag , 'f'}}; ...
    'GSE' , {{@j2000_twixt_gei, 'f'}, {@gei_twixt_gse, 'f'}}
    'GSM' , {{@j2000_twixt_gei, 'f'}, {@gei_twixt_gse, 'f'}, ...
      {@gse_twixt_gsm , 'f'}}; ...
    'SM'  , {{@j2000_twixt_gei, 'f'}, {@gei_twixt_gse, 'f'} , ...
      {@gse_twixt_gsm , 'f'}, {@gsm_twixt_sm , 'f'}}; ...
    'RTN' , {{@j2000_twixt_gei, 'f'}, {@gei_twixt_gse, 'f'}, ...
      {@gse_twixt_rtn , 'f'}}; ...
    'GSEQ', {{@j2000_twixt_gei, 'f'}, {@gei_twixt_gse, 'f'}, ...
      {@gse_twixt_gseq, 'f'}}; ...
    'HEE' , {{@j2000_twixt_gei, 'f'}, {@gei_twixt_gse, 'f'}, ...
      {@gse_twixt_hee , 'f'}}; ...
    'HAE' , {{@j2000_twixt_gei, 'f'}, {@gei_twixt_gse, 'f'},  ...
      {@gse_twixt_hee , 'f'}, {@hae_twixt_hee, 'b'}}; ...
    'HEEQ', {{@j2000_twixt_gei, 'f'}, {@gei_twixt_gse, 'f'},  ...
      {@gse_twixt_hee , 'f'}, {@hae_twixt_hee, 'b'}, {@hae_twixt_heeq, 'f'}}};
  % If you need to check the whole transformation table, uncomment below!  
  % checkSourceConvCell('J2000', '', sourceConv, 'table')
  if (strcmpi(sourceFrame, 'J2000'))
    convCell = sourceConv{ismember(sourceConv(:, 1), destFrame), 2};
  else
    % Inverse transformation
    convCell = sourceConv{ismember(sourceConv(:, 1), sourceFrame), 2};
    invFlag  = 1;
  end

elseif (strcmpi(sourceFrame, 'GEO') || strcmpi(destFrame, 'GEO'))
  sourceConv = {...
    'MAG' , {{@geo_twixt_mag, 'f'}}; ...
    'GSE' , {{@gei_twixt_geo, 'b'}, {@gei_twixt_gse, 'f'}}
    'GSM' , {{@gei_twixt_geo, 'b'}, {@gei_twixt_gse, 'f'}, ...
      {@gse_twixt_gsm , 'f'}}; ...
    'SM'  , {{@gei_twixt_geo, 'b'}, {@gei_twixt_gse, 'f'} , ...
      {@gse_twixt_gsm , 'f'}, {@gsm_twixt_sm , 'f'}}; ...
    'RTN' , {{@gei_twixt_geo, 'b'}, {@gei_twixt_gse, 'f'}, ...
      {@gse_twixt_rtn , 'f'}}; ...
    'GSEQ', {{@gei_twixt_geo, 'b'}, {@gei_twixt_gse, 'f'}, ...
      {@gse_twixt_gseq, 'f'}}; ...
    'HEE' , {{@gei_twixt_geo, 'b'}, {@gei_twixt_gse, 'f'}, ...
      {@gse_twixt_hee , 'f'}}; ...
    'HAE' , {{@gei_twixt_geo, 'b'}, {@gei_twixt_gse, 'f'},  ...
      {@gse_twixt_hee , 'f'}, {@hae_twixt_hee, 'b'}}; ...
    'HEEQ', {{@gei_twixt_geo, 'b'}, {@gei_twixt_gse, 'f'},  ...
      {@gse_twixt_hee , 'f'}, {@hae_twixt_hee, 'b'}, {@hae_twixt_heeq, 'f'}}};
  % If you need to check the whole transformation table, uncomment below!  
  % checkSourceConvCell('GEO', '', sourceConv, 'table')
  if (strcmpi(sourceFrame, 'GEO'))
    convCell = sourceConv{ismember(sourceConv(:, 1), destFrame), 2};
  else
    % Inverse transformation
    convCell = sourceConv{ismember(sourceConv(:, 1), sourceFrame), 2};
    invFlag  = 1;
  end

elseif (strcmpi(sourceFrame, 'MAG') || strcmpi(destFrame, 'MAG'))
  sourceConv = {...
    'GSE' , {{@geo_twixt_mag, 'b'}, {@gei_twixt_geo, 'b'}, ...
      {@gei_twixt_gse, 'f'}}
    'GSM' , {{@geo_twixt_mag, 'b'}, {@gei_twixt_geo, 'b'}, ...
      {@gei_twixt_gse , 'f'}, {@gse_twixt_gsm , 'f'}}; ...
    'SM'  , {{@geo_twixt_mag, 'b'}, {@gei_twixt_geo, 'b'} , ...
      {@gei_twixt_gse , 'f'}, {@gse_twixt_gsm, 'f'}, {@gsm_twixt_sm, 'f'}}; ...
    'RTN' , {{@geo_twixt_mag, 'b'}, {@gei_twixt_geo, 'b'}, ...
      {@gei_twixt_gse , 'f'}, {@gse_twixt_rtn , 'f'}}; ...
    'GSEQ', {{@geo_twixt_mag, 'b'}, {@gei_twixt_geo, 'b'}, ...
      {@gei_twixt_gse, 'f'}, {@gse_twixt_gseq, 'f'}}; ...
    'HEE' , {{@geo_twixt_mag, 'b'}, {@gei_twixt_geo, 'b'}, ...
      {@gei_twixt_gse , 'f'}, {@gse_twixt_hee , 'f'}}; ...
    'HAE' , {{@geo_twixt_mag, 'b'}, {@gei_twixt_geo, 'b'},  ...
      {@gei_twixt_gse , 'f'}, {@gse_twixt_hee, 'f'}, {@hae_twixt_hee, 'b'}}; ...
    'HEEQ', {{@geo_twixt_mag, 'b'}, {@gei_twixt_geo, 'b'},  ...
      {@gei_twixt_gse , 'f'}, {@gse_twixt_hee, 'f'}, ...
      {@hae_twixt_hee, 'b'}, {@hae_twixt_heeq, 'f'}}};
  % If you need to check the whole transformation table, uncomment below!  
  % checkSourceConvCell('MAG', '', sourceConv, 'table')
  if (strcmpi(sourceFrame, 'MAG'))
    convCell = sourceConv{ismember(sourceConv(:, 1), destFrame), 2};
  else
    % Inverse transformation
    convCell = sourceConv{ismember(sourceConv(:, 1), sourceFrame), 2};
    invFlag  = 1;
  end

elseif (strcmpi(sourceFrame, 'GSE') || strcmpi(destFrame, 'GSE'))
  sourceConv = {...
    'GSM' , {{@gse_twixt_gsm , 'f'}}; ...
    'SM'  , {{@gse_twixt_gsm , 'f'}, {@gsm_twixt_sm , 'f'}}; ...
    'RTN' , {{@gse_twixt_rtn , 'f'}}; ...
    'GSEQ', {{@gse_twixt_gseq, 'f'}}; ...
    'HEE' , {{@gse_twixt_hee , 'f'}}; ...
    'HAE' , {{@gse_twixt_hee , 'f'}, {@hae_twixt_hee, 'b'}}; ...
    'HEEQ', {{@gse_twixt_hee , 'f'}, {@hae_twixt_hee, 'b'},  ...
      {@hae_twixt_heeq , 'f'}}};
  % If you need to check the whole transformation table, uncomment below!  
  % checkSourceConvCell('GSE', '', sourceConv, 'table')
  if (strcmpi(sourceFrame, 'GSE'))
    convCell = sourceConv{ismember(sourceConv(:, 1), destFrame), 2};
  else
    % Inverse transformation
    convCell = sourceConv{ismember(sourceConv(:, 1), sourceFrame), 2};
    invFlag  = 1;
  end
  
elseif (strcmpi(sourceFrame, 'GSM') || strcmpi(destFrame, 'GSM'))
  sourceConv = {...
    'SM'  , {{@gsm_twixt_sm , 'f'}}; ...
    'RTN' , {{@gse_twixt_gsm, 'b'}, {@gse_twixt_rtn , 'f'}}; ...
    'GSEQ', {{@gse_twixt_gsm, 'b'}, {@gse_twixt_gseq, 'f'}}; ...
    'HEE' , {{@gse_twixt_gsm, 'b'}, {@gse_twixt_hee , 'f'}}; ...
    'HAE' , {{@gse_twixt_gsm, 'b'}, {@gse_twixt_hee , 'f'}, ...
       {@hae_twixt_hee , 'b'}}; ...
    'HEEQ', {{@gse_twixt_gsm , 'b'}, {@gse_twixt_hee , 'f'},  ...
      {@hae_twixt_hee , 'b'}, {@hae_twixt_heeq , 'f'}}};
  % If you need to check the whole transformation table, uncomment below!  
  % checkSourceConvCell('GSM', '', sourceConv, 'table')
  if (strcmpi(sourceFrame, 'GSM'))
    convCell = sourceConv{ismember(sourceConv(:, 1), destFrame), 2};
  else
    % Inverse transformation
    convCell = sourceConv{ismember(sourceConv(:, 1), sourceFrame), 2};
    invFlag  = 1;
  end
  
elseif (strcmpi(sourceFrame, 'SM') || strcmpi(destFrame, 'SM'))
  sourceConv = {...
    'RTN' , {{@gsm_twixt_sm , 'b'}, {@gse_twixt_gsm , 'b'}, ...
      {@gse_twixt_rtn , 'f'}}; ...
    'GSEQ', {{@gsm_twixt_sm , 'b'}, {@gse_twixt_gsm , 'b'}, ...
      {@gse_twixt_gseq , 'f'}}; ...
    'HEE' , {{@gsm_twixt_sm , 'b'}, {@gse_twixt_gsm , 'b'}, ...
      {@gse_twixt_hee , 'f'}}; ...
    'HAE' , {{@gsm_twixt_sm , 'b'}, {@gse_twixt_gsm , 'b'}, ...
      {@gse_twixt_hee , 'f'}, {@hae_twixt_hee , 'b'}}; ...
    'HEEQ', {{@gsm_twixt_sm , 'b'}, {@gse_twixt_gsm , 'b'}, ...
      {@gse_twixt_hee , 'f'}, {@hae_twixt_hee , 'b'}, {@hae_twixt_heeq, 'f'}}};
  % If you need to check the whole transformation table, uncomment below!  
  % checkSourceConvCell('SM', '', sourceConv, 'table')
  if (strcmpi(sourceFrame, 'SM'))
    convCell = sourceConv{ismember(sourceConv(:, 1), destFrame), 2};
  else
    % Inverse transformation
    convCell = sourceConv{ismember(sourceConv(:, 1), sourceFrame), 2};
    invFlag  = 1;
  end
  
elseif (strcmpi(sourceFrame, 'RTN') || strcmpi(destFrame, 'RTN'))
  sourceConv = {...
    'GSEQ', {{@gse_twixt_rtn, 'b'}, {@gse_twixt_gseq, 'f'}}; ...
    'HEE' , {{@gse_twixt_rtn, 'b'}, {@gse_twixt_hee , 'f'}}; ...
    'HAE' , {{@gse_twixt_rtn, 'b'}, {@gse_twixt_hee , 'f'}, ...
      {@hae_twixt_hee , 'b'}}; ...
    'HEEQ', {{@gse_twixt_rtn, 'b'}, {@gse_twixt_hee , 'f'}, ...
      {@hae_twixt_hee , 'b'}, {@hae_twixt_heeq, 'f'}}};
  % If you need to check the whole transformation table, uncomment below!  
  % checkSourceConvCell('RTN', '', sourceConv, 'table')
  if (strcmpi(sourceFrame, 'RTN'))
    convCell = sourceConv{ismember(sourceConv(:, 1), destFrame), 2};
  else
    % Inverse transformation
    convCell = sourceConv{ismember(sourceConv(:, 1), sourceFrame), 2};
    invFlag  = 1;
  end
  
elseif (strcmpi(sourceFrame, 'GSEQ') || strcmpi(destFrame, 'GSEQ'))
  sourceConv = {...
    'HEE' , {{@gse_twixt_gseq, 'b'}, {@gse_twixt_hee , 'f'}}; ...
    'HAE' , {{@gse_twixt_gseq, 'b'}, {@gse_twixt_hee , 'f'}, ...
      {@hae_twixt_hee , 'b'}}; ...
    'HEEQ', {{@gse_twixt_gseq, 'b'}, {@gse_twixt_hee , 'f'}, ...
      {@hae_twixt_hee , 'b'}, {@hae_twixt_heeq, 'f'}}};
  % If you need to check the whole transformation table, uncomment below!  
  % checkSourceConvCell('GSEQ', '', sourceConv, 'table')
  if (strcmpi(sourceFrame, 'GSEQ'))
    convCell = sourceConv{ismember(sourceConv(:, 1), destFrame), 2};
  else
    % Inverse transformation
    convCell = sourceConv{ismember(sourceConv(:, 1), sourceFrame), 2};
    invFlag  = 1;
  end
  
elseif (strcmpi(sourceFrame, 'HEE') || strcmpi(destFrame, 'HEE'))
  sourceConv = {...
    'HAE' , {{@hae_twixt_hee, 'b'}}; ...
    'HEEQ', {{@hae_twixt_hee, 'b'}, {@hae_twixt_heeq , 'f'}}};
  % If you need to check the whole transformation table, uncomment below!  
  % checkSourceConvCell('HEE', '', sourceConv, 'table')
  if (strcmpi(sourceFrame, 'HEE'))
    convCell = sourceConv{ismember(sourceConv(:, 1), destFrame), 2};
  else
    % Inverse transformation
    convCell = sourceConv{ismember(sourceConv(:, 1), sourceFrame), 2};
    invFlag  = 1;
  end
  
elseif (strcmpi(sourceFrame, 'HAE') || strcmpi(destFrame, 'HAE'))
  sourceConv = {...
    'HEEQ', {{@hae_twixt_heeq, 'f'}}};
  % If you need to check the whole transformation table, uncomment below!  
  % checkSourceConvCell('HAE', '', sourceConv, 'table')
  if (strcmpi(sourceFrame, 'HAE'))
    convCell = sourceConv{ismember(sourceConv(:, 1), destFrame), 2};
  else
    % Inverse transformation
    convCell = sourceConv{ismember(sourceConv(:, 1), sourceFrame), 2};
    invFlag  = 1;
  end
  
else
  % This should not happen!
  error('Failed to get conversion specification!')
end

if (invFlag)
  % Inverse transformation
  convCell = convCell(end:-1:1);
  for i = 1 : length(convCell)
    if (convCell{i}{2} == 'f')
      convCell{i}{2} = 'b';
    else
      convCell{i}{2} = 'f';
    end
  end
end

% Check conversion specification
checkSourceConvCell(sourceFrame, destFrame, convCell, 'single')

end



function checkSourceConvCell(sourceFrame, destFrame, sourceConv, mode)
% Help function to check if the conversion cells are correct!

frames = {'J2000', 'GEI', 'GEO', 'MAG', 'GSE', 'GSM', 'SM', 'RTN', 'GSEQ', ...
  'HEE', 'HAE', 'HEEQ'};
errID  = 'xcframe:conversionSpecification';

switch mode
  
  case 'single'
    currFrame = sourceFrame;
    nSteps    = length(sourceConv);
    for stepNo = 1 : nSteps
      errStr    = ['Conversion cell from ', sourceFrame, ' to ', ...
        destFrame, ' is incorrect (step ', num2str(stepNo), ')! Please ', ...
        'check and correct the code'];
      stepCell  = sourceConv{stepNo};
      try
        stepStr = func2str(stepCell{1});
      catch
        error(errID, errStr)
      end
      stepStart = regexp(stepStr, '^([a-zA-Z0-9])*?_', 'tokens');
      if (isempty(stepStart)), error(errID, errStr), end
      stepStart = stepStart{1}{1};
      stepStop  = regexp(stepStr, '_([a-zA-Z0-9])*?$', 'tokens');
      if (isempty(stepStop)), error(errID, errStr), end
      stepStop  = stepStop{1}{1};
      if (~ismember(upper(stepStart), frames) || ...
          ~ismember(upper(stepStop ), frames))
        error(errID, errStr)
      end
      dirStr = strtrim(stepCell{2});
      if (~strcmpi(dirStr, 'f') && ~strcmpi(dirStr, 'b'))
        error(errID, errStr)
      end
      if (strcmpi(dirStr, 'b'))
        tmpStr    = stepStart;
        stepStart = stepStop;
        stepStop  = tmpStr;
      end
      if (~strcmpi(stepStart, currFrame))
        error(errID, errStr)
      end
      currFrame = stepStop;
    end
    if (~strcmpi(currFrame, destFrame))
      error(errID, errStr)
    end
    
  case 'table'
    nConv = size(sourceConv, 1);
    for convNo = 1 : nConv
      destFrame = sourceConv{convNo, 1};
      convRow   = sourceConv{convNo, 2};
      currFrame = sourceFrame;
      nSteps    = length(convRow);
      for stepNo = 1 : nSteps
        errStr    = ['Conversion cell from ', sourceFrame, ' to ', ...
          destFrame, ' is incorrect (step ', num2str(stepNo), ')! Please ', ...
          'check and correct the code'];
        stepCell  = convRow{stepNo};
        try
          stepStr = func2str(stepCell{1});
        catch
          error(errID, errStr)
        end
        stepStart = regexp(stepStr, '^([a-zA-Z0-9])*?_', 'tokens');
        if (isempty(stepStart)), error(errID, errStr), end
        stepStart = stepStart{1}{1};
        stepStop  = regexp(stepStr, '_([a-zA-Z0-9])*?$', 'tokens');
        if (isempty(stepStop)), error(errID, errStr), end
        stepStop  = stepStop{1}{1};
        if (~ismember(upper(stepStart), frames) || ...
            ~ismember(upper(stepStop ), frames))
          error(errID, errStr)
        end
        dirStr = strtrim(stepCell{2});
        if (~strcmpi(dirStr, 'f') && ~strcmpi(dirStr, 'b'))
          error(errID, errStr)
        end
        if (strcmpi(dirStr, 'b'))
          tmpStr    = stepStart;
          stepStart = stepStop;
          stepStop  = tmpStr;
        end
        if (~strcmpi(stepStart, currFrame))
          error(errID, errStr)
        end
        currFrame = stepStop;
      end
      if (~strcmpi(currFrame, destFrame))
        error(errID, errStr)
      end
    end
    
end

end



%--------------------BELOW IS A MATLAB IMPLEMENTATION OF CXFORM-----------------
% MAIN CHANGES FROM CXFORM
% * Replaced "date2es" with call to MATLAB function "juliandate"
% * New code for calculation of Greenwich mean sidereal time, see above.
% * Added  12th generation IGRF coefficients and extrapolation is done after
%   2015 based on predicted change. A warning is issued for extrapolation
%   before/after current IGRF time range and if the IGRF constants should be
%   updated

% cxform-manual.c  --  manually coded functions for coordinate transforms
%
% This file is part of Ed's "cxform" Coordinate Transformation package.
% It contains the hand-coded functions for converting between some
% coordinate systems.  
%
% All code in here is derived from Mike Hapgood <M.Hapgood@rl.ac.uk>'s
% excellent introduction to coordinate transformations, along with
% Markus Fraenz' and Christopher Russell's work:
%
%	http://sspg1.bnsc.rl.ac.uk/Share/Coordinates/ct_home.htm
%      http://www.space-plasma.qmul.ac.uk/heliocoords/
%      http://www-ssc.igpp.ucla.edu/personnel/russell/papers/gct1.html/
%
% Converted to C by Ed Santiago <esm@pobox.com>
% Modified & maintained by Ryan Boller <Ryan.A.Boller@nasa.gov>
%
% Version history:
%   2000/06/21  v0.2  Ed Santiago:  Last released version from Ed.
%   2003/09/12  v0.3  Ryan Boller:  First modified version from Ryan.  
%                     Added RTN and GSEQ systems, IGRF model, and slightly
%                     different time manipulation.
%   2004/03/19  v0.4  Ryan Boller: Updated Makefile to auto-detect platform and
%                     to build under Mac OS 0.  HEEQ system now implemented by
%                     Kristi Keller.
%   2004/05/21  v0.5  Ryan Boller: Corrected small discrepancy when calculating
%				       T0 and lambda0.  Results now match SSCWeb's output very
%					   closely when magnetic lat & lon are forced to their
%                     values (they haven't updated their IGRF coefficients).
%   2004/11/25  v0.6  Ryan Boller: Moved utility functions from main.c into 
%                     shared library for convenience;  updated IGRF coefficients
%                     to the 9th generation which adds definitive values for
%                     1995 & 2000
%	 2006/09/24  v0.7  Ryan Boller: Updated IGRF coefficients to 10th generation


% /%%%%%%%%%*\
% |* hapgood_matrix  *|  defines a rotation matrix for a given angle & axis
% \%%%%%%%%%*|
%  *
%  * Rotation matrices are a special case.  They can be defined by two 
%  * parameters: an axis about which to rotate (0, 1, 2) and an angle.
%  * Given those two, we can fill in all nine elements of a 3x3 matrix.
%  *
%  * See http://sspg1.bnsc.rl.ac.uk/Share/Coordinates/matrix.htm
%  %
function mat = hapgood_matrix(theta, axis)

% 1.calculate sin(zeta) and cos(zeta)
sin_theta = sind(theta);
cos_theta = cosd(theta);

% compute the indices for the other two axes (e.g., "0,2" for 1) %
t1 = mod((axis + 1), 3);
t2 = mod((axis + 2), 3);
if (t1 > t2)
  tmp = t1;
  t1  = t2;
  t2  = tmp;
end

% 4.set the remaining off-diagonal terms to zero.
mat = zeros(3, 3);

% 2.determine the matrix diagonal:
%   1.put 1 in the Nth term, where N=1 if the rotation axis is 0, etc
%   2.put cos(zeta) in the other two terms
mat(axis+1, axis+1) = 1;
mat(t1+1, t1+1)     = cos_theta;
mat(t2+1, t2+1)     = cos_theta;

% 3.locate the two off-diagonal terms in the same columns and rows as
%   the cos(zeta) terms - put sin(zeta) in the term above the diagonal
%   and -sin(zeta) in the term below,
mat(t1+1, t2+1)     =  sin_theta;
mat(t2+1, t1+1)     = -sin_theta;

end


% /%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*\
% |*                                                                          *|
% |*                          ELEMENTAL FUNCTIONS                             *|
% |*                                                                          *|
% |*  This section provides functions required by the actual transformations  *|
% |*                                                                          *|
% \%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% MJD
%
% Modified Julian Date
%
function jd = MJD(et)

  jd = (et/86400) - 0.5;

  jd = jd + 51545;
end


%
% T0
%
% Julian Centuries from a certain time to 1 Jan 2000 12:00
%
% Note that this version takes into account hours, minutes, and seconds
%
function jc = T0(et)

% %  Old method
% jd = (et/86400) - 0.5;
% jd = fix(jd) - 0.5;
% jc = jd/36525;

%  Seconds --> days --> centuries  %
jc = (et/86400)/36525;

end


%
% lambda0
%
% The Sun's ecliptic longitude (lambdaO) can be calculated using the 
% series of formulae:
%
%     M = 357.528 + 35999.050T0 + 0.04107H degrees 
%     Lambda = 280.460 + 36000.772T0 + 0.04107H degrees 
%     lambdaO = Lambda + (1.915 - 0.0048T0) sinM + 0.020 sin2M
%
% where T0 is the time in Julian centuries from 12:00 UT on 1 January 2000
% to the midnight Universal Time (UT) preceding the time of interest and 
% H is the time in hours since that preceding UT midnight. Formulae 
% derived from the Almanac for Computers. In the intermediate formulae, 
% M is the Sun's mean anomaly and Lambda its mean longitude.
%
function eclng = lambda0(et)

M      = 357.528 + 35999.050*T0(et);  % + 0.04107 * H(et); %
lambda = 280.460 + 36000.772*T0(et);  % + 0.04107 * H(et); %

eclng  = lambda + (1.915 - 0.0048*T0(et))*sind(M) + 0.020*sind(2*M);

end


%
% The obliquity of the ecliptic (epsilon) can be calculated using the formula:
%
%	epsilon = 23.439 - 0.013T0 degrees
%
% where T0 is the time in Julian centuries from 12:00 UT on 1 January 2000
% to the midnight (UT) preceding the time of interest. Formula derived from
% the Almanac for Computers.
%
function ecob = epsilon(et)
ecob = 23.439 - 0.013 * T0(et);
end


%
% The following functions calculate the first three coefficients of the IGRF,
% used in determining the location of the magnetic dipole.  The values
% are interpolated to the day of interest.
% 
% Inputs:
%   fracYearIndex:  The double-precision "array index" of the appropriate
%                   IGRF year [0, 23], corresponding to
%                   the years 1900-2005.  This can be found by converting et to
%                   seconds past 1900, converting to years, and dividing by 5
%                   to get the appropriate array index.  I.e.,
%                   fracYearIndex = (et+3155803200.0) / 157788000.0
%   fracYear:  The remainder of the fracYearIndex [0, 1] for use during
%              interpolation.  
%
% The IGRF/DGRF coefficients are taken from the IGRF2000 (8th gen) model:
%   http://www.ngdc.noaa.gov/IAGA/wg8/igrf2000.html
%
% Note: the coefficients are now updated to the IGRF 9th generation model:
%   http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
%
% Note: IGRF coefficients are now updated to 10th generation
%
% Note: IGRF coefficients are now updated to 12th generation
%

function [G01, warnFlag] = calcG01(fracYearIndex, fracYear)

warnFlag = 0;
g01 = ...
  [-31543, -31464, -31354, -31212, -31060, -30926, -30805, -30715, ...
  -30654, -30594, -30554, -30500, -30421, -30334, -30220, -30100, ...
  -29992, -29873, -29775, -29692, -29619.4, -29554.63, -29496.57, -29442.0];
vs  = 10.3;

% START ADDED CODE
if (floor(fracYearIndex)+1 < 1)
  % Just use the first available value
  G01      = g01(1);
  warnFlag = -1;     % Extraploation done before
elseif (ceil(fracYearIndex)+1 > length(g01))
  % Add predicted change to last value
  G01 = g01(end) + (fracYearIndex + 1 - length(g01))*vs;
  if (fracYearIndex + 1 - length(g01) > 1.1)
    warnFlag = 1;    % Extraploation done more than 5 years after
  end
  % END ADDED CODE
else
  G01 = g01(floor(fracYearIndex)+1)*(1-fracYear) + ...
    g01(ceil(fracYearIndex)+1)*fracYear;
end

end



function G11 = calcG11(fracYearIndex, fracYear)

g11 = ...
  [-2298, -2298, -2297, -2306, -2317, -2318, -2316, -2306, -2292, -2285, ...
  -2250, -2215, -2169, -2119, -2068, -2013, -1956, -1905, -1848, -1784, ...
  -1728.2, -1669.05, -1586.42, -1501.0];
vs  = 18.1;

% START ADDED CODE
if (floor(fracYearIndex)+1 < 1)
  % Just use the first available value
  G11 = g11(1);
elseif (ceil(fracYearIndex)+1 > length(g11))
  % Add predicted change to last value
  G11 = g11(end) + (fracYearIndex + 1 - length(g11))*vs;
  % END ADDED CODE
else
  G11 = g11(floor(fracYearIndex)+1)*(1-fracYear) + ...
    g11(ceil(fracYearIndex)+1)*fracYear;
end

end



function H11 = calcH11(fracYearIndex, fracYear)

h11 = ...
  [5922, 5909, 5898, 5875, 5845, 5817, 5808, 5812, 5821, 5810, 5815, ...
  5820, 5791, 5776, 5737, 5675, 5604, 5500, 5406, 5306, 5186.1, 5077.99, ...
  4944.26, 4797.1];
vs = -26.6;

% START ADDED CODE
if (floor(fracYearIndex)+1 < 1)
  % Just use the first available value
  H11 = h11(1);
elseif (ceil(fracYearIndex)+1 > length(h11))
  % Add predicted change to last value
  H11 = h11(end) + (fracYearIndex + 1 - length(h11))*vs;
  % END ADDED CODE
else
  H11 = h11(floor(fracYearIndex)+1)*(1-fracYear) + ...
    h11(ceil(fracYearIndex)+1)*fracYear;
end

end


%
% Latitude and longitude of Earth's magnetic pole using above coefficients
%
function lambda0 = mag_lon(et)

fracYearIndex = (et + 3155803200)/157788000;
fracYear      = rem(fracYearIndex, 1);

% if (fracYearIndex + 1 >= 24)
%   error('Specified year is greater than IGRF implementation.');
% end

g11 = calcG11(fracYearIndex, fracYear);
h11 = calcH11(fracYearIndex, fracYear);

% lambda0 / longitude  %
lambda0 = atan2(h11, g11) + pi;
% lambda0 = (360.0-71.115)*(pi/180.0);  // SSC year 1990 value  %
% lambda0 = (360.0-71.381)*(pi/180.0);  // SSC year 1995 value  %
% lambda0 = (360.0-71.647)*(pi/180.0);  // SSC year 2000 value  %

end



function phi0 = mag_lat(et)

fracYearIndex = (et + 3155803200)/157788000;
fracYear      = rem(fracYearIndex, 1);

% if (fracYearIndex + 1 >= 24)
%   error('Specified year is greater than IGRF implementation.');
% end

g01 = calcG01(fracYearIndex, fracYear);
g11 = calcG11(fracYearIndex, fracYear);
h11 = calcH11(fracYearIndex, fracYear);
lambda = mag_lon(et);

%  phi0 / latitude  %
phi0 = (pi/2) - atan((g11*cos(lambda) + h11*sin(lambda))/g01);
% phi0 = (90.0-10.872)*(pi/180.0);   // SSC year 1990 value  %
% phi0 = (90.0-10.730)*(pi/180.0);   // SSC year 1995 value  %
% phi0 = (90.0-10.588)*(pi/180.0);   // SSC year 2000 value  %

end

% /%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*\
% |*                                                                          *|
% |*                    TRANSFORMATION MATRICES                               *|
% |*                                                                          *|
% |* Hapgood defines all his transformations in terms of various matrices T1, *|
% |* T2, ... Tn.  Since they're used in various places throughout the paper,  *|
% |* it makes sense to define them as subroutines.                            *|
% |*                                                                          *|
% |* Comments preceding each function are extracted verbatim from Hapgood.    *|
% |*                                                                          *|
% \%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% The GEI2000 to GEI transformation is given by the matrix 
%
%	P = <-zA,2>*<thetaA,1>*<-zetaA,2>
%
% where the rotation angles zA, thetaA and zetaA are the precession 
% angles. This transformation is a precession correction as described
% by Hapgood (1995).
%
%     zA     = 0.64062 T0 + 0.00030 T0^2 degrees 
%     thetaA = 0.55675 T0 - 0.00012 T0^2 degrees 
%     zetaA  = 0.64062 T0 + 0.00008 T0^2 degrees
%
function mat = mat_P(et)
t0      = T0(et);

mat     = hapgood_matrix(-1*(0.64062*t0 + 0.00030*t0*t0), 2);

mat_tmp = hapgood_matrix(0.55675*t0 - 0.00012*t0*t0, 1);
mat     = mat*mat_tmp;

mat_tmp = hapgood_matrix(-1*(0.64062*t0 + 0.00008*t0*t0), 2);
mat     = mat*mat_tmp;
end


%
% The GEI to GEO transformation is given by the matrix
%
%	T1 = <theta,2>
%
% where the rotation angle theta is the Greenwich mean sidereal time. This
% transformation is a rotation in the plane of the Earth's equator from 
% the First Point of Aries to the Greenwich meridian. 
%---
% The Greenwich Mean Sidereal Time (theta) can be calculated using the formula:
%
%  theta = 280.46061837 + 360.98564736629*d0 + 0.0003875*T0^2 - 
%          0.0000000258333*T0^3
%
%  where d0 is Julian days from J2000 and T0 is Julian centuries from J2000.
%  The forumla can be found in:
%    http://www.space-plasma.qmul.ac.uk/heliocoords/systems2art/
%    (Heliospheric Coordinate Systems, M. Franz, D. Harper)
%
%
% GST was formerly calculated as: 
% 
%	theta = 100.461 + 360000.770T0 + 15.04107H degrees
%
% where T0 is the time in Julian centuries from 12:00 UT on 1 January 2000
% to the midnight Universal Time (UT) preceding the time of interest and 
% H is the time in hours since that preceding UT midnight. Formula derived
% from the Almanac for Computers.
%

function mat = mat_T1(et)
% Original code
% theta = 100.461 + 36000.770*T0(et) + 360*(H(et)/24);  % + 15.04107 * H(et); %
% 
% % theta = rem(theta, 360);
% % if (theta < 0)
% %   theta = theta + 360;
% % end
theta = gmst(et);

mat   = hapgood_matrix(theta, 2);
end


%
% The GEI to GSE transformation is given by the matrix 
%
%	T2 = <lambdaO,2>*<epsilon,0>
%
% where the rotation angle lambdaO is the Sun's ecliptic longitude and 
% the angle epsilon is the obliquity of the ecliptic. This transformation 
% is a rotation from the Earth's equator to the plane of the ecliptic
% followed by a rotation in the plane of the ecliptic from the First Point 
% of Aries to the Earth-Sun direction. 
%
function  mat = mat_T2(et)
mat     = hapgood_matrix(lambda0(et), 2);

mat_tmp = hapgood_matrix(epsilon(et), 0);
mat     = mat*mat_tmp;
end


%
% vec_Qe
%
% don't ask.
%
function Qe = vec_Qe(et)

lat = mag_lat(et);
lon = mag_lon(et);

cos_lat = cos(lat);
sin_lat = sin(lat);
cos_lon = cos(lon);
sin_lon = sin(lon);

Qg(1) = cos_lat*cos_lon;
Qg(2) = cos_lat*sin_lon;
Qg(3) = sin_lat;

mat     = mat_T2(et);
mat_tmp = mat_T1(et);
mat_tmp = transpose(mat_tmp);
mat     = mat*mat_tmp;
Qe      = mat*Qg(:);

end


%
% The GSE to GSM transformation is given by the matrix
%
%	T3 = <-psi,0>
%
% where the rotation angle psi is the the GSE-GSM angle. This
% transformation is a rotation in the GSE YZ plane from the GSE 2 axis
% to the GSM 2 axis. 
%
function mat = mat_T3(et)
Qe  = vec_Qe(et);

psi = atan2d(Qe(2), Qe(3));

mat = hapgood_matrix(-psi, 0);
end


%
% The GSM to SM transformation is given by the matrix
%	T4 = <- mu,1>
%
% where the rotation angle mu is the dipole tilt. This transformation 
% is a rotation in the GSM XZ plane from the GSM 2 axis to the 
% geomagnetic dipole axis. 
%
function mat = mat_T4(et)
Qe = vec_Qe(et);

mu = atan2(Qe(1), sqrt(Qe(2)*Qe(2) + Qe(3)*Qe(3)))/(pi/180);

mat = hapgood_matrix(-mu, 1);
end


%
% The GEO to MAG transformation is given by the matrix 
%
%	T5 = <lat-90,1>*<long,2>
%
% where the rotation angle lat is the latitude and angle long is the 
% longitude of the geomagnetic pole (as defined by the axis of the 
% dipole component of the geomagnetic field). This transformation is 
% a rotation in the plane of the Earth's equator from the Greenwich 
% meridian to the meridian containing the dipole axis, followed by 
% a rotation in that meridian from the rotation axis to the dipole axis. 
%
function mat = mat_T5(et)
mat     = hapgood_matrix((mag_lat(et)*(180/pi)) - 90, 1);

mat_tmp = hapgood_matrix((mag_lon(et)*(180/pi))     , 2);

mat     = mat*mat_tmp;
end

%
% The GSE to GSEQ transformation is given by the matrix
%
%  T6 = <theta, 0>
%
% where theta is the angle between the 1-axes in the two systems. A full
% description can be found at
%  http://www-ssc.igpp.ucla.edu/personnel/russell/papers/gct1.html/#s3.5
%  (Geophysical Coordinate Transformations, C. T. Russell 1971)
%
function mat = mat_T6(et)
%  Get Earth-Sun vector in GEI  %

GSE_ES(1) = 1;
GSE_ES(2) = 0;
GSE_ES(3) = 0;

%  Convert GSE --> GEI  %
matT2  = mat_T2(et);
matT2  = transpose(matT2);
GEI_ES = matT2*GSE_ES(:);

%  Rotation axis of the Sun (GEI): (1.217,- 0.424, 0.897)  %
thetaN    = GEI_ES(1)*(-0.032) + GEI_ES(2)*(-0.112) + GEI_ES(3)*(-0.048);
thetaD(1) = (-0.424)*GEI_ES(3) - 0.897*GEI_ES(2);
thetaD(2) = 0.897*GEI_ES(1) - 0.1217*GEI_ES(3);
thetaD(3) = 0.1217*GEI_ES(2) - (-0.424)*GEI_ES(1);
magThetaD = sqrt(power(thetaD(1), 2) + power(thetaD(2), 2) + ...
  power(thetaD(3), 2));

theta = asin(thetaN/magThetaD);

mat = hapgood_matrix((theta*(180/pi)), 0);

%  TODO: Unknown why transpose is necessary to match previous results  %
mat = transpose(mat);

end


% /%%%%%%%  Heliocentric transformations %%%%%%/

%
% The HAE to HEE transformation is given by the matrix
%
%	S1 = <lambdaO + 180,2>
%
% where the rotation angle lambdaO is the Sun's ecliptic longitude.  
% This transformation is a rotation in the plane of the ecliptic from 
% the First Point of Aries to the Sun-Earth direction. 
%
function mat = mat_S1(et)
mat = hapgood_matrix(lambda0(et) + 180, 2);
end


%
% The HAE to HEEQ transformation is given by the matrix
%
% 	S2 = <theta0,2>*<i,0>*<Omega,2>
%
% where the rotation angle theta0 is the the longitude of the Sun's
% central meridian, i is the the inclination of the Sun's equator and
% Omega is the the ecliptic longitude of the ascending node of the Sun's
% equator. This transformation comprises a rotation in the plane of the
% ecliptic from the First Point of Aries to the ascending node of the
% solar equator, then a rotation from the plane of the ecliptic to the 
% plane of the equator and finally a rotation in the plane of the solar
% equator from the ascending node to the central meridian.
%
% Implemented by Kristi Keller on 2/2/2004
%
function mat = mat_S2(et)

Omega   = 73.6667+0.013958*(MJD(et)+3242)/365.25;
theta0  = atand(cosd(7.25)*tand(lambda0(et)-Omega));
angle_1 = lambda0(et) - Omega;
angle_1 = rem(angle_1, 360);
if (angle_1 < 0), angle_1 = angle_1 + 360; end

theta0 = rem(theta0, 360);
if (theta0 < 0),     theta0 = theta0 + 360; end
if (angle_1 < 180)
  if (theta0 < 180), theta0 = theta0 + 180; end
end
if (angle_1 > 180)
  if (theta0 > 180), theta0 = theta0 - 180; end
end
mat = hapgood_matrix(theta0, 2);

mat_tmp = hapgood_matrix(7.25, 0);
mat     = mat*mat_tmp;

mat_tmp = hapgood_matrix(Omega, 2);
mat     = mat*mat_tmp;

end



% /%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*\
% |*                                                                          *|
% |*                    TRANSFORMATIONS BEGIN HERE                            *|
% |*                                                                          *|
% |*                                                                          *|
% \%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%*\
%  simple_rotation  %  utility function used by all the "twixt" functions
%%%%%%%%%%%
%
% This is basically what all the "twixt" functions do:
%
%    1) define a rotation matrix
%    2) If doing an inverse transformation, transpose that matrix.
%    3) multiply the rotation matrix by the input vector.
%
% To save all that work in the "twixt" functions, they just call this
% function, passing us a pointer to a function that defines the matrix.
%
function v_out = simple_rotation(et, v_in, d, tmatFun)

%
% Call the user-specified function to get a rotation matrix.
%
mat = tmatFun(et);

%
% To do the inverse transformation, we use the transposition of the matrix
%
if (strcmpi(d, 'b'))
  mat = transpose(mat);
end

%
% Multiply the rotation matrix by the input vector, and that's it!
%
v_out = mat*v_in(:);

end


%  Bo-ring.  Just call simple_rotation() with the appropriate matrix. %
function v_out = j2000_twixt_gei(et, v_in, direction)
v_out = simple_rotation(et, v_in, direction, @mat_P);
end

function v_out = gei_twixt_geo(et, v_in, direction)
v_out = simple_rotation(et, v_in, direction, @mat_T1);
end

function v_out = geo_twixt_mag(et, v_in, direction)
v_out = simple_rotation(et, v_in, direction, @mat_T5);
end

function v_out = gei_twixt_gse(et, v_in, direction)
v_out = simple_rotation(et, v_in, direction, @mat_T2);
end

function v_out = gse_twixt_gsm(et, v_in, direction)
v_out = simple_rotation(et, v_in, direction, @mat_T3);
end

function v_out = gsm_twixt_sm(et, v_in, direction)
v_out = simple_rotation(et, v_in, direction, @mat_T4);
end


%
% Define conversion from GSE to Earth-centered RTN
%
% NOTE:  This coordinate system has not been tested and its results should be
% treated as such
%
function v_out = gse_twixt_rtn(et, v_in, ~)

%  Convert time  %
dYear = (et/31557600) + 2000.0013689254;
iYear = floor(dYear);
doy   = (rem(dYear, 1)*365.25) + 1;

%  Longitude of the ascending node for the year  %
alpha = 74.3666667 + 0.014*(iYear - 1900);

%  Angle between fall equinox and alpha  %
phi = (pi/180)*(alpha - (360*(doy-266)/365.25));

%  Rotation angle around GSE 0 / RTN R to change coordinates:  %
delta = atan(tan((pi/180)*7.25)*cos(phi));

%  Rotation is the same, regardless of direction  %
v_out(1) = -v_in(1);
v_out(2) = (-v_in(2))*cos(delta) + v_in(3)*sin(delta);
v_out(3) = v_in(2)*sin(delta) + v_in(3)*cos(delta);

end


function v_out = gse_twixt_gseq(et, v_in, direction)
v_out = simple_rotation(et, v_in, direction, @mat_T6);
end



% /%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*\
% |*              HELIOCENTRIC COORDINATE SYSTEMS                             *|
% \%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Hapgood defines a transformation between GSE and HEE in his 1992
% paper (section 6), but this part isn't online.  
%
% The gist of it is, we rotate 180 degrees about 2, and then translate
% along 0.
%
% But we also need to add "R", a constant vector defined by
%
%      R = [ Rsun, 0, 0 ]
%
% where 
%
%             r0 (1 - e^2)
%    Rsun =   ------------
%             1 + e cos(v)
%
%   r0 = 1.495985E8 km        	mean distance of the Sun from Earth.
%
%    e = 0.016709 - 0.0000418T0	eccentricity of the Sun's apparent
%					orbit around the Earth.
%
%    w = 282.94 + 1.72 T0		longitude of perigee of that orbit
%
%    v = lambda0 - w			(see lambda0 above)
%
%
% Implemented by Ed Santiago, Updated by Kristi Keller
%
function v_out = gse_twixt_hee(et, v_in, ~)

mat = hapgood_matrix(180, 2);

%
% Note that there's no transposition here if the direction is "back";
% the operation works identically in both directions.
%
v_out = mat*v_in(:);

% Translate the 0 axis about the earth-sun distance %
r0   = 1.495985e8;
e    = 0.016709 - 0.0000418*T0(et);
w    = 282.94 + 1.72*T0(et);
v    = lambda0(et) - w;
Rsun = r0*(1 - e*e)/(1 + e*cos(v*(pi/180)));
%  v_out[0] += (double)1.5e8;  %

v_out(1) = v_out(1) + Rsun;

end


function v_out = hae_twixt_hee(et, v_in, direction)
v_out = simple_rotation(et, v_in, direction, @mat_S1);
end

function v_out = hae_twixt_heeq(et, v_in, direction)
v_out = simple_rotation(et, v_in, direction, @mat_S2);
end
%------------------------------END CXFORM CODE----------------------------------