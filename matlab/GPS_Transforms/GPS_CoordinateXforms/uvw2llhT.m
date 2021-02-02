function llh = uvw2llhT ( uvw, origin )
%
% Description:
%   This function will convert a UVW coordinates vector to geodetic
%   LLH coordinates (Longitude, Latitude, Height)
%
%   Given a three-dimensional vector, uvw, represented in the Joint Stars
%   Geocentric system described below, uvw2llh calculates the corresponding
%   longitude, latitude (radians), and the height (meters).
%
% NOTE: This routine uses an iterative method; to compute the exact solution
%   for the spherical earth model, see the routine suvw2ecf.
%
% Usage:
%   llh=uvw2llhT(uvw,origin)
%
% Inputs:
%   origin - 3xN or 1xN vector [lat; lon; hgt] (radians,radians,meters)
%   uvw    - 3xN Vector, [x; y; z]
%
% Outputs:
%   llh   - 1xN Structure such that
%             llh.lon - longitude  (rad)
%             llh.lat - latitude   (rad)
%             llh.hgt - height above ellipsoid (meters)
%

%***********************************************************************
%* Convert Joint STARS geocentric to geodetic (radians)
%***********************************************************************
%*    Joint STARS UVW coordinate system:
%*      + W axis is coincident with ellipsoid axis of rotation.
%*      + W axis exits ellipsoid at 90 degrees north latitude
%*        (north pole).
%*      + U axis exits ellipsoid at the equator.
%*      + U axis exits ellipsoid at meridian of topocentric site center.
%*      + V axis exits ellipsoid at the equator.
%*      + V axis is 90 degrees east of U axis in longitude.
%*      + The axes form a right-handed (UVW order) cartesian system.
%*      + Values are specified in meters.
%*      + This is similar to a Universal Space Rectangular coordinate
%*        system.
%*
%*    Joint STARS generalized topocentric coordinate system XYZ:
%*      + XYZ axes from a right-handed (XYZ order) cartesian system.
%*      + Z axis is perpendicular to a plane tangent to the ellipsoid
%*        at the topocentric center point. It corresponds to the
%*        geodetic latitude vector at that point.
%*      + Y axis lies in the plane of the meridian of the topocentric
%*        site center.
%*
%*    The plane formed by the XY axes is usually taken to be tanget to
%*    to the ellipsoid at the topocentric site center.  It can be
%*    formed at any distance from the ellipsoid center along the z axis.
%*
%*    The Z axis is positive away from the spheroid.  Values are
%*    specified in meters.  This is local space rectangular coordinate
%*    system and is obtained by translation and rotation from uvw.
%*
%*    Joint STARS latitude, longitude, and elevation are defined in the
%*    usual manner with respect to the ellipsoid (i.e. geodetic).
%*    Thus, elevation is the distance above the ellisoid along the
%*    normal to a plane tangent to the ellipsoid at a point.
%*    Latitude is positive in the northern hemisphere.
%*    Longitude is positive in the eastern hemisphere.
%*
%***********************************************************************

if isvector(uvw)
    uvw = uvw(:);
end
if isvector(origin)
    origin = origin(:);
end

%  Set up WGS-84 constants.
[alpha,f] = EarthModel;

% eccentricity squared for WGS84.
ecc_sq = (2.0 - f) * f;

llhorigin.lat = origin(1,:);
llhorigin.lon = origin(2,:);
llhorigin.hgt = origin(3,:);
lorigin=length(llhorigin.lat);

% Radius of curvature of ellipsoid in a plane perpendicular to
% a meridian and perpendicular to a plane tangent to the surface
% The value N here is for the origin and is used as the initial
% value in the iterations in the geodetic to uvw transformation.

denom = 1.0 - ecc_sq*sin(llhorigin.lat).^2;
N = alpha./sqrt( denom );

% Compute the offset of the geodetic and geocentric centers - a magic
% number first guess.
esqNsin = ecc_sq* N.*sin(llhorigin.lat);

% Compute derivative of N with latitude as help for first guess
dNdlat = esqNsin.*cos(llhorigin.lat)./denom;

tmp1 = sqrt( uvw(1,:).^2 + uvw(2,:).^2 );

for n=1:length(uvw(1,:))
    if ( tmp1(n) == 0.0 )  % At North or South Pole.
        if lorigin==1
            llh.lon(n) = llhorigin.lon;
        else
            llh.lon(n) = llhorigin.lon(n);
        end
        if ( uvw(3,n) > 0.0 )
            llh.hgt(n) = uvw(3,n) - (alpha/sqrt(1.0 - ecc_sq));
            llh.lat(n) = asin(1.0);
        else
            llh.hgt(n) = -uvw(3,n) - (alpha/sqrt(1.0 - ecc_sq));
            llh.lat(n) = asin(-1.0);
        end

    else  % Position is NOT at the Pole.
        if lorigin==1
            llh.lon(n) = llhorigin.lon + atan2( uvw(2,n), uvw(1,n) );

            % Take initial guess at latitude and effective radius, then iterate.
            lat = atan2( uvw(3,n)+esqNsin,tmp1(n));
            % radius of earth in meters
            re = N + dNdlat*(lat-llhorigin.lat);

            dlat = 1.0;
            % Go until roughly half meter error on surface i.e. 1e-7 * 6.38e6
            while ( dlat > 1.0e-7 )
                olatsav = lat;
                tmp2 = uvw(3,n) + ecc_sq*re*sin(lat);
                lat = atan2( tmp2, tmp1(n) );
                re = alpha / sqrt(1.0 - ecc_sq * sin(lat)^2);
                dlat = abs(lat - olatsav);
            end
        else
            llh.lon(n) = llhorigin.lon(n) + atan2( uvw(2,n), uvw(1,n) );
            
            lat = atan2( uvw(3,n)+esqNsin(n),tmp1(n));
            re = N(n)+dNdlat(n)*(lat-llhorigin.lat(n));
            
            dlat = 1.0;
            % Go until roughly half meter error on surface i.e. 1e-7 * 6.38e6
            while ( dlat > 1.0e-7 )
                olatsav = lat;
                tmp2 = uvw(3,n) + ecc_sq*re*sin(lat);
                lat = atan2( tmp2, tmp1(n) );
                re = alpha / sqrt(1.0 - ecc_sq * sin(lat)^2);
                dlat = abs(lat - olatsav);
            end

        end

        llh.hgt(n) = tmp1(n) / cos(lat)  - re;		% height in meters
        llh.lat(n) = lat;

    end
end

return