                function ned = wgsxyz2ned2(ecef,ecef_ref)
%----------------------------------------------------------------------
%               function ned = wgsxyz2ned2(ecef,ecef_ref)
%
%   Converts a vector (other than a position vector) given in ECEF 
%   coordinates to a vector in North East Down coordinates centered
%   at the coordinates given by ecef_ref.
%
%   Programmer: Demoz Gebre-Egziabher 
%   Created: 12/31/98
%---------------------------------------------------------------------

[lat, lon, ~] = wgsxyz2lla(ecef_ref);

enu(3,1)= cosd(lat)*cosd(lon)*ecef(1)+cosd(lat)*sind(lon)*ecef(2)+sind(lat)*ecef(3);
enu(1,1)=-sind(lon)*ecef(1) + cosd(lon)*ecef(2);
enu(2,1)=-sind(lat)*cosd(lon)*ecef(1)-sind(lat)*sind(lon)*ecef(2)+cosd(lat)*ecef(3);

C = eul2Cbn([pi 0 pi/2])';
ned = C*enu;
