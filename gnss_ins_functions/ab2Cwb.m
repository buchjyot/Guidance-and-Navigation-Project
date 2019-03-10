        	    function Cwb = ab2Cwb(aoa,ssa)
% %===========================================================%
% %             function Cwb = ab2Cwb(aoa,ssa)                %
% %                                                           %
% %   This functions determines the transformation matrix     %
% %   Cwb that maps vectors from the wind frame to the        %
% %   body frame.  Thus, the input aoa is angle of attack     %
% %   in radians and ssa is the slide slip angle also in      %
% %   radians.  Note that aoa does not follow the standard    %
% %   sign convention (right handed coordinate system) for    %
% %   angles.
% %                                                           %
% %   Programmer:     Demoz Gebre-Egziabher                   %
% %   Created:        July 2, 1998                            %
% %   Last Modified:  March 26, 2009                          %
% %                                                           %
% %===========================================================%

PSI_nb = [0 -aoa ssb];

Cwb = eul2Cbn(PSI_nb)';
%===========================================================%
