function DC = Get_DirCos_ForwardT (A, MatrixFlavor);
%
% Description:
%	Fills a direction cosine matrix defined by positive right-hand rule Euler
%	angles that transforms from an INS type basis to a body type basis.    
%
% Usage:
%  DC = Get_DirCos_ForwardT (A, MatrixFlavor);
%
% Inputs:
%  A            - Angle in Radians
%  MatrixFlavor - Axis: ROLL_TYPE, PITCH_TYPE, YAW_TYPE
%
% Outputs:
%  DC           - Direction Cosine Matrix
%
% Modified by CLR to accept N vector of angles & return (3x3xN) rotation matrices
%

YAW_TYPE = 1;
PITCH_TYPE = 2;
ROLL_TYPE = 3;

cosA =  cos(A(:));
sinA =  sin(A(:));
DC   =  zeros(3,3,length(A));
   
switch (MatrixFlavor) 

	case YAW_TYPE,
	        DC(1,1,:) = cosA;
	        DC(1,2,:) = sinA;
	        %DC(1,3) =  0;
	        DC(2,1,:) = -sinA;
	        DC(2,2,:) = cosA;
	        %DC(2,3) =  0;
	        %DC(3,1) =  0;
	        %DC(3,2) =  0;
	        DC(3,3,:) =  1;

	case PITCH_TYPE,
           DC(1,1,:) = cosA;
           %DC(1,2) =  0;
           DC(1,3,:) = -sinA;
           %DC(2,1) =  0;
           DC(2,2,:) =  1;
           %DC(2,3) =  0;
           DC(3,1,:) = sinA;
           %DC(3,2) =  0;
           DC(3,3,:) = cosA;

	case ROLL_TYPE,
	        DC(1,1,:) =  1;
	        %DC(1,2) =  0;
	        %DC(1,3) =  0;
	        %DC(2,1) =  0;
	        DC(2,2,:) = cosA;
	        DC(2,3,:) = sinA;
	        %DC(3,1) =  0;
	        DC(3,2,:) = -sinA;
	        DC(3,3,:) = cosA;
end

return