%just checking ray corner positions and things

rayN = 55;

i = 1;
distance1 = sqrt( (stf(1).ray(rayN).targetPoint(1) - stf(1).ray(rayN).rayCorners_SCD(i,1))^2 + (stf(1).ray(rayN).targetPoint(2) - stf(1).ray(rayN).rayCorners_SCD(i,2))^2 + (stf(1).ray(rayN).targetPoint(3) - stf(1).ray(rayN).rayCorners_SCD(i,3))^2 )
i = 2;
distance2 = sqrt( (stf(1).ray(rayN).targetPoint(1) - stf(1).ray(rayN).rayCorners_SCD(i,1))^2 + (stf(1).ray(rayN).targetPoint(2) - stf(1).ray(rayN).rayCorners_SCD(i,2))^2 + (stf(1).ray(rayN).targetPoint(3) - stf(1).ray(rayN).rayCorners_SCD(i,3))^2 )
i = 3;
distance3 = sqrt( (stf(1).ray(rayN).targetPoint(1) - stf(1).ray(rayN).rayCorners_SCD(i,1))^2 + (stf(1).ray(rayN).targetPoint(2) - stf(1).ray(rayN).rayCorners_SCD(i,2))^2 + (stf(1).ray(rayN).targetPoint(3) - stf(1).ray(rayN).rayCorners_SCD(i,3))^2 )
i = 4;
distance4 = sqrt( (stf(1).ray(rayN).targetPoint(1) - stf(1).ray(rayN).rayCorners_SCD(i,1))^2 + (stf(1).ray(rayN).targetPoint(2) - stf(1).ray(rayN).rayCorners_SCD(i,2))^2 + (stf(1).ray(rayN).targetPoint(3) - stf(1).ray(rayN).rayCorners_SCD(i,3))^2 )


%get distance from point Q to plane:
v1 = [ (stf(1).ray(rayN).rayCorners_SCD(1,1) - stf(1).ray(rayN).rayPos(1)) (stf(1).ray(rayN).rayCorners_SCD(1,2) - stf(1).ray(rayN).rayPos(2)) (stf(1).ray(rayN).rayCorners_SCD(1,3) - stf(1).ray(rayN).rayPos(3)) ];
v2 = [ (stf(1).ray(rayN).rayCorners_SCD(2,1) - stf(1).ray(rayN).rayPos(1)) (stf(1).ray(rayN).rayCorners_SCD(2,2) - stf(1).ray(rayN).rayPos(2)) (stf(1).ray(rayN).rayCorners_SCD(2,3) - stf(1).ray(rayN).rayPos(3)) ];

v1x2 = cross(v1, v2);
N    = v1x2 / norm(v1x2);

% Vector from P (point defines plane) to Q:
P = stf(1).ray(rayN).rayCorners_SCD(4,:);
Q = stf(1).ray(rayN).rayPos;
PQ = Q - P;

% Dot product between line from Q to P1 and normal of the plane:
Dist = dot(PQ, N)


%get distance from point Q2 (target point) to plane:
v3 = [ (stf(1).ray(rayN).rayCorners_SCD(1,1) - stf(1).ray(rayN).targetPoint(1)) (stf(1).ray(rayN).rayCorners_SCD(1,2) - stf(1).ray(rayN).targetPoint(2)) (stf(1).ray(rayN).rayCorners_SCD(1,3) - stf(1).ray(rayN).targetPoint(3)) ];
v4 = [ (stf(1).ray(rayN).rayCorners_SCD(2,1) - stf(1).ray(rayN).targetPoint(1)) (stf(1).ray(rayN).rayCorners_SCD(2,2) - stf(1).ray(rayN).targetPoint(2)) (stf(1).ray(rayN).rayCorners_SCD(2,3) - stf(1).ray(rayN).targetPoint(3)) ];

v3x4 = cross(v3, v4);
N2    = v3x4 / norm(v3x4);

% Vector from P (point defines plane) to Q:
P2 = stf(1).ray(rayN).rayCorners_SCD(4,:);
Q2 = stf(1).ray(rayN).targetPoint;
PQ2 = Q2 - P2;

% Dot product between line from Q to P1 and normal of the plane:
Dist2 = dot(PQ2, N2)

% for rayN = 1:length(stf(1).ray)
%     
%     v1 = [ (stf(1).ray(rayN).rayCorners_SCD(1,1) - stf(1).ray(rayN).rayPos(1)) (stf(1).ray(rayN).rayCorners_SCD(1,2) - stf(1).ray(rayN).rayPos(2)) (stf(1).ray(rayN).rayCorners_SCD(1,3) - stf(1).ray(rayN).rayPos(3)) ];
%     v2 = [ (stf(1).ray(rayN).rayCorners_SCD(2,1) - stf(1).ray(rayN).rayPos(1)) (stf(1).ray(rayN).rayCorners_SCD(2,2) - stf(1).ray(rayN).rayPos(2)) (stf(1).ray(rayN).rayCorners_SCD(2,3) - stf(1).ray(rayN).rayPos(3)) ];
% 
%     v1x2 = cross(v1, v2);
%     N    = v1x2 / norm(v1x2);
% 
%     % Vector from P (point defines plane) to Q:
%     P = stf(1).ray(rayN).rayCorners_SCD(4,:);
%     Q = stf(1).ray(rayN).rayPos;
%     PQ = Q - P;
% 
%     % Dot product between line from Q to P1 and normal of the plane:
%     Dist = dot(PQ, N);
%     fprintf("%f    ",Dist);
% 
%     %get distance from point Q2 (target point) to plane:
%     v3 = [ (stf(1).ray(rayN).rayCorners_SCD(1,1) - stf(1).ray(rayN).targetPoint(1)) (stf(1).ray(rayN).rayCorners_SCD(1,2) - stf(1).ray(rayN).targetPoint(2)) (stf(1).ray(rayN).rayCorners_SCD(1,3) - stf(1).ray(rayN).targetPoint(3)) ];
%     v4 = [ (stf(1).ray(rayN).rayCorners_SCD(2,1) - stf(1).ray(rayN).targetPoint(1)) (stf(1).ray(rayN).rayCorners_SCD(2,2) - stf(1).ray(rayN).targetPoint(2)) (stf(1).ray(rayN).rayCorners_SCD(2,3) - stf(1).ray(rayN).targetPoint(3)) ];
% 
%     v3x4 = cross(v3, v4);
%     N2    = v3x4 / norm(v3x4);
% 
%     % Vector from P (point defines plane) to Q:
%     P2 = stf(1).ray(rayN).rayCorners_SCD(4,:);
%     Q2 = stf(1).ray(rayN).targetPoint;
%     PQ2 = Q2 - P2;
% 
%     % Dot product between line from Q to P1 and normal of the plane:
%     Dist2 = dot(PQ2, N2);
%     fprintf("%f\n",Dist2);
% end


%let's look at beam angles:
refX = [1 0 0];
refY = [0 1 0];

for rayN = 1:length(stf(1).ray)
    ray = stf(1).ray(rayN).targetPoint - stf(1).ray(rayN).rayPos;
    thetaX = acosd( dot(refX,ray)/norm(ray));
    fprintf("%f, ",thetaX);
    thetaY = acosd( dot(refY,ray)/norm(ray));
    fprintf("%f\n",thetaY);
end


% rayposbev, targetpointbev,raypos,targetpoint = {[-40,0,-40],[-80,1000,-80],[-40,0,-40],[-80,1000,-80]}
% raycorners = [-18.75,-500,-18.75;-21.25,-500,-18.75;-21.25,-500,-21.25;-18.75,-500,-21.25]