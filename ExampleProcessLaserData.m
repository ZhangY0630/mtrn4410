
% Example program, for off-line processing of saved laser scans. 
% Example source code, useful for Parst A and B of project 01.. 
% AAS - 2020.T1
% Jose Guivant.

% Note: Plotting of results is implemented, here, via a brute force approach.
% See other provided examples, for better implementation (dynamic updates using "set()").
% I expect you to avoid "brute force" implementations, in your projects. 






function MyProgram(DataFileName)

clc(); close all;

% In case the caller does not specify the input argument, we propose a default one.
if ~exist('DataFileName','var'), DataFileName ='Laser__2.mat'; end;

% load data file.
load(DataFileName); 
% The variable that contains the data structure is named "dataL"
% (because it was created under that name and saved)  

N = dataL.N;                         % number of scans in this squence.

disp('Brute force plotting');
disp('There is no GUI. Use "Control-C" to break the loop.');

for i=1:1:N,                        % in this example, I skip some of them..
    scan_i = dataL.Scans(:,i);
    ProcessScan(scan_i);
    
    s=sprintf('(Brute force plotting...)\nShowing scan #[%d]/[%d]\r',i,N);
    title(s);
    
    pause(0.01) ;                   % wait for ~10ms
end;
disp('Done. Bye.');

return;
end


%.............................
function ProcessScan(scan)

% Extract range and intensity information, from raw measurements.
% Each "pixel" is represented by its range and intensity of reflection.
% It is a 16 bits number whose bits 0-12 define the distance (i.e. the range)
% in cm (a number 0<=r<2^13), and bits 13-15 indicate the intensity 
%( a number 0<=i<8 ).

% We extract range and intensity, here.
%useful masks, for dealing with bits.
mask1FFF = uint16(2^13-1);
maskE000 = bitshift(uint16(7),13)  ;

intensities = bitand(scan,maskE000);

ranges    = single(bitand(scan,mask1FFF))*0.01; 
% Ranges expressed in meters, and unsing floating point format (e.g. single).

% 2D points, expressed in Cartesian. From the sensor's perpective.
angles = [0:360]'*0.5* pi/180 ;         % associated angle, for each individual range in a scan
X = cos(angles).*ranges;
Y = sin(angles).*ranges;    

% Plot. "BRUTE FORCE" plot (see note 1).
figure(1) ; clf(); hold on;
plot(X,Y,'b.');                     % all points
ii = find(intensities~=0);          % find those "pixels" that had intense reflection (>0) (aka: Highly Reflective pixels, HR)
plot(X(ii),Y(ii),'+r');             % plot highly reflective ones
axis([-10,10,0,20]);                % focuses plot on this region ( of interest in L220)

% To be done (by you)
OOIs = ExtractOOIs(ranges,intensities) ;
PlotOOIs(OOIs);


return;
end
function r = ExtractOOIs(ranges,intensities)
    
    temp = [];%for index
    temp1 =[];%for corresponding range value
    num = 1;
    r.N = 1;
    r.Centers= [];
    r.Sizes=[];
    scanNum = length(intensities);
    for n = 1:scanNum
        if intensities(n)>0 
            temp(num)= n;
            temp1(num)=ranges(n);
            num = num+1;
        end
        if (intensities(n)==0 && ~isempty(temp))
            %r.Centers = objectDetail(temp,temp1);
            if(length(temp)>=2)
%                 r.Centers,r.Sizes = ObjectDetail(temp,temp1);
                info = objDetail(temp,temp1);
                r.Centers(r.N).x = info(1);
                r.Centers(r.N).y = info(2);
                r.Sizes(r.N) = info(3);
                r.N = r.N+1;

            end
% initialise the count
            temp=[];
            temp1= [];
            num = 1;
        end
        
    end
    r.N = r.N-1;
    disp(r.N);
return;
end
function  info = objDetail(points,pointsRange)
coor = [];
angles = (points-1).*0.5.* pi./180 ;
X = cos(angles).*pointsRange;
Y = sin(angles).*pointsRange;


coor = CenterFind(X,Y);
cx = coor(1);
cy = coor(2);

dia = sqrt(max((X-cx).^2+(Y-cy).^2));
info = [cx,cy,dia];

end



%Randy bullock's paper
function center = CenterFind(X,Y)
    x = mean(X);
    y = mean(Y);

     center = [x;y];
    return ;
end


function PlotOOIs(OOIs)
    if OOIs.N<1, return ; end
    % your part....
    for n = 1:OOIs.N
        c = OOIs.Centers(n);
       plot(c.x,c.y,"g*")
    end
return;
end






    
% --------------------------------
% note 1: for a more efficient way of dynamically plotting, see example
% "ExampleUseLaserData.m", where we used graphical objects handles.
% --------------------------------
% Questions?  Ask the lecturer, j.guivant@unsw.edu.au
% --------------------------------

