
function MyProgram(DataFileName)

clc(); close all;

% In case the caller does not specify the input argument, we propose a default one.
if ~exist('DataFileName','var'), DataFileName ='Laser__2.mat'; end;

% load data     file.
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
MyGUIHandles.bluePlot = plot(X,Y,'b.');
               
ii = find(intensities~=0);          % find those "pixels" that had intense reflection (>0) (aka: Highly Reflective pixels, HR)
MyGUIHandles.RedPlot = plot(X(ii),Y(ii),'+r');             % plot highly reflective ones
axis([-10,10,0,20]);                % focuses plot on this region ( of interest in L220)
MyGUIHandles.GreenPlot = plot(0,0,"g*");
% To be done (by you)
OOIs = ExtractOOIs(ranges,intensities) ;
PlotOOIs(OOIs,MyGUIHandles);


return;
end
function r = ExtractOOIs(ranges,intensities)
    
    global temp 
    temp = [];%for index
    global temp1 
    temp1 =[];%for corresponding range value
    global num;
    num = 1;
    r.N = 1;
    r.Centers= [];
    r.Sizes=[];
    global scanNum;
    scanNum = length(intensities);
    for n = 1:scanNum
        if intensities(n)>0 
            temp(num)= n;
            temp1(num)=ranges(n);
            
            num = num+1;
            findNearByponits(n,ranges,num,intensities);
            
            
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
return;
end
function  info = objDetail(points,pointsRange)

angles = (points-1).*0.5.* pi./180 ;
X = cos(angles).*pointsRange;
Y = sin(angles).*pointsRange;


coor = CenterFind(X,Y);
cx = coor(1);
cy = coor(2);

dia = sqrt(max((X-cx).^2+(Y-cy).^2));
info = [cx,cy,dia];

end
function findNearByponits(num,r,num1,i)
    global temp;
    global temp1;
    global scanNum;
    constraint = 0.3;
    currentAngle = num;
    currentLength = r(num);
    for s =num-1:-1:1
       
        if(EulDistance(s,r(s),currentAngle,currentLength)<=constraint && i(s)==0 )
            temp(num1) = s;
            temp1(num1) = r(s);
        elseif(i(s)>0)
            continue;
        else
            break;
        end
    end
    
    for s =num+1:1:scanNum
        if(EulDistance(s,r(s),currentAngle,currentLength)<=constraint&& i(s)==0 )
            temp(num1) = s;
            temp1(num1) = r(s);
        elseif(i(s)>0)
            continue;
        else
            break;
        end
    end
end
function d =  EulDistance(a1,l1,a2,l2)
angle1 = a1*0.5* pi/180;
angle2 = a1*0.5* pi/180;

X1 = cos(angle1).*l1;
Y1 = sin(angle1).*l1;

X2 = cos(angle2).*l2;
Y2 = sin(angle2).*l2;

d = sqrt((X1-X2)^2+(Y1-Y2)^2);
disp(d);

end

%Randy bullock's paper
function center = CenterFind(X,Y)
%     x = mean(X);
%     y = mean(Y);
      x = (max(X)+min(X))/2;
      y = (max(Y)+min(Y))/2;

     center = [x;y];
    return ;
end


function PlotOOIs(OOIs,handler)
    if OOIs.N<1, return ; end
    % your part....
%     for n = 1:OOIs.N
%         c = OOIs.Centers(n);
% %        MyGUIHandles.GreenPlot = plot(c.x,c.y,"g*")
%        set(handler.GreenPlot,'xdata',c.x,'ydata',c.y);
%     end
    table = struct2table(OOIs.Centers())
    x  = table.x;
    y = table.y;
    set(handler.GreenPlot,'xdata',x,'ydata',y);
  
return;
end





