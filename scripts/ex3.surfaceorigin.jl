#= %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Example 3: Find the surface origin of water for some interior box %
 %                                                                   %
 % This is equivalent to solving a sensitivity problem:              %
 % The total volume is V = v^T c , where v is the volume of a
 % given interior box,
 % and c is the fraction of volume from a given source which         %
 % satisfies the equation A c = d.                                   %
 % Next, dV/d(d) = A^(-T) v, and dV/d(d) is exactly the volume       %
 % originating from each source.      
 % Very similar mathematically to example 2.                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
=#

 nz = maximum(kt);
 ny = maximum(jt);
 nx = maximum(it);
 nfield = size(A,1);
 
# % choose an interior location X (Xlon[lon], Xlat [lat], Xdepth [m depth]).
#% -7.38, 115.26E
 xlon = 125.26; #% deg E.
 xlat = -6.38; 
 xdepth = 500; #% meters.
 
 % Find the coordinate on the grid by linear interpolation. 
 LONtmp = [LON(1)-4 ; LON; LON(end)+4];  %% watch out for the
                                         %wraparound of the grid.
 iX = interp1(LONtmp,0:length(LON)+1,Xlon,'linear')
 jX = interp1(LAT,1:length(LAT),Xlat,'linear')
 kX = interp1(DEPTH,1:length(DEPTH),Xdepth,'linear')
 kX(Xdepth>5500) = 33; % if deeper than the bottom level, put it on
                       % the bottom level.

 % watch out for wraparound again.
 i2 = i;
 i2( i-iX > NX/2) = i2( i-iX > NX/2) - NX;
 i2( iX-i > NX/2) = i2( iX-i > NX/2) + NX;
  
 % find the closest gridpoints to the location of interest.
 [dist,loc] = sort((i2-iX).^2 + (j-jX).^2 + (k-kX).^2);
 dist = dist+0.1; % to eliminate singularity if point matches up
                  % exactly with the grid. 
 
 mask = zeros(Nfield,1); 
 mask(loc(1:6)) =(1./dist(1:6))./sum(1./dist(1:6)); % mask adds up
                                                    % to 1.
 
 vtot = R' \ (P' * (L' \ (U' \ (Q' * mask)))); 
 Vtot = vector_to_field(vtot,i,j,k);
 
 Vtot = squeeze(Vtot(1,:,:)); % just keep the surface as that is
                              % the ultimate source region.

 contourf(LON,LAT,Vtot,[0 .001 .005 .01 .05 .1 .5]) % one way to
                                                    % plot it.
