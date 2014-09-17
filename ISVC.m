function [tspan, isvc, h, ai] = ISVC(force, w,eta,tfinal, th, hinf, A, plot_ai, plot_movie,noplot, a)
%ISVC this function computes the inter strand void content evolution during
%compression forming of Randomly Oriented Strands of Thermoplastic
%Composite.
%     [tspan, isvc] = ISVC(force, w,eta,tfinal, th, hinf, A, a, plot_ai,
%     plot_movie,noplot) returns two vectors tspan for the time stamps and
%     isvc for the inter strand void content. force is the closing force, w
%     the width of the strands, eta the viscosity of the melt, tfinal the
%     final time computed for, th the thickness of one strand, hinf the
%     expected final nominal thickness of the part, A the area of the mold,
%     a the initial column distribution (vector). plot_ai a boolean to
%     express whether the column distribution evolution should be plot and
%     plot_movie a boolean to express wether the movie should be plot.If
%     noplot=1, nothing is plot. tspan is the time stamp vector at which
%     the isvc is returned.
%
%     ISVC is the same as ISVC(5000, 3.18e-3, 4.01e4, 5*60, 0.146e-3, 2.67e-3,
%     (4*25.4e-3)^2, ai, 0, 0) where ai is the distribution representative
%     of a perfect randomly positionning of the strands (Poisson based).
%     This is the reference case of the article.
%
%     [tspan, isvc, h] = ISVC(...) returns, in addition the vector h
%     giving the gap between plattens at each time
%
%     [tspan, isvc, h, ai] = ISVC(...) returns, in addition the column
%     representation ai at each time. ai is a 2D array.
%      
%     ISVC(force, w,eta,tfinal, th,hinf, A) is the same as ISVC(force,
%     w,eta,tfinal, th, hinf, A, ai) where ai is the distribution
%     representative of a perfect randomly positionning of the strands
%     (Poisson based).
%
%see also: GUI_ISVC



%% initialization
if(nargin>9 && noplot==1)
	plot_ai = 0;
	plot_movie = 0;
else
	noplot=0;
end
	
if (nargin<8)
	%flags not provided:
	plot_ai = 0;
	plot_movie =0;
	
	if (nargin<7) %nothing provided
		
		force = 5000; %force applied [N]
		w=3.18e-3; %width of a strand [m]
		eta = 4.01e4;%newtonian viscosity [Pa.s]
		tfinal = 5*60; %final time [s]
		th =1.46e-4; %thickness of a strand [m]
		hinf = 2.67e-3; %nominal final thickness expected
		A = (4*25.4e-3)^2; %sample area [m2]
	end
	

	
end

%% %%%%% distribution %%%%%%%%
if (nargin~=11) %a not provided
	a = poisson_distrib(hinf, th, A); %initial a(i) distrib.
end

a_init = a; %store initial value

%max number of strands stacked
n = length(a)-1;
% *note that a(1) is the column area with no strand called $a_0$ in the doc.

%initial thickness
h = th*n;

%final thickness
hfinal = get_final_h( a_init,th); %final minimum thickness (for isvc=0)

%% Unknown: concatenation of h and a.
%my unknown initial value
X0 = [h; a];

%% solving
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-17);

[tspan, X] = ode45(...
	@(t,X_) residual(X_,force/eta,w, th,a_init, hfinal),...
	[0,tfinal],...
	X0,options);

%% output:
h = X(:,1); %thickness
isvc = (X(:,1)-hfinal)./X(:,1); % interstrand void content
ai = X(:,2:end); % column representation

%% plotting
if(~noplot)
	%thickness
	figure(1);
	hold all
	plot(tspan, h);
	xlabel('Time (s)');
	ylabel('Thickness (m)');
	plot([0,tspan(end)],[hfinal,hfinal], ':r');
	
	%void content
	figure (2);
	hold all
	plot(tspan, isvc)
	xlabel('Time (s)');
	ylabel('Inter-Strand Void Content');
	
	%evolution of ai
	if(plot_ai)
		figure (3);
		plot(tspan, X(:,2:end));
		legend(num2str((0:n)'));
		xlabel('Time (s)');
		ylabel('a(i)');
	end
	
	%movie rectangles:
	if(plot_movie)
		figure(4)
		axis auto;
		set(gca,'NextPlot','replacechildren');
		F(fix(size(X,1)/20)) = struct('cdata',[],'colormap',[]);
		for J=1:20:size(X,1);
			plot(0,0);
			a = X(J,2:end);
			thicknesses = 0:th:th*(n);
			thicknesses(thicknesses>X(J,1))=X(J,1);
			x=0;
			for k=1:n+1
				rectangle('Position',[x,0,a(k)+1e-10,thicknesses(k)+eps],...
					'LineStyle','-','EdgeColor','k')
				x=x+a(k);
			end
			
			if J == 1
				drawnow;
				axis manual;
			end
			
			legend(['time: ' num2str(tspan(J)), ' s']);
			F((J-1)/20 + 1) = getframe;
		end
		
		movie(F,1);
		movie2avi(F,'myavifile.avi');
	end
	
end
end


% given a current configuration (X) returns the ode RHS Xdot
function Xdot = residual(X,F_eta,w, th, a_init, hfinal )

h = X(1);
if (h>hfinal) %checks that we did not fill up the voids yet.
	a = X(2:end);
	n = length(a)-1;%max number of strands stacked
	A = sum(a_init);
	
	%column that will be filled (ie lowest non zero column)
	p = find(a>0,1,'first');
	
	%smallest stack number that will be squeezed 
	j = ceil (h/th) +1;
	%ie all columns between j and n are squeezed
	
	%thickness of interest, ie thickness of the highest column minus
	%thickness of the lowest. This is the height of the "stairs". also
	%called "deformable height"
	h_interest = h - (p-1)*th;
	
	%initialize
	adot=zeros(length(X)-1,1);
	
	%% squeeze of a n-j columns:
	%number of strands per column (in the height)
	k = (1:n+1)'-p;
	
	% the s parameter:
	s_i = k.^5.*a_init/A;
	s = sum(s_i(j:n+1));
	
	%squeezing using the s parameter
	if(s == 0) %no contact yet
		hdot = -1e6;
	else
		hdot = - F_eta * h_interest^6 / (w^2 * th^3 * A * s);
	end
	
	%% evolutions according to hdot:
	%increasing of columns spreaded (incompressibility):
	adot(j:n+1) = -hdot*a(j:n+1) / h_interest;

	%reducing of this column p due to spreading:
	adot(p) = -sum(adot(j:n+1));%matter ./th_lowest_col;
	
	%% Concatenate
	Xdot = [hdot; adot];
else % we filled the void : no more evolution:
	Xdot =zeros(length(X),1);
end

end


% this returns the final thickness that would correspond to zero
% interstrand void content
function h_final = get_final_h( a_init,th)

	%total area
	A = sum(a_init);
	
	n = length(a_init)-1;
	
	% total volume of chips (sum for each column):
	V = sum(a_init.*(0:n)'*th);
	
	%final h for which void content = 0
	h_final = V / A;
end

%this returns the ideal poisson distribution of a_i
function a = poisson_distrib(hinf, th, A)
	mean_str = hinf/th;%mean number of strands in thickness
	max_str = 4*mean_str; %max number of strands modeled
	poisson = exp(-mean_str)*mean_str.^(0:max_str)...
		./factorial(fix(0:max_str));%this is my poisson distrib.
	%discard negligible colums
	max_str = find(poisson>1e-5,1,'last'); %the last non negligible column
	poisson((max_str+1):end) = [];%discard largers
	a = poisson'*A;
end
