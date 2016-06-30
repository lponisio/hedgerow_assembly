function [hnam,hleaf,hint]=hrgplot(f1, varargin)
% HRGPLOT reads files written by fitHRG (C++ program) and renders the
%    corresonding hierarchical random graph as a radial dendrogram.
%    Source: http://www.santafe.edu/~aaronc/randomgraphs/
% 
%    HRGPLOT returns handles to certain parts of the figure, for further
%    customization. hnam is a vector of handles to the text labels on the
%    outer edge (h1), nleaf is a vector of handles to the leaf shapes (h2), 
%    and hint is a vector of handles to the internal vertices (h3).
%    
%    The optional argument 'scale' allows the user to specify as scalar 
%    value by which to rescale these aspects on the initial rendering. The
%    optional argument 'groups' causes hrgplot to import a .groups file
%    that gives the numerical group assigments of each input vertex.
%    
%    Note: In the example below, hrgplot also imports a file named
%    karate-names.lut and, if 'groups' is called, a file named
%    karate.groups. Both names are derived from the input argument f1.
% 
%    Example:
%       filename   = '/Users/myname/Documents/karate-dendro.hrg';
%       h = hrgplot(filename);
%       [h1,h2,h3] = hrgplot(filename,'scale',2);
%       [h1,h2,h3] = hrgplot(filename,'groups');
%
%    For more information, try 'type hrgplot'

% Version 1.0   (2008 April)
% Copyright (C) 2008 Aaron Clauset (Santa Fe Institute)
% Distributed under GPL 2.0
% http://www.gnu.org/copyleft/gpl.html
% HRGPLOT comes with ABSOLUTELY NO WARRANTY
% 
% Notes:
% 
% 1. Verbose mode, which outputs the text of the imported .hrg file, can be
%    activated like so
%    
%       h = hrgplot(filename,'v');

verbose  = false;
lim      = 1.0;
f_groups = false;

% parse command-line parameters; trap for bad input
i=1; 
while i<=length(varargin), 
  argok = 1; 
  if ischar(varargin{i}), 
    switch varargin{i},
        case 'scale',        lim      = varargin{i+1}; i = i + 1;
        case 'groups',       f_groups = true;          i = i + 1;
        case 'v',            verbose  = true;          i = i + 1;
        otherwise, argok=0; 
    end
  end
  if ~argok, 
    disp(['(HRGPLOT) Ignoring invalid argument #' num2str(i+1)]); 
  end
  i = i+1; 
end
if ~isempty(lim) && (~isscalar(lim) || lim<0)
	fprintf('(HRGPLOT) Error: ''scale'' argument must be a scalar value near 1');
    fprintf('; using default.\n');
    lim = 1.0;
end;

% parse filename
k = findstr(f1,'_best-dendro.hrg');
if isempty(k),
	fprintf('Error: %s does not appear to be a valid HRG file.\n',f1);
    return;
end;
f2 = strcat(f1(1:k-1),'-names.lut');
if f_groups, f3 = strcat(f1(1:k-1),'.groups'); else f3 = ''; end;

% --- Import HRG data---
% import the -dendro.hrg file
try
	[I,L,Ls,R,Rs,p,e,n] = textread(f1,'[ %u ] L= %u%s R= %u%s p= %n e=%u n=%u');
catch
	fprintf('Error: failed to read %s\n',f1);
	return;
end;

% if necessary, report the contents of the file
if verbose
	for i=1:length(I)
		fprintf('[ %i ] L= %i %s R= %i %s p= %9.7f e= %i n= %i\n', ... 
			I(i),L(i),Ls{i},R(i),Rs{i},p(i),e(i),n(i));
	end;
end;
I = I+1; L = L+1; R = R+1;	% adjust the indexing so that it starts at 1 not 0

% import the -names.lut file. this is a table that maps the node names in the
% .hrg file onto the node names in the .pairs file
try
	names = dlmread(f2,'\t',1,0);
	names(:,1) = names(:,1) + 1;
catch
	fprintf('Warning: failed to read %s. Using defaults.\n',f2);
	names = [(1:length(I)+1)' (1:length(I)+1)'];
end;

% import the .groups file. this maps node names in the .pairs file onto their
% group indices
if f_groups
    try
        ss = textread(f3,'%s');
        gg = sortrows([str2num(char(ss(1:2:end))) str2num(char(ss(2:2:end)))]);
        groups = zeros(size(names));
        for i=1:size(groups,1)
            groups(i,:) = [names(i,1) gg(gg(:,1)==names(i,2),2)];
        end;
    catch
        fprintf('Warning: failed to read %s. Using defaults.\n',f3);
        groups = [names(:,1) ones(size(names(:,1)))];
    end;
else
    groups = [names(:,1) ones(size(names(:,1)))];
end;

% --- Parse HRG data ---
% do DFS traversal of tree structure (starting at the root, which is always 
% the first entry L(1), and create list of leaf vertices in the order they 
% were visited. simultaneously, build a matrix structure that one could pass to
% the dendrogram() function; we also build a variation on that structure where
% instead of labeling internal vertices with new numbers, we give them the
% minimum leaf index in their subtree (this is Q)

leaf  = [];
D     = []; Q = [];
curr  = I(1);
marks = zeros(size(I)); marks(curr) = 1;
pars  = -1*ones(size(I));
depth = zeros(size(I));
dep   = 0;
H     = [I zeros(size(I))];
con   = max([L; R]);
N     = [I I];
while (true)
	if (marks(curr) == 1)
		marks(curr) = 2;
		if (strcmp(Ls{curr}, '(D)'))
			pars(L(curr)) = curr;
			curr = L(curr);
			dep = dep + 1;
			marks(curr) = 1;
		else % Ls{curr} == '(G)'
			leaf = [leaf; L(curr)];
			depth(L(curr)) = dep + 1;
		end;
	end;
	if (marks(curr) == 2)
		marks(curr) = 3;
		if (strcmp(Rs{curr},'(D)'))
			pars(R(curr)) = curr;
			curr = R(curr);
			dep = dep + 1;
			marks(curr) = 1;
		else % Rs{curr} == '(G)'
			leaf = [leaf; R(curr)];
			depth(R(curr)) = dep + 1;
		end;
	end;
	if (marks(curr) == 3)
		marks(curr) = 4;
		con         = con + 1;
		N(curr,2)   = con;
		if (strcmp(Ls{curr},'(G)'))
			left  = L(curr);
			hleft = 0;
			qleft = left;
		else
			left  = N(L(curr),2);
			hleft = H(L(curr),2);
			qleft = N(L(curr),3);
		end;
		if (strcmp(Rs{curr},'(G)'))
			right  = R(curr);
			hright = 0;
			qright = right;
		else
			right  = N(R(curr),2);
			hright = H(R(curr),2);
			qright = N(R(curr),3);
		end;
		N(curr,3) = min(qleft,qright);
		H(curr,2) = 1+max(hleft,hright);
		D = [D; [left  right  H(curr,2)]];
		Q = [Q; [qleft qright H(curr,2)]];
		curr = pars(curr);
		dep  = dep - 1;
		if (curr == -1)
			break;		% root has no parent
		end;
	end;
end;
% build map of row number of D matrix to dendrogram node index from the input
% file, and the branch probability that lives there.
pmap = [N(:,1) N(:,2)-max([L;R])];
[u,b] = sort(pmap(:,2));
pmap = [u pmap(b,1) p(b)];
% cleanup memory
clear marks H con dep pars N;

% --- Draw radial HRG ---
m    = size(Q,1)+1;
cmap = [0 0 0];
T    = zeros(2*(m-1),1);
R    = zeros(2*(m-1),2);
Y    = ones(m,1);
W    = Q;

% get list of leafs in the order they first appear in Q. this will be used
% to label the leafs in the figure. then create X vector, which maps a
% leaf index in Q to its clockwise-ordered position in the leaf ordering
S     = D(:,1:2)';
perm  = S(S<(m+1));
label = num2str(names(perm,2));
glabs = groups(perm,2);
[u,X] = sort(perm);

dtheta = 2*pi/m;
angles = dtheta*[0:(m-1)];
		
% create end-points for posts (we'll do the crossbars later), in polar
% coordinates. the posts are pairs of lines, where T(i) stores the angle of
% the ith post, and R(i,:) stores the radii of the outer and inner
% endpoints of the post.
h      = zeros(m-1,1);						% list of handles to posts
W(:,3) = (1 - (Q(:,3)./max(Q(:,3))));		% rescale distances
A      = angles(X)';						% angles of leaf nodes
for n=1:(m-1)
    [i,j,w]    = deal(W(n,1),W(n,2),W(n,3));% 
    T(2*n-1)   = A(i);						% angle of left-post
    R(2*n-1,:) = [Y(i) w];					% radii of left-post endpoints
    T(2*n)     = A(j);						% angle of right-post
    R(2*n,:)   = [Y(j) w];					% radii of right-post endpoints
    A(i)       = (A(i)+A(j))/2;				% angle to center of crossbar
    Y(i)       = w;                         % new depth of subtree
    C(n,:)     = [A(i) w];					% angle and radius of branch point 
end
	
% create arcs for crossbars. these are just portions of a circle that span
% the inner endpoints of the two posts, converted to cartesian coordinates
% for plotting.
for n=1:(m-2)
    t      = linspace(T(2*n-1),T(2*n));
    cb{n}  = [R(2*n,2).*sin(t)' R(2*n,2).*cos(t)'];
    ext(n) = abs(T(2*n-1)-T(2*n));
end;

% draw posts
figure;
set(gcf,'Position', [0, 0, 500, 500]);
for n=1:length(T)
    x = [R(n,1)*sin(T(n)) R(n,2)*sin(T(n))];	% convert to cartesian coords.
    y = [R(n,1)*cos(T(n)) R(n,2)*cos(T(n))];	%
    h(n) = line(x,y,'color',cmap); hold on;		% plot line
    set(h(n),'LineWidth',lim*2);
end
% draw crossbars
for n=1:(m-2)
    plot(cb{n}(:,1),cb{n}(:,2),'k','LineWidth',2*lim);
end;
% draw branch probabilities
msize = floor(lim*9);
hint  = zeros(m-1,1);
for i=1:(m-1)
    hint(i) = plot(C(i,2)*sin(C(i,1)),C(i,2)*cos(C(i,1)),'o','MarkerSize', ...
        msize, 'MarkerFaceColor', pmap(i,3)*ones(1,3), 'MarkerEdgeColor', ...
        [0 0 0]);
end;
	
% pretty-up the figure: lay out the leaf nodes
x = sin(angles)';
y = cos(angles)';
mar = {'s','v','^','o','s','v','^','o','s','v','^','diamond'};
if f_groups,
    gn  = length(unique(groups(:,2)));
    mc  = linspace(0,1,ceil(log(gn)/log(3))); mcn = length(mc);
    mfc = zeros(gn,3);
    for i=1:gn
        mfc(i,:) = [mc(mod(floor((i-1)/mcn^2),mcn)+1) mc(mod(floor((i-1)/mcn),mcn)+1) mc(mod(i-1,mcn)+1)];
    end;
else
    mfc = [1 1 1];
end;

hleaf = zeros(m,1);
for i=1:m
    hleaf(i)=plot(x(i),y(i),mar{glabs(i)},'MarkerSize',msize,'MarkerFaceColor', ...
        mfc(glabs(i),:), 'MarkerEdgeColor',[0 0 0]); hold on;
    set(hleaf(i),'Color',cmap);
end
axis square;
axis([-1 1 -1 1]);
axis off;

% pretty-up the figure: lay out their labels
fac  = 1.05;
hnam = zeros(m,1);
for i=1:m
    if (i <= ceil(m/2))
        hnam(i) = text(fac*x(i),fac*y(i),label(i,:));
        % right-side of circle
        set(hnam(i),'Rotation', 90-180*angles(i)/pi, 'FontWeight','bold', ... 
            'HorizontalAlignment','left', 'VerticalAlignment', ... 
            'Middle','Color',cmap,'FontSize',10*lim);
    else
        hnam(i) = text(fac*x(i),fac*y(i),label(i,:));
        % left-side of circle
        set(hnam(i),'Rotation',-90-180*angles(i)/pi, 'FontWeight','bold', ... 
            'HorizontalAlignment','right', 'VerticalAlignment', ... 
            'Middle','Color',cmap,'FontSize',10*lim);
    end;
end;
hold off;
