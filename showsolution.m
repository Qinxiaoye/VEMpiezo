function h=showsolution(node,elem,u,varargin)
%Showsolution displays the solution corresponding to a mesh given by [node,elem] in 2-D.
%
% Copyright (C) Terence Yu. Modified by Bingbing Xu

data = node;
if ~iscell(elem)
    patch('Faces', elem,...
        'Vertices', data,...
        'FaceColor', 'interp',...
        'CData', u);
else
    max_n_vertices = max(cellfun(@length, elem));
    padding_func = @(vertex_ind) [vertex_ind,...
        NaN(1,max_n_vertices-length(vertex_ind))];  % function to pad the vacancies
    tpad = cellfun(padding_func, elem, 'UniformOutput', false);
    tpad = vertcat(tpad{:});
    h=patch('Faces', tpad,...
        'Vertices', data,'EdgeColor','k',...
        'FaceColor', 'interp',...
        'CData', u);
end
axis equal; 
% sh = 0.015;
% xlim([min(node(:,1)) - sh, max(node(:,1)) + sh])
% ylim([min(node(:,2)) - sh, max(node(:,2)) + sh])
% zlim([min(u) - sh, max(u) + sh])
% xlabel('x'); ylabel('y'); zlabel('u');

view(2); axis off;

colorbar;
end
