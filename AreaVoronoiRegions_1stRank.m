function [Area,idx,arearank,neigh] = AreaVoronoiRegions_1stRank(V,C,dt)
    % This code computes the voronoi areas and the 1st rank voronoi areas:
    % The indices for the localizations:
    idx = (1:length(C))';
    % For better performance, I already allocate a matrix for the Voronoi
    % areas:
    Area = ones(length(C),1);
    % I here take the Voronoi neighbors in order to compute the 1st Rank
    % Voronoi areas:
    Ed = edges(dt);
    % In order to compute the 1st rank voronoi areas, I will calculate it with
    % a multiplication of matrices: arearank = M*Area:
    K = [idx idx; Ed;fliplr(Ed)];
    M = sparse(K(:,1),K(:,2),1);
    % Computing the voronoi areas for each Voronoi cell:
    parfor i = 1:length(C)
        Area(i) = Get_polyarea(V,C,i);
    
    end
    % Calculate the 1st rank voronoi areas and the number of neighbors+1:
    % Don't take into account the NaN (open areas):
    Area_nNaN = Area; Area_nNaN(isnan(Area))=0;
    arearank = M*Area_nNaN;
    % Those that had open areas, its area1strank continues to be NaN:
    arearank(isnan(Area))=NaN;
    neigh=full(sum(M,2))-1;
