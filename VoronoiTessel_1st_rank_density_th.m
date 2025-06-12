function [Data_voronoi, Data_voronoi_high, threshold] = VoronoiTessel_1st_rank_density_th(XY, thL)
%==========================================================================
% VoronoiTessel_1st_rank_density_th
%
% Author:        Nicolas Mateos
% Affiliation:   ICFO - The Institute of Photonic Sciences
% Email:         nicolas.mateos@icfo.eu
% GitHub:        https://github.com/nmateosHub
%
% Created:       2020-03-20
% Last Updated:  2025-06-12
%
%==========================================================================
% Description:
% This function computes the first-rank density of 2D points using Voronoi
% tessellation, normalizes it against a random distribution, and identifies
% high-density regions based on a statistical threshold.
%
% INPUTS:
%   XY   - Nx2 or Nx3 matrix with 2D coordinates (and optional intensity)
%   thL  - Threshold level (quantile) for identifying high-density points
%
% OUTPUTS:
%   Data_voronoi       - Voronoi-based metrics for input points:
%                        [X Y (optional Z) density norm. rank score]
%   Data_voronoi_high  - Subset of Data_voronoi above the density threshold
%   threshold          - Density threshold (in log10 scale) for filtering
%
% NOTE:
%   This function requires the auxiliary function:
%   AreaVoronoiRegions_1stRank(V, C, dt)
% 
% Handle edge cases: empty input or too few points
%
% License:
%
% This code is open-source and distributed under the MIT License.


if isempty(XY)
    Data_voronoi = [];
elseif size(XY,1) <= 2
    Data_voronoi = [];
else
    % Remove duplicate coordinates
    [~, ia, ~] = unique(XY(:,1:2), 'rows');
    XY = XY(ia, :);
    
    % Estimate point density
    Mx = max(XY(:,1)) - min(XY(:,1));
    My = max(XY(:,2)) - min(XY(:,2));
    delta = size(XY,1) / (Mx * My);  % Average density per unit area

    % Voronoi tessellation
    dt = delaunayTriangulation(XY(:,1), XY(:,2));
    [V, C] = voronoiDiagram(dt);

    % Compute Voronoi region areas and neighborhood metrics
    [Area, indxs, arearank, neigh] = AreaVoronoiRegions_1stRank(V, C, dt);
    
    % Calculate density and normalized score
    Data_voronoi = [XY(indxs,1:3), (1./Area)/delta, neigh./arearank, (neigh./arearank)/delta];
    
    % Remove NaN entries
    Data_voronoi = Data_voronoi(~isnan(Data_voronoi(:,6)), :);

    % Generate random data for comparison
    if ~isempty(Data_voronoi)
        xyrand = [Mx * rand(length(Data_voronoi), 1), My * rand(length(Data_voronoi), 1)];
        rdt = delaunayTriangulation(xyrand(:,1), xyrand(:,2));
        [rV, rC] = voronoiDiagram(rdt);
        [rArea, rindxs, rarearank, rneigh] = AreaVoronoiRegions_1stRank(rV, rC, rdt);

        % Compute metrics for random data
        randData_voronoi = [xyrand(rindxs,1:2), rArea, (1./rArea)/delta, rneigh./rarearank, (rneigh./rarearank)/delta];
        randData_voronoi = randData_voronoi(~isnan(randData_voronoi(:,6)), :);

        % Determine density threshold from empirical CDF
        if ~isempty(randData_voronoi)
            rlogdata = log10(randData_voronoi(:,6));
            [f, rx] = ecdf(rlogdata);           % Empirical cumulative distribution
            ax = find(f >= thL);                % Find quantile index
            threshold = rx(ax(1));              % Corresponding threshold value
        else
            threshold = 0;
        end
    else
        threshold = 0;
    end

    % Filter high-density points
    Data_voronoi_high = Data_voronoi;
    Data_voronoi_high(:,6) = log10(Data_voronoi_high(:,6));
    Data_voronoi_high = Data_voronoi_high(Data_voronoi_high(:,6) >= threshold, :);
end