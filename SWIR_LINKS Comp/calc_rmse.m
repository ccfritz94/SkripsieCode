function Rmse = calc_rmse( ActualY, EstimatedY )
%CALC_RMSE Calculate RMSE 

Rmse = sqrt( sum( ( repmat( ActualY, 1, size(EstimatedY,2) ) - EstimatedY ).^2 ) / length(ActualY) );

end

