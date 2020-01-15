function r2 = calc_r2( ActualY, EstimatedY )
%CALC_R2 Calculate r2 

Rss = sum( ( repmat( ActualY, 1, size(EstimatedY,2) ) - EstimatedY ).^2 );


r2 = 1 - Rss / sum((ActualY-mean(ActualY)).^2);


end

