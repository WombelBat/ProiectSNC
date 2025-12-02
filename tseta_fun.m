function [tseta] = tseta_fun(delta)
    if delta == 0
     tseta = 1;
    else
        temp = log(delta / 100)^2;
        tseta = sqrt(temp / (pi^2 + temp));
    end
end