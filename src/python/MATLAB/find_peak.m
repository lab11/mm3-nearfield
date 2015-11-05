function [peak, m] = find_peak(C)
clear pos;
i = 1;
last_max = 0;
peak(1) = 0;
m = [];


    for n = 1: length(C)
        if(C(n) > 0.45)
            m = [m n];
            if(C(n) > last_max)
                last_max = C(n);
                peak(i) = n;
                V(i) = last_max;
            else
                if((n - peak(i)) >= 2000)
                    last_max = 0;
                    i = i+1;
                end
            end 
        else
            if(last_max > 0)
                if((n - peak(i)) > 2000)
                    last_max = 0;
                    i = i+1;
                end
            end
        end
    end
end


