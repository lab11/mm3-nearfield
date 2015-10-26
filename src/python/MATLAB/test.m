
i = 1;
last_max = 0;
value(1) = 0;


for n = 1: length(C)
    if(C(n) > 1.5)
        if(C(n) > last_max)
            last_max = C(n);
        else
            if(i == 1) 
                value(i) = n;
                V(i) = last_max;
                last_max = 0;
                i = i+1;
            else
                if((n - value(i-1)) > 1000)
                    value(i) = n;
                    V(i) = last_max;
                    last_max = 0;
                    i = i+1;
                end
            end
        end 
    end
end
