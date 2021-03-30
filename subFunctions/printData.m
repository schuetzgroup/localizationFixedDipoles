function [] = printData(x)

def_grid = length(x);

for k=0:def_grid-1
    disp(['(',num2str(1000/(def_grid-1)*k-500), ',',num2str(x(k+1)),')' ])
end

end