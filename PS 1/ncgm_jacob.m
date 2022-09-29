function jacob = ncgm_jacob(x0, ncgm_par)

jacob = zeros(length(x0), length(x0));

for i = 1:length(x0)
    dx = zeros(length(x0),1);
    dx(i) = x0(i)*1e-3;
    x1 = x0 + dx;
    jacob(:,i) = (rbc_obj_start(x1, ncgm_par)- rbc_obj_start(x0, ncgm_par))/ dx(i);
end
