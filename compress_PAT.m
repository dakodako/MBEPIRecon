new_image_data_PAT_meas1_coil1 = zeros(98,43);
n = 1;
for i = 1:size(image_data_PAT_meas1_coil1,2)
 
    if(mod(i,2) ~= 1)
        new_image_data_PAT_meas1_coil1(:,n) = image_data_PAT_meas1_coil1(:,i);
        n = n + 1;
    end
    
end