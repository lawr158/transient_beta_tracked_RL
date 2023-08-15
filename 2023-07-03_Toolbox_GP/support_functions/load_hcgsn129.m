function [net_layout,net_1020] = load_hcgsn129()
    load('egi_hcgsn129.mat','sensor_layout')
    net_layout = sensor_layout;
    try
        load('10_20_cap19.mat','sensor_layout')
        net_1020 = sensor_layout;
    catch
        net_1020 = [];
    end
end