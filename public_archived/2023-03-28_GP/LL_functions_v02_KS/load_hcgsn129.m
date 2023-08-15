function net_layout = load_hcgsn129()
    load('egi_hcgsn129.mat','sensor_layout')
    net_layout = sensor_layout;
end