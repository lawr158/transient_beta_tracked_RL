The only difference between original GP toolbox and the KS copy is found in 
support_functions > batch_convert_beapp_data

line 72 (original) : F.id = str2double(F.name(1:4));

Commented out, added : 

    id = strrep(F.name,'-','_');
    id = strrep(id,' ','_');
    id = erase(id,'.mat');
    F.id = str2double(id);