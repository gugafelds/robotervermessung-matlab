function [transformed_data] = coord_transformation(data_ist, trafo_rot, trafo_trans)
    if istable(data_ist)
        pos_ist = table2array(data_ist(1,2:4));
        pos_ist = [pos_ist{1,1} pos_ist{1,2} pos_ist{1,3}];
        pos_ist_trafo = pos_ist * trafo_rot + trafo_trans;
        data_ist(:,2:4) = array2table([{pos_ist_trafo(:,1)} {pos_ist_trafo(:,2)} {pos_ist_trafo(:,3)}]);
        transformed_data = data_ist;
    elseif size(data_ist,2) == 3
        transformed_data = data_ist * trafo_rot + trafo_trans;
    else
        error('Unerwartetes Datenformat für Koordinatentransformation');
    end
end