function data_ist_trafo = coordTransformation(data_ist,trafo_rot, trafo_trans)

if istable(data_ist)
    pos_ist = table2array(data_ist(1,2:4));
    pos_ist = [pos_ist{1,1} pos_ist{1,2} pos_ist{1,3}];
    pos_ist_trafo = pos_ist * trafo_rot + trafo_trans;
    data_ist(:,2:4) = array2table([{pos_ist_trafo(:,1)} {pos_ist_trafo(:,2)} {pos_ist_trafo(:,3)} ]);
    assignin("base","pos_ist_trafo",data_ist)
elseif size(data_ist,2) == 3 
    data_ist = data_ist * trafo_rot + trafo_trans;
    data_ist_trafo = data_ist;
end
