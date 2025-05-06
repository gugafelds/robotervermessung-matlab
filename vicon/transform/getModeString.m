function mode_str = getModeString(transform_only, upload_single, upload_all)
    if transform_only
        mode_str = 'Nur Transformation';
    elseif upload_single
        mode_str = 'Upload einzelne Bahn';
    elseif upload_all
        mode_str = 'Upload alle Bahnen';
    else
        mode_str = 'Kein Modus ausgewählt';
    end
end

