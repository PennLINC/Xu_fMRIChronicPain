% Read in coordinates
coordinates=readtable('coordinates_20200131.csv');

% Setup new columns for MNI space XYZ
coordinates.mni_x=NaN([149 1]);
coordinates.mni_y=NaN([149 1]);
coordinates.mni_z=NaN([149 1]);

% Convert anything Tal to MNI
for i=1:149
    if strcmp(coordinates.space(i), 'MNI')
        coordinates.mni_x(i)=coordinates.x(i);
        coordinates.mni_y(i)=coordinates.y(i);
        coordinates.mni_z(i)=coordinates.z(i);
    else
        % Temporarily save Tal coordinates and transformed MNI coordinates
        tal_coords = ([coordinates.x(i) coordinates.y(i) coordinates.z(i)]);
        mni_coords = my_tal2icbm_spm(tal_coords);
        % Fill in coordinates with transformed ones
        coordinates.mni_x(i) = mni_coords(1);
        coordinates.mni_y(i) = mni_coords(2);
        coordinates.mni_z(i) = mni_coords(3);
    end
end 

% Export
writetable(coordinates, 'coordinatesMNI_20200131.csv', 'Delimiter', ',');