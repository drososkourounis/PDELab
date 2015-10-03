function writeMeshToVTKFile(mesh, vtkfile)

% 2. read the .node file
% %%%%%%%%%%%%%%%%%%%%%%

% this opens a file in text 't mode for read access 'r'
filename = strcat(vtkfile, '.vtk');

fprintf('writting mesh file %s\n', filename);

fout = fopen(filename, 'w');
fprintf(fout,'# vtk DataFile Version 5.10\n');
fprintf(fout,'Hexahedral mesh with data\n');
fprintf(fout,'ASCII\n');
fprintf(fout,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fout,'POINTS %d float\n', mesh.N);


% now write the PointList
% -----------------------


for i = 1:mesh.N
    x_i = mesh.Points(i,:);
    fprintf(fout,'%25.16e %25.16e %25.16e\n', x_i(1), x_i(2), 0.0);
end

fprintf(fout,'\n');
entries = (mesh.N_v+1)*mesh.N_e;
fprintf(fout,'CELLS %d %d\n', mesh.N_e, entries);
first_number = 1;

  for e = 1:mesh.N_e
      
      v_e = mesh.Elements(e, :);
      v_e = v_e - first_number;
      
      fprintf(fout,'%d ', mesh.N_v);
      for i=1:mesh.N_v
          fprintf(fout,'%d ', v_e(i));
      end
      fprintf(fout, '\n');
  end


  fprintf(fout,'\n');
  fprintf(fout,'CELL_TYPES %d\n', mesh.N_e);

  type = 0;
  if (mesh.grid_type == 'triangles')
     type = 5;
  elseif (mesh.grid_type == 'quadrilaterals')
     type = 9;
  end

  
  for e = 1:mesh.N_e
      fprintf(fout,'%d\n', type);
  end

  fprintf(fout,'\n');
  fprintf(fout,'POINT_DATA %d\n', mesh.N);
  fprintf(fout,'SCALARS RealPressureField float 1\n');
  fprintf(fout,'LOOKUP_TABLE default\n');
  
   for e = 1:mesh.N_e
      fprintf(fout,'%25.16e\n', 1); 
   end

  fclose(fout);
end


