% user-defined parameters
max_mode = 2;
dataset_name = 'charge';
phi = 0;
n_files = 2500;
step_size = 10;
src_dir = './Charge';
des_dir = './Charge_combine';

% set hdf5 file handle for m=0 mode.
hfile_0.name = [src_dir, '/Re0/', dataset_name, '_(*).h5'];
hfile_0.from = 0;
hfile_0.to	 = n_files;
hfile_0.padwidth = 8;
hfile_0.step = step_size;

hfile_re = cell( max_mode );
hfile_im = cell( max_mode );

% set hdf5 file handle for m>0 mode.
for ii = 1:max_mode
  hfile_re{ii}.name = [src_dir, '/Re', num2str(ii), '/', dataset_name, '_(*).h5']; % user-defined parameters
  hfile_re{ii}.from = hfile_0.from;
  hfile_re{ii}.to   = hfile_0.to;
  hfile_re{ii}.padwidth = 8;
  hfile_re{ii}.step = step_size;

  hfile_im{ii}.name = [src_dir, '/Im', num2str(ii), '/', dataset_name, '_(*).h5']; % user-defined parameters
  hfile_im{ii}.from = hfile_0.from;
  hfile_im{ii}.to   = hfile_0.to;
  hfile_im{ii}.padwidth = 8;
  hfile_im{ii}.step = step_size;
end

hfile_out.name = [des_dir, '/', dataset_name, '_(*).h5']; % user-defined parameters
hfile_out.padwidth = 8;
hfile_out.step = step_size;
hfile_out.from = hfile_0.from;
hfile_out.to = hfile_0.to;

if max_mode > 0
  file_name_re = cell( max_mode );
  file_name_im = cell( max_mode );
end

file_name_0 = genFileName( hfile_0, 1 );

for ii = hfile_0.from:hfile_0.step:hfile_0.to

  disp( ['step = ', num2str(ii)] )

  file_name_0 = genFileName( hfile_0, ii );
  finfo = h5info( file_name_0 );
  natt = size( finfo.Attributes, 1 );
  natt_dataset = size( finfo.Datasets(1).Attributes, 1 );
  f_half1 = h5read( file_name_0, ['/',dataset_name] );
  f_half2 = f_half1;
  r = h5read( file_name_0, '/AXIS/AXIS1' );
  z = h5read( file_name_0, '/AXIS/AXIS2' );

  if max_mode > 0
    for m = 1:max_mode
      file_name_re{m} = genFileName( hfile_re{m}, ii );
      file_name_im{m} = genFileName( hfile_im{m}, ii );
      f_re = h5read( file_name_re{m}, ['/',dataset_name] );
      f_im = h5read( file_name_im{m}, ['/',dataset_name] );
      f_half1 = f_half1 + 2 * ( f_re * cos(m*phi) - f_im * sin(m*phi) );
      f_half2 = f_half2 + 2 * ( f_re * cos(m*(phi+pi)) - f_im * sin(m*(phi+pi)) );
    end
  end

  f_comb = [ flipud(f_half2(2:end,:)); f_half1 ];

  % get hdf5 file info
  file_name_out = genFileName( hfile_out, ii );
  h5create( file_name_out, ['/',dataset_name], size(f_comb) );
  h5write( file_name_out, ['/',dataset_name], f_comb );

  % write attribute to root
  for k = 1:natt
    val = finfo.Attributes(k).Value;
    if iscell( val )
      h5writeatt( file_name_out, '/', ...
       finfo.Attributes(k).Name, strtrim(finfo.Attributes(k).Value{:}) );
    else
      h5writeatt( file_name_out, '/', ...
       finfo.Attributes(k).Name, finfo.Attributes(k).Value );
    end
  end

  % write dataset attributes
  for k = 1:natt_dataset
    val = finfo.Datasets(1).Attributes(k).Value;
    if iscell( val )
      h5writeatt( file_name_out, ['/',dataset_name], ...
       finfo.Datasets(1).Attributes(k).Name, ...
       strtrim(finfo.Datasets(1).Attributes(k).Value{:}) );
    else
      h5writeatt( file_name_out, ['/',dataset_name], ...
       finfo.Datasets(1).Attributes(k).Name, ...
       finfo.Datasets(1).Attributes(k).Value );
    end
  end

  h5create( file_name_out, '/AXIS/AXIS1', 2 );
  h5create( file_name_out, '/AXIS/AXIS2', 2 );
  h5write( file_name_out, '/AXIS/AXIS1', [-r(2);r(2)] );
  h5write( file_name_out, '/AXIS/AXIS2', z );
  h5writeatt( file_name_out, '/AXIS/AXIS1', 'TYPE', 'linear' );
  h5writeatt( file_name_out, '/AXIS/AXIS1', 'UNITS', 'c / \omega_p' );
  h5writeatt( file_name_out, '/AXIS/AXIS1', 'NAME', 'x' );
  h5writeatt( file_name_out, '/AXIS/AXIS1', 'LONG_NAME', 'x' );
  h5writeatt( file_name_out, '/AXIS/AXIS2', 'TYPE', 'linear' );
  h5writeatt( file_name_out, '/AXIS/AXIS2', 'UNITS', 'c / \omega_p' );
  h5writeatt( file_name_out, '/AXIS/AXIS2', 'NAME', '\xi' );
  h5writeatt( file_name_out, '/AXIS/AXIS2', 'LONG_NAME', '\xi' );

end