% Version 2.2.3

clear;

% 0. declare reference variables
HYDRO = 1;
ELBDM = 3;


% 1. input parameters and set the default values
Model        = input( 'Model (HYDRO=1, ELBDM=3) [1] ? ' );
if isempty( Model )
Model        = 1;
end

MultiFiles   = input( 'Process multiple files (0/1->No/Yes) [0] ? ' );
if isempty( MultiFiles )
MultiFiles   = 0;
end

if ( MultiFiles )
FilePrefix   = input( 'File name prefix ? ', 's' );
FileID_Start = input( 'File start-ID ? ' );
FileID_Stop  = input( 'File stop-ID ? ' );
FileID_Delta = input( 'File delta-ID [1] ? ' );
FileSuffix   = input( 'File name suffix ? ', 's' );

if isempty( FileID_Delta )
FileID_Delta = 1;
end

else % if ( MultiFiles )

FileName     = input( 'FileName ? ', 's' );
FileID_Start = 0;
FileID_Stop  = 0;
FileID_Delta = 1;
end % if ( MultiFiles ) ... else ...

LoadBinary   = input( 'Load binary data (0/1/2->No/float4/float8) [0] ? ' );
if isempty( LoadBinary )
LoadBinary   = 0;
end

XYZ          = input( '(x[1] / y[2] / z[3]) plane [3] ? ' );
if isempty( XYZ )
XYZ          = 3;
end

if ( LoadBinary )
switch XYZ
   case 1
   Size1     = input( 'NY ? ' );
   Size2     = input( 'NZ ? ' );
   case 2
   Size1     = input( 'NZ ? ' );
   Size2     = input( 'NX ? ' );
   case 3
   Size1     = input( 'NX ? ' );
   Size2     = input( 'NY ? ' );
   otherwise
   error( 'ERROR : please provide the correct targeted plane [1/2/3] !!' );
end

else % if ( LoadBinary )

NHeader      = input( 'The number of header lines [1] ? ' );
if isempty( NHeader )
NHeader      = 1; 
end

WithCoord    = input( 'New output format with physical coordinates (0/1->No/Yes) [1] ? ' );
if isempty( WithCoord )
WithCoord    = 1;
end
end % if ( LoadBinary ) ... else ...

if ( LoadBinary  ||  ~WithCoord )
dh           = input( 'Cell physical size [cell scale in ASCII input / 1 in binary input] ? ' );
end
if ( LoadBinary && isempty(dh)  )
dh           = 1;
end

SaveFig      = input( 'Save figures (0/1/2->No/fig/png) [0] ? ' );
if isempty( SaveFig )
SaveFig      = 0;
end

if ( MultiFiles )
SaveMov      = input( 'Save an AVI movie (0/1->No/Yes) [0] ? ' );
else
SaveMov      = 0;
end
if isempty( SaveMov )
SaveMov      = 0;
end

Log10Plot    = input( 'Plot in log10 scale (0/1->No/Yes) [0] ? ' );
if isempty( Log10Plot )
Log10Plot    = 0;
end

if ( LoadBinary )
CMax         = input( 'Upper limit of colorbar [no limit] ? ' );
if ( ~isempty(CMax) )
CMin         = input( 'Lower limit of colorbar ? ' );
end
else
CMax         = [];   % set it as a empty array
end

Title0       = input( 'Title ? ', 's' );


% check the input parameters
if ( Model ~= HYDRO  &&  Model ~= ELBDM )
   error( 'ERROR : please provide the correct model [1/3] !!' );
end

if ( MultiFiles )

   if isempty( FilePrefix )
      error( 'ERROR : please provide the file name prefix !!' );
   end

   if isempty( FileID_Start )
      error( 'ERROR : please provide the starting file ID !!' );
   end

   if isempty( FileID_Stop )
      error( 'ERROR : please provide the stoping file ID !!' );
   end

else % if ( MultiFiles )

   if isempty( FileName )
      error( 'ERROR : please provide the file name !!' );
   end

end % if ( MultiFiles ) ... else ...

if ( LoadBinary < 0  ||  LoadBinary > 2 )
   error( 'ERROR : please provide the correct option for loading binary data [0/1/2] !!' );
end

if ( LoadBinary )
   if isempty( Size1 )
      error( 'ERROR : please provide the horizontal size of the input binary file !!' );
   end

   if isempty( Size2 )
      error( 'ERROR : please provide the vertical size of the input binary file !!' );
   end
else
   if ( NHeader < 0 )
      error( 'ERROR : please provide the correct number of header lines [>=0] !!' );
   end
end

if ( XYZ < 1  ||  XYZ > 3 )
   error( 'ERROR : please provide the correct targeted plane [1/2/3] !!' );
end

if ( exist( 'dh' )  &&  ~isempty(dh)  &&  dh <= 0.0 )
   error( 'ERROR : cell physical size (%f) <= 0.0 !!', dh );
end

if ( SaveFig < 0  ||  SaveFig > 2 )
   error( 'ERROR : please provide the correct option for saving figures [0/1/2] !!' );
end

if ( ~isempty(CMax)  &&  isempty(CMin) )
   error( 'ERROR : please provide the lower limit of colorbar !!' );
end


% create the avi object
if ( SaveMov )
   MovName = [ FilePrefix FileSuffix ];
   Mov = avifile( MovName, 'fps', 1 );
end


% loop over all input files
for FileID=FileID_Start: FileID_Delta: FileID_Stop

%  2. load data
%  2-1. set the file name and title
   if ( MultiFiles )

      fprintf( 1, '---------------------------\n' );
      fprintf( 1, 'Processing file ID %6d\n', FileID );
      fprintf( 1, '---------------------------\n' );

      ID0=int2str(  floor( FileID/100000 )             );
      ID1=int2str(  floor( mod(FileID,100000)/10000 )  );
      ID2=int2str(  floor( mod(FileID,10000)/1000 )    );
      ID3=int2str(  floor( mod(FileID,1000)/100 )      );
      ID4=int2str(  floor( mod(FileID,100)/10 )        );
      ID5=int2str(  floor( mod(FileID,10) )            );

      FileName = [ FilePrefix '_' ID0 ID1 ID2 ID3 ID4 ID5 FileSuffix ];

      Title    = [ Title0 ' ( FileID ' ID0 ID1 ID2 ID3 ID4 ID5 ' )' ];

   else
      Title = Title0;
   end


%  2-2. begin to load data
   fprintf( 1, 'loading data from %s\n', FileName );

   FID = fopen( FileName, 'r' );
   
   if ( FID == -1 )
      error( 'ERROR : file %s is not found !!', FileName );
   end
   
   if ( LoadBinary )
      if ( LoadBinary == 1 ) % single precision
         format = 'float32';
      else                   % double precision
         format = 'float64';
      end

      if ( XYZ ~= 2 )
         Data = fread( FID, [Size1 Size2], format );
         Data = Data'; % transpose the data of xy and yz planes to give consistent figure produced by ASCII input
      else
         Data = fread( FID, [Size2 Size1], format );
      end

      clear format;

   else % if ( LoadBinary )

      if ( Model == HYDRO )
         Data_tp = textscan( FID, '%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Headerlines', NHeader, 'CollectOutput', 1 );
      else 
         Data_tp = textscan( FID, '%f%f%f%f%f%f%f%f%f%f',       'Headerlines', NHeader, 'CollectOutput', 1 );
      end
      
      Data = Data_tp{1};
      clear Data_tp;
   end % if ( LoadBinary ) ... else ...

   fclose( FID );
   
   

%  3. set all parameters
   fprintf( 1, 'setting parameters\n' );
   
   if ( LoadBinary )

%     set plot range and axis label for creating figures
      Range1 = [ 0.5, Size1-0.5 ] * dh;
      Range2 = [ 0.5, Size2-0.5 ] * dh;
      Label1 = char( 120+mod(XYZ  ,3) );
      Label2 = char( 120+mod(XYZ+1,3) );

   else % if ( LoadBinary )

%     3-0. set variable indices
      if ( WithCoord )  dIdx = 3;
      else              dIdx = 0;
      end
      
      if ( Model == HYDRO )
         DENS =  4 + dIdx;
         MOMX =  5 + dIdx;
         MOMY =  6 + dIdx;
         MOMZ =  7 + dIdx;
         ENGY =  8 + dIdx;
         PRES =  9 + dIdx;
         POTE = 10 + dIdx;
      else
         DENS =  4 + dIdx;
         REAL =  5 + dIdx;
         IMAG =  6 + dIdx;
         POTE =  7 + dIdx;
      end
      

%     3-1. set data range
      DataSize = size( Data, 1 );
      Sort_ijk = sort( Data(:,1:3) );
      Min_ijk  = Sort_ijk(        1, 1:3 );
      Max_ijk  = Sort_ijk( DataSize, 1:3 );
      TID1     = mod(  XYZ, 3 ) + 1;
      TID2     = mod( TID1, 3 ) + 1;
      

%     3-2. get spatial interval
      dScale = 0;
      for t=1:DataSize
         if ( Sort_ijk(t,TID1) > Min_ijk(TID1) )
            dScale = Sort_ijk(t,TID1) - Min_ijk(TID1);
            break;
         end
      end
      
      if ( WithCoord )

         Sort_xyz = sort( Data(:,4:6) );
         Min_xyz  = Sort_xyz( 1, 1:3 );
      
         dh = 0;
         for t=1:DataSize
            if ( Sort_xyz(t,TID1) > Min_xyz(TID1) )
               dh = Sort_xyz(t,TID1) - Min_xyz(TID1);
               break;
            end
         end

      else % if ( WithCoord )
      
         if ( isempty( dh ) )
            dh = dScale;
         end

      end % if ( WithCoord ) ... else ...
      

%     3-3. set array size
      Size1 = ( Max_ijk(TID1) - Min_ijk(TID1) ) / dScale + 1;
      Size2 = ( Max_ijk(TID2) - Min_ijk(TID2) ) / dScale + 1;
      Dim1  = ( Data(:,TID1)  - Min_ijk(TID1) ) / dScale + 1;
      Dim2  = ( Data(:,TID2)  - Min_ijk(TID2) ) / dScale + 1;
      

%     3-4. set plot range and axis label for creating figures
      Range1 = [ Min_ijk(TID1)+0.5*dScale, Max_ijk(TID1)+0.5*dScale ] / dScale * dh;
      Range2 = [ Min_ijk(TID2)+0.5*dScale, Max_ijk(TID2)+0.5*dScale ] / dScale * dh;
      Label1 = char( 120+mod(XYZ  ,3) );
      Label2 = char( 120+mod(XYZ+1,3) );

      
%     3-5. set variables
      fprintf( 1, 'setting variables\n' );
      
      if ( Model == HYDRO )
         Dens = zeros(Size2, Size1);
         VelX = zeros(Size2, Size1);
         VelY = zeros(Size2, Size1);
         VelZ = zeros(Size2, Size1);
         Engy = zeros(Size2, Size1);
         Pres = zeros(Size2, Size1);
         Pote = zeros(Size2, Size1);
      
         for t=1:DataSize
            Dens( Dim2(t), Dim1(t) ) = Data(t,DENS);
            VelX( Dim2(t), Dim1(t) ) = Data(t,MOMX)/Data(t,DENS);
            VelY( Dim2(t), Dim1(t) ) = Data(t,MOMY)/Data(t,DENS);
            VelZ( Dim2(t), Dim1(t) ) = Data(t,MOMZ)/Data(t,DENS);
            Engy( Dim2(t), Dim1(t) ) = Data(t,ENGY);
            Pres( Dim2(t), Dim1(t) ) = Data(t,PRES);
            Pote( Dim2(t), Dim1(t) ) = Data(t,POTE);
         end
      else
         Dens = zeros(Size2, Size1);
         Real = zeros(Size2, Size1);
         Imag = zeros(Size2, Size1);
         Pote = zeros(Size2, Size1);
      
         for t=1:DataSize
            Dens( Dim2(t), Dim1(t) ) = Data(t,DENS);
            Real( Dim2(t), Dim1(t) ) = Data(t,REAL);
            Imag( Dim2(t), Dim1(t) ) = Data(t,IMAG);
            Pote( Dim2(t), Dim1(t) ) = Data(t,POTE);
         end
      end

   end % if ( LoadBinary ) ... else ...



%  4. plot
   fprintf( 1, 'plotting\n' );
   
%  get the screen size for full-screen display
   ScreenSize = get( 0, 'ScreenSize' );

   if ( LoadBinary )

      figure( 'OuterPosition', ScreenSize );

      if isempty( CMax )
         if ( Log10Plot )
            imagesc( Range1, Range2, log10(Data) );
         else
            imagesc( Range1, Range2,       Data );
         end
      else
         if ( Log10Plot )
            imagesc( Range1, Range2, log10(Data), [CMin CMax] );
         else
            imagesc( Range1, Range2,       Data,  [CMin CMax] );
         end
      end

      if ~isempty( Title )    title( Title );   end
      axis equal; colorbar; set( gca, 'YDir', 'normal' ); xlabel( Label1 ); ylabel( Label2, 'Rotation', 0 );

      if ( SaveFig )
         if     ( SaveFig == 1 )    
            FigName = [ FileName '.fig' ];
            saveas( gcf, FigName, 'fig' );
         elseif ( SaveFig == 2 )    
            FigName = [ FileName '.png' ];
            saveas( gcf, FigName, 'png' );
         end
      end

      if ( SaveMov )    Mov = addframe( Mov, getframe(gcf) );  end

   else % if ( LoadBinary )

      if ( Model == HYDRO )
      
%        all variables
         figure( 'OuterPosition', ScreenSize );    if ~isempty( Title )    suptitle( Title );   end
         
         subplot(2,3,1);   imagesc( Range1, Range2, Dens );    title( 'Density' );
         axis equal; colorbar; set( gca, 'YDir', 'normal' ); xlabel( Label1 ); ylabel( Label2, 'Rotation', 0 );
         
         subplot(2,3,2);   imagesc( Range1, Range2, Engy );    title( 'Energy' );
         axis equal; colorbar; set( gca, 'YDir', 'normal' ); xlabel( Label1 ); ylabel( Label2, 'Rotation', 0 );
         
         subplot(2,3,3);   imagesc( Range1, Range2, Pres );    title( 'Pressure' );
         axis equal; colorbar; set( gca, 'YDir', 'normal' ); xlabel( Label1 ); ylabel( Label2, 'Rotation', 0 );
         
         subplot(2,3,4);   imagesc( Range1, Range2, VelX );    title( 'Vx' );
         axis equal; colorbar; set( gca, 'YDir', 'normal' ); xlabel( Label1 ); ylabel( Label2, 'Rotation', 0 );
         
         subplot(2,3,5);   imagesc( Range1, Range2, VelY );    title( 'Vy' );
         axis equal; colorbar; set( gca, 'YDir', 'normal' ); xlabel( Label1 ); ylabel( Label2, 'Rotation', 0 );
         
         subplot(2,3,6);   imagesc( Range1, Range2, VelZ );    title( 'Vz' );
         axis equal; colorbar; set( gca, 'YDir', 'normal' ); xlabel( Label1 ); ylabel( Label2, 'Rotation', 0 );
         

%        density
         figure( 'OuterPosition', ScreenSize );    if ~isempty( Title )    suptitle( Title );   end
         if ( Log10Plot )
            imagesc( Range1, Range2, log10(Dens) );
         else
            imagesc( Range1, Range2, Dens );
         end
         title( 'Density' );
         axis equal; colorbar; set( gca, 'YDir', 'normal' ); xlabel( Label1 ); ylabel( Label2, 'Rotation', 0 );

         if ( SaveFig )
            if     ( SaveFig == 1 )    
               FigName = [ FileName '_Dens' '.fig' ];
               saveas( gcf, FigName, 'fig' );
            elseif ( SaveFig == 2 )    
               FigName = [ FileName '_Dens' '.png' ];
               saveas( gcf, FigName, 'png' );
            end
         end

         if ( SaveMov )    Mov = addframe( Mov, getframe(gcf) );  end

         
%        density contour
         figure( 'OuterPosition', ScreenSize );    if ~isempty( Title )    suptitle( Title );   end
         x = Range1(1):dh:Range1(2);
         y = Range2(1):dh:Range2(2);
         if ( Log10Plot )
            contour( x, y, log10(Dens) );
         else
            contour( x, y, Dens );
         end
         title( 'Density contour' );
         axis equal; colorbar; set( gca, 'YDir', 'normal' ); xlabel( Label1 ); ylabel( Label2, 'Rotation', 0 );
         

%        potential
%        figure( 'OuterPosition', ScreenSize );    if ~isempty( Title )    suptitle( Title );   end
%        imagesc( Range1, Range2, Pote );
%        title( 'Potential' );
%        axis equal; colorbar; set( gca, 'YDir', 'normal' ); xlabel( Label1 ); ylabel( Label2, 'Rotation', 0 );
      

      else % ELBDM
      

%        density
         figure( 'OuterPosition', ScreenSize );    if ~isempty( Title )    suptitle( Title );   end
         if ( Log10Plot )
            imagesc( Range1, Range2, log10(Dens) );
         else
            imagesc( Range1, Range2, Dens );
         end
         title( 'Density' );
         axis equal; colorbar; set( gca, 'YDir', 'normal' ); xlabel( Label1 ); ylabel( Label2, 'Rotation', 0 );

         if ( SaveFig )
            if     ( SaveFig == 1 )    
               FigName = [ FileName '_Dens' '.fig' ];
               saveas( gcf, FigName, 'fig' );
            elseif ( SaveFig == 2 )    
               FigName = [ FileName '_Dens' '.png' ];
               saveas( gcf, FigName, 'png' );
            end
         end

         if ( SaveMov )    Mov = addframe( Mov, getframe(gcf) );  end
      

%        real part
%        figure( 'OuterPosition', ScreenSize );    if ~isempty( Title )    suptitle( Title );   end
%        imagesc( Range1, Range2, Real );
%        title( 'Real' );
%        axis equal; colorbar; set( gca, 'YDir', 'normal' ); xlabel( Label1 ); ylabel( Label2, 'Rotation', 0 );
      

%        imag part
%        figure( 'OuterPosition', ScreenSize );    if ~isempty( Title )    suptitle( Title );   end
%        imagesc( Range1, Range2, Imag );
%        title( 'Imaginary' );
%        axis equal; colorbar; set( gca, 'YDir', 'normal' ); xlabel( Label1 ); ylabel( Label2, 'Rotation', 0 );
      
      end % if ( Model == HYDRO ) ... else ...

   end % if ( LoadBinary ) ... else ...

end % for FileID=FileID_Start: FileID_Delta: FileID_Stop


% close the avi object
if ( SaveMov )    Mov = close( Mov );  end


fprintf( 1, '\n~~ program finished successfully ~~\n\n' );

