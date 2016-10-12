%Creat atom postion data file

clear all



%water
a_w = 3.111582;     %T = 300K

nx_w = 10;
lx_w = nx_w*a_w;
ny_w = 10;
ly_w = ny_w*a_w;
nz_w = 10;
lz_w = nz_w*a_w;

%graphene
lx0 = 150;
ly0 = 150;

a_g = 1.42;

nx_g = ceil(lx0/3/a_g);
lx = nx_g*3*a_g;
ny_g = ceil(ly0/sqrt(3)/a_g);
ly = ny_g*sqrt(3)*a_g;

natom = 3*nx_w*ny_w*nz_w + 2*nx_g*2*ny_g;
nbond = 2*nx_w*ny_w*nz_w;
nangle = nx_w*ny_w*nz_w;

shiftx = (lx - lx_w)/2;
shifty = (ly - ly_w)/2;
shiftz = 0.5;

n = 1;
m = 1;
l = 1;


%water
for i = 1:nx_w
     for j = 1:ny_w
           for k = 1:nz_w
               
               
            Bond1(m) = n;         %Bonds
            Bond2(m) = n+1; 
            m = m + 1;
            Bond1(m) = n;
            Bond2(m) = n+2;
            m = m + 1;
            
            
            Angle1(l) = n+1;      %Angles
            Angle2(l) = n;
            Angle3(l) = n+2;
            l = l + 1;
            
        
            x(n) = shiftx + 1.555791 + (i-1)*a_w;    %Oxygen atom
            y(n) = shifty + 1.555791 + (j-1)*a_w;
            z(n) = shiftz + 1.5 + (k-1)*a_w;
             
            Atomn(n) = n;
            Moleculart(n) = 1;
            Atomt(n) = 1;
            Charge(n) = -1.1128;
            n = n + 1; 
           
            x(n) = shiftx + 1.555791 + (i-1)*a_w;    %hydrogen atom
            y(n) = shifty + 2.312741 + (j-1)*a_w;
            z(n) = shiftz + 2.085882 + (k-1)*a_w;
            
            Atomn(n) = n;
            Moleculart(n) = 1;
            Atomt(n) = 2;
            Charge(n) = 0.5564;
            n = n + 1; 

            x(n) = shiftx + 1.555791 + (i-1)*a_w;    %hydrogen atom
            y(n) = shifty + 0.798841 + (j-1)*a_w;
            z(n) = shiftz + 2.085882 + (k-1)*a_w;
            
            Atomn(n) = n;
            Moleculart(n) = 1;
            Atomt(n) = 2;
            Charge(n) = 0.5564;
            n = n + 1; 
            
            

           end
    end
end

%graphene
for j = 1:2*ny_g
     for i = 1:2*nx_g

        
        if(mod(j,2)==1)
            x(n) = (i - 1)*3/2*a_g + a_g/2 - (mod(i,2)==0)*a_g/2;
        else
            x(n) = (i - 1)*3/2*a_g + (mod(i,2)==0)*a_g/2;
        end
        
        y(n) = sqrt(3)/2*a_g*j;
        
        z(n) = 0;
        
        Atomn(n) = n;
        Moleculart(n) = 2;
        Atomt(n) = 3;
        Charge(n) = 0;
        n = n + 1;   
    end
end



%Output
fid=fopen('data.txt', 'wt'); 
fprintf(fid, '#water_graphen data ');
fprintf(fid,'%d',lx0); 
fprintf(fid, ' * ');
fprintf(fid,'%d',ly0);
fprintf(fid, ' * ');
fprintf(fid,'%d',50);
fprintf(fid, '\n');
fprintf(fid, '#water density =9.93000e+02 ');
fprintf(fid, '\n');
fprintf(fid, '\n');

fprintf(fid,'%d',natom);
fprintf(fid, ' atoms');
fprintf(fid, '\n');

fprintf(fid,'%d',nbond);
fprintf(fid, ' bonds');
fprintf(fid, '\n');

fprintf(fid,'%d',nangle);
fprintf(fid, ' angles');
fprintf(fid, '\n');

fprintf(fid,'%d',3);
fprintf(fid, ' atom types');
fprintf(fid, '\n');

fprintf(fid,'%d',1);
fprintf(fid, ' bond types');
fprintf(fid, '\n');

fprintf(fid,'%d',1);
fprintf(fid, ' angle types');
fprintf(fid, '\n');
fprintf(fid, '\n');

fprintf(fid,'%d',0);
fprintf(fid, '   ');
fprintf(fid,'%10.6f',lx);
fprintf(fid, '  xlo xhi');
fprintf(fid, '\n');

fprintf(fid,'%d',0);
fprintf(fid, '   ');
fprintf(fid,'%10.6f',ly);
fprintf(fid, '  ylo yhi');
fprintf(fid, '\n');

fprintf(fid,'%d',-10);
fprintf(fid, '   ');
fprintf(fid,'%10.6f',70);
fprintf(fid, '  zlo zhi');
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '\n');

fprintf(fid, 'Masses');
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '1   15.9994 \n');
fprintf(fid, '2   1.0080 \n');
fprintf(fid, '3   12');
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '\n');

fprintf(fid, 'Atoms');
fprintf(fid, '\n');
fprintf(fid, '\n');

for n = 1:natom
    
    fprintf(fid,'%d',Atomn(n));
    fprintf(fid, '   ');
    fprintf(fid,'%d',Moleculart(n));
    fprintf(fid, '   ');
    fprintf(fid,'%d',Atomt(n));
    fprintf(fid, '   ');
    fprintf(fid,'%d',Charge(n));
    fprintf(fid, '   ');
    fprintf(fid,'%10.6e',x(n));
    fprintf(fid, '   ');
    fprintf(fid,'%10.6e',y(n));
    fprintf(fid, '   ');
    fprintf(fid,'%10.6e',z(n));
    fprintf(fid, '\n');    
 
end

fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, 'Bonds');
fprintf(fid, '\n');
fprintf(fid, '\n');

for n = 1:nbond
    
    fprintf(fid,'%d',n);
    fprintf(fid, '   ');
    fprintf(fid,'%d',Moleculart(n));
    fprintf(fid, '   ');
    fprintf(fid,'%d',Bond1(n));
    fprintf(fid, '   ');
    fprintf(fid,'%d',Bond2(n));
    fprintf(fid, '\n'); 
    
end


fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, 'Angles');
fprintf(fid, '\n');
fprintf(fid, '\n');


for n = 1:nangle
    
    fprintf(fid,'%d',n);
    fprintf(fid, '   ');
    fprintf(fid,'%d',Moleculart(n));
    fprintf(fid, '   ');
    fprintf(fid,'%d',Angle1(n));
    fprintf(fid, '   ');
    fprintf(fid,'%d',Angle2(n));
    fprintf(fid, '   '); 
    fprintf(fid,'%d',Angle3(n));
    fprintf(fid, '\n'); 
    
end



