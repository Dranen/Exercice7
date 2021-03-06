function Simulation(name, Ninter, xL, xR, equation, question, u_l, u_r, hocean, xocean, hplage, CFL, tfinal, bord_l, bord_r, bc_l, bc_r, A, omega, mode, nscan, omega_stop, CFL_stop, Ninter_stop, ech_t)

workingfolder = './';
binfilename = 'Exercice7';
input_file = [name '_inp.dat'];

%create the input data file
fid = fopen( [ workingfolder, input_file], 'wt' ); %create or overwrite (empty file, text mode)

%fill the file
fprintf( fid, name );
fprintf( fid, '\n');
fprintf( fid, num2str(Ninter));
fprintf( fid, '\n');
fprintf( fid, num2str(xL));
fprintf( fid, '\n');
fprintf( fid, num2str(xR));
fprintf( fid, '\n');
fprintf( fid, num2str(equation));
fprintf( fid, '\n');
fprintf( fid, num2str(question));
fprintf( fid, '\n');

if question == 0
    fprintf( fid, num2str(u_l));
    fprintf( fid, '\n');
else if question == 1
    fprintf( fid, num2str(u_l));
    fprintf( fid, '\n');
    fprintf( fid, num2str(u_r));
    fprintf( fid, '\n');
else if question == 2
    fprintf( fid, num2str(hocean));
    fprintf( fid, '\n');
    fprintf( fid, num2str(xocean));
    fprintf( fid, '\n');    
    fprintf( fid, num2str(hplage));
    fprintf( fid, '\n');
    end
    end
end

fprintf( fid, num2str(CFL));
fprintf( fid, '\n');
fprintf( fid, num2str(tfinal));
fprintf( fid, '\n');

fprintf( fid, num2str(bc_l));
fprintf( fid, '\n');
if bc_l == 0
    fprintf( fid, num2str(bord_l));
    fprintf( fid, '\n');
else if bc_l == 2
    fprintf( fid, num2str(A));
    fprintf( fid, '\n');
    fprintf( fid, num2str(omega));
    fprintf( fid, '\n');
    end
end

fprintf( fid, num2str(bc_r));
fprintf( fid, '\n');
if bc_r == 0
    fprintf( fid, num2str(bord_r));
    fprintf( fid, '\n');
else if bc_r == 2
    fprintf( fid, num2str(A));
    fprintf( fid, '\n');
    fprintf( fid, num2str(omega));
    fprintf( fid, '\n');
    end
end

fprintf( fid, num2str(mode));
fprintf( fid, '\n');
if mode == 0;
fprintf( fid, num2str(ech_t));
fprintf( fid, '\n');    
else if mode == 1
    fprintf( fid, num2str(nscan));
    fprintf( fid, '\n');
    if nscan > 1
        fprintf( fid, num2str(omega_stop));
        fprintf( fid, '\n');
    end
else if mode == 2
    fprintf( fid, num2str(nscan));
    fprintf( fid, '\n');
    if nscan > 1
        fprintf( fid, num2str(CFL_stop));
        fprintf( fid, '\n');
    end
else if mode == 3
    fprintf( fid, num2str(nscan));
    fprintf( fid, '\n');
    if nscan > 1
        fprintf( fid, num2str(Ninter_stop));
        fprintf( fid, '\n');
    end
    end
    end
    end
end


fclose( fid );


eval( [ '!', workingfolder, binfilename, ' < ', input_file] );

fprintf('\n');
delete(input_file);

end

