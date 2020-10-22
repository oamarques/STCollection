function plot_eig(plot_eig_list)
%***********************************************************************
%*                                                                     *
%* plot_eig reads file names from the file plot_eig_list, reads N      *
%* ======== and eig_1, eig_2, ... eig_N from each file, and plots      *
%*          them in linear and logarithmic scales.                     *
%*                                                                     *
%* Usage:                                                              *
%*    plot_eig(plot_eig_list)                                          *
%* Input arguments:                                                    *
%*    plot_eig_list : file defining a list of files from which N       *
%*                    and EIG will be read from. The end of the list   *
%*                    should be indicated with a STOP or stop. The     *
%*                    entries of the list should look like             *
%*                       T_0005.eig                                    *
%*                       ...                                           *
%*                    so as plot_eig will look for a file T_0005.eig   *
%*                    which should contain                             *
%*                       N                                             *
%*                       EIG(1)                                        *
%*                       ...                                           *
%*                       EIG(N)                                        *
%*                                                                     *
%***********************************************************************

dev = '-djpeg';  % device driver for 'print': -djpeg or -dpsc2
ext = 'jpeg';    % extension for output: jpeg or ps

% open plot_eig_list ...................................................

plot_eig_in = fopen(plot_eig_list,'r');

if ( plot_eig_in == -1 )
   fprintf('File %s cannot be opened!\n',plot_eig_list)
   return
end

for k = 1:1000

%.. read an entry from plot_eig_list ...................................

    [file_eig,count] = fscanf(plot_eig_in,'%s*');

    if ( count==0 | file_eig(1:4)=='STOP' | file_eig(1:4)=='stop' )
       break
    end

%.. read data from file ................................................

    fprintf('.. Processing file %s\n',file_eig)
    x = load(file_eig); n = x(1); eig = x(2:n+1); eig = sort(eig);
    i = regexp(file_eig,'.eig');

%.. plot eigenvalue distribution .......................................

    figure(1)
    plot([1:n],eig,'r*');
    string = (['Eigenvalue distribution: ',file_eig, ...
                sprintf(' (n=%i)',n)]);
    xlabel('index','FontSize',9); ylabel('eigenvalue','FontSize',9); 
    title(string,'Interpreter','none','FontWeight','normal','FontSize',9);
    file_out = [file_eig(1:i-1) '.' ext];
    eval(['print ' dev ' ' file_out]);
 
    eig(eig<eps^2) = eps^2

    figure(2)
    [km] = find(eig<0); 
    [kp] = find(eig>0); 
    [kz] = find(eig==0);
    semilogy(km,abs(eig(km)),'m*',kp,eig(kp),'b*'); 
    if ( size(kz,1)>0 ) 
       hold on 
       ylim = get(gca,'Ylim'); eig(kz)=ylim(1);
       semilogy(kz,eig(kz),'r*')
       hold off
    end
    string = (['Eigenvalue distribution (log scale): ',file_eig, ...
                sprintf(' (n=%i)',n)]);
    xlabel('index','FontSize',9); ylabel('log(abs(eigenvalue))','FontSize',9); 
    title(string,'Interpreter','none','FontWeight','normal','FontSize',9);
    file_out = [file_eig(1:i-1) '_log.' ext];
    eval(['print ' dev ' ' file_out]);

    pause(1)

end

status = fclose('all');

return
